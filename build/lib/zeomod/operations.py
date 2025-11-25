import numpy as np
import random
from itertools import combinations
from ase import Atoms
from ase.neighborlist import build_neighbor_list
from matscipy.neighbours import neighbour_list

from .utils import find_best_atom_position

# --- 1. 智能超胞创建 (Proactive V3) ---
def create_supercell_automatically(atoms: Atoms, config) -> Atoms:
    """根据孔道形状(球/方/椭/柱)和缓冲要求，智能计算并创建超晶胞"""
    print(f"正在智能计算超胞倍率 (Proactive V3 - Shape: {config.pore_shape})...")

    cell_vectors = atoms.get_cell()
    volume = atoms.get_volume()
    
    # 1. 计算垂直宽度
    cross_products = [
        np.linalg.norm(np.cross(cell_vectors[1], cell_vectors[2])),
        np.linalg.norm(np.cross(cell_vectors[0], cell_vectors[2])),
        np.linalg.norm(np.cross(cell_vectors[0], cell_vectors[1]))
    ]
    perp_widths = np.array([volume / area if area > 1e-9 else 0.0 for area in cross_products])
    
    # 2. 解析孔道核心尺寸
    shape = config.pore_shape.lower()
    pore_dims = np.zeros(3)
    
    # 尝试读取自定义 dimensions
    custom_dims = config.pore_dimensions

    if shape in ['sphere', 'cube']:
        d = config.pore_diameter
        pore_dims = np.array([d, d, d])
    elif shape in ['box', 'ellipsoid']:
        if custom_dims is not None and len(custom_dims) == 3:
            pore_dims = np.array(custom_dims)
        else:
            d = config.pore_diameter
            pore_dims = np.array([d, d, d])
    elif shape == 'cylinder':
        d = config.pore_diameter
        pore_dims = np.array([d, d, d])
    else:
        raise ValueError(f"不支持的形状: {shape}")

    # 3. 缓冲掩码与轴向修正
    buffer_mask = np.array([1.0, 1.0, 1.0])

    if shape == 'cylinder':
        axis_map = {'x': 0, 'y': 1, 'z': 2}
        axis_idx = axis_map.get(config.cylinder_axis.lower(), 2)
        print(f"  - 圆柱轴向: {config.cylinder_axis} (idx: {axis_idx})")
        # 轴向无缓冲，孔径也不占空间(因为是贯穿的)
        buffer_mask[axis_idx] = 0.0
        pore_dims[axis_idx] = 0.0

    # 4. 计算倍率
    required_total_size = pore_dims + (2 * config.pore_buffer * buffer_mask)
    repeat_factors = np.ceil(np.maximum(required_total_size, 0) / (perp_widths + 1e-9)).astype(int)
    np.maximum(repeat_factors, 1, out=repeat_factors)

    if config.supercell_dims is not None:
        repeat_factors = np.maximum(repeat_factors, np.array(config.supercell_dims)).astype(int)

    print(f"  - 最终扩胞倍率: {repeat_factors}")
    return atoms.repeat(repeat_factors)


# --- 2. 形状感知切孔与钝化 ---
def create_passivated_mesopore(atoms: Atoms, config) -> Atoms:
    """
    全形状支持的介孔切割与钝化。
    """
    structure = atoms.copy()
    print(f"步骤 0: 正在切割介孔 (Shape: {config.pore_shape})...")
    
    center = structure.cell.dot(np.array(config.pore_center_frac))
    positions = structure.get_positions()
    
    # 形状感知切割逻辑
    shape = config.pore_shape.lower()
    mask = np.zeros(len(structure), dtype=bool) # True表示保留，False表示删除(稍后取反)
    
    if shape == 'sphere':
        dists = np.linalg.norm(positions - center, axis=1)
        mask = dists < (config.pore_diameter / 2.0)
        
    elif shape == 'cylinder':
        axis_map = {'x': 0, 'y': 1, 'z': 2}
        axis_idx = axis_map.get(config.cylinder_axis.lower(), 2)
        # 计算到轴线的投影距离
        # 移除轴向坐标后计算模长
        proj_vecs = positions - center
        proj_vecs[:, axis_idx] = 0 
        dists = np.linalg.norm(proj_vecs, axis=1)
        mask = dists < (config.pore_diameter / 2.0)
        
    elif shape in ['box', 'cube']:
        # 简单的 AABB 包围盒检测
        half_dims = np.array([config.pore_diameter]*3) / 2.0
        if shape == 'box' and config.pore_dimensions:
             half_dims = np.array(config.pore_dimensions) / 2.0
        
        rel_pos = np.abs(positions - center)
        # 在所有三个方向都小于半长
        mask = np.all(rel_pos < half_dims, axis=1)

    elif shape == 'ellipsoid':
        radii = np.array([config.pore_diameter]*3) / 2.0
        if config.pore_dimensions:
            radii = np.array(config.pore_dimensions) / 2.0
        
        # 椭球方程: (x/a)^2 + (y/b)^2 + (z/c)^2 < 1
        rel_pos = positions - center
        normalized = (rel_pos / radii) ** 2
        mask = np.sum(normalized, axis=1) < 1.0

    # 执行删除
    del structure[np.where(mask)[0]]
    print(f"切割完成，剩余原子: {len(structure)}")

    # --- 迭代清理 (使用 Matscipy 加速) ---
    print("步骤 1: 迭代清理不饱和位点...")
    for i in range(20):
        # 快速构建邻居
        i_idx, j_idx = neighbour_list('ij', structure, cutoff=2.0)
        coordinations = np.bincount(i_idx, minlength=len(structure))
        symbols = np.array(structure.get_chemical_symbols())
        
        to_del_mask = ((symbols == 'Si') & (coordinations < 4)) | (coordinations == 0)
        to_del_indices = np.where(to_del_mask)[0]

        if len(to_del_indices) == 0:
            print(f"-> 骨架在第 {i+1} 轮稳定。")
            break
        del structure[to_del_indices]
    
    # --- 智能钝化 ---
    print("步骤 2: 智能钝化 (Smart H Placement)...")
    # 重新计算邻居
    i_idx, j_idx = neighbour_list('ij', structure, cutoff=2.0)
    adj = [[] for _ in range(len(structure))]
    for i, j in zip(i_idx, j_idx): adj[i].append(j)
    
    symbols = structure.get_chemical_symbols()
    h_positions = []
    
    for idx, symb in enumerate(symbols):
        if symb == 'O' and len(adj[idx]) == 1:
            # 悬挂氧
            t_idx = adj[idx][0]
            # 寻找 T 的另一个邻居作为平面参考 (plane_defining_point)
            t_neighbors = adj[t_idx]
            other_n = [n for n in t_neighbors if n != idx]
            if not other_n: continue # 极端情况
            plane_ref = other_n[0]
            
            # 调用核心智能搜索
            pos = find_best_atom_position(
                base_atoms=structure,
                anchor_info=idx,
                primary_axis_point_info=t_idx,
                plane_defining_point_info=plane_ref,
                new_atom_symbol='H',
                config=config,
                bond_angle=config.t_o_h_angle_pas,
                bond_length_override=config.o_h_bond_length_pas
            )
            h_positions.append(pos)
            
    if h_positions:
        structure.extend(Atoms('H'*len(h_positions), positions=h_positions))
        print(f"钝化完成，添加了 {len(h_positions)} 个 H 原子。")
        
    return structure


# --- 3. 严格的 Löwenstein 掺杂 ---
def _select_doping_sites_proportional(candidate_pool, num_needed, config, graph):
    """根据 Campaign 比例选择位点"""
    final_sites = set()
    campaign = sorted(config.doping_campaign, key=lambda x: x.get('ratio') == 'remainder')
    
    remaining = num_needed
    
    for job in campaign:
        if remaining <= 0: break
        mode = job['mode']
        ratio = job.get('ratio')
        
        if ratio == 'remainder':
            count = remaining
        else:
            count = int(round(num_needed * ratio))
            
        if count <= 0: continue
        
        available = list(set(candidate_pool) - final_sites)
        
        if mode == 'pair':
            # 寻找配对
            pair_type = job.get('pair_type', 'NNN')
            target_dist = {'NNN': 2, 'NNNN': 3}.get(pair_type, 2)
            
            pairs_to_find = count // 2
            # 简化逻辑：随机采样并验证
            random.shuffle(available)
            # 注意：全组合搜索在大体系太慢，这里采用随机尝试策略
            # 实际应用中可能需要更复杂的图搜索，这里遵循原有逻辑
            for p1, p2 in combinations(available[:min(len(available), 500)], 2): # 限制搜索空间
                if len(final_sites) >= num_needed: break
                if graph.get_distance(p1, p2) == target_dist:
                    # 检查是否与已选位点冲突 (Al-O-Al)
                    if not any(graph.get_distance(p1, x) <= 1 for x in final_sites) and \
                       not any(graph.get_distance(p2, x) <= 1 for x in final_sites):
                        final_sites.add(p1)
                        final_sites.add(p2)
                        
        elif mode == 'auto_single':
            random.shuffle(available)
            added = 0
            for site in available:
                if added >= count: break
                if not any(graph.get_distance(site, x) <= 1 for x in final_sites):
                    final_sites.add(site)
                    added += 1
        
        remaining = num_needed - len(final_sites)
        
    return list(final_sites)

def add_bronsted_sites(atoms: Atoms, config, graph) -> Atoms:
    """创建 B 酸位点，支持 Ratio 模式和 Fixed Count 模式"""
    structure = atoms.copy()
    
    # 1. 计算目标 Al 总数量 (Total Al = Framework Al + EFAl)
    t_sites = [a for a in structure if a.symbol in ('Si', 'Al')]
    
    # === [核心修改逻辑开始] ===
    total_al = 0
    
    # 模式 A: 显式指定数量 (优先级高)
    if config.num_al_atoms is not None and config.num_al_atoms > 0:
        total_al = config.num_al_atoms
        print(f"\n--- 掺杂模式: 固定数量 (Fixed Count) ---")
        print(f"  - 指定总 Al 数量: {total_al}")
        
    # 模式 B: 使用硅铝比计算
    elif config.si_al_ratio > 0:
        # Formula: N_Al = N_Total / (Ratio + 1)
        total_al = int(round(len(t_sites) / (config.si_al_ratio + 1)))
        print(f"\n--- 掺杂模式: 硅铝比 (Si/Al Ratio) ---")
        print(f"  - 设定 Ratio: {config.si_al_ratio}")
        print(f"  - 计算总 Al 数量: {total_al}")
    
    if total_al <= 0:
        print("  - 警告: 计算出的 Al 数量为 0，跳过掺杂。")
        return structure
    # === [核心修改逻辑结束] ===

    # 计算骨架 Al (B酸) 和 非骨架 Al (EFAl) 的分配
    num_efal = int(round(total_al * config.efal_ratio))
    bronsted_al = total_al - num_efal
    
    print(f"  - 骨架 Al (B酸) 目标: {bronsted_al}")
    print(f"  - 非骨架 Al (EFAl) 目标: {num_efal} (暂未实现EFAl结构生成)") # 你的代码暂时主要处理B酸
    
    if bronsted_al <= 0: return structure
    print(f"\n--- 掺杂 B 酸位点 (目标: {bronsted_al}) ---")
    
    # 2. 筛选表面 Si 候选池 (连接有 H 的)
    nl = build_neighbor_list(structure, bothways=True, self_interaction=False)
    symbols = np.array(structure.get_chemical_symbols())
    si_indices = [a.index for a in structure if a.symbol == 'Si']
    
    surface_si = []
    for si in si_indices:
        # 检查邻居O是否有H邻居
        is_surface = False
        for o_idx in nl.get_neighbors(si)[0]:
            if symbols[o_idx] == 'O':
                for h_idx in nl.get_neighbors(o_idx)[0]:
                    if symbols[h_idx] == 'H':
                        is_surface = True; break
            if is_surface: break
        if is_surface: surface_si.append(si)
        
    pool = surface_si if len(surface_si) >= bronsted_al else si_indices
    
    # 3. 循环重试机制 (Löwenstein Rule)
    sites_to_replace = []
    for trial in range(config.max_doping_trials):
        candidates = _select_doping_sites_proportional(pool, bronsted_al, config, graph)
        
        # 最终验证
        valid = True
        if len(candidates) < bronsted_al:
            valid = False
        else:
            for s1, s2 in combinations(candidates, 2):
                if graph.get_distance(s1, s2) <= 1:
                    valid = False; break
        
        if valid:
            sites_to_replace = candidates
            print(f"  -> 尝试 {trial+1}: 成功找到合法位点。")
            break
    else:
        raise RuntimeError("达到最大重试次数，无法满足 Löwenstein 规则。")
        
    # 4. 执行替换与加质子
    for idx in sites_to_replace:
        structure[idx].symbol = 'Al'
        
    # 添加质子
    h_to_add = []
    # 重新构建 neighbor list 因为符号变了
    nl = build_neighbor_list(structure, bothways=True, self_interaction=False)
    
    for al_idx in sites_to_replace:
        # 找桥氧
        candidates = []
        for o_idx in nl.get_neighbors(al_idx)[0]:
            if structure[o_idx].symbol == 'O':
                neighbors = nl.get_neighbors(o_idx)[0]
                # 必须连接两个 T 原子 (Al-O-Si)
                t_neighbors = [n for n in neighbors if structure[n].symbol in ('Si', 'Al')]
                if len(t_neighbors) == 2:
                    candidates.append(o_idx)
        
        if not candidates: continue
        target_o = random.choice(candidates)
        # 锚点计算
        t_neighbor = [n for n in nl.get_neighbors(target_o)[0] if n != al_idx][0]
        
        pos = find_best_atom_position(
            base_atoms=structure,
            anchor_info=target_o,
            primary_axis_point_info=al_idx,
            plane_defining_point_info=t_neighbor,
            new_atom_symbol='H',
            config=config,
            bond_angle=109.5, # 或者是 config.t_o_h_angle_pas
            bond_length_override=0.97
        )
        h_to_add.append(pos)

    if h_to_add:
        structure.extend(Atoms('H'*len(h_to_add), positions=h_to_add))
        
    return structure
