import numpy as np
from ase import Atoms
from ase.data import covalent_radii
from typing import List, Union, Tuple

def get_rotation_matrix(axis: np.ndarray, theta: float) -> np.ndarray:
    """计算绕任意轴旋转的旋转矩阵 (Rodrigues)"""
    axis = axis / np.linalg.norm(axis)
    a = np.cos(theta / 2.0)
    b, c, d = -axis * np.sin(theta / 2.0)
    aa, bb, cc, dd = a * a, b * b, c * c, d * d
    bc, ad, ac, ab, bd, cd = b * c, a * d, a * c, a * b, b * d, c * d
    return np.array([
        [aa + bb - cc - dd, 2 * (bc + ad), 2 * (bd - ac)],
        [2 * (bc - ad), aa + cc - bb - dd, 2 * (cd + ab)],
        [2 * (bd + ac), 2 * (cd - ab), aa + dd - bb - cc]
    ])

def get_bond_length(symbol1: str, symbol2: str) -> float:
    """根据共价半径估算键长"""
    r1 = covalent_radii[Atoms(symbol1).numbers[0]]
    r2 = covalent_radii[Atoms(symbol2).numbers[0]]
    return 1.2 * (r1 + r2)

def find_best_atom_position(
    base_atoms: Atoms,
    anchor_info: Union[int, np.integer, Tuple[np.ndarray, str]], # 目标原子 (O)
    primary_axis_point_info: Union[int, np.ndarray],             # 轴向参考点 (T site)
    plane_defining_point_info: Union[int, np.ndarray],           # 平面参考点 (Neighbor T)
    new_atom_symbol: str,                                        # 'H'
    config,                                                      # Config对象
    bond_angle: float,                                           # 角度
    bond_length_override: float = None,                          # 强制键长
    exclude_indices: List[int] = None
) -> np.ndarray:
    """
    全功能的智能原子位置搜索算法。
    通过绕轴旋转，寻找不与现有骨架发生碰撞的最佳位置。
    """
    # 1. 解析输入坐标
    if isinstance(anchor_info, (int, np.integer)):
        anchor_idx = anchor_info
        anchor_pos = base_atoms.positions[anchor_info]
        anchor_symbol = base_atoms[anchor_info].symbol
    else:
        anchor_idx = -1
        anchor_pos, anchor_symbol = anchor_info

    p1_pos = base_atoms.positions[primary_axis_point_info] if isinstance(primary_axis_point_info, (int, np.integer)) else primary_axis_point_info
    p2_pos = base_atoms.positions[plane_defining_point_info] if isinstance(plane_defining_point_info, (int, np.integer)) else plane_defining_point_info

    # 2. 构建几何向量
    # 主轴向量 (从 T 指向 O)
    primary_axis = anchor_pos - p1_pos
    primary_axis_norm = primary_axis / np.linalg.norm(primary_axis)

    # 定义平面的向量 (T1 -> T2)
    plane_defining_vec = p2_pos - p1_pos
    
    # 计算垂直向量 (构建局部坐标系)
    # 这里的逻辑严格遵循原始代码：先尝试与 T1-T2 叉乘，如果不稳健则与固定向量叉乘
    perp_vec = np.cross(primary_axis_norm, plane_defining_vec)
    # 如果平行，使用备用向量 (0.577, 0.577, 0.577) 逻辑
    if np.linalg.norm(perp_vec) < 1e-6:
        # 你的原始代码逻辑：
        perp_vec = np.cross(primary_axis_norm, np.array([0.577, 0.577, 0.577]))
        if np.linalg.norm(perp_vec) < 1e-6:
            perp_vec = np.cross(primary_axis_norm, np.array([-0.577, 0.577, 0.577]))
            
    perp_vec_norm = perp_vec / np.linalg.norm(perp_vec)

    # 3. 初始猜测向量
    dist = bond_length_override if bond_length_override else get_bond_length(anchor_symbol, new_atom_symbol)
    angle_rad = np.deg2rad(180.0 - bond_angle)
    
    # 在平面内计算初始 OH 向量
    initial_vec = (primary_axis_norm * np.cos(angle_rad) + perp_vec_norm * np.sin(angle_rad)) * dist
    
    # 4. 旋转搜索最佳位置
    all_positions = base_atoms.positions
    best_pos = None
    max_min_dist = -1.0
    
    threshold = config.collision_threshold_h if new_atom_symbol == 'H' else config.collision_threshold

    for step in range(config.rotation_steps_h_placement):
        theta = np.deg2rad(step * (360.0 / config.rotation_steps_h_placement))
        rot_matrix = get_rotation_matrix(primary_axis, theta)
        
        candidate_pos = anchor_pos + np.dot(rot_matrix, initial_vec)
        
        # 碰撞检测
        distances = np.linalg.norm(all_positions - candidate_pos, axis=1)
        
        # 排除自身和指定的原子
        indices_to_ignore = [idx for idx in (([anchor_idx] if anchor_idx != -1 else []) + (exclude_indices or [])) if idx < len(distances)]
        if indices_to_ignore:
            distances[indices_to_ignore] = np.inf
            
        min_dist = np.min(distances)
        
        if min_dist > threshold and min_dist > max_min_dist:
            max_min_dist = min_dist
            best_pos = candidate_pos
            
    if best_pos is None:
        print(f"  Warning: 无法找到无碰撞位置 (Anchor: {anchor_idx}), 使用默认位置。")
        return anchor_pos + initial_vec
        
    return best_pos