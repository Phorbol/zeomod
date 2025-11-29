import numpy as np
from ase import Atoms
from ase.io import read
from scipy.spatial import cKDTree
from tqdm import tqdm
from typing import List, Union, Optional, Dict, Any

# ==========================================
# 1. 数学工具: 完美均匀随机旋转
# ==========================================
def get_uniform_random_rotation_matrix() -> np.ndarray:
    """
    生成一个严格符合 Haar Measure 的均匀随机旋转矩阵 (3x3)。
    算法来源: Fast Random Rotation Matrices (James Arvo, 1992)
    这比欧拉角采样更均匀，且不会产生万向节死锁或极点聚集。
    """
    x1, x2, x3 = np.random.random(3)
    
    # 旋转矩阵 R (绕 Z 轴)
    R = np.array([
        [np.cos(2 * np.pi * x1), np.sin(2 * np.pi * x1), 0],
        [-np.sin(2 * np.pi * x1), np.cos(2 * np.pi * x1), 0],
        [0, 0, 1]
    ])
    
    # Householder 矩阵 V
    sqrt_x3 = np.sqrt(x3)
    v = np.array([
        np.cos(2 * np.pi * x2) * sqrt_x3,
        np.sin(2 * np.pi * x2) * sqrt_x3,
        np.sqrt(1 - x3)
    ])
    H = np.eye(3) - 2 * np.outer(v, v)
    
    # 最终矩阵 M = -H @ R
    return -np.dot(H, R)

def apply_random_rotation(atoms: Atoms):
    """对分子应用均匀随机旋转"""
    matrix = get_uniform_random_rotation_matrix()
    com = atoms.get_center_of_mass()
    atoms.translate(-com) # 归零
    # 坐标变换: new = old @ matrix.T
    atoms.set_positions(np.dot(atoms.get_positions(), matrix.T))
    atoms.translate(com)  # 复位

# ==========================================
# 2. 物理引擎: 空间碰撞检测
# ==========================================
class SmartCollider:
    """基于 cKDTree 的高性能碰撞检测器 (支持 PBC)"""
    def __init__(self, host_atoms: Atoms, cutoff: float = 1.6):
        self.cell = host_atoms.get_cell()
        self.pbc = host_atoms.get_pbc()
        self.cutoff = cutoff
        self.host_pos = host_atoms.get_positions()
        self.tree = self._build_pbc_tree(self.host_pos)

    def _build_pbc_tree(self, positions):
        """构建 3x3x3 超胞树以处理周期性边界"""
        if not np.any(self.pbc):
            return cKDTree(positions)
        
        # 简单的全镜像构建 (对于复杂晶胞最稳健)
        images = []
        for x in [-1, 0, 1]:
            for y in [-1, 0, 1]:
                for z in [-1, 0, 1]:
                    if x==0 and y==0 and z==0:
                        images.append(positions)
                    else:
                        disp = np.dot([x, y, z], self.cell)
                        images.append(positions + disp)
        return cKDTree(np.vstack(images))

    def is_colliding(self, guest_atoms: Atoms, strictness: float = 1.0) -> bool:
        eff_cut = self.cutoff * strictness
        g_pos = guest_atoms.get_positions()
        dists, _ = self.tree.query(g_pos, k=1)
        if np.min(dists) < eff_cut:
            return True
        return False

    def update_host(self, new_atoms: Atoms):
        """插入成功后更新环境"""
        new_pos = new_atoms.get_positions()
        self.host_pos = np.vstack([self.host_pos, new_pos])
        self.tree = self._build_pbc_tree(self.host_pos)

# ==========================================
# 3. 业务逻辑: 吸附引擎
# ==========================================
class AdsorptionEngine:
    def __init__(self, host_atoms: Atoms, seed: int = None):
        self.host_atoms = host_atoms.copy()
        self.collider = SmartCollider(self.host_atoms)
        self.cell = host_atoms.get_cell()
        self.inv_cell = np.linalg.inv(self.cell)
        self.rng = np.random.RandomState(seed)

    def _resolve_targets(self, target_input) -> Optional[np.ndarray]:
        """解析目标：支持 元素(str), 索引(int list), 坐标(float list)"""
        if target_input is None: return None
        
        # 1. 元素符号 "Al"
        if isinstance(target_input, str):
            indices = [a.index for a in self.host_atoms if a.symbol == target_input]
            if not indices:
                print(f"    ⚠️ Warning: Target element '{target_input}' not found. Switching to random.")
                return None
            return self.host_atoms.positions[indices]

        # 2. 列表/数组
        if isinstance(target_input, (list, tuple, np.ndarray)):
            arr = np.array(target_input)
            
            # 单个坐标 [10.5, 10.5, 5.0] (float)
            if arr.ndim == 1 and len(arr) == 3 and np.issubdtype(arr.dtype, np.floating):
                return np.array([arr])
            
            # 多个坐标 [[1.1, 1.1, 1.1], [2.2, 2.2, 2.2]]
            if arr.ndim == 2 and arr.shape[1] == 3:
                return arr
                
            # 索引列表 [1, 5, 10] (int)
            if arr.ndim == 1 and np.issubdtype(arr.dtype, np.integer):
                # 过滤越界索引
                valid_idx = arr[arr < len(self.host_atoms)]
                return self.host_atoms.positions[valid_idx]

        return None

    def _get_valid_seed(self, target_pos_list: np.ndarray, radius: float) -> Optional[np.ndarray]:
        """在目标周围采样合法的空隙中心"""
        for _ in range(100): # 种子采样最大尝试次数
            if target_pos_list is not None and len(target_pos_list) > 0:
                # 定向模式：随机选一个目标点 -> 球壳偏移
                center = target_pos_list[self.rng.randint(len(target_pos_list))]
                vec = self.rng.normal(0, 1, 3)
                vec /= np.linalg.norm(vec)
                r = self.rng.uniform(2.0, radius)
                cart_pos = center + vec * r
            else:
                # 全局随机模式
                frac = self.rng.uniform(0, 1, 3)
                cart_pos = np.dot(frac, self.cell)

            # 快速预检: 种子点是否直接在墙里? (dist < 1.8A)
            d, _ = self.collider.tree.query(cart_pos, k=1)
            if d > 1.8:
                # Wrap to unit cell
                f = np.dot(cart_pos, self.inv_cell)
                f -= np.floor(f)
                return np.dot(f, self.cell)
        return None

    def insert_batch(self, guest_atoms: Atoms, count: int, 
                    target_spec: Any, target_radius: float, name: str) -> int:
        """执行单一批次的插入"""
        g_temp = guest_atoms.copy()
        g_temp.center()
        g_temp.translate(-g_temp.get_center_of_mass()) # 归零
        
        targets = self._resolve_targets(target_spec)
        target_str = str(target_spec) if not isinstance(target_spec, (list, np.ndarray)) else "Coords/Indices"
        desc_str = f"Target: {target_str}" if targets is not None else "Target: Random"

        success = 0
        pbar = tqdm(total=count, desc=f"  -> {name} ({desc_str})", unit="mol", leave=False)
        
        for _ in range(count):
            inserted = False
            for _ in range(500): # 单个分子的最大尝试次数
                seed = self._get_valid_seed(targets, target_radius)
                if seed is None: continue
                
                g_trial = g_temp.copy()
                apply_random_rotation(g_trial) # 使用完美随机旋转
                g_trial.translate(seed)
                
                if not self.collider.is_colliding(g_trial):
                    self.host_atoms.extend(g_trial)
                    self.collider.update_host(g_trial) # 更新环境
                    success += 1
                    inserted = True
                    pbar.update(1)
                    break
        pbar.close()
        return success

def perform_adsorption(host_atoms: Atoms, config) -> Atoms:
    """ZeoMod 集成入口"""
    if not config.adsorption_specs: return host_atoms
    
    print(f"\n>>> [5/5] Adsorption Module (KDTree Accelerated)...")
    engine = AdsorptionEngine(host_atoms)
    
    total = 0
    for i, spec in enumerate(config.adsorption_specs):
        fpath = spec.get('file')
        count = spec.get('count', 1)
        target = spec.get('target', None)
        radius = spec.get('radius', 5.0)
        name = spec.get('name', f"Batch_{i+1}")
        
        try:
            guest = read(fpath)
            n = engine.insert_batch(guest, count, target, radius, name)
            total += n
        except Exception as e:
            print(f"  ❌ Error loading/inserting {fpath}: {e}")
            
    print(f"✅ Adsorption complete. Total inserted: {total}")
    return engine.host_atoms