# zeomod/config.py

import yaml
from dataclasses import dataclass, field, asdict
from typing import List, Dict, Any, Tuple, Optional
from pathlib import Path

@dataclass
class BatchZeoliteConfig:
    """批量生成沸石结构的配置类"""
    
    # --- I/O (可通过命令行覆盖) ---
    input_cif: str = "MFI.cif"
    output_dir: str = "output"

    # --- 几何与超胞 ---
    pore_shape: str = 'cylinder'   # sphere, cylinder, box, ellipsoid
    pore_diameter: float = 10.0
    pore_buffer: float = 3.0
    cylinder_axis: str = 'z'
    pore_dimensions: Optional[List[float]] = None
    pore_center_frac: List[float] = field(default_factory=lambda: [0.5, 0.5, 0.5])
    supercell_dims: Optional[List[int]] = None  # YAML中读取通常是 List

    # --- 化学配比 ---
    si_al_ratio: float = 30.0
    efal_ratio: float = 0.0
    num_al_atoms: Optional[int] = None
    # --- 掺杂策略 ---
    doping_campaign: List[Dict[str, Any]] = field(default_factory=lambda: [
        {'mode': 'pair', 'pair_type': 'NNNN', 'ratio': 0.4}, 
        {'mode': 'auto_single', 'ratio': 'remainder'} 
    ])
    max_doping_trials: int = 100

    # --- 物理参数 ---
    collision_threshold: float = 1.2
    o_h_bond_length_pas: float = 0.97
    t_o_h_angle_pas: float = 109.5
    collision_threshold_h: float = 0.8
    rotation_steps_h_placement: int = 72


    
    cluster: bool = False
    skip_supercell: bool = False  # 跳过扩胞 (假设输入已经是超胞)
    skip_cut: bool = False        # 跳过切孔/清洗/钝化 (假设输入已经有孔或不需要孔)
    skip_doping: bool = False     # 跳过掺杂 (只生成纯硅骨架)
    adsorption_specs: List[Dict[str, Any]] = field(default_factory=list)

    
    # [新增] 全局吸附碰撞阈值 (单位: Å)
    # 默认 1.6 Å，控制吸附质与骨架原子的最小允许距离
    adsorption_cutoff: float = 1.6
    adsorbate_file: Optional[str] = None
    adsorbate_count: int = 0
    target_element: Optional[str] = None  # 对应 --ad-target
    target_radius: float = 3.0            # 对应 --ad-radius [新增]

    @classmethod
    def from_yaml(cls, yaml_path: str):
        """从 YAML 文件加载配置"""
        path = Path(yaml_path)
        if not path.exists():
            raise FileNotFoundError(f"配置文件未找到: {yaml_path}")
        
        with open(path, 'r', encoding='utf-8') as f:
            data = yaml.safe_load(f)
        
        # 过滤掉不在 dataclass 定义中的多余键值，防止报错
        valid_keys = cls.__annotations__.keys()
        filtered_data = {k: v for k, v in data.items() if k in valid_keys}
        
        return cls(**filtered_data)

    def dump_yaml(self, save_path: str):
        """将当前配置保存为 YAML 模板"""
        with open(save_path, 'w', encoding='utf-8') as f:
            yaml.dump(asdict(self), f, sort_keys=False, allow_unicode=True)
