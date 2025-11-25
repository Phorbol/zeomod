# tests/test_geometry.py
import pytest
import numpy as np
from zeomod.operations import create_supercell_automatically, create_passivated_mesopore

# 移除了 Cluster 相关的参数化，只测试不同的形状
@pytest.mark.parametrize("shape", ["sphere", "cylinder", "box"])
def test_pore_cutting_modes(dummy_zeolite, base_config, shape):
    """
    测试不同形状下，切介孔 (Periodic Mesopore) 的逻辑是否正确。
    仅验证周期性结构及其切割效果。
    """
    # 1. 配置参数
    base_config.pore_shape = shape
    base_config.cluster = False  # 明确指定非 Cluster 模式
    base_config.pore_diameter = 10.0
    
    # 为了测试 cylinder，指定一个轴向
    if shape == 'cylinder':
        base_config.cylinder_axis = 'z'
    
    # 2. 执行切割
    # 注意：我们先不做扩胞，直接在 dummy_zeolite 上切，为了测试掩码逻辑
    processed_atoms = create_passivated_mesopore(dummy_zeolite, base_config)
    
    original_count = len(dummy_zeolite)
    new_count = len(processed_atoms)
    
    # 3. 验证原子数量确实减少了 (发生了切割)
    assert new_count < original_count, f"Shape {shape} did not cut any atoms!"
    
    # 4. 验证介孔模式的核心属性：周期性必须保持
    # 如果变成了 Cluster，这里会是 [False, False, False]
    assert all(processed_atoms.pbc), "Mesopore mode must maintain Periodic Boundary Conditions (PBC)."

def test_supercell_calculation(dummy_zeolite, base_config):
    """测试 Proactive V3 智能扩胞"""
    base_config.pore_diameter = 20.0 # 设置一个很大的孔
    base_config.pore_buffer = 5.0
    base_config.pore_shape = 'sphere'
    
    # 原始尺寸约 21Å (3x7.16)
    # 需求: 20 + 10 = 30Å
    # 预期扩胞: 至少需要 2x2x2 的倍率来包裹
    
    supercell = create_supercell_automatically(dummy_zeolite, base_config)
    
    vol_orig = dummy_zeolite.get_volume()
    vol_super = supercell.get_volume()
    
    assert vol_super > vol_orig * 1.5 # 肯定变大了
    assert len(supercell) > len(dummy_zeolite)
