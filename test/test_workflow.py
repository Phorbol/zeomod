# tests/test_workflow.py
import pytest
from zeomod.operations import create_supercell_automatically, create_passivated_mesopore, add_bronsted_sites
from zeomod.topology import GraphProcessor

def test_workflow_skip_logic(dummy_zeolite, base_config):
    """
    模拟 main.py 的流程，测试 skip 参数是否有效。
    测试场景: 跳过扩胞，跳过掺杂，只做切孔。
    """
    # 1. 设置 Flags
    base_config.skip_supercell = True
    base_config.skip_cut = False
    base_config.skip_doping = True
    
    current_atoms = dummy_zeolite.copy()
    original_len = len(current_atoms)
    
    # 2. 模拟 Pipeline
    
    # Stage 1: Supercell
    if not base_config.skip_supercell:
        current_atoms = create_supercell_automatically(current_atoms, base_config)
    # Assert: 应该没变
    assert len(current_atoms) == original_len
    
    # Stage 2: Cut
    if not base_config.skip_cut:
        current_atoms = create_passivated_mesopore(current_atoms, base_config)
    # Assert: 应该变小了 (被切了)
    assert len(current_atoms) < original_len
    len_after_cut = len(current_atoms)
    
    # Stage 3: Doping
    if not base_config.skip_doping:
        graph = GraphProcessor(current_atoms)
        current_atoms = add_bronsted_sites(current_atoms, base_config, graph)
    
    # Assert: 应该没有 Al (跳过了掺杂)
    symbols = current_atoms.get_chemical_symbols()
    assert 'Al' not in symbols
    # 原子数应该只增加了钝化的 H，或者保持不变(如果没切出悬挂键)
    # 这里主要验证没有 Al 产生