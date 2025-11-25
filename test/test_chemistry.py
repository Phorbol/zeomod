# tests/test_chemistry.py
import pytest
from collections import Counter
from zeomod.operations import add_bronsted_sites
from zeomod.topology import GraphProcessor

def count_elements(atoms):
    return Counter(atoms.get_chemical_symbols())

# 参数化测试：测试两种掺杂模式
@pytest.mark.parametrize("mode", ["ratio", "fixed"])
def test_doping_modes(dummy_zeolite, base_config, mode):
    """测试 Si/Al 比模式和固定数量模式"""
    
    # 初始化图处理器 (在未掺杂结构上)
    graph = GraphProcessor(dummy_zeolite)
    
    target_al_count = 0
    
    if mode == "ratio":
        base_config.si_al_ratio = 10.0 # 约 10% 的 T 位点
        base_config.num_al_atoms = None
        # 估算预期数量: Total_T / (10+1)
        n_t = len([a for a in dummy_zeolite if a.symbol == 'Si'])
        expected_al = int(round(n_t / 11.0))
        target_al_count = expected_al
        
    elif mode == "fixed":
        base_config.num_al_atoms = 5
        base_config.si_al_ratio = 1000 # 这个应该被忽略
        target_al_count = 5

    # 执行掺杂
    doped_atoms = add_bronsted_sites(dummy_zeolite, base_config, graph)
    counts = count_elements(doped_atoms)
    
    # 验证 Al 数量
    assert counts['Al'] <= target_al_count + 1
    
    # 验证电荷平衡 (Al数量应等于引入的 H 数量，假设无 EFAl)
    assert counts['Al'] == counts['H']

def test_lowenstein_rule(dummy_zeolite, base_config):
    """验证生成的结构不存在 Al-O-Al"""
    base_config.num_al_atoms = 10 # 强制塞入较多 Al 增加冲突概率
    base_config.max_doping_trials = 50
    
    graph = GraphProcessor(dummy_zeolite)
    doped_atoms = add_bronsted_sites(dummy_zeolite, base_config, graph)
    
    # 检查结果
    al_indices = [a.index for a in doped_atoms if a.symbol == 'Al']
    
    # 遍历所有 Al 对
    for i in al_indices:
        for j in al_indices:
            if i >= j: continue
            # 拓扑距离
            dist = graph.get_distance(i, j)
            # 距离 1 代表 T-O-T 直接连接
            assert dist > 1, f"Found Al-O-Al connection between {i} and {j}"

def test_impossible_doping_loading(dummy_zeolite, base_config):
    """
    鲁棒性测试：尝试进行物理上不可能的高密度掺杂 (Si/Al < 1)。
    验证程序是否能正确识别并抛出 RuntimeError (Löwenstein规则限制)，而不是死循环或崩溃。
    """
    graph = GraphProcessor(dummy_zeolite)
    
    # 设置极端的掺杂要求
    # dummy体系约72个T位点，要求掺杂60个Al，几乎肯定违反Löwenstein规则
    base_config.num_al_atoms = 60 
    base_config.max_doping_trials = 10 # 限制重试次数，防止测试卡死

    # 断言：这里必须抛出 RuntimeError，且错误信息包含特定关键词
    with pytest.raises(RuntimeError, match=r".*Löwenstein.*"):
        add_bronsted_sites(dummy_zeolite, base_config, graph)