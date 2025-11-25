# tests/conftest.py
import pytest
import numpy as np
from ase.build import bulk
from zeomod.config import BatchZeoliteConfig
from ase.io import read
@pytest.fixture
def dummy_zeolite():
    """
    创建一个简单的周期性 SiO2 结构作为测试对象。
    使用方石英 (Cristobalite) 扩胞，模拟一个 mini 沸石。
    """
    # 创建 SiO2 晶胞
    #atoms = bulk('SiO2' , cubic=True)
    atoms = read('./test/MFI.cif')
    # 扩胞以容纳操作 (4x4x4 约 28Å 边长，足够切孔)
    #atoms = atoms.repeat((3, 3, 3)) 
    return atoms

@pytest.fixture
def base_config():
    """返回一个默认配置对象"""
    return BatchZeoliteConfig(
        input_cif="dummy.cif",
        output_dir="test_output"
    )
