# tests/test_config.py
import pytest
import os
from zeomod.config import BatchZeoliteConfig
from zeomod.main import parse_override_value

def test_parse_override_values():
    """测试命令行参数的智能解析逻辑"""
    assert parse_override_value("10") == 10
    assert parse_override_value("10.5") == 10.5
    assert parse_override_value("True") is True
    assert parse_override_value("[1, 2, 3]") == [1, 2, 3]
    assert parse_override_value("cylinder") == "cylinder"

def test_config_yaml_loading(tmp_path):
    """测试 YAML 读取与默认值回退"""
    # 创建临时 YAML
    yaml_file = tmp_path / "test.yaml"
    yaml_content = """
    pore_diameter: 15.0
    doping_campaign:
      - mode: auto_single
        ratio: remainder
    """
    yaml_file.write_text(yaml_content, encoding='utf-8')
    
    config = BatchZeoliteConfig.from_yaml(str(yaml_file))
    assert config.pore_diameter == 15.0
    assert config.si_al_ratio == 30.0 #(检查默认值保持不变)

def test_invalid_yaml_key(tmp_path):
    """测试 YAML 中包含无效键时是否健壮"""
    yaml_file = tmp_path / "bad.yaml"
    yaml_file.write_text("invalid_key: 123\npore_diameter: 8.0", encoding='utf-8')
    
    # 应该忽略 invalid_key 而不报错
    config = BatchZeoliteConfig.from_yaml(str(yaml_file))
    assert config.pore_diameter == 8.0
    assert not hasattr(config, 'invalid_key')
