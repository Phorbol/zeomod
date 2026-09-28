import argparse
import sys
import ast
from pathlib import Path
from ase.io import read, write

# 导入模块
from .config import BatchZeoliteConfig
from .topology import GraphProcessor
from .operations import (
    create_supercell_automatically, 
    create_passivated_mesopore, 
    add_bronsted_sites
)
from .validators import validate_pure_silica, validate_bronsted

def parse_override_value(value_str):
    """
    智能解析命令行传入的字符串值
    '10' -> 10 (int)
    '10.5' -> 10.5 (float)
    'True' -> True (bool)
    '[1,2,3]' -> [1, 2, 3] (list)
    'cylinder' -> 'cylinder' (str)
    """
    try:
        # 尝试作为 Python 字面量解析 (处理数字, bool, list, dict)
        return ast.literal_eval(value_str)
    except (ValueError, SyntaxError):
        # 如果解析失败，说明它只是一个普通字符串
        return value_str

def main():
    parser = argparse.ArgumentParser(description="ZeoMod CLI: 分步执行与流程控制")
    
    # --- 1. 基础配置 ---
    parser.add_argument('config_file', type=str, nargs='?', default='input.yaml',
                        help='YAML 配置文件路径')
    parser.add_argument('--generate-template', action='store_true',
                        help='生成默认配置文件模板并退出')

    # --- 2. 高频参数覆盖 (Explicit Flags) ---
    group_io = parser.add_argument_group('I/O Overrides')
    group_io.add_argument('--input', '-i', type=str, help='覆盖输入 CIF 文件')
    group_io.add_argument('--output', '-o', type=str, help='覆盖输出目录')

    group_geom = parser.add_argument_group('Geometry Overrides')
    group_geom.add_argument('--shape', '-s', type=str, choices=['sphere', 'cylinder', 'box', 'ellipsoid'],
                            help='孔道形状')
    group_geom.add_argument('--diameter', '-d', type=float, help='孔道主要直径 (Å)')
    group_geom.add_argument('--buffer', '-b', type=float, help='扩胞缓冲距离 (Å)')
    group_geom.add_argument('--axis', type=str, choices=['x', 'y', 'z'], help='圆柱轴向')

    group_chem = parser.add_argument_group('Chemistry Overrides')
    group_chem.add_argument('--si-al', '-r', type=float, help='硅铝比 (Si/Al Ratio)')
    group_chem.add_argument('--num-al', '-n', type=int, help='直接指定 Al 原子的总数量 (覆盖硅铝比)')  # [新增]
    group_chem.add_argument(
        '--doping-sites', nargs='+', type=int, metavar='ATOM_INDEX',
        help='指定要替换为 Al 的 Si 原子索引（从 0 开始）'
    )
    group_chem.add_argument('--efal', type=float, help='非骨架铝比例 (EFAl Ratio)')

    # --- 3. 万能覆盖参数 (Dynamic Overrides) ---
    parser.add_argument('--set', nargs='*', metavar='KEY=VALUE',
                        help='高级覆盖: 修改任意配置参数 (例如: --set collision_threshold=1.5 max_doping_trials=20)')

    group_flow = parser.add_argument_group('Workflow Control')
    group_flow.add_argument('--skip-supercell', action='store_true', help='跳过自动扩胞')
    group_flow.add_argument('--skip-cut', action='store_true', help='跳过切孔与钝化')
    group_flow.add_argument('--skip-doping', action='store_true', help='跳过Al掺杂')

    args = parser.parse_args()

    # --- 生成模板 ---
    if args.generate_template:
        BatchZeoliteConfig().dump_yaml('template_config.yaml')
        print("✅ 已生成默认模板: template_config.yaml")
        sys.exit(0)

    # --- 步骤 1: 加载 YAML ---
    print(f"Reading configuration from: {args.config_file}")
    try:
        if Path(args.config_file).exists():
            config = BatchZeoliteConfig.from_yaml(args.config_file)
        else:
            print(f"⚠️ Warning: 配置文件 {args.config_file} 不存在，将使用默认参数。")
            config = BatchZeoliteConfig()
    except Exception as e:
        print(f"❌ YAML 加载失败: {e}")
        sys.exit(1)

    # --- 步骤 2: 应用显式 CLI 参数 (Explicit Overrides) ---
    # 仅当用户在命令行提供了该参数时才覆盖
    if args.input: config.input_cif = args.input
    if args.output: config.output_dir = args.output
    if args.shape: config.pore_shape = args.shape
    if args.diameter: config.pore_diameter = args.diameter
    if args.buffer: config.pore_buffer = args.buffer
    if args.axis: config.cylinder_axis = args.axis
    if args.si_al: config.si_al_ratio = args.si_al
    if args.efal: config.efal_ratio = args.efal
    if args.num_al: config.num_al_atoms = args.num_al  # [新增] 同步到 config
    if args.doping_sites is not None: config.doping_sites = args.doping_sites

    # --- 步骤 3: 应用万能覆盖参数 (Generic Overrides) ---
    # 处理 --set key=value 列表
    if args.set:
        print("Processing generic overrides:")
        for item in args.set:
            if '=' not in item:
                print(f"  ⚠️ 忽略无效格式: {item} (应为 key=value)")
                continue
            
            key, val_str = item.split('=', 1)
            
            # 检查 key 是否有效
            if not hasattr(config, key):
                print(f"  ⚠️ 警告: 配置中不存在参数 '{key}'，将被忽略。")
                continue
                
            # 智能转换类型
            val = parse_override_value(val_str)
            
            # 执行覆盖
            setattr(config, key, val)
            print(f"  -> Set '{key}' = {val} ({type(val).__name__})")

    if args.skip_supercell: config.skip_supercell = True
    if args.skip_cut: config.skip_cut = True
    if args.skip_doping: config.skip_doping = True

    # --- 准备工作目录 ---
    Path(config.output_dir).mkdir(parents=True, exist_ok=True)

    # --- 步骤 4: 执行核心流程 (与之前相同) ---
    try:
        # 1. Loading
        print(f"\n>>> [1/4] Loading Structure: {config.input_cif}")
        atoms = read(config.input_cif)
        current_atoms = atoms # 使用一个变量跟踪当前结构

        # 2. Supercell
        if not config.skip_supercell:
            print(f"\n>>> [2/4] Creating Supercell (Proactive V3)...")
            current_atoms = create_supercell_automatically(current_atoms, config)
        else:
            print(f"\n>>> [2/4] Skipping Supercell (User Request).")

        # 3. Pore Cutting
        if not config.skip_cut:
            print(f"\n>>> [3/4] Cutting Pore & Passivating...")
            current_atoms = create_passivated_mesopore(current_atoms, config)
            validate_pure_silica(current_atoms)
        else:
            print(f"\n>>> [3/4] Skipping Cutting & Passivation (User Request).")

        # 4. Doping
        if not config.skip_doping:
            # 只有当 Si/Al 比合理时才执行，且用户没跳过
            need_doping = config.doping_sites is not None or \
                      (config.num_al_atoms is not None and config.num_al_atoms > 0) or \
                      (config.si_al_ratio > 0 and config.si_al_ratio < 10000)

            if need_doping:
                print(f"\n>>> [4/4] Doping Al & Adding Protons...")
                graph = GraphProcessor(current_atoms)
                current_atoms = add_bronsted_sites(current_atoms, config, graph)
                validate_bronsted(current_atoms, graph)
            else:
                print("    Skipping doping (No Al requested via Ratio or Num-Al).")
        else:
            print(f"\n>>> [4/4] Skipping Doping (User Request).")
        # 5. Save
        # 自动生成一个带有状态标识的文件名
        status_tag = ""
        if config.skip_cut: status_tag += "_NoCut"
        if config.skip_doping: status_tag += "_PureSi"
        
        save_name = f"result_{config.pore_shape}{status_tag}.cif"
        save_path = Path(config.output_dir) / save_name
        write(save_path, current_atoms)
        
        config.dump_yaml(Path(config.output_dir) / "run_config_snapshot.yaml")
        print(f"\n✅ Success! Structure saved to: {save_path}")

    except Exception as e:
        import traceback
        traceback.print_exc()
        print(f"\n❌ Error: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
