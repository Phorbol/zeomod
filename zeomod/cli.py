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
    add_bronsted_sites,
)
from .validators import validate_pure_silica, validate_bronsted
from .adsorption import perform_adsorption
def parse_override_value(value_str):
    try:
        return ast.literal_eval(value_str)
    except (ValueError, SyntaxError):
        return value_str

def main():
    parser = argparse.ArgumentParser(description="ZeoMod CLI: 分步执行与流程控制")
    
    # --- 1. 基础配置 ---
    parser.add_argument('config_file', type=str, nargs='?', default='input.yaml',
                        help='YAML 配置文件路径')
    parser.add_argument('--generate-template', action='store_true',
                        help='生成默认配置文件模板并退出')

    # --- 2. 高频参数覆盖 ---
    group_io = parser.add_argument_group('I/O Overrides')
    group_io.add_argument('--input', '-i', type=str, help='覆盖输入 CIF 文件')
    group_io.add_argument('--output', '-o', type=str, help='覆盖输出目录')

    group_geom = parser.add_argument_group('Geometry Overrides')
    group_geom.add_argument('--shape', '-s', type=str, choices=['sphere', 'cylinder', 'box', 'ellipsoid'], help='孔道形状')
    group_geom.add_argument('--diameter', '-d', type=float, help='孔道主要直径 (Å)')
    group_geom.add_argument('--buffer', '-b', type=float, help='扩胞缓冲距离 (Å)')
    group_geom.add_argument('--axis', type=str, choices=['x', 'y', 'z'], help='圆柱轴向')

    group_chem = parser.add_argument_group('Chemistry Overrides')
    group_chem.add_argument('--si-al', '-r', type=float, help='硅铝比')
    group_chem.add_argument('--num-al', '-n', type=int, help='固定 Al 原子数量')
    group_chem.add_argument(
        '--doping-sites', nargs='+', type=int, metavar='ATOM_INDEX',
        help='指定要替换为 Al 的 Si 原子索引（从 0 开始）'
    )
    group_chem.add_argument('--efal', type=float, help='非骨架铝比例')

    # --- 3. 吸附参数 (Adsorption Overrides) ---
    group_ad = parser.add_argument_group('Adsorption Overrides')
    group_ad.add_argument('--ad-file', action='append', type=str, help='吸附分子文件 (可多次指定)')
    group_ad.add_argument('--ad-count', action='append', type=int, help='对应吸附数量')
    group_ad.add_argument('--ad-target', action='append', type=str, help='对应目标元素')
    group_ad.add_argument('--ad-radius', action='append', type=float, help='对应吸附半径')

    # --- 4. 流程控制 ---
    parser.add_argument('--set', nargs='*', metavar='KEY=VALUE', help='高级参数覆盖')
    group_flow = parser.add_argument_group('Workflow Control')
    group_flow.add_argument('--skip-supercell', action='store_true', help='跳过自动扩胞')
    group_flow.add_argument('--skip-cut', action='store_true', help='跳过切孔')
    group_flow.add_argument('--skip-doping', action='store_true', help='跳过掺杂')

    args = parser.parse_args()

    # --- 模板生成 ---
    if args.generate_template:
        BatchZeoliteConfig().dump_yaml('template_config.yaml')
        print("✅ 已生成默认模板: template_config.yaml")
        sys.exit(0)

    # --- 加载配置 ---
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

    # --- 应用参数覆盖 ---
    if args.input: config.input_cif = args.input
    if args.output: config.output_dir = args.output
    if args.shape: config.pore_shape = args.shape
    if args.diameter: config.pore_diameter = args.diameter
    if args.buffer: config.pore_buffer = args.buffer
    if args.axis: config.cylinder_axis = args.axis
    if args.si_al: config.si_al_ratio = args.si_al
    if args.efal: config.efal_ratio = args.efal
    if args.num_al: config.num_al_atoms = args.num_al
    if args.doping_sites is not None: config.doping_sites = args.doping_sites

    if args.set:
        print("Processing generic overrides:")
        for item in args.set:
            if '=' not in item: continue
            key, val_str = item.split('=', 1)
            if hasattr(config, key):
                setattr(config, key, parse_override_value(val_str))

    if args.skip_supercell: config.skip_supercell = True
    if args.skip_cut: config.skip_cut = True
    if args.skip_doping: config.skip_doping = True

    # --- 准备目录 ---
    Path(config.output_dir).mkdir(parents=True, exist_ok=True)

    # --- 核心流程 ---
    try:
        # [0] 初始化 status_tag (防止 UnboundLocalError)
        status_tag = ""
        if config.skip_cut: status_tag += "_NoCut"
        if config.skip_doping: status_tag += "_PureSi"

        # [1] Load
        print(f"\n>>> [1/5] Loading Structure: {config.input_cif}")
        atoms = read(config.input_cif)
        current_atoms = atoms

        # [2] Supercell
        if not config.skip_supercell:
            print(f"\n>>> [2/5] Creating Supercell...")
            current_atoms = create_supercell_automatically(current_atoms, config)
        else:
            print(f"\n>>> [2/5] Skipping Supercell.")

        # [3] Cut
        if not config.skip_cut:
            print(f"\n>>> [3/5] Cutting Pore...")
            current_atoms = create_passivated_mesopore(current_atoms, config)
            validate_pure_silica(current_atoms)
        else:
            print(f"\n>>> [3/5] Skipping Cutting.")

        # [4] Doping
        if not config.skip_doping:
            need_doping = config.doping_sites is not None or \
                          (config.num_al_atoms and config.num_al_atoms > 0) or \
                          (config.si_al_ratio > 0 and config.si_al_ratio < 10000)
            if need_doping:
                print(f"\n>>> [4/5] Doping Al...")
                graph = GraphProcessor(current_atoms)
                current_atoms = add_bronsted_sites(current_atoms, config, graph)
                validate_bronsted(current_atoms, graph)
            else:
                print("    Skipping doping (No Al requested).")
        else:
            print(f"\n>>> [4/5] Skipping Doping.")

        # [5] Adsorption (修复版逻辑)
        # 优先构建任务清单
        if args.ad_file:
            specs = []
            files = args.ad_file
            n = len(files)
            # 补全列表
            counts = args.ad_count if args.ad_count else []
            if len(counts) < n: counts.extend([1] * (n - len(counts)))
            
            targets = args.ad_target if args.ad_target else []
            if len(targets) < n: targets.extend([None] * (n - len(targets)))
            
            radii = args.ad_radius if args.ad_radius else []
            if len(radii) < n: radii.extend([5.0] * (n - len(radii)))
            
            for i in range(n):
                specs.append({
                    'file': files[i], 
                    'count': counts[i], 
                    'target': targets[i], 
                    'radius': radii[i]
                })
            config.adsorption_specs = specs

        # 执行任务
        if config.adsorption_specs:
            print(f"\n>>> [5/5] Adsorption Module (Running {len(config.adsorption_specs)} tasks)...")
            total_inserted = 0
            
            for i, spec in enumerate(config.adsorption_specs):
                f_path = spec.get('file')
                count = spec.get('count', 1)
                
                if not f_path: continue

                # [关键修复]：确保 f_path 是字符串，不是列表！
                if isinstance(f_path, list):
                    f_path = str(f_path[0]) # 兜底逻辑
                else:
                    f_path = str(f_path)

                print(f"  -> Task {i+1}: {Path(f_path).name} (n={count})")
                
                # 注入运行时参数
                config.adsorbate_file = f_path  # 此时必为 str
                config.adsorbate_count = int(count)
                config.target_element = spec.get('target')
                config.target_radius = float(spec.get('radius', 5.0))
                
                try:
                    current_atoms = perform_adsorption(current_atoms, config)
                    total_inserted += count
                except Exception as e:
                    print(f"     ❌ Task Failed for {f_path}: {e}")
                    import traceback
                    traceback.print_exc()

            print(f"✅ Adsorption batch complete.")
            status_tag += f"_AdsMix{total_inserted}"
        else:
            print(f"\n>>> [5/5] Skipping Adsorption.")

        # [6] Save
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
