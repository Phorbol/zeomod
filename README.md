# ZeoMod: Intelligent Mesoporous Zeolite Modeling Framework





![image-20251125185928111](ZeoMod.png)

**ZeoMod** is a high-performance Python toolkit designed for computational chemists to automate the construction of complex mesoporous zeolite models. It aims to bridge the gap between ideal crystal structures (CIF) and realistic, functionalized catalyst models.↳

By combining geometric analysis, graph-theoretic algorithms, and automated workflows, ZeoMod efficiently generates doped structures with specific pore shapes, rational surface passivation, and strict adherence to physicochemical rules (such as Löwenstein's rule).↳

------



## 🌟 Key Features



- **🧠 Intelligent Supercell Construction (Proactive Supercell V3)**:↳
  - Automatically calculates the minimal required supercell repetition factors based on the input unit cell geometry and target pore shape, minimizing computational cost while ensuring the pore fits perfectly.↳
- **📐 Shape-Aware Sculpting**:
  - Supports cutting various parameterized geometric shapes: **Sphere**, **Cylinder**, **Box/Cube**, and **Ellipsoid**.
  - Supports generating periodic mesoporous structures or non-periodic **Cluster** models.
- **🛡️ Smart Cleaning & Passivation**:
  - Iteratively cleans unstable, low-coordination atoms resulting from cutting.
  - **Collision-Free H Placement Algorithm**: Uses a rotation-matrix-based search to find optimal positions for passivating protons on dangling oxygens, ensuring reasonable bond lengths, angles, and no steric clashes.
- **⚛️ Topology-Guided Doping**:
  - **Strict Enforcement of Löwenstein's Rule**: Utilizes `matscipy` and sparse matrix calculations to determine graph topological distances, ensuring absolutely no Al-O-Al linkages in large, disordered supercells.
  - **Flexible Doping Modes**: Supports doping by target **Si/Al Ratio** or explicitly specifying a **Fixed Count** of Al atoms.
  - Automatically adds charge-balancing protons (Brønsted acid sites) to introduced framework aluminum atoms.
- **⚙️ Flexible Workflow Control**:
  - Complete pipeline design allows skipping specific stages via command-line flags (e.g., cut pores only without doping, or dope an existing mesoporous structure).
  - Supports a mixed configuration mode of "YAML config file + CLI argument overrides" for batch jobs and rapid debugging.

------



## 📦 Installation





### Prerequisites



- Python 3.8 or higher
- Installation within a separate virtual environment (conda or venv) is recommended.



### Install from Source (Recommended)



Bash

```
# 1. Clone the repository
git clone https://github.com/Phorbol/zeomod.git
cd zeomod

# 2. Install in editable mode (allows code updates without re-installation)
pip install -e .
```

Once installed, the `zeomod` command will be registered in your system path and can be called from any directory.

------



## 🚀 Quick Start





### 1. Prepare Input File



Prepare a CIF file of a pure silica zeolite (e.g., `MFI.cif`) in your working directory.

> **⚠️ Important Scientific Note**: To ensure charge balance and chemical rationality, it is strongly recommended to **always start building from a Pure Silica framework**. ZeoMod's designed workflow is to first cut pores and passivate on pure silica, and then perform controlled doping.



### 2. Generate Configuration Template



If you are unsure how to write the configuration, let ZeoMod generate a template for you:

Bash

```
zeomod --generate-template
# This generates a 'template_config.yaml' file in the current directory.
# You can rename it to 'input.yaml' and modify it.
mv template_config.yaml input.yaml
```



### 3. Run Standard Workflow



After editing `input.yaml`, run:

Bash

```
zeomod input.yaml
```

The program will automatically execute the supercell creation, cutting, passivation, doping, and validation processes. Results are saved in the `output/` directory.



### 4. CLI Rapid Override



Try different parameter combinations directly via command-line arguments without modifying the YAML file:

Bash

```
# Quickly test cylindrical pores with different diameters
zeomod input.yaml --shape cylinder --axis z --diameter 12.0 --output run_cyl_12

# Force doping of exactly 4 Al atoms, ignoring the ratio in the config file
zeomod input.yaml --num-al 4 --output run_fixed_al4

# Replace specific Si atoms (zero-based atom indices) with Al
zeomod input.yaml --doping-sites 12 28 34
```

------



## 📖 Configuration Details



The configuration file (`input.yaml`) supports rich parameter settings.



### Basic I/O & Geometry



| **Parameter**   | **Type** | **Default** | **Description**                                              |
| --------------- | -------- | ----------- | ------------------------------------------------------------ |
| `input_cif`     | str      | "MFI.cif"   | Path to the input crystal structure file.                    |
| `output_dir`    | str      | "output"    | Directory for output results.                                |
| `pore_shape`    | str      | "cylinder"  | Pore shape: `sphere`, `cylinder`, `box`, `ellipsoid`.        |
| `pore_diameter` | float    | 10.0        | Main pore diameter dimension (Å).                            |
| `cylinder_axis` | str      | "z"         | Axis for cylinder alignment: `x`, `y`, or `z`.               |
| `pore_buffer`   | float    | 3.0         | Buffer distance from the pore edge to the unit cell boundary during smart supercell scaling (Å). |
| `cluster`       | bool     | False       | Whether to generate a non-periodic cluster model (True) instead of a periodic mesoporous structure (False). |



### Chemical Doping



| **Parameter**       | **Type** | **Default** | **Description**                                              |
| ------------------- | -------- | ----------- | ------------------------------------------------------------ |
| `si_al_ratio`       | float    | 30.0        | Target Si/Al ratio.                                          |
| `num_al_atoms`      | int      | None        | **[Priority]** Explicitly specify the total number of Al atoms to dope. If set and >0, ignores `si_al_ratio`. |
| `doping_sites`      | list[int] | None       | Exact zero-based atom indices of Si atoms to replace with framework Al. When set, this list determines the number and locations of dopants and takes precedence over `num_al_atoms`, `si_al_ratio`, and `doping_campaign`. Indices refer to the structure entering the doping stage, after any supercell creation and pore cutting. |
| `doping_campaign`   | list     | [...]       | (Advanced) List defining site selection strategies. See template comments for details. |
| `max_doping_trials` | int      | 100         | Maximum attempts to find sites satisfying Löwenstein's rule. |

Explicit sites must identify distinct Si atoms. ZeoMod rejects invalid indices and any selected set that would violate Löwenstein's rule. This mode places framework Al only; `efal_ratio` is not applied.



### Advanced Algorithm Parameters (Usually no need to change)



Can be overridden via CLI using `--set key=value`.

| **Parameter**                | **Description**                                              |
| ---------------------------- | ------------------------------------------------------------ |
| `collision_threshold`        | Threshold for detecting collisions between framework atoms (Å). |
| `collision_threshold_h`      | Threshold for detecting collisions during H atom placement (Å). |
| `rotation_steps_h_placement` | Rotation search precision (number of steps) for smart H placement. |

------



## 🛠️ Advanced CLI Usage



The zeomod command supports three levels of configuration precedence:

CLI --set > CLI Explicit Flags (--diameter, etc.) > YAML Config File > Code Defaults



### Workflow Control Flags



Used to skip specific stages for pipelined operations:

- `--skip-supercell`: Skip automatic supercell creation (assumes input is already a supercell).
- `--skip-cut`: Skip pore cutting and passivation (assumes input already has pores).
- `--skip-doping`: Skip Al doping (generates pure silica mesopore only).

**Example: Generate Pure Silica Mesoporous Structure Only**

Bash

```
zeomod input.yaml --skip-doping --output pure_silica_meso
```



### Universal Override Argument `--set`



Used to dynamically modify any configuration item, including lists or booleans.

Bash

```
# Force specific supercell dimensions and increase Löwenstein retry attempts
zeomod input.yaml --set supercell_dims=[3,3,3] max_doping_trials=500
```

------



## 📄 License



This project is licensed under the MIT License - see the [LICENSE](https://www.google.com/search?q=LICENSE) file for details.

------



## 👥 Authors & Acknowledgments



- Developed by [Geng Jianrui/XMCAO Group]
- Core topology analysis powered by `matscipy` and `scipy.sparse`.↳
- Structure manipulation built upon the `ASE` (Atomic Simulation Environment) library.
