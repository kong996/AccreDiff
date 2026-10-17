# Fig. 2 前期模块开发方案

## 目标与边界

目标是为 Zheng、Lin 与 Kouwenhoven (2016, arXiv:1610.09670v1) 的 Fig. 2 建立可复用的气盘与 10 km 小天体气体阻力计算基础。第一阶段只实现论文专用盘结构、气体流速处方、阻力加速度及其轨道平均诊断；它不能单独重现图中的共振激发和大幅内移。

**代码位置：**第一阶段所有新增 Python 模块、notebook 和测试均放在 `examples/Expansion_module/GasDrag_Type1/`，不修改或新增 `src/accrediff/` 文件。此阶段只读取现有 `accrediff.gas_migration` 接口。代码和数值结果经审核通过后，再单独决定哪些稳定接口迁入 AccreDiff 核心目录；迁入不属于第一阶段工作。

本文是开发设计，不是对论文数值结果的验证声明。论文的正文 §3.1–3.2、式 (1)–(6)、§4.1 和 Fig. 2 是参数及物理范围的依据。

## Fig. 2 的目标数据

- 质量由**物理半径**决定：四颗粒子的半径均为 10 km，起始半长轴分别为 3.5、3.2、2.7、2.3 AU。
- 输出为随时间变化的半长轴 `a`、近日点 `q=a(1-e)`、远日点 `Q=a(1+e)`；论文横轴为 `t/T_dep`，约从 0.1 到 10，默认 `T_dep=1 Myr`。
- 黑、灰曲线是 `ν5`、`ν6` 世俗共振位置，黄色水平线是木星的平均运动共振位置。这些曲线及轨道上的共振尖峰需要后续的行星引力、气盘引力和共振诊断模块。
- 论文的 10 km 天体主要受气体阻力。Type-I 机制不纳入本阶段的 10 km 基准计算。

## 第一阶段文件布局

```text
GasDrag_Type1/
  FIG2_DEVELOPMENT_PLAN.md
  paper_fig2_disk.py          # 论文气盘和现有盘适配器
  paper_fig2_drag.py          # 粒子、气体流和瞬时阻力
  paper_fig2_orbit_average.py # 对同一瞬时力做轨道平均，仅供诊断
  run_paper_fig2_drag.ipynb   # 参数、基准值、两种盘结构比较
  tests/
    test_paper_fig2_disk.py
    test_paper_fig2_drag.py
```

上述文件已按此布局创建。notebook 只负责组装参数、调用模块和画图；物理公式放在 `.py` 文件中，避免 notebook 与今后的轨道积分器出现两套实现。使用目录内导入，不在 `src/accrediff/__init__.py` 增加导出。

## 统一的气盘接口

阻力模块只依赖下面的物理单位接口，不直接访问具体盘对象的 `cfg` 或其时间单位：

```python
class DiskProfile(Protocol):
    def sigma_cgs(self, r_au: ArrayLike, t_myr: float) -> ArrayLike: ...
    def scale_height_au(self, r_au: ArrayLike, t_myr: float) -> ArrayLike: ...
```

约定 `r_au > 0`、`t_myr >= 0`，输入与输出保持标量或 NumPy 数组形状；表面密度为 `g/cm²`，尺度高为 `AU`。气盘对象不负责阻力系数或积分。对于论文平面模型，体密度归一化和气体速度由独立的 `GasFlowPolicy` 读取 `DiskProfile` 后计算。这样更换盘结构不会暗中改变阻力物理处方。

### `PaperFig2Disk`

默认参数不可变、可显式覆盖：

| 参数 | 默认值 | 论文依据 |
| --- | ---: | --- |
| `sigma0_cgs` | 1700 g/cm² | 式 (1)，1 AU |
| `surface_density_exponent` | 1.5 | 式 (1) 的 `k` |
| `h1_au` | 0.025 AU | 式 (2)，1 AU |
| `height_exponent` | 1.25 | 式 (2) |
| `tau_decay_myr` | 1 Myr | 默认模型 |

实现公式：

```text
Sigma(r,t) = sigma0 * (r/AU)^(-k) * exp(-t/tau_decay)
H(r)       = h1 * (r/AU)^(height_exponent)
```

本阶段的 `PaperFig2Disk` 只描述无间隙的局部表面密度与尺度高。论文式 (3)–(4) 的盘引力和木星附近间隙属于后续独立的 `DiskGravity`，不要把其力学公式藏进阻力盘接口。

### `GasMigrationDiskAdapter`

持有现有 `Gas_DiskModel` 和 `Gas_UnitSystem`。`sigma_cgs(r_au,t_myr)` 内部使用 `units.myr_to_code(t_myr)`；`scale_height_au` 调用现有 `disk.scale_height(r_au)`。适配器仅完成接口和单位转换，不重算盘公式。

现有 `Gas_DiskModel` 的径向式为 `Sigma ∝ exp[-(a_alpha1 ln r + a_alpha2) ln r]`、`H = def_h_1*r^(1+beta)`。因此 **`a_alpha1=0`、`a_alpha2=1.5`、`def_h_1=0.025`、`beta=0.25`、`Sigma0_cgs=1700`、`tau_decay_myr=1` 时，它可以表示与论文相同的无间隙 `Sigma` 和 `H`**。这要成为两个盘提供者的等价性测试。气体体密度及速度处方仍需分别配置，不能因为 `Sigma`、`H` 相同就视为完整物理模型相同。

## 粒子与气体流

`Planetesimal` 以 `radius_km` 和 `solid_density_cgs` 为主输入，质量由球体体积派生。论文基准：`radius_km=10`、`solid_density_cgs=1`，对应质量约 `4.19e18 g`，即约 `7.0e-10 M_earth`。这与当前质量扫描 notebook 的最小质量不是同一对象。参数校验要求半径、密度为正。

`GasFlowPolicy` 必须显式选择以下约定：

| 处方 | 体密度 | 方位气体速度 |
| --- | --- | --- |
| `paper_2016` | 平面中 `rho_g=Sigma/H` | `v_g=v_K(1-eta)`；论文式 (5) 给出 `eta≈(c_s/(4v_K))²`。若取 `c_s=H*Omega_K`，须在结果中注明这一闭合关系 |
| `genga_style` | `rho_g=Sigma/(sqrt(2π)H)`，离盘面时再乘垂向因子 | 按 `gas.cu` 的局部压力梯度和 `eta` 公式计算 |

`PaperFig2Disk` 与 `GasFlowPolicy` 是两个独立配置。不能将论文的 `rho_g=Sigma/H` 与现有 notebook 的高斯中面归一化混用。力模型统一用当地 `v_K=r*Omega_K`，不得把 `Omega_K*1 AU` 当作任意半径处的开普勒线速度。

## 阻力加速度与轨道平均

基础接口使用位置和速度，返回同一坐标系下的加速度：

```python
drag = AerodynamicDrag(planetesimal, flow_policy, law="paper_2016")
accel = drag.acceleration(position_au, velocity_au_per_code_time,
                          t_myr, disk_profile)
```

统一定义相对速度 `u = v_particle - v_gas`。论文式 (5)–(6) 写作 `a_drag=-u/T_damp`、`T_damp^(-1)=C_D*rho_g*|u|/(rho_p*R)`；按其字面形式，`a_drag=-[C_D*rho_g/(rho_p*R)]*|u|*u`，默认 `C_D=0.165`。`gas.cu` 使用 `pi*C_D*R²*rho_g/(2m)` 的系数；对球体等于 `3*C_D*rho_g/(8*rho_p*R)`。两种归一化不相等，必须分别实现为 `paper_2016`、`genga_style`，在代码和图例中注明选用者；审核前不要通过随意调整 `C_D` 吸收差异。论文关于与 Whipple (1972) 的 `3/8` 系数差别也要在公式注释中标明。

所有内部计算应明确选择一种单位体系。建议力模型内部使用现有 GENGA-like 代码单位，入口将半径、密度和时间统一换算，出口明确标注加速度单位。每个公开函数在 docstring 中写清单位；避免在积分回调中重复构造转换常数。

轨道平均模块应在指定 `a,e` 上采样偏近点角 `E`，按 `dM=(1-e*cos E)dE` 加权，利用同一瞬时阻力力由轨道能量和偏心率矢量计算 `<da/dt>`、`<de/dt>`。它用于与已有 `gas_drag_dispersal_tracks.ipynb` 比较及无共振基准测试。后续包含木星、土星共振的积分器应直接调用**瞬时力**，不能用轨道平均结果代替共振动力学。

## 第一阶段交付与验收

1. 配置和两个盘提供者可在不改阻力计算代码的情况下互换；论文参数与现有盘的等价配置在多个半径、时间点给出相同 `Sigma` 和 `H`。
2. 数值基准覆盖 `t=0`、`t=T_dep`，以及 1、2.3、3.5 AU；表面密度满足 `Sigma(r,T_dep)/Sigma(r,0)=e^-1`。在 1 AU、`t=0`，论文模型应给出 `Sigma=1700 g/cm²`、`H=0.025 AU`。
3. 10 km、`1 g/cm³` 的派生质量、两种体密度约定、开普勒速度及两种阻力系数均有独立的数值核对；零相对速度时加速度为零，阻力功率 `u·a_drag` 非正。
4. 轨道平均结果随采样数加密而收敛；瞬时阻力与平均导数使用同一个力实现。
5. notebook 展示论文盘和现有盘适配器的剖面比较、基准阻力时间尺度和无共振轨道演化。结果应明确标为“阻力基准”，不称为 Fig. 2 复现。
6. 第一阶段不运行或修改 GENGA、不计算 `ν5/ν6` 或 MMR，也不更改 `src/accrediff/`。审核时检查物理公式、单位、数值基准、接口稳定性及 notebook 输出，再决定核心模块的迁移方案。

## 审核后迁入核心的条件

只有在论文与现有盘的等价性、阻力系数约定和数值收敛均通过审核后，才评估将 `DiskProfile`、适配器及瞬时阻力实现迁入 `src/accrediff/`。论文 Fig. 2 的特定默认值与绘图脚本仍可留在示例目录。核心迁移需要单独检查公开 API、向后兼容和文档，不在本阶段预先改动核心模块。

## 实施记录（待审核）

第一阶段的三个 Python 模块、示例 notebook 和目录内测试已建立。`PaperFig2Disk` 与等价参数的 `GasMigrationDiskAdapter` 在多个半径、时间点的 `Sigma` 和 `H` 数值一致。10 km、`1 g/cm³` 天体的派生质量为 `4.18879e18 g`，约 `7.01364e-10 M_earth`。默认 `e0=0.01` 的纯阻力示例在 10 Myr 内只产生约 `0.0015–0.0019 AU` 的内移；此结果是无共振基准，不用于声称复现 Fig. 2。测试命令：

```bash
python -m pytest -q examples/Expansion_module/GasDrag_Type1/tests
```

额外核对的原始研究：[Adachi、Hayashi 与 Nakazawa (1976)](https://academic.oup.com/ptp/article/56/6/1756/1828048) 研究了气体阻力下椭圆轨道的半长轴、偏心率和倾角变化；[Beaugé 等 (2010)](https://academic.oup.com/mnras/article/408/1/503/1059401) 明确讨论了 km 级天体阻力与相对速度平方的关系。它们支持使用同一个瞬时二次阻力力计算轨道平均变化，但本项目的论文特定系数仍以 Zheng 等 (2016) 的式 (5)–(6) 为准。
