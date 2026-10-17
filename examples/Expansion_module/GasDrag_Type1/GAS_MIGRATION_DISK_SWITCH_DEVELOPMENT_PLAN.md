# Fig. 2 动力学实验切换至 `gas_migration.py` 气盘：开发方案

## 目标与边界

在现有太阳、木星、土星、四颗 10 km 测试粒子的直接积分框架中，将气盘结构提供者从 `PaperFig2Disk` 切换到 `accrediff.gas_migration.Gas_DiskModel`。气盘结构指 `Σ(r,t)`、`H(r)` 及其消散时间；10 km 天体仍使用论文的气体阻力处方，木星和土星仍以 N-body 引力激发共振。Type-I 力矩不参与本实验。

所有新实现、测试及 notebook 改动留在 `examples/Expansion_module/GasDrag_Type1/`。本阶段只读取 `src/accrediff/gas_migration.py`，不改动核心模块。输出称为“使用 AccreDiff 气盘结构的 Fig. 2 式实验”；只有等价参数对照组仍对应论文气盘的无间隙局部剖面。当前盘引力是软化环积分近似，不是论文式 (3)–(4) 的 Ward 级数。

## 现状与关键差异

| 部件 | 现有能力 | 本次工作 |
| --- | --- | --- |
| `GasMigrationDiskAdapter` | 已把 `Gas_DiskModel` 转为 `sigma_cgs(r_au,t_myr)`、`scale_height_au(r_au,t_myr)`、`sigma_slope_p(r_au)` | 沿用并补齐验证；避免在积分代码中直接读取 `Gas_ModelConfig` |
| `AerodynamicDrag` | Python 路径只依赖气盘接口 | 固定 `paper_2016` 流速、体密度和阻力系数；两个盘共用同一物理处方 |
| `PaperGapDiskGravity` | 构造时读取任意盘的初始 `Σ,H`，计算带间隙环势；通过统一消散因子缩放 | 为每个新盘对象重新建表，验证可分离消散、边界和表格精度 |
| `paper_fig2_resonances.py` | 使用盘引力给出的力及导数计算候选 ν5/ν6 | 对每种盘重算共振位置；不能沿用论文盘的旧曲线 |
| `paper_fig2_compiled.py` / `.c` | 当前拒绝非 `PaperFig2Disk`，C 中写死幂律 `Σ,H` | 扩展参数化 C 力；仍与 Python 版计算同一瞬时阻力 |
| `run_paper_fig2_resonances.ipynb` | 默认短跑，显示一份已保存的论文盘 10 Myr 轨迹 | 用显式实验配置选择盘；图和数据携带盘标识，防止旧结果被误认为新结果 |

`Gas_ModelConfig` 的**模块默认值**（例如 `Sigma0_cgs=1600`、`tau_decay_myr=2`）与现有 `gas_profile.ipynb` 的实际配置不同。后者明确设置 `Sigma0_cgs=1700`、`tau_decay_myr=1`、`a_alpha1=1`、`a_alpha2=0.5`、`alpha_visc=1e-4`；`beta=0.25`、`def_h_1=0.0335743` 沿用模块默认。第二阶段使用这组 notebook 参数作为目标盘基准，并在输出元数据中逐项记录，不以隐式默认值代替。

## 两阶段实验配置

### A. 等价替换：只检验盘提供者

建立 `Gas_DiskModel` 与 `GasMigrationDiskAdapter`，设：

```python
Gas_ModelConfig(
    Sigma0_cgs=1700.0,
    a_alpha1=0.0,
    a_alpha2=1.5,
    def_h_1=0.025,
    beta=0.25,
    tau_decay_myr=1.0,
)
```

此时 `Σ=1700 r^-1.5 exp(-t/1 Myr)`、`H=0.025 r^1.25 AU`，与 `PaperFig2Disk` 的无间隙局部剖面代数等价。使用同一 `GasFlowPolicy("paper_2016")`、`AerodynamicDrag`、`GapGravityConfig`、行星与粒子初始条件、积分器步长及输出时刻。除气盘实现类外，不改变其他变量。这一阶段的目标是证明适配与编译路径正确，而不是发现新的物理差别。

### B. 目标盘：`gas_profile.ipynb` 所用结构

将盘结构换成该 notebook 的实际剖面：

```text
Sigma0_cgs = 1700 g/cm²
a_alpha1 = 1.0, a_alpha2 = 0.5
def_h_1 = 0.0335743, beta = 0.25
tau_decay_myr = 1.0 Myr
```

其 `Σ(r,t)=Σ0 exp[-(a_alpha1 ln r+a_alpha2)ln r]exp(-t/τ)`，局部对数斜率 `p(r)=2 a_alpha1 ln r+a_alpha2`，不再是单一幂律。`alpha_visc`、`prandtl`、`gamma_eff` 和 `use_saturation` 属于现有 Type-I 力矩模型；本实验不调用 `Gas_TorqueModel`，这些字段可记录以重建配置，但不得解释为对 10 km 阻力施加 Type-I 效应。`particle_mass_mearth`、`a0`、`t_max_myr`、`n_steps` 也不控制这里的 N-body 粒子或积分器；分别由 `Planetesimal`、`Fig2SystemConfig`、`IntegrationConfig` 决定。

## 物理边界与间隙策略

第一轮固定 `GapGravityConfig` 的间隙 4.5–11 AU、内边界 0.05 AU、外边界 50 AU、软化 `0.6H`；盘引力的环质量及软化高度均从当前选中的盘重建。局部气体阻力仍读取无间隙的 `Gas_DiskModel` 剖面，延续当前 Fig. 2 实验的“内盘测试粒子阻力、带间隙盘引力”分工。若轨道进入间隙或越出力表域，须明确报错或记录事件，不得默默套用无间隙阻力。

这项固定间隙选择是**控制变量假设**，并不宣称新盘能自洽地形成同一个间隙。完成目标盘基准后，另做间隙边界、外盘截断和软化系数敏感性试验，分别报告对行星进动率、ν5/ν6 位置和粒子轨迹的影响。盘密度在径向有弯曲时，不能只用某一半径的比值代替全部盘面；`PaperGapDiskGravity` 当前通过多个半径检查时间衰减可分离性，实施时应把验证扩展到完整环半径网格和多个时间点。

气体流速与阻力归一化保持论文约定：`ρ_g=Σ/H`、`η=(H/4r)^2`、`C_D=0.165`，并继续使用论文的二次阻力系数。切换 `Σ,H` 时这三者的数值自然变化，但不混入 `genga_style` 的体密度或阻力系数。若以后要比较 GENGA 处方，建立独立实验组。

## 编译回调的设计

优先把现有 C 回调扩展为 `Gas_DiskModel` 的**解析径向参数化**，而非每个固定步长从 C 调 Python：

```text
ln Σ(r,t) = ln Σ0 - a_alpha1(ln r)^2 - a_alpha2 ln r - t/τ
H(r)      = def_h_1 r^(1+beta)
η(r)      = [H(r)/(4r)]²
```

保留既有径向盘引力表及 Hermite 插值，C 只需新增盘结构参数和模型分支；`PaperFig2Disk` 继续走原参数分支，`GasMigrationDiskAdapter` 走上述解析分支。所有 C 参数由 Python bridge 从所持 `Gas_DiskModel.cfg` 生成，并校验物理单位、`G=1`、恒星质量和气体流速闭合。`σ`、`H` 与阻力的 Python 和 C 实现必须在多个半径、时间与非圆轨道状态逐点比较，不能仅以最终 `a` 接近作为验收。若未来盘结构无法写成此解析式，再考虑带误差控制的气体性质查表接口。

C 回调目前使用进程全局参数，不能在同一进程并行运行不同配置的模拟；本阶段维持串行运行并在文档中说明。编译库缓存应依据 C 源码和 REBOUND 头文件版本更新。缺少编译器时允许 Python 短段诊断，但不得将 Python 回调作为默认 10 Myr 路径。

## notebook 与输出组织

在 `run_paper_fig2_resonances.ipynb` 增加清晰的盘选择，例如 `paper_reference`、`gas_migration_matched`、`gas_profile_baseline`。配置单元只组装盘、阻力、盘引力和系统；已有计算模块不复制到 notebook。短段运行先验收；完整运行须显式开启。显示已保存数据时，用当前配置哈希或盘标签检查元数据，不匹配则标记为旧结果，不显示成当前配置的预测。

输出文件名包含实验标签，且元数据保存：全部 `Gas_ModelConfig` 字段、适配器类型、阻力与流速处方、间隙参数、行星及粒子初始角、REBOUND 版本、C/Python 后端、步长、输出采样、配置与源码哈希。保留 `t,a,e,ϖ,q,Q` 及笛卡尔状态。图题和图注注明“paper reference / matched adapter / gas profile baseline”及软化环势近似。默认不覆盖现有 `fig2_exploratory_10myr.npz`；用户重新运行时给出明确的新路径。

## 开发顺序与验收门槛

1. **参数和单位审计：**确认 `gas_profile.ipynb` 的有效配置与 `Gas_ModelConfig` 默认值；验证 Myr↔代码时间、g/cm²↔代码表面密度、`H` 的 AU 单位。检查 `Gas_ModelConfig.Mstar/G` 与 N-body 系统一致。
2. **等价盘静态测试：**在 0.5、1、2.3、3.5、5.2、9.58、12 AU 及 0、1、5 Myr 比较 `Σ,H,p`、`ρ_g,η,v_g`、瞬时阻力和盘引力/导数；以相对误差为主，接近零的量设绝对误差。严格等价处误差目标约 `1e-10`，但盘力数值积分另按表格收敛评估。
3. **C/Python 力一致性：**在两种盘及多组位置、速度、时刻比对盘力与阻力；100 yr 轨迹的 `a,e` 差达到现有编译路径测试量级后，才开展长时计算。额外覆盖径向域外和无效参数。
4. **短段动力学对照：**同一初始条件下比较论文盘与等价适配盘的 500 yr 轨迹、行星进动和 ν5/ν6 根。两者在既定数值容差内重合；目标盘则允许偏离，但必须能由 `Σ,H` 和盘进动差解释。
5. **长期数值检查：**至少比较两个步长的 0.1 Myr 与具有代表性的共振窗口，验证 `IAS15` 短窗口、盘力表分辨率、输出采样对结果的影响。当前 10 Myr 论文盘运行四轨道趋于约 1.69 AU，其成因尚未验证，不能以此作为新盘正确性的目标值。
6. **10 Myr 目标盘运行：**保存可重绘数据与图，报告各粒子的 `a,e,q,Q`、候选 ν5/ν6 和 MMR 位置、越界/近距离事件；并列展示论文盘、等价适配盘和目标盘。轨道形态按“物理配置变化、共振响应、数值误差”三层解释，不宣称复现论文原图的单次轨迹。

完成第 1–4 项才进入长期生产运行。时间步收敛、间隙敏感性和初始相位敏感性应作为长期结果的可信度门槛，而不是仅检查程序运行结束。

## 实施记录与当前限制

已增加 `paper_fig2_experiment.py` 的三种命名盘配置，扩展 `paper_fig2_compiled.py` 和 `paper_fig2_force.c` 以解析计算 `Gas_DiskModel` 的弯曲 `Σ(r,t)`，并更新 `run_paper_fig2_resonances.ipynb` 的盘选择、结果路径及配置匹配检查。目标盘默认采用 `gas_profile.ipynb` 的参数。论文盘、等价适配盘和目标盘的短段积分、Python/C 力对照及目录内测试已通过；最终测试数量和结果以本次开发的验证记录为准。

目标盘的 10 km、10 Myr 探索性轨迹已保存于 `outputs/gas_profile_baseline_10myr.npz`，图为同名 `.png`。使用 `WHFast dt=0.025 yr`，约 4 亿固定步，耗时约 364 s；四颗粒子的末态半长轴约为 `2.231214、2.229208、2.233043、2.216353 AU`。**这不是收敛解：**在 0.1 Myr 处，步长 0.025 与 0.0125 yr 的首颗粒子半长轴相差约 0.0149 AU；0.0125 与 0.00625 yr 的差约 0.0137 AU。该粒子靠近强共振，对步长和相位敏感。长期轨迹仅展示模型可运行及定性共振响应；需局部高精度积分、事件窗口采样及相位/软化/间隙敏感性分析后，才能作定量比较。

旧文件 `outputs/fig2_exploratory_10myr.npz` 经元数据核对实际使用 **100 km** 天体。notebook 的配置匹配检查会拒绝在当前 10 km 设置下显示该旧图；此前开发文档中将其称为 10 km 的记录已更正。它不可作为本次论文盘等价实验的长期基准。论文盘与等价适配盘的 **10 km 长期对照仍待完成**。

### 自定义盘越界诊断（2026-09-19）

对 notebook 中保存的 `Sigma0_cgs=2000`、`a_alpha1=1`、`a_alpha2=0.5`、100 km、`rho_p=3 g/cm³` 配置，编译回调现记录首次越界的天体、积分时刻、半径和太阳二体能量。默认 0.3–13 AU 表下，三个步长均由 `particle_0` 向外越界：`dt=0.025 yr` 为 `0.625618 Myr`，`0.0125 yr` 为 `0.485371 Myr`，`0.00625 yr` 为 `0.517781 Myr`。事件身份和方向稳定，但发生时间不收敛，符合强共振/巨行星散射的混沌敏感性。

最小步长在首次穿过 13 AU 时的太阳二体能量为 `-0.0198741`，仍受束缚，所以 13 AU 只是原插值表边界。将力表扩至盘外缘 50 AU、同时把表点增至 3072 后，`particle_0` 在 `0.493545 Myr` 穿过 50 AU，二体能量为 `+0.627636`，已经非束缚。由此，当前实验应把该轨迹分类为强散射/逃逸并终止或移出内盘样本，而不应仅扩大力表来生成一条连续 10 Myr 内盘曲线。具体散射时刻对步长敏感，不能作为收敛预测；若研究目标转为逃逸轨迹，则需改用适合近距离遭遇的积分方案和明确的逃逸判据。
