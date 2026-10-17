# 含硫分异求解中的无精确根事件

> 历史结果说明（2026-10-17）：本文的事件数、无根统计及残差范围来自 2026-09-16 的旧公式/旧组成口径计算。Boujibar FeO 单位和氧化物归一化修复后，必须重新计算完整演化及独立审计才能更新这些结论；本文数字不能当作修复后的验证结果。求解策略属于历史工作流记录，后续调整需单独确认。文中 `notebooks/` 为被忽略的本地研究目录。

## 文档目的

本文档记录 Fe–Ni–Si–O–S 多阶段分异计算中出现的“目标氧 KD 无精确物理根”问题、当前采用的通用处理方案，以及所有 N-body 模型运行完成后需要开展的系统性研究。

本文档是问题记录和后续研究入口，不把当前单个 N-body 模型的统计结果解释为普遍规律。

记录日期：2026-09-16。

## 问题描述

后期碰撞事件使用 Boujibar S7 硫分配模型，并在给定压力、温度和物质组成下联合求解 Fe、Ni、Si、O 和 S 的金属–硅酸盐分配。

部分事件可以同时满足：

- Ni、Si 和 S 的分配平衡；
- 元素守恒；
- 各相含量非负；
- 目标氧分配系数，即 `KD_O_model = KD_O_target`。

另一些事件可以满足前三项，却无法达到目标氧 KD。表现为物理解分支存在，但目标氧 KD 位于连续可行分支之外或分支间隙中。

这类情况不能简单解释为普通数值误差。可能原因包括：

- 当前 Ni、Si、O 参数化在实验范围以外的压力或温度上被外推；
- 含硫活动度修正改变了可行物理解分支；
- 当前模型尚未加入硫对氧分配的额外修正；
- 给定组成和经验 KD 参数在当前方程组中不完全相容。

## 当前样例的证据

在 `note_02` 工作流中共有 264 个后期碰撞事件：

- 188 个事件由有界联合求解得到完整精确解；
- 76 个事件的快速搜索没有找到精确氧 KD 根；
- 对这 76 个事件执行 thorough 连通根复核后，仍有 76 个事件未找到精确根；
- 未发现“fast 搜索漏掉、thorough 搜索找到”的事件；
- 76 个事件中包含 73 个 `small` 和 3 个 `global` 事件；
- 56 个事件的压力低于 5 GPa；
- thorough 最近物理解的氧 KD 对数残差绝对值为 0.149–1.318 dex。

这里的“未找到精确根”是当前模型、输入状态和既定搜索域内的数值结论，不是对所有热力学模型的数学无解证明。

这些数字只用于验证和记录当前流程，不能在其他 N-body 模型运行完成前解释为无根事件的一般发生率。

## 当前采用的通用求解方案

当前生产流程按以下顺序处理每个事件。

### 1. 有界联合求解

优先同时求解 Fe–Ni–Si–O–S 方程。结果必须通过 Ni、Si、S 平衡残差、氧 KD 残差、非负性、物理状态和独立元素守恒检查。

通过后记录：

```text
solution_strategy = simultaneous
oxygen_converged = True
accepted_approximate = False
```

### 2. 连通物理根搜索

联合求解失败时，搜索相互连接的物理解区间。找到满足完整方程的根后记录：

```text
solution_strategy = connected_root
oxygen_converged = True
accepted_approximate = False
```

### 3. 最近物理解回退

没有精确根时，允许使用满足 Ni、Si、S 平衡、元素守恒和物理条件的最近状态继续计算。

当前工作流设置的最大氧 KD 偏差为：

```text
abs(log10(KD_O_model / KD_O_target)) <= 1.5 dex
```

这种状态必须记录为：

```text
solution_strategy = nearest_physical
oxygen_converged = False
overall_converged = False
accepted_approximate = True
```

1.5 dex 是为了让当前外推工作流能够继续运行的操作上限，不是实验误差，也不是已经获得文献支持的物理容差。此类结果不能与完整平衡解合并称为“全部收敛”。

## 独立外推情景

对 thorough 复核后仍未找到精确根的事件，保留两个独立端元。

### `thorough_nearest_physical`

- 使用 thorough 搜索得到的最近物理解；
- Ni、Si、S 平衡及物理、守恒审计通过；
- 目标氧 KD 不满足；
- 用于模型约束下的外推敏感性分支。

### `phase_preserving_no_repartition`

- 不强制发生新的化学重新分配；
- 再平衡库存保持进入事件前的硅酸盐/金属相分配；
- 不声称满足新的 Ni、Si、S、O 平衡；
- 用作“不发生重新平衡”的端元。

两个情景必须单独传播和报告，不能与精确解混合为一条没有标记的结果。

## 必须保存的诊断字段

每个事件至少保存：

```text
N-body model
planet
event_number
Time
event
target_id
impactor_id
m_target
m_impactor
Impact ratio
P_equil
T_equil
k_mantle
S
S_sil
S_met
KD_O_target
KD_O_model
O_KD_log_residual
D_S_mole
sulfur_residual
solution_strategy
oxygen_converged
overall_converged
accepted_approximate
max_conservation_error
```

后续统计不能只保留成功结果，否则会低估无根事件比例。

## 全部 N-body 模型完成后的研究计划

### 发生比例

- 无精确根事件数及占全部事件的比例；
- 按 N-body 模型和最终行星分别统计；
- 按 `small`、`global`、`partial` 分组；
- 区分联合求解失败、fast 搜索失败和 thorough 复核失败。

### 参数范围

比较精确解和无根事件的：

- 压力和温度；
- 目标氧 KD；
- 撞击质量、撞击比和事件类型；
- `k_mantle`；
- 总硫以及 `S_sil/S_met`；
- Fe、Ni、Si 和 O 的总库存及相比例；
- Ni、S 活度系数；
- `D_S_mole`；
- 最近物理解的氧 KD 残差及残差符号。

### 分支结构

- 检查残差是否聚集在若干稳定分支；
- 判断目标 KD 位于可行域外部还是不连续分支间隙；
- 对代表性事件绘制 `KD_O_model` 随 FeO、硫分配和 Ni 状态的变化；
- 检查 fast 与 thorough 是否选择不同分支。

### 模型敏感性

- 分别传播精确解、thorough 最近物理解和相态保持情景；
- 比较最终行星 FeO、核质量、`S_sil`、`S_met` 和总硫；
- 测试氧 KD 参数化和低压/高压分支；
- 评估硫对氧分配修正是否能减少无根事件；
- 依据实验回归不确定度重新确定可接受的外推范围。

## 当前相关文件

- 主流程：`notebooks/work_flow_S/D_01_differentiation_s.ipynb`
- 当前基准结果：`notebooks/work_flow_S/Tables/06_differentiation_S.xlsx`
- thorough 复核脚本：`notebooks/work_flow_S/audit_06_thorough.py`
- thorough 复核结果：`notebooks/work_flow_S/Tables/06_thorough_root_audit.xlsx`
- 外推情景构建脚本：`notebooks/work_flow_S/build_06_extrapolation_scenarios.py`
- 外推情景结果：`notebooks/work_flow_S/Tables/06_extrapolation_scenarios.xlsx`

## 当前决策

在全部 N-body 数据完成前：

1. 保留当前三层通用求解方案，使批量计算可以继续；
2. 精确解与最近物理解必须保留不同标记；
3. 不把 1.5 dex 解释为物理误差范围；
4. 不根据 `note_02` 单独推断无根事件的一般发生率；
5. 保留完整事件诊断，待全部模型完成后统一研究发生比例、参数范围和模型敏感性。
