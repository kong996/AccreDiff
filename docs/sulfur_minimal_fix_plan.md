# 含硫模块第一步审查：最小修改方案

> 归档说明（2026-10-17）：本方案中的两处公式修复及回归测试已独立提交为 `a0a2896`；文中 52 项测试记录针对当时完整工作区，公式修复独立提交快照为 42 项通过。后续求解器增强提交 `c200cf0` 的独立快照为 52 项通过。

状态：已按方案完成最小修改并验证（2026-10-17）。

执行结果：

- 修正 Boujibar FeO wt% 对数项和硫容量 Al2O3 归一化分母；公开接口未变。
- 修正旧参考式，补充固定参考值及含铝/无铝两种完整相组成回归。
- 使用现有 `/opt/anaconda3/envs/data-analysis/bin/python` 环境，无需安装依赖。
- Boujibar 定向测试：8 passed；项目默认测试集：52 passed。
- S8 和 S7 示例全部代码单元在内存执行通过，原有断言通过；未保存 notebook 或覆盖已有结果。绘图采用 Agg 后端。
- `git diff --check` 通过。其余原有工作区改动保留。

## 目标与边界

修正 Boujibar 路径的 FeO 单位和硫容量氧化物归一化基准，保持公开函数签名、参数含义、返回结构、模型名称、事件选择和调用顺序不变。现有 notebook 无需修改调用语句。

接口兼容不等于结果保持不变：Boujibar 结果将发生科学上预期的变化，耦合求解的根、可达性和耗时也可能变化。不得通过放宽容差、自动接受近似解或修改失败处理来掩盖变化。

## 修改 1：FeO 的 wt% 对数项

文件：`src/accrediff/differentiation_sulfur.py`。

函数：`_log_D_S_model_boujibar_molar`。

- 保留 `w_feo` 为质量分数，保留现有物理域校验。
- 将 `math.log10(w_feo)` 改为 `math.log10(w_feo) + 2.0`。
- 注释说明：文献直接项用 FeO wt%，因此 log10(100*w_feo) = log10(w_feo)+2；其他金属组成项仍使用质量分数。
- 同步修正文档字符串和“包括 FeO 均为质量分数”的误导性注释。
- 不在最终返回值处整体乘 100，不改质量—摩尔分配系数转换。

采用对数加 2 的写法可保留现有对数域计算，也避免额外乘法。

## 修改 2：仅调整硫容量的氧化物分母

函数：`ForwardKDOSolverSulfur._phase_properties`。

将：

```python
oxide_sum = sol["x_"] + sol["y_"] + sol["z_"] + p.u + p.m + p.n
```

改为：

```python
# p.m counts AlO1.5 units; capacity uses Al2O3 formula units.
oxide_sum = sol["x_"] + sol["y_"] + sol["z_"] + p.u + 0.5 * p.m + p.n
```

依据：`sulfur_params_from_composition` 将 Al2O3 摩尔量乘 2 后存入 `p.m`；硫容量的常规氧化物摩尔分数需要将其还原。当前源码中 `oxide_sum` 仅被 Boujibar 硫容量路径消费。

保留返回字典中的 `oxide_sum` 键，不新增参数或分支。保持 `sil_sum`、`metal_sum`、相质量和平均摩尔质量的现有定义；不修改 `p.m` 的存储约定和 AlO1.5 摩尔质量。

## 测试修改

文件：`tests/test_differentiation_sulfur.py`。

1. 修正现有 `test_boujibar_function_returns_molar_coefficient` 的参考式：明确形成 FeO wt% 后再取对数，避免参考式继续复制实现错误。现有输入的独立参考结果为 D_mole = 42.23846633617064。
2. 增加组成转换回归：由完整氧化物质量表独立计算各氧化物摩尔量和硫容量；使用非零 Al2O3，核对 `_phase_properties` 与传入 Boujibar 函数的 FeO/CaO/MgO 摩尔分数。测试同时核对相质量、`sil_sum` 和平均摩尔质量仍满足原定义。
3. 增加组合路径回归：固定完整两相库存，独立从质量和氧化物化学式单位计算文献 D_weight，再核对 `_sulfur_coefficients` 的 D_weight、D_mole 和 ln(D_mole)。包含无铝和含铝状态，隔离两项修正的影响。
4. 复用现有硫守恒、分配恒等式、库存缩放不变性、S8/S7 模型选择、物理状态审计及不可达异常测试；出现失败时先判断是旧预期值还是实际回归，禁止直接更新全部预期或放宽断言。

上述库存是公式验证状态，无需声称它们同时满足所有主元素平衡方程。耦合求解验证由现有求解器测试承担。

## 验证顺序与验收

1. 在具有 numpy/scipy/pandas/pytest 的环境中运行定向公式与组成转换测试。
2. 运行 `tests/test_differentiation_sulfur.py` 全部测试，再运行项目默认测试集。
3. 用现有 S8 初始事件和 S7 后续事件示例进行隔离 smoke check，输出至临时目录；不覆盖已有科研结果，不直接运行会覆盖结果的批处理入口。
4. 检查调用仍使用原签名、字段和模型选择。S8 公式路径不应因本次补丁改变；S7 以新文献参考值和原有守恒/残差阈值验收。

执行时记录实际通过数、失败项和环境限制。此前隔离审查所用环境缺少 pytest/scipy；正式验证已找到并使用具有完整依赖的现有 data-analysis 环境。

## 本次不修改

- Rose–Weston 公式、Ma 活度和 KD 参数。
- 硫容量未跟踪的 TiO2/Na2O/K2O 和碳模型。
- 氧 KD、IW 定义、求解策略、容差和近似解规则。
- 事件再平衡比例、物质混合、CMF 口径、notebook 编排和批处理结果管理。
- 对外接口、`__init__.py` 导出、结果字段和文件格式。

## 科研结果衔接

修复后旧 Boujibar 结果不代表新公式结果，不可与重算结果直接混合解释。完整多阶段影响需后续另行重算评估，不属于本次最小代码补丁。不增加保留错误公式的兼容开关。

## 预计变更文件

计算代码 1 个文件，测试 1 个文件；本方案文档记录修改范围。无需修改其他调用代码。
