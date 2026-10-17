# 含硫 N 体模型批量工作流

日常入口是 [run_batch.ipynb](run_batch.ipynb)。在参数单元填写 `MODEL_LIST`、`DATA_ROOT`、`OUTPUT_ROOT`、运行终点与 A–E 的完整参数；先运行计划预览，再运行正式执行单元。`LAST_STAGE` 可截断于任一阶段，默认运行到 E_01。每次正式运行都从 A_01 重新计算全部启用阶段，不使用缓存。Earth/Mars 比较值来自 `E_00_CMF.ipynb`，不是拟合参数。

`OUTPUT_ROOT` 是所有产物的指定根目录，默认 `kong_2027_backup/work_flow_S_batch/results`。每个模型的最新**成功**结果保存在 `OUTPUT_ROOT/models/<model_id>/`：`tables/`、`stage_workspaces/`、`parameters.json`、`manifest.json`、参数化 notebook 和执行 notebook 一起对应同一次计算。同一模型再次成功运行会整体覆盖旧结果；其他模型不受影响。`scenario_id` 记录在参数快照中，不另建历史方案目录。保存副本的首单元显示最终模型目录路径；notebook metadata 记录执行时使用的暂存路径。

`OUTPUT_ROOT/latest_run/` 只保留最近一次批处理尝试的 `batch_status.csv`、`attempt.json`、`logs/<model_id>/` 和本次成功模型的汇总表。下次运行会覆盖它。失败模型停止下游，原有完整模型结果保持不变；状态表明确说明本次失败及旧结果是否仍在。汇总表只纳入本次成功模型，不混入其他模型的旧结果。`accepted_approximate` 和氧 KD 残差保留原值，不能把近似解视作精确根。

执行器在 `OUTPUT_ROOT/.staging/` 隔离运行、验证并按模型提交，完成后清理暂存文件。A_01 在同一批次算一次并复制给各模型；原始 A–E notebook 与 N 体数据保持只读。此前 `runs/<run_id>/...` 历史目录不会自动清理或迁移。本目录仍属于被 Git 忽略的试验区。输出规则见[最新结果设计](../../docs/sulfur_batch_latest_results_design.md)，参数设计见[批量流程设计](../../docs/sulfur_batch_workflow_design.md)。

## 2026-10-17 路径与组成口径修订

- 批处理固定读取本备份目录相邻的 `../work_flow_S/` 模板；入口显式加载本目录执行器，默认结果写入备份目录，N-body 输入仍可引用主仓库数据。
- 手动执行的 A–E notebook 从当前目录的祖先定位仓库；默认使用备份工作流目录，并尊重显式 `WORKFLOW_DIR`，不再覆盖隔离工作目录。
- A_01 新增 `Meteorites_inventory_per_100g.csv` 和 `Meteorites_mass_coverage.csv`；批处理参数、文件复制及输出校验同步。
- B_01 按 N-body 来源质量比例混合每 100 g 完整端元的摩尔库存；D_01 传递未建模质量，E_01 输出条件 CMF、质量覆盖率及完整物体 CMF 上下界。
- 参数快照以 `composition_basis=mol_per_100g_whole_body` 标记新口径；旧表格和已执行 notebook 保留，不能与新结果混合解释。重新运行应从 A_01 开始。
- D_01 的同时求解、fast 搜索和最近物理解回退逻辑未调整；`ROOT_SEARCH='thorough'` 不代表所有回退都执行 thorough 搜索，求解策略待另行确认。
- 本轮不运行完整批处理、不刷新 notebook 保存输出，也不改动 `results/`、`runs/` 和原 `Tables/`。
- 验证：本目录批处理测试 9 项、核心默认测试 52 项全部通过；新增 A→B 实际 S8 混合和 D→E 三类碰撞质量记账检查，后者隔离化学求解，不代表完整 N-body 重算。
