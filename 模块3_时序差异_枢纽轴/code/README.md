# 模块3·前1/3 防御性脚本

> 对应技术路线文档「模块3 分期差异分析与hub轴锁定」步骤3.1–3.2
> 编制日期：2026-05-30

## 一、覆盖范围

| 步骤 | 内容 | 文件 |
|------|------|------|
| 3.1 | limma三组对比（hyperacute/acute/subacute vs control，全转录组） | `02_limma_three_stage.R` |
| 3.2 | 与1059全集 / 92高置信集双轨取交集 → stage DE-FRG | `03_intersect_ferroptosis.R` |

3.3 Mfuzz / 3.4 WGCNA / 3.5 sensitivity / 3.6 三重交集 hub锁定 → 后2/3 实现。

## 二、防御机制矩阵

### 2.1 数据偏移（drift）

| 风险 | 防御 | 文件 |
|------|------|------|
| 模块1输出缺失或schema不符 | 9条fail-fast硬断言 | `01_load_inputs.R` |
| 分期标签命名漂移 | 大小写归一化 + 白名单过滤 | `01_load_inputs.R` |
| 每分期样本数过少 | `MIN_N_PER_STAGE=3` 硬断言（control强制，其他分期缺失则跳过对应contrast） | `01_load_inputs.R` |
| 隐性混杂变量（批次/采样时间/性别） | `sva::num.sv` + `sva::sva()` 自动检测并入设计矩阵 | `02_limma_three_stage.R` |
| 异常样本拉偏DE | `limma::arrayWeights()` 自动降权，|z|>2样本入DRIFT日志 | `02_limma_three_stage.R` |
| 离群值膨胀t统计 | `lmFit(method="robust")` + `eBayes(robust=TRUE)` | `02_limma_three_stage.R` |
| 基因池版本依赖（hub基因池一变就换） | 1059全集 vs 92高置信双轨，Jaccard<0.5 报DRIFT | `03_intersect_ferroptosis.R` |

### 2.2 过拟合（实际是「signature列表膨胀」防御）

模块3阶段无模型训练，但**DE基因列表会喂给模块7建模**，因此防止列表被噪声/少数样本主导：

| 风险 | 防御 |
|------|------|
| 假阳性膨胀 | BH-FDR adjusted p-value（adj.P.Val<0.05） |
| 单样本主导 | arrayWeights自动降权 + 报告低权样本 |
| p-value分布异常 | 每contrast输出p-value直方图，<6%或>40%触发DRIFT告警 |
| 基因池选择bias | 双轨敏感性分析，仅保留两池都入选的「pool-stable」基因（保存到`_pool_stable_genes.rds`） |
| 完整过拟合防御 | 留给模块5（bootstrap）+ 模块7（LASSO 10折CV + one-SE规则 + 置换检验） |

## 三、前置依赖（**未满足前不能运行**）

1. **模块1未跑** → 缺 `模块1_预检/output/expr_mrna.csv`、`metadata.csv`
   - 我交付的模块1脚本只覆盖步骤1.1–1.3（下载+断言+漂移初筛），**未实现ComBat、CPM归一化、metadata.csv产物**
   - 需先补足模块1后段（log2(TPM+1)、ComBat、低表达过滤、stage标签对齐）并执行
2. **模块2已跑** ✅ → `模块2_铁死亡基因集/output/ferroptosis_geneset.csv` 与 `ferroptosis_geneset_high_confidence.csv` 已就绪

`01_load_inputs.R` 会在缺失任一依赖时fail-fast退出，不会带病前行。

## 四、执行命令

```bash
cd "/Volumes/拓展盘/论文编写/期刊论文/医学领域/个跨队列、跨平台验证的分期 signature/ 交付资料 /实验/模块3_时序差异_枢纽轴"
export PROJ_ROOT="$(pwd)"
Rscript code/run_module3.R
```

依赖：R + CRAN(`ggplot2`) + Bioconductor(`limma`, `sva`)

## 五、产物清单

| 路径 | 内容 |
|------|------|
| `output/_inputs_cache.rds` | 输入缓存（快速重跑） |
| `output/_limma_cache.rds` | 三contrast的limma fit缓存 |
| `output/stage_vs_control_DE_full.csv` | 全基因×3 contrast的DE结果长表 |
| `output/stage_DE_FRG_summary.csv` | 每contrast×每基因池的命中数+富集率 |
| `output/stage_DE_FRG_hits.csv` | 命中的具体基因+logFC+FDR |
| `output/_pool_stable_genes.rds` | 双池都入选的稳定基因列表 |
| `log/module3_<runid>.log` | 全量结构化日志（含ASSERT/DRIFT/STAB标签） |
| `log/pHist_<contrast>.png` | 每contrast的p-value直方图 |

## 六、调整决策树

| 触发条件 | 编号 | 调整 |
|---------|-----|------|
| 模块1输出schema不符 | D1 | 修补模块1后段后重跑 |
| 某分期n<3 | D2 | 该contrast被跳过；考虑分期重新分组（如hyperacute+acute合并） |
| `n.sv`>=设计矩阵秩 | D3 | 跳过SVA，先核查批次结构是否与stage混杂 |
| arrayWeights低权样本>样本总数20% | D4 | 怀疑全局质量问题，回到模块1重做QC |
| p-value直方图非均匀（左峰高>40%） | D5 | 检查是否有未建模的强混杂变量 |
| 全集vs高置信集Jaccard<0.5 | D6 | hub筛选必须用pool-stable基因（已自动保存） |

## 七、与原始技术路线的差异

| 点 | 原技术路线 | 本次实现 |
|---|----------|---------|
| 对比设计 | 三对相邻分期两两比较 | **stage vs control 三对**（与用户v1.2清单一致） |
| 阈值 | \|logFC\|>0.585, adj.P<0.05 | 同 |
| 拓展 | Mfuzz + WGCNA | 留给后2/3 |
