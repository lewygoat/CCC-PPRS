# 模块1 完整流程脚本（步骤1.1–1.5）

> 对应技术路线「阶段1 模块1：数据预处理」全部步骤
> 编制日期：2026-05-30

## 一、覆盖范围

| 步骤 | 内容 | 文件 |
|------|------|------|
| 1.1 | 三套GEO数据下载+缓存 | `01_preflight_download.R` |
| 1.1 | 结构断言（样本数、双矩阵、时间字段、GPL） | `02_structure_assert.R` |
| 1.1 | 漂移初筛（PCA离群、缺失率、自动log判定） | `03_drift_scan.R` |
| 1.4 | 分期标签提取+样本配对表生成 | `04_metadata_align.R` |
| 1.2/1.3/1.5 | mRNA/miRNA归一化 + ComBat去批次 + 低表达过滤 | `05_normalize_combat.R` |

## 二、执行命令

```bash
cd "/Volumes/拓展盘/论文编写/期刊论文/医学领域/个跨队列、跨平台验证的分期 signature/ 交付资料 /实验/模块1_预检"
export PROJ_ROOT="$(pwd)"
Rscript code/run_module1.R
```

依赖：
- CRAN：`ggplot2`、`matrixStats`
- Bioconductor：`GEOquery`、`Biobase`、`sva`、`edgeR`、`limma`

## 三、产出物

| 路径 | 内容 |
|------|------|
| `data_raw/<GSE>.rds` | getGEO缓存 |
| `data_raw/_eset_list.rds` | 三套合并 |
| `data_raw/_meta_all.rds` | 每slot的metadata |
| `output/metadata.csv` | **训练集（GSE296792）配对表** ← 模块3硬依赖 |
| `output/metadata_external.csv` | 外验集（GSE125512）配对表 |
| `output/expr_mrna.csv` | **训练集mRNA归一化+ComBat后矩阵** ← 模块3硬依赖 |
| `output/expr_mirna.csv` | **训练集miRNA归一化+ComBat后矩阵** ← 模块4硬依赖 |
| `output/expr_mrna_external.csv` | 外验集mRNA矩阵 ← 模块5使用 |
| `pre_check/<GSE>_slot<i>_stage_mapping.csv` | **分期标签推断结果**（人工核对用） |
| `log/preflight_*.log` | 全量日志 |
| `log/preflight_figs/PCA_*.png` | 归一化前后PCA图 |

## 四、防御机制

### 4.1 数据偏移（drift）

| 阶段 | 机制 |
|------|------|
| 1.1下载 | 三次指数退避重试 + 落盘缓存 |
| 1.1结构 | 9条fail-fast断言（样本数、双矩阵、GPL、时间字段） |
| 1.1漂移初筛 | top-var 2000基因 PCA + 多维距离z-score（|z|>3报DRIFT） |
| 1.4分期 | 4套正则规则推断stage（hyperacute/acute/subacute/control）；推断失败入WARN日志，要求人工review |
| 1.2归一化 | 自动dtype检测（counts/intensity/log_intensity），按类型分发normalize方法 |
| 1.3miRNA | CPM>1在≥30%样本的低表达过滤（防止低counts噪声主导） |
| 1.2ComBat | stage作为协变量保留生物信号；批次<2时自动跳过 |
| 1.2归一化后 | 再做一次PCA，PC1>60%入WARN（残余批次怀疑） |

### 4.2 stage标签覆盖机制

`04_metadata_align.R`正则推断不准时：

1. 检查`pre_check/<GSE>_slot<i>_stage_mapping.csv`找出inferred_stage=NA的样本
2. 在`pre_check/`下手写`stage_override.csv`，两列：`sample_id`、`stage`
3. 重跑`05_normalize_combat.R`（或整个`run_module1.R`），覆盖会自动加载

### 4.3 「过拟合」相关

模块1**无模型训练**，过拟合不适用。完整防御推后到模块7（LASSO 10折CV + one-SE规则 + 置换检验）。

## 五、与v1.2阶段1清单对照

| v1.2步骤 | 实现 |
|---------|------|
| 1.1 GSE296792下载、boxplot+PCA | ✅ PCA（boxplot留作可选） |
| 1.2 mRNA log2 + ComBat | ✅ 自动dtype检测+TMM+log2(CPM+1)+ComBat |
| 1.3 miRNA CPM + 低表达过滤 | ✅ |
| 1.4 stage标签对齐、metadata.csv | ✅ |
| 1.5 GSE125512同流程独立处理 | ✅ 自动分流到`expr_mrna_external.csv` |
| 阶段0 SampleSheet预检 | ✅ 合并到`02_structure_assert.R` |

## 六、调整决策树

| 触发 | 编号 | 动作 |
|------|-----|------|
| 任一GSE下载失败 | D1 | 检查网络，重跑（缓存自动复用） |
| 双矩阵断言FAIL | D2 | 手动从GEO suppl下载miRNA Agilent txt，改用`limma::read.maimages()` |
| stage推断NA>20% | D3 | 写`pre_check/stage_override.csv`，重跑 |
| ComBat后PC1>60% | D4 | 检查批次变量是否选对，可能要加入采样时间或患者ID |
| miRNA低表达过滤后<200 features | D5 | 放宽到CPM>0.5或SAMPLE_FRAC=0.2 |
| 不同GSE的样本ID冲突 | D6 | 在合并步加`gse`前缀（脚本已隔离不同GSE的处理） |

## 七、运行时间估计（在干净R环境）

| 阶段 | 时间 |
|------|------|
| Bioconductor首装（GEOquery+sva+edgeR+limma） | 20–40 min |
| GEO三套数据下载 | 10–30 min（取决于网速） |
| 结构断言+漂移初筛 | <1 min |
| metadata对齐 | <1 min |
| 归一化+ComBat | 2–5 min |
| **总计（含首装）** | **约30–80 min** |

重跑时缓存命中，主要时间在ComBat（2–5 min）。

---

*README更新，覆盖完整模块1流程*
