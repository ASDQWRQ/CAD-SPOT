# 冻结评分生成链只读审计

审计日期：2026-11-15。

范围：记录此前已经完成的本地只读核查，不扩展调查、不重新运行长分析。原软件、文稿和数据均未修改。本记录单独保存于新的 follow-up 目录。

## 核心结论

1. 当前冻结输入与旧工程坐标表逐字节一致，可以确认表级来源。
2. 可以从冻结表中两个已存模块精确恢复目标评分的最后一步，但未找到这些模块及外部参照从 counts 生成的原始代码。
3. 声明的目标/参照基因集为 23/21；当前公开 filtered 矩阵只含其中 19/20。声明集合大小不能作为旧评分实际使用全部基因的证据。
4. 冻结评分的调整后 rho = 0.199254 与重新计算评分的 rho = 0.324403 不是同一评分输入的两个估计。差异不能仅归因于缺失基因，也不能仅归因于最终 Z 标准化。
5. 现有材料支持从冻结表复现 0.199，不支持从公开 counts 精确重建 0.199。新评分审计仍沿用冻结组成模块，因此也不是完整的 counts 端到端重建。

## 1. 冻结表来源与 Hash

当前冻结输入：

`D:\CADSPOT_REVISION_20261115\Source_Data\archived_coordinate_table.tsv.gz`

旧工程来源表：

`D:\BMC_Bioinformatics_CADSPOT_REBUILD_20260730\Source_Data\cadspot_candidate_coordinate_table.tsv.gz`

两份压缩文件的 SHA-256 均为：

```text
ef2155d233e21cd12e0bf8dc3b03c7a061e00895f781a967a309ab6d818f0a12
```

解压后内容的 SHA-256 也一致，均为：

```text
7b19a5f0166632528098ac0b3b00a238809e9e8b261af90fff8f891a32412f10
```

冻结表包含 23,013 个 spot、5 张切片、79 列。Hash 证明文件身份，不证明上游评分生成过程正确或已被重建。

## 2. 已核实的关键公式

对全部 23,013 行进行了短时数值核对：

```text
target_score = MMS_z
reference_score = External_HG_zero_overlap_z

u_i = (Metabolic_stress_i + Mitochondrial_injury_i) / 2
target_score_i = (u_i - mean_all_spots(u)) / SD_all_spots_ddof0(u)
```

上述两个列别名逐值一致，最大绝对差为 0。目标评分公式的最大绝对误差为 `8.881784197001252e-16`。

辅助核对：`mean(u) = 1.312217714761407e-18`，`SD_ddof0(u) = 0.2818418301001304`。采用样本标准差 `ddof=1` 时不再精确相等。

另一个容易混淆的已存列满足：

```text
HG_score = (Hypoxia + Glycolysis) / 2
```

该公式的最大绝对误差为 `2.220446049250313e-16`。`HG_score` 不能据此当作本次主分析实际使用的 `reference_score`；主分析使用的是 `External_HG_zero_overlap_z`。

重要边界：以上是从已存数值恢复的代数关系，不是找到了原始生成脚本。`counts -> Metabolic_stress / Mitochondrial_injury / External_HG_zero_overlap_z` 仍有来源缺口。

## 3. 基因集与缺失基因

旧工程的声明基因集位于：

`D:\BMC_Bioinformatics_CADSPOT_REBUILD_20260730\Tables\Supplementary_Tables_S1-S10.xlsx`

工作表：`S1_gene_sets`。其 23/21 个基因与新版 `Analysis\prepare_pdac.py` 第 15-16 行的 `TARGET`、`REFERENCE` 完全一致。

| 集合 | 声明数量 | 当前矩阵可用数量 | 缺失基因 |
|---|---:|---:|---|
| Target | 23 | 19 | LDHA; ALDOA; PKM; TPI1 |
| Reference | 21 | 20 | GAPDH |

直接检查了以下五份 H5 的 feature 名称及对应 ID；五份文件均缺少上述基因条目：

```text
D:\CADSPOT_REVISION_20261115\Data\GSE274103\raw\GSM8443449_PDAC-p1_filtered_feature_bc_matrix.h5
D:\CADSPOT_REVISION_20261115\Data\GSE274103\raw\GSM8443450_PDAC-p2_filtered_feature_bc_matrix.h5
D:\CADSPOT_REVISION_20261115\Data\GSE274103\raw\GSM8443451_PDAC-p3_filtered_feature_bc_matrix.h5
D:\CADSPOT_REVISION_20261115\Data\GSE274103\raw\GSM8443452_PDAC-p4_filtered_feature_bc_matrix.h5
D:\CADSPOT_REVISION_20261115\Data\GSE274103\raw\GSM8443453_PDAC-p5_filtered_feature_bc_matrix.h5
```

这些是 GEO 公开的 filtered 计数矩阵，而不是原始测序 reads。缺少 feature 行不能解释为该基因表达为零。也不能仅凭旧基因集列出了 23/21 个基因，就断定旧评分实际使用了全部基因；旧流程是否已跳过缺失基因，尚未核实。

公开来源：NCBI Gene Expression Omnibus，GSE274103，样本 GSM8443449-GSM8443453。

本地来源记录：`D:\CADSPOT_REVISION_20261115\Data\GSE274103_family.soft.gz`。

下载入口记录：`D:\CADSPOT_REVISION_20261115\Analysis\download_data.py`。

## 4. 新评分的确定计算链

代码：`D:\CADSPOT_REVISION_20261115\Analysis\prepare_pdac.py`，特别是第 88-103 行。

```text
公开 filtered counts
-> 按基因符号合并重复行，并匹配冻结表的 sample_id/barcode
-> 每个 spot 归一化到 10,000 counts
-> log1p
-> 每个基因跨全部 23,013 个分析 spot 做 Z 标准化，ddof=0
-> 对明确可用的 19 个目标基因或 20 个参照基因等权平均
```

不对缺失基因补零。校准实际调用的共享评分函数为：

`D:\CADSPOT_REVISION_20261115\Software\cadspot_null.py`，`score_gene_set()`，第 106 行起。

该函数接收已经处理的表达矩阵，不在内部隐式归一化。实际调用采用 `mean_gene_z`。

新旧评分一致性已有文件记录：

`D:\CADSPOT_REVISION_20261115\Source_Data\score_reconstruction_audit.tsv`

| 评分 | 冻结与新重建 Spearman | 冻结与新重建 Pearson |
|---|---:|---:|
| Target | 0.5720166297912258 | 0.5652122468026664 |
| Reference | 0.5217542775554649 | 0.5035070416585818 |

这些结果表明两者并不相同。冻结目标的已核实定义是两个已存模块等权合成；新审计直接对可用基因等权平均。旧模块本身的归一化、基因权重、实际使用基因和标准化范围仍未查明。

## 5. Cross-score 对照：rho = 0.199 与 0.324

本次只读核查使用已有文件中的新旧评分列，保持相同的 23,013 个 spot、同一协变量表、切片固定效应与残差算法，进行了四组短时计算。未重新运行匹配随机集或 bootstrap。

| Target 列 | Reference 列 | 原始 Spearman rho | 组分调整后 rho |
|---|---|---:|---:|
| target_score | reference_score | 0.5395119057478159 | 0.1992540732014225 |
| fresh_target_score | reference_score | 0.4954705390847205 | 0.3607343771945077 |
| target_score | fresh_reference_score | 0.3558736011789769 | 0.2010467312782161 |
| fresh_target_score | fresh_reference_score | 0.5536148928489854 | 0.324402948638828 |

四行对照的数值输入来自：

`D:\CADSPOT_REVISION_20261115\Source_Data\raw_matrix_pair_scores.tsv`

四行汇总是在只读核查时计算的，原文件并不预存这张四行汇总表。本记录保存该核查结果。

技术协变量：`log_nCount`、`log_nFeature`。

组成协变量：`Epithelial_malignant_proxy`、`ECM_Stroma`、`CAF`、`Macrophage`、`T_cell`。

模型：每个未排名评分分别对标准化协变量及切片固定效应作带截距 OLS，然后计算两组残差的 Spearman 相关。

实现：`D:\CADSPOT_REVISION_20261115\Software\cadspot.py`，`design()` 第 68 行起、`fit_pair()` 第 115 行起。

独立的已有主结果记录：

- `D:\CADSPOT_REVISION_20261115\Source_Data\fresh_associations.tsv`：冻结评分的 0.199254。
- `D:\CADSPOT_REVISION_20261115\Source_Data\fresh_matched_null_500_summary.tsv`：新评分的 0.324403。
- `D:\CADSPOT_REVISION_20261115\Source_Data\fresh_run_manifest.json`：冻结输入 Hash、实际参数和协变量。

解释边界：在这组固定设计的对照中，差异来自评分输入变化，不能归咎于更换了协变量或残差算法。但四组对照不是对缺失基因、模块权重或归一化差异的独立因果分解，不能据此给它们分别分配贡献。最终整体 Z 标准化是正向仿射变换，本身不改变这里的相关结果。

## 6. 已找到与未找到的代码

已找到的下游入口：

- `D:\BMC_Bioinformatics_CADSPOT_REBUILD_20260730\Software\cadspot.py` 第 105 行读取评分表，第 107-114 行计算原始和调整后相关；不从 counts 生成原评分。
- `D:\CADSPOT_REVISION_20261115\Analysis\run_reanalysis.py` 第 32 行起读取冻结评分表；第 51-64 行定义独立表达矩阵审计。
- `D:\CAD-SPOT_GitHub` 的可用 Git 历史中，初始版本已包含处理后评分表；所检查的历史中未找到更早的原评分生成代码。

本次检查范围限于指定主工程、Git 仓库及 D 盘此前 CADSPOT/BMC 工程和相关已有压缩包，没有全盘大数据扫描。E 盘未挂载，不能据此断言未挂载盘上不存在原始代码。

仍缺少的原始上游材料：

- `Metabolic_stress`、`Mitochondrial_injury` 的 counts 生成脚本及实际使用基因记录。
- `External_HG_zero_overlap_z` 的生成脚本、具体基因处理和标准化范围。
- 五个组成模块的生成脚本、基因定义和标准化记录。
- 原始处理对象或等价的中间矩阵，以及能连接这些对象与冻结表的执行记录。

旧稿文字称采用 z-mean 原则，但仅有文字不足以确认实际执行过程。可读旧稿位于：

`D:\BMC_Bioinformatics_CADSPOT_REBUILD_20260730\work\old_manuscript_recovered\recovered_text.txt`，尤其第 53、61、63 行。

## 7. 后续使用边界

目前可以准确表述为“从发布的冻结评分表复现主分析”，不能表述为“从 GEO counts 精确重建原评分及主分析”。

此外，新评分的 0.324403 审计仍读取冻结表中的五个组成模块，见 `Analysis\run_reanalysis.py` 第 62 行及 `Software\cadspot_null.py` 第 200 行起。因此，目标和参照评分虽可从公开矩阵重新计算，整个调整流程的上游组成评分仍未独立重建。

若找不到原始生成链，合理选择是另行建立并明确命名一个可完整重建的主示例；不能把新定义的 19/20 基因评分追认为旧评分的精确复现。本记录不执行该变更。
