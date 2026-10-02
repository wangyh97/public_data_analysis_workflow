# Fig1：肿瘤表达与 MHC-I 相关基因/步骤

目的：检查 METTL3、STRAP、PTBP1 的肿瘤细胞表达是否与 MHC-I 结构表达及各加工呈递步骤相关；不预设方向或显著性。主队列 GSE236581 Pre（19 人），验证队列 GSE205506 Pre（8 人）。

| 热图行 | 基因 |
|---|---|
| MHC-I structure | HLA-A、HLA-B、HLA-C、B2M |
| Immunoproteasome | PSMB8、PSMB9、PSMB10 |
| Peptide transport into ER | TAP1、TAP2 |
| ER peptide trimming | ERAP1、ERAP2 |
| Folding and peptide loading | TAPBP、CALR、CANX、PDIA3 |
| Transcriptional regulators | NLRC5、STAT1、IRF1、NFKB1、NFKB2、RELA、RELB、REL |

保留用户给出的全部基因。原 `peptide_loading_and_ER_processing` 拆为转运、修剪、折叠装载三步。NF-κB 按用户确认加入五个家族成员；转录调控这一行是这些基因的等权表达分数，不是 NF-κB 活性。

`figure/Fig1_MHC_I_module_heatmap.pdf/png`：六模块热图。

`figure/Fig1_scatter_<目标基因>_<结构基因>.pdf/png`：12 张独立散点图，每张都有两个队列。散点纵轴是对应结构基因的表达，横轴是目标基因表达；两者都来自恶性标签细胞。每点为患者，灰线为描述性线性拟合；Spearman rho、p、q、n 位于图上方，避免遮挡数据。

必要作图数据为 `source/patient_values.csv` 和 `source/correlation_statistics.csv`；`gene_sets.csv` 冻结本次实际分组。`sample_compartment_expression.csv`、`input_sample_audit.csv`、`gene_coverage.csv` 提供清洗审计。统计表 CI 为 2,000 次患者 bootstrap 的 percentile 95% CI。

`script/config.json` 改基因集/阈值，`script/data_processing.R` 重建 source，`script/plot.R` 的 STYLE 区改外观。完整运行说明和源数据路径见上级 README。

检验组为每队列 18 次模块相关、12 次结构单基因相关，各自 BH；本轮均无 q<0.05。保留所有非显著结果。两个队列肿瘤注释来源不同，GSE205506 malignant 标签没有在本轮做独立 CNV 验证。
