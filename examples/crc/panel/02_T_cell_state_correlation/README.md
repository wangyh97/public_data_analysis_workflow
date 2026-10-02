# Fig2：肿瘤表达与同样本 T 细胞状态

目的：关联恶性标签细胞中的 METTL3、STRAP、PTBP1 表达与同一份肿瘤样本中的 T 细胞状态。主队列 GSE236581 Pre（19 人），验证队列 GSE205506 Pre（8 人）。

| 状态模块 | 基因 | 热图细胞群 |
|---|---|---|
| Cytotoxic | PRF1、GZMA、GZMB、GZMH、NKG7、CCL5 | CD8 T 与全部 T 分别计算 |
| Dysfunction | PDCD1、LAG3、HAVCR2、TIGIT、TOX、ENTPD1 | CD8 T 与全部 T 分别计算 |
| Stem/memory | TCF7、CCR7、IL7R、LEF1 | CD8 T 与全部 T 分别计算 |
| Treg program | FOXP3、IL2RA、CTLA4 | 全部 T 细胞 |

遵循原配置不增删 T 状态基因。All-T Treg program 是全部 T 细胞的表达分数，不是已分离 Treg 亚群内的功能强度，也不是 Treg 比例本身。Dysfunction 分数包含可随激活上调的分子，不能单靠该分数确诊功能耗竭。

GEO：全部 T=`MajorCellType == T`；CD8 T 进一步要求作者 SubCellType 包含 CD8。scCT：按数据库 T/CD4/CD8/Treg 等标签识别全部 T，再按 CD8 标签识别 CD8 亚群。细节在共享 process.R；没有重新做 T 细胞亚群注释。

`figure/Fig2_T_cell_state_heatmap.pdf/png` 包含 7 行，前 3 行是 CD8，后 4 行是全部 T。目标表达始终来自肿瘤细胞，纵向状态始终来自行名指定的 T 细胞群，不能把该图理解为肿瘤细胞表达 T 细胞状态基因。

`source/patient_values.csv` 记录逐患者 x、y、匹配样本和两侧细胞数；`correlation_statistics.csv` 为热图数据；`gene_sets.csv` 冻结实际 gene set。数据处理先按同一肿瘤样本连接，随后患者层面等权平均。具体运行方式见上级 README。

每队列 21 次相关分别做 BH。GSE205506 的 PTBP1–All-T Treg program rho=0.9286，精确 p=0.002232，q=0.046875；主队列未复现，其余 q≥0.05。CD8 和 All-T 共用细胞，不是额外独立队列。本轮也未使用外周血 T 细胞，因此不能推广到患者全身体内所有 T 细胞。
