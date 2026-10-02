# CellChat 扩展：设计、配置与复现

这是前一轮报告中未执行的 CellChat 扩展。使用真实 CellChat 推断；前一轮患者级 MHC/APM 和疗效比较不因本次扩展而改变。结果及实际完成情况见 `CellChat补充报告.md`。

## 固定设计

输入为完整 GEO GSE236581 的 19 个治疗前肿瘤样本，以及 scCT-DB GSE205506 的 8 个治疗前肿瘤样本，每样本一位患者。平台同源数据不重复作为独立队列。本次不包括外周血队列和仅 1 位患者的队列。

直接沿用 `tables/patient_high_low_groups.csv` 中三个基因的患者分组：恶性细胞全部合格细胞的平均 log1p(CP10K) 高于研究内中位数为 High，否则 Low；分组不由通讯结果决定。用于分组的表达值不依赖 CellChat 抽样。

按作者注释定义互不重叠的 Malignant、CD8_T、CD4_T、Treg、Other_T、Myeloid、B_plasma、Endothelial、Fibroblast、NK_ILC。不明确恶性状态的其他上皮群不被标为 Malignant。每样本至少 30 个恶性细胞，其他群至少 20 个细胞。每个合格细胞群固定随机抽取最多 300 个细胞；seed=20260905，共 46,015 个细胞。完整名单保存在 `data/CellChat/selected_cells.csv`，群体资格见 `audit/cellchat_cell_group_eligibility.csv`。

重新从完整矩阵提取上述细胞的**全部基因**，使用全基因文库大小进行 CP10K、log1p 规范化。每样本逐细胞核对文库大小以及 METTL3/STRAP/PTBP1 计数，要求与前一轮独立提取完全一致。之后才调用 CellChat `subsetData()`，保留数据库涉及的信号、复合物和辅因子基因。

## CellChat 参数

使用 [jinworks/CellChat 官方源码](https://github.com/jinworks/CellChat/tree/75253cd0c9e68410e6e721a6d3a0419a1d7e358f)，固定 commit `75253cd0c9e68410e6e721a6d3a0419a1d7e358f`，版本 2.2.0.9001。安装包由 [r-universe](https://blaserlab.r-universe.dev/CellChat) 从该上游源码构建。官方方法见 [CellChat Nature Protocols](https://doi.org/10.1038/s41596-024-01045-4)。

数据库为随包的 CellChatDB.human，仅保留 Secreted Signaling、ECM-Receptor、Cell-Cell Contact，排除 Non-protein Signaling；具体条目和冻结 RDS 均保存。使用 `identifyOverExpressedGenes(do.fast=FALSE)`、`identifyOverExpressedInteractions()`、`computeCommunProb(type='triMean', raw.use=TRUE, population.size=FALSE, nboot=100, seed.use=20260905)`，随后 `filterCommunication(min.cells=20)`、`computeCommunProbPathway()`、`aggregateNet()`。其余值为该固定源码的默认值，实际函数和参数保存在 audit/ 与结果 RDS 中。未使用 PPI 投影、空间距离或 UMAP。

`raw.use=TRUE` 在这里指 CellChat 的未投影规范化表达，不是把原始 counts 直接输入模型。关闭 population.size，避免把通讯强度直接乘以抽样后的群体比例。CellChat 仍会按样本最大信号表达值缩放，因此跨患者强度是相对的模型分数。

## 患者级比较

- 主要展示恶性细胞与 CD8_T、CD4_T、Treg 的两个方向，保留其他肿瘤相关网络用于描述。
- 每个样本单独推断。通讯强度为通过 CellChat 内部 P≤0.05 门槛的配体–受体概率之和，即标准 pathway 输出。群体存在但未检出通讯记为 0；群体缺失或不达细胞数门槛记为不可用，绝不补零。上游 `aggregateNet()` 的描述性总网络另用严格 P<0.05；本次保留原函数行为并明确此边界差异，患者通路检验统一使用 pathway 的 ≤0.05。
- High/Low 检验至少每组 3 位患者。统计量为两组平均秩之差；若所有分组组合数 ≤10,000，枚举全部组合做双侧精确检验（8 位患者 4/4 分组共 70 种）。否则按患者标签置换 10,000 次，经验 P 的分子/分母均加 1。实际模式和次数保存在结果表。连续表达另作患者级 Spearman：n≤9 且两变量均无 ties 时用精确 P，否则用渐近 P；常量分数的 rho 不可估计，记 NA，P=1 保留在校正 family 中；热图用短横线标明，不把不可估计相关显示为 rho=0。
- 每研究将 3 基因 × 全部通路 × 6 个方向作为探索性 BH family。原始强度和除以患者全网络总强度的相对强度分别校正，后者仅作尺度敏感性分析。
- 预先单列 HLA-A/B/C→CD8A/B 六种经典边在 Malignant→CD8_T 方向的强度；三个基因为一个 focused BH family，不与全通路探索 family 混读。

数据库完整 “MHC-I” 通路包含非经典 HLA、NK 受体及 MICA/MICB/ULBP 等压力配体。经典边的单列有助于解释，但它们同样不能测量抗原肽、TCR 特异性或实际细胞接触。CellChat 内部细胞置换 P 不作为患者间 High/Low 的 P 值。

## 环境与命令

依赖安装在本项目 `runtime/R-library`；没有更新系统 R 或原用户 R 库。使用现有 R 4.4.2，CellChat 的预编译包构建于 R 4.5.3，R 会显示构建版本警告。本次必须以实际加载、矩阵核对、完整推断和概率范围检查确认运行情况；不能只凭安装成功声称兼容。最终版本见 `audit/CellChat_final_sessionInfo.txt`，安装清单与上游来源见 `audit/cellchat_install_plan.txt`、`cellchat_binary_provenance.json`。

```powershell
# 已安装依赖且数据仍在原输出目录时：
.\run_cellchat.ps1 -Mode all -Workers 4
# 只重算患者比较和图表：
.\run_cellchat.ps1 -Mode compare
# 只继续尚未完成的逐样本推断：
.\run_cellchat.ps1 -Mode infer -Workers 4
# 新环境需要安装时（仅项目库）：
.\run_cellchat.ps1 -Mode install
```

脚本 18 固定样本与抽样，19/20 提取全基因稀疏矩阵，21 调用 CellChat，22 作患者检验与 R 绘图，23 为限并发逐样本调度。默认 4 个 R 进程；内存不足时 Workers=1。每样本使用约 1–2 GB 内存，耗时随候选 LR 和细胞数变化。全矩阵提取需额外读取原始 GEO 的 13.1 亿条非零记录。

中间数据在 `data/CellChat/`，紧凑但保留准确 net/netP/LR/options 的结果在 `derived/CellChat/`。已完成结果可复用；**改变抽样或参数后必须使用新的结果目录或先备份再移走旧检查点**，否则会复用旧推断。分享 ZIP 不包含 data、derived、runtime；仅解压分享包不能直接执行 compare。

本扩展沿用作者 QC，不重跑从未过滤 droplets 开始的 doublet/ambient RNA 处理。脚本 10 的慢速全量读取仍不是必需步骤，其原本用途为独立读取器交叉核对。
