# Fig3 | Tumor gene expression before and during/after treatment

**统计更新：当前图使用双侧配对 t 检验，新增 4 张配对箱线+散点图。完整说明见 [配对t检验说明](配对t检验说明.md)。下方旧 Wilcoxon 说明仅对应保留的敏感性表。**

## 本次实际绘制

每位患者一条轨迹，每个点为一个真实治疗节点；responder 红色、non-responder 蓝色，另用圆点/三角形辅助辨认。分别提供两队列三基因总览，以及每个基因的逐患者小图，避免多条红线交叉后无法识别患者。现共 12 张图（含新增 4 张箱线图），PDF + 600 dpi PNG。

| 研究及矩阵来源 | 原始节点 | 进入纵向图的患者 | 各节点患者数 |
|---|---|---|---|
| GSE236581，GEO 矩阵/作者注释 | I=Pre；II；III；IV | 13：12 R、1 NR | Pre 13；II 12；III 7；IV 2 |
| GSE205506，scCT-DB 矩阵/注释 | Baseline=Pre；Post-treatment=Post | 8：7 pCR、1 non-pCR | Pre 8；Post 8 |

**GSE205506 没有已确认的 II、III、IV 取样序列。** 本轮另外下载 GEO 官方 family SOFT 并逐项核对：211/213、241/243 等是样本标识，同一患者对应治疗后/未治疗样本，不能把末位 1/3 当作阶段 I/III。16 个用于绘图的 GSM 的患者和 Pre/Post 标签全部与 GEO 匹配。完整官方记录和精简样本表保存在 source 中。GSE236581 的 II/III/IV 保留作者原始序数，不等同于第 2/3/4 周；本图没有假定均匀的实际时间间隔。

## 细胞范围和患者筛选

- 延续 Fig1/2 的范围：**肿瘤组织内的恶性标签细胞**，不是整份肿瘤组织的所有细胞混合 pseudobulk。
- GSE236581 限作者 `c91_Epi_Tumor`；GSE205506 限 scCT-DB `Malignant epithelial cell`。后者没有在本轮独立做 CNV 验证；尤其 pCR 后标为 malignant 的细胞不能据此认定是存活残余癌细胞。
- 仅 CRC 肿瘤组织；正常组织、血液、未映射 CRC 患者不纳入。相同研究内重复样本仅保留首个版本，不把两数据库重复数据当作独立患者。
- 每个患者节点合计至少 30 个恶性标签细胞、全基因文库大于 0 才合格。与此前“每个样本先达到 30”不同，本轮按用户指定的患者节点合计后筛选，阈值适用于一个 pseudobulk 点。
- 纵向图要求一个合格 Pre 加至少一个合格后续节点。不要求四节点齐全；没有补零、插值或借用其他患者表达。
- 缺失中间节点时，两个真实观察点之间使用虚线；缺失节点本身无点。缺失可能来自未取样或未获得足够恶性细胞，不能解释为表达值为 0。
- `patient_pseudobulk_all.csv` 保留全部患者节点；未入图者、细胞不足节点及无恶性细胞样本均可在 source 审计表追溯。

## Pseudobulk 的明确计算

对于患者 p、节点 t、基因 g：

```
count[p,t,g] = 该患者该节点所有纳入恶性标签细胞的基因 g 原始 counts 之和
library[p,t] = 同一批细胞全部基因原始 counts 之和
CPM[p,t,g] = count[p,t,g] / library[p,t] * 1,000,000
图上表达 = log2(CPM + 1)
```

同一患者同一节点有多份肿瘤样本时，先相加 counts 和 library，再计算 CPM；不是平均样本 logCPM，也不是 Fig1/2 的细胞平均 log1pCP10K。多样本合并会受到各部位回收细胞量的影响，因此样本 ID、样本数和每份样本计数均保留。

使用的是简单 library-size CPM，无 TMM、无全基因 edgeR 模型；这里用于目标基因轨迹的可复现展示。相对表达下降不能直接解读为每个细胞绝对转录本减少或肿瘤负荷下降。

## 数据、脚本、运行

`source/trajectory_plot_data.csv` 是全部作图必要数据，包含 cohort、patient、response、stage、gene、raw_count、full_gene_library、CPM、expression、恶性细胞数及原始样本。其它关键表：

- `sample_pseudobulk_counts.csv`：每份样本的原始基因计数及全基因文库。
- `patient_pseudobulk_all.csv`、`patient_eligibility.csv`：所有患者节点、纵向纳入情况。
- `timepoint_availability.csv`：各阶段/疗效组实际人数。
- `change_from_baseline.csv`：每个可配对后续节点相对该患者 Pre 的差值。
- `paired_change_statistics.csv`：探索性配对变化检验；不在折线图上堆叠显著性标注。
- `input_sample_audit.csv`、`gene_coverage.csv`、`input_manifest.csv`：筛选、基因覆盖、原始输入与 SHA256。
- `GSE205506_GEO_sample_catalog.csv`、`timepoint_mapping.csv`：真实节点依据。

从此条目目录运行 PowerShell：

```powershell
$rExe = 'Rscript'  # 或指定本机 Rscript 可执行文件
# 从已有完整计数提取缓存清洗、汇总、分析
& $rExe ./script/data_processing.R
# 从 source CSV 重画（修改外观时只需要这一步）
& $rExe ./script/plot.R
# 如需重新读取原始提交矩阵：
& $rExe ./script/data_processing.R --rebuild-raw
```

脚本可从任意工作目录用绝对路径调用。`data_processing.R` 已完整写出筛选、汇总和分析；原始矩阵读取复用上级 `_shared/extract_geo.mjs` 与 `_shared/extract_scct.py`。原始输入沿用项目 `data/GEO`、`data/scCT-DB` 与两个公开样本注释表，具体路径见上级 README 和 input_manifest。需要移动完整分析目录时保留这些相对位置，并在项目 `analysis_config.json` 修改本机 R/Python/Node 路径。

R 4.4.2，依赖 data.table、jsonlite、ggplot2；原始流式读取另需已有 Python/numpy、Node 运行时。包精确版本见 source 与 figure 中的 sessionInfo。PDF 用 R Cairo；PNG 用 R Cairo 600 dpi。SVG 为可选接口，默认不需要 svglite。

`script/plot.R` 开头 STYLE 区控制红蓝颜色、线宽/透明度、点大小/形状、图尺寸、逐患者列数、坐标是否统一及是否画缺失节点虚线。默认逐患者小图共享 y 轴，不采用平滑曲线或组均值替代患者轨迹。基因和阈值配置在 `script/config.json`。

## 结果及解释范围

GSE236581 的 Pre→II 有 12 对患者，Pre→III 有 7 对，Pre→IV 只有 2 对；GSE205506 有 8 对 Pre→Post。GSE236581 纵向总人数是 13，而不是 12，因为有人仅有合格的较晚节点。阶段越晚的患者集合不同，不能通过不同阶段人群均值推断同一人持续下降。

总体配对差值的中位数：

| 队列/比较 | METTL3 | STRAP | PTBP1 |
|---|---:|---:|---:|
| GSE236581 Pre→II | −0.049 | −0.225 | −0.193 |
| GSE236581 Pre→III | −0.231 | −0.345 | −0.062 |
| GSE205506 Pre→Post | −0.098 | −1.666 | −0.314 |

单位为 `log2(CPM+1)` 的差值。患者方向并不一致，不能把这些中位数称为所有患者一致下降。总体配对比较均未达到 BH q<0.05。

探索性双侧 Wilcoxon signed-rank 使用患者差值，exact=FALSE，默认连续性校正；至少 3 对才检验。BH 按队列和分析组（全部患者/R/NR）分别校正该组所有实际可检验的基因×节点比较。2 人的 IV 与 1 人 NR 不给出可靠群体推断。

两队列纵向样本都仅有 1 位 NR，无法据此可靠检验“治疗变化与疗效存在关联”。红蓝用于标识患者的最终疗效，不能代替 time×response 交互作用检验。GSE236581 使用 RECIST R/NR，GSE205506 使用 pCR/non-pCR，两终点不合并。

本轮图像生成前有两个可见的代码校验中断（数据库标签 `Non_response` 拼写、data.table 逻辑列筛选），均修正后完整重跑成功；未使用失败中断的中间结果绘图。

## 来源

- [GSE236581 官方数据](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE236581)
- [GSE205506 官方数据](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE205506)
- [GSE205506 原始研究](https://doi.org/10.1016/j.ccell.2023.04.011)
- [scCT-DB](https://scctdb.ncpsb.org.cn/)
