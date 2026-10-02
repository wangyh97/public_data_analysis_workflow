# Fig3 配对 t 检验更新

当前图上的统计标注均改用双侧配对 t 检验。旧的 `paired_change_statistics.csv` 为 Wilcoxon 敏感性结果单独保留，不是当前图中的 p/q 来源。现在共 12 张图：2 张带统计标注的队列轨迹总览、6 张逐患者图、4 张配对箱线+散点图；PDF + 600 dpi PNG。

## 配对与统计

`script/paired_t_tests.R` 按 cohort、patient、gene 显式连接 Pre 与各后续节点，执行 `t.test(post, pre, paired=TRUE, alternative='two.sided')`。图上是**全部配对患者**的检验，红蓝用于标识疗效，不代表分别检验两组。单个患者图无可用于患者层面 t 检验的独立重复，因此不标个体 p 值。

检验目标是患者差值的均值是否为 0。`source/paired_t_test_statistics.csv` 提供均值差、差值 SD、t、df、95% CI、p、q；`source/paired_t_test_pairs.csv` 保留逐患者两端值。BH 校正按队列覆盖 3 基因×后续节点：GSE236581 9 项，GSE205506 3 项。另在表中提供 R/NR 分层结果，按 cohort×group 独立校正；NR 只有 1 人，不可检验。未计算 time×response 交互作用。

IV 的 2 对患者可在数学上计算 t 检验（df=1），本轮列出但标为高度不稳定。小样本配对 t 检验依赖患者差值近似正态且不被异常差值主导，2 人无法检查该假设。所有主比较均 p≥0.05，BH 后亦无显著变化；不显著不等于证明前后完全相同。

## 箱线图为何单独按比较绘制

箱线图可补充分布信息，配对点和连线仍是主要证据。每个后续节点单独绘制一张 `Fig3_<GSE>_Pre_vs_<节点>_paired_boxplot`，两侧严格使用相同患者。因此不同图的 Pre 箱体可能不同，这是匹配人群不同造成的。

箱体为中位数/IQR；须延伸至 1.5 IQR 内最远观察值。显示全部观察点，不另外重复标出离群点。横向 jitter 为 seed=20260910 的固定患者偏移，同一患者两端偏移相同；纵轴不扰动。IV 仅 2 人，应重点查看点和连线，不把箱体解释为稳定的总体分布。

R=红色圆点，NR=蓝色三角形；两个队列的临床终点仍分别是 RECIST 与 pCR，不能合并解释。GSE205506 Post 的 malignant 标签仍保留此前的注释限制。

## 复现与样式

只需运行 `script/plot.R`：它先从已有 patient CSV 重新执行配对 t 检验，再绘制总览、逐患者图和箱线图，无需读取原始单细胞矩阵。`data_processing.R` 完整重处理时也会生成新检验表。

STYLE 新接口：`show_t_tests`、`stat_font_size`、`show_paired_boxes`、`box_width_mm`、`box_height_mm`、`box_width`、`jitter_width`、`jitter_seed`。原有红蓝颜色、点形、字号和导出尺寸接口继续有效。

Python 从配对差值独立复算 t，并通过 Student-t 密度的数值积分复核 p，同时验证 BH；不需额外安装 scipy。PDF 文字按实际变换矩阵核验，结果保存在 source/QA.json。详见同目录 `配对t检验结果.md` 的具体 p/q 表。
