"""Refresh prior deliverables only after the completed CellChat report exists."""
from pathlib import Path
import pandas as pd
ROOT=Path(__file__).resolve().parents[1]
assert (ROOT/'CellChat补充报告.md').exists()
assert len(pd.read_csv(ROOT/'audit/cellchat_completed_samples.csv'))==27
def edit(name,old,new):
 p=ROOT/name;s=p.read_text(encoding='utf-8')
 if old in s:p.write_text(s.replace(old,new),encoding='utf-8')
edit('README.md','CellChat 未执行，已明确列为后续扩展。','CellChat 扩展现已执行；结果单独见 [CellChat 补充报告](CellChat补充报告.md)，方法和运行入口见 [CellChat 运行与设计](CellChat运行与设计.md)。')
edit('结果报告.md','4. **CellChat 放在扩展层**：本次没有运行 CellChat。已准备患者 high/low 分组及明确的后续设计，但现环境缺少该依赖，而且它不能直接测量肽–MHC/TCR 特异性。不要拿 79 个功能基因的表达乘积代替 CellChat 结果。',
     '4. **CellChat 扩展已补执行**：27 个治疗前肿瘤样本已从完整基因矩阵分别完成真实 CellChat 推断，并进行了患者级 high/low 比较。见 [CellChat 补充报告](CellChat补充报告.md) 和 [配置与复现](CellChat运行与设计.md)。该扩展仍不能直接测量肽–MHC/TCR 特异性，前述主分析结果不变。')
edit('QA与文件索引.md','CellChat 为后续设计，本次没有执行。','CellChat 扩展现已执行，另见 CellChat补充报告.md 和 CellChat运行与设计.md。')
edit('QA与文件索引.md','**不含大型 data/ 和 derived/**','**不含大型 data/、derived/ 和 runtime/**')
edit('QA与文件索引.md','本次没有安装 svglite，也未输出 SVG/TIFF','首轮没有安装 svglite；CellChat 扩展已将它作为依赖装入项目库，仍未输出 SVG/TIFF')
edit('运行说明.md','没有修改系统 R/Python，也没有安装新包。','首轮没有修改系统 R/Python，也没有安装新包；后续 CellChat 扩展在项目 runtime/R-library 中安装了缺失依赖，详见 CellChat运行与设计.md。')
p=ROOT/'运行说明.md';s=p.read_text(encoding='utf-8')
old=s.split('## CellChat 与原始 QC 的扩展边界')[1] if '## CellChat 与原始 QC 的扩展边界' in s else None
if old is not None:
 s=s.split('## CellChat 与原始 QC 的扩展边界')[0]+'''## CellChat 扩展和原始 QC 边界

此前未执行的 CellChat 扩展现已完成。见 `CellChat补充报告.md`、`CellChat运行与设计.md`，配置为 `cellchat_config.json`，入口为 `run_cellchat.ps1`。快速重画使用 `./run_cellchat.ps1 -Mode compare`；完整扩展流程使用 `./run_cellchat.ps1 -Mode all -Workers 4`。原先 `run_analysis.ps1` 保持原主分析范围。

本扩展重新提取完整基因稀疏矩阵，按患者建模；没有用 79 个功能基因冒充全表达矩阵。也没有重跑从未过滤 droplets 开始的 doublet 或 ambient RNA 处理，这类重新 QC 需要相应原始输入。本次继续使用作者 QC 和细胞注释。
'''
 p.write_text(s,encoding='utf-8')
print('Prior reports linked to completed CellChat extension')
