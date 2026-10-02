# CRC 公共单细胞分析交付

目标基因：**METTL3、STRAP、PTBP1**。

- [结果报告](结果报告.md)：完整性核查、补下载、主要结果和解释。
- [设计与复现说明](设计与复现说明.md)：分析优先级、统计单位、阈值和局限。
- [运行说明](运行说明.md)：实际运行时、脚本顺序和复现命令。
- `analysis_config.json`：统一配置。
- `run_analysis.ps1`：一键复现入口；默认重跑分析/图表，`-Mode full` 重建完整流程。
- `figures/`：可编辑 PDF 和 300-dpi PNG。
- `tables/`：全部检验、患者/样本源数据与敏感性分析。
- `audit/`：官方文件核对、SHA256、矩阵验证、版本和执行日志。
- `data/GEO/`：补齐的完整公开数据；原有 OneDrive 文件未覆盖。

**主要结论：尚无经过多重校正、跨独立研究稳健支持的恶性细胞 MHC-I/APM 关联，完整 GSE236581 中三个基因的基线 R/NR 差异也不显著。** 小样本 PTBP1–Treg 程序关联属于探索发现，尚未独立验证。CellChat 扩展现已执行；结果单独见 [CellChat 补充报告](CellChat补充报告.md)，方法和运行入口见 [CellChat 运行与设计](CellChat运行与设计.md)。
