# public_data_analysis

可扩展、模块化、配置驱动、可复现的公共生物医学数据分析框架。
当前有可运行的 TCGA adapter、通用 correlation/survival 模块，以及完整版本化的 CRC 单细胞研究案例（含旧版图、CellChat 和新版 Fig1–Fig3）。后者位于 `examples/crc/`；GEO、CellResDB 和 scCT-DB 是数据来源，CRC 和 GSE 编号是研究选择，不作为通用 adapter 名称。
以下第 1–11 节记录原有 TCGA 工作流的接口和示例；CRC 案例使用独立的版本化入口。GEO、CellResDB、scCT-DB 的来源级通用 adapter 尚待实现，不能仅替换 GSE 编号就复用 CRC 案例脚本。

先读 [文档导航](docs/README.md) 和 [CRC 案例运行说明](examples/crc/README.md)。历史运行可用 `python workflow/replay.py --manifest <run_manifest.json>` 重放。

- 语言：R（base R 优先，兼容 HPC `module load R/4.2.0-container`）
- 依赖：`yaml` 必需；`data.table` 强烈建议；`ggplot2` 仅画图需要
- 工作流：Snakemake-ready（PoC 已附），但不强依赖 Snakemake

---

## 先选择运行入口

| 需求 | 入口 | 说明 |
|---|---|---|
| TCGA 相关性或生存分析 | `datasets/tcga/` → `analyses/` → `plotting/` | 使用标准化 bulk 表达/临床数据契约 |
| 完整复现本次 CRC 研究全部图表 | [`examples/crc/run.py`](examples/crc/README.md) | 版本化的研究案例；原始矩阵在仓库外 |
| 一键重放某次历史运行 | `python workflow/replay.py --manifest <run_manifest.json>` | 按该次 Git 版本及参数重跑，先检查输入身份 |
| 增加新研究 | [`docs/how_to/run_new_study.md`](docs/how_to/run_new_study.md) | 选择来源和 accession，复用已有模块 |
| 增加来源、分析或图型 | [`docs/how_to/`](docs/README.md) | 按数据契约实现相应模块 |

日常路径通过命令行传入；相对稳定的科学规则留在配置文件。每次运行产生
`run_manifest.json`，记录实际使用的路径、规则、代码版本、步骤和日志。
`examples/crc/` 的脚本固定解释这次研究所用的已发表注释；未来的 GEO 数据
应由按来源/文件格式开发的通用 adapter 处理，不能只替换一个 GSE 编号。

---

## 1. 架构与数据生命周期

```text
Download ──► Raw data ──► Prepare/Clean ──► Standardized processed ──► Analysis ──► Tables ──► Plot ──► Figures
 (per dataset adapter)      (per adapter)      (统一 contract)          (通用模块)                (通用模块)
```

四层严格解耦：

| 层 | 位置 | 输入 | 输出 |
|---|---|---|---|
| Download | `datasets/<src>/download.R` | dataset YAML | `raw/...` + `manifest.yaml` |
| Prepare | `datasets/<src>/prepare.R` | raw + dataset YAML | `processed/...` 统一 contract + `manifest.yaml`(fingerprint) |
| Analysis | `analyses/<type>/run.R` | processed contract + analysis YAML | `results/<name>/tables/*.tsv` + metadata |
| Plot | `plotting/<type>/plot.R` | result tables + plotting YAML | `results/<name>/figures/*` |

保证：改绘图参数不重算统计；改分析参数不重下载；raw 存在即复用；
processed 参数指纹不变即复用；任何步骤可单独执行。

## 2. 目录结构

```text
config/
  dataset_registry.yaml        # 数据源注册表（唯一事实来源）
  datasets/tcga.yaml           # TCGA：数据是什么、怎么获取/准备
  analyses/correlation_example.yaml
  plotting/default.yaml
  gene_sets/antigen_presentation.txt
datasets/tcga/                 # TCGA adapter：download.R / prepare.R / validate.R
processing/                    # 模块间数据契约及通用处理规则
analyses/correlation/          # 相关分析模块（dataset-agnostic）
analyses/survival/             # 生存分析模块（log-rank + Cox per SD，dataset-agnostic）
plotting/correlation/          # 热图模块（只读结果表，不重算）
plotting/survival/             # KM 曲线 + 森林图模块（只读结果表）
examples/crc/                  # 版本化 CRC 研究案例及全部图表生成代码
docs/                          # 长期通用的架构、数据契约、扩展指南
utils/                         # cli / logging / io / config / validation / gene_id
workflow/Snakefile             # PoC：仅依赖编排
tests/                         # 冒烟测试（暂缓，见 §9）
```

## 3. 数据目录与代码分离（强制）

代码仓库不含大矩阵。TCGA 现有模块支持下列环境变量；新入口优先通过命令行接收路径，并把生效值写入 run manifest：

```bash
export PDA_DATA_ROOT=/data/public_data     # raw/ processed/ cache/ 的根
export PDA_RESULTS_ROOT=/project/results   # 分析结果根
```

布局：

```text
$PDA_DATA_ROOT/raw/TCGA/…                  # 原始文件，只增不改
$PDA_DATA_ROOT/processed/TCGA/<cohort>/    # 标准化产物（contract 见下）
$PDA_DATA_ROOT/cache/
$PDA_RESULTS_ROOT/<analysis_name>/{tables,figures,logs,metadata,sessionInfo.txt}
```

配置里的路径支持 `${VAR}` 展开，也可用 CLI flag 覆盖。

## 4. Processed-data contract（所有 adapter 必须满足）

| 文件 | 内容 |
|---|---|
| `expression.rds` | numeric matrix，**feature × sample**，值域/单位写入 manifest（TCGA 为 log2(TPM+0.001) 透传） |
| `sample_metadata.rds/.tsv` | `sample_id` 主键 + patient/cohort/sample_type/sex 等 |
| `feature_metadata.rds/.tsv` | symbol、来源 Ensembl、折叠策略审计列 |
| `clinical.rds/.tsv` | 项目级标准字段（OS/PFS/DSS/DFI_time+event、age、sex、stage、response、treatment…）；缺失即 NA，绝不伪造 |
| `manifest.yaml` | 来源、参数 fingerprint（sha256 of 预处理参数+输入校验和）、计数、字段可得性 |

ID 政策显式化：Ensembl 版本剥离、无映射 ID 回退保留并标记、重复 symbol 按
`max_iqr|max_variance|first_lexical` 折叠并留 audit。

## 5. 快速开始（最小运行案例）

```bash
# HPC
module load R/4.2.0-container
export PDA_DATA_ROOT=/data/public_data
export PDA_RESULTS_ROOT=/project/results
Rscript -e 'install.packages(c("yaml","data.table","ggplot2"))'

cd public_data_analysis   # 项目根（脚本对 cwd 不敏感，但从根运行最直观）

# 1) download（幂等；已有且校验和匹配则跳过）
Rscript datasets/tcga/download.R --config config/datasets/tcga.yaml

# 2) prepare（fingerprint 未变则复用）
Rscript datasets/tcga/prepare.R --config config/datasets/tcga.yaml --cohorts SKCM

# 3) validate
Rscript datasets/tcga/validate.R --config config/datasets/tcga.yaml --cohort SKCM

# 4) analyze
Rscript analyses/correlation/run.R --config config/analyses/correlation_example.yaml

# 5) plot（只读第4步的结果表）
Rscript plotting/correlation/plot.R \
  --input $PDA_RESULTS_ROOT/mettl3_mhci_tcga_skcm/tables/correlation_results.tsv
```

或用 PoC 工作流：

```bash
snakemake -n --cores 1 && snakemake --cores 1
```

## 6. Config 系统

三层解析：**内置 default < YAML 配置 < CLI flag**；路径字符串支持 `${ENV}`。

* **Dataset config**（`config/datasets/*.yaml`）：描述"数据是什么/如何获取与准备"。
  TCGA-specific 字段（cohorts、sample_types、Xena URL）只出现在这里。
* **Analysis config**（`config/analyses/*.yaml`）：描述"对已备好的数据做什么"。
  **gene list 只在这里**（§8：下载永远拿全矩阵）。
  支持 inline `genes:` 与 `geneset_file:`，二者可组合。
* **Plotting config**（`config/plotting/default.yaml`）：纯外观参数。
* required/optional/default 在各模块 README 中列表化，缺 required 直接 fail。

## 7. 如何新增一个 dataset（如 CPTAC）

1. 新建 `datasets/cptac/{download.R, prepare.R, validate.R, README.md}`
   （复制 tcga 的骨架，替换获取逻辑与解析逻辑）；
2. 在 `config/dataset_registry.yaml` 把对应条目 `status: planned → active`；
3. 满足 §4 contract。
完成。correlation/survival/plotting 等**零改动**即可消费新数据源。

## 8. 如何新增一个 analysis 模块（如 survival）

1. 新建 `analyses/survival/run.R`：只 import `utils/`，只读 processed contract；
2. 输出固定到 `results/<analysis.name>/tables/<module>_results.tsv` + `metadata/`；
3. 需要图则新建 `plotting/survival/plot.R`（只读表，不重算）。
CLI 规范：非交互、cwd 无关、`--help`、exit 0/1/2、fail loudly、原子写、写 provenance。

## 9. 测试状态（当前阶段）

TCGA 模块尚未在本次重构中重新运行。CRC 案例已在现有数据上重绘和校验，执行范围详见 `examples/crc/README.md`。建议下一步在 HPC 上对 TCGA 做冒烟：
小 cohort（如 SKCM）走完 download→prepare→validate→analyze→plot，
再验证：缺基因仅告警；改 plot config 只有 plot 重跑；改 analysis config 重跑 analyze+plot；
Snakemake `-n` DAG 符合预期。

## 10. Snakemake / HPC 兼容性

* 每个脚本 = 独立 CLI 单元，输入输出路径确定 ⇒ 可直接包成 rule；
* Snakefile 只做依赖/参数/资源编排，逻辑永远在 R 脚本内；
* R 脚本内部绝无 `#SBATCH`、绝无 `module load`——资源请求留给 cluster profile；
* exit code：0 成功 / 1 运行失败 / 2 用法错误。

## 11. Provenance

每次分析产出：`metadata/config_used.yaml`（生效配置回显）、
`metadata/analysis_metadata.yaml`（输入、软件版本、时间戳、git commit 若可用）、
`sessionInfo.txt`、`logs/run.log`（含 warnings 计数与 STATUS 行）。
