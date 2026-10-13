# Codex Omics Skills

44 个可复用 Codex skills，覆盖单细胞 RNA-seq、Visium、Visium HD、Xenium、
CODEX/IMC、DotPlot，以及环境、运行时和发表级质量控制。发布副本按本地最新
方法与证据契约更新；项目专用路径、队列示例和运行结果不随包发布。

Reusable skills for single-cell and spatial omics, with portable scripts,
shared resources, explicit evidence contracts and synthetic regression tests.

## 发布范围

| 分组 | 数量 | 内容 |
|---|---:|---|
| `_meta` | 4 | taxonomy、环境部署、运行时预检、conda/pip 环境策略 |
| `single-cell` | 12 | QC、整合聚类、marker 注释、样本比较、通信、子聚类、轨迹、regulon、VDJ、Monocle2、R 桥接、workflow |
| `spatial` | 25 | shared、Visium spot/HD、Xenium、CODEX/IMC，含空间零模型和重分割输入审计 |
| 发表级分析 | 1 | 证据、统计与审稿检查 |
| DotPlot | 2 | scRNA 与 Xenium，共享 native Seurat 数值和渲染引擎 |
| **总计** | **44** | [`catalog.json`](catalog.json) |

同时包含 `_shared/dotplot`、`_shared/spatial`、各 skill 的脚本、参考文档和
合成数据测试。不包含患者数据、项目 annotation registry、结果、日志、缓存、
服务器绝对路径或项目专用绘图 skill。历史 Git 提交保留。

路由入口为 [skills-taxonomy](skills/_meta/skills-taxonomy/SKILL.md)。
选择最小相关 skill，并在实际需要时加载其资源；安装依赖不意味着自动运行分析。

## 安装

```bash
git clone https://github.com/ning2978/codex-omics-skills.git
cd codex-omics-skills
python scripts/install_skills.py --dry-run
python scripts/install_skills.py
```

默认目标为 `${CODEX_HOME}/skills`；未设置时为 `~/.codex/skills`。
用 `--target /path/to/skills` 指定其他位置。保留原目录层级，不要把叶节点扁平化。

```bash
python scripts/install_skills.py --name render-scrna-dotplot
python scripts/install_skills.py --name sc-subcluster-analysis
python scripts/install_skills.py --category spatial
python scripts/install_skills.py --list
python scripts/install_skills.py --force
```

按名称、分类或全量安装均递归包含 catalog 声明的 skill 依赖和共享资源。
[`dependencies.json`](dependencies.json) 声明资源所有者、代码依赖和必要路由叶节点；
条件性方法建议不一律视为依赖。taxonomy 安装包含完整路由集合。`--name` 可重复，
与 `--category` 同时使用时先取交集，再展开依赖。

目标目录已存在时默认跳过并返回 1，包括共享资源目录；不会悄悄覆盖本机修改。
更新时显式 `--force` 覆盖同名文件，但不删除无关文件。安装到已有旧版目录时，
历史额外文件也保留；需要精确发布文件集时安装到新目录。`--dry-run` 不写入文件；
名称或分类无匹配时返回 2。

## 新手入口：先部署分析环境

首次使用可让 Codex「使用 `$omics-environment-setup`，检查我的电脑并部署
Python 单细胞和 R/Seurat 基础环境」。支持 Linux x86_64 和 Apple Silicon
macOS；Windows 请在 WSL2 内操作。macOS 原生验收目前待运行，不能视为 PASS。

```bash
python3 skills/_meta/omics-environment-setup/scripts/setup.py plan \
  --profile both --report-dir "$HOME/omics-reports/plan-01"
python3 skills/_meta/omics-environment-setup/scripts/setup.py create \
  --profile both --bootstrap --install-dir "$HOME/omics-environments" \
  --report-dir "$HOME/omics-reports/create-01"
```

需要 Python 3.10+ 启动部署脚本；无需先安装 conda。`inspect`/`plan` 只检查并写
报告，`create` 可引导安装固定版本 Miniforge，并从平台锁定文件创建独立环境。
不会修改 shell 配置或升级原环境。使用报告中的准确解释器路径运行后续任务。
详见 [部署 skill](skills/_meta/omics-environment-setup/SKILL.md) 和
[本次验收](ENVIRONMENT_VALIDATION.md)。

**环境可用不等于已有完整分析 runner。** 基础环境覆盖 Scanpy、Seurat、桥接、
DotPlot 和 edgeR；没有新增完整 Seurat 分析 skill，也不保证全部 44 skills
的方法运行条件。scVI、cell2location、Monocle2 等仍需对应方法环境与 runtime-preflight。

## 运行依赖与可移植性

安装器只安装 skills 文件，不安装分析软件。使用实际任务已验证的 conda/pip
环境。Python 脚本按功能需要 NumPy、pandas、SciPy、AnnData/Scanpy、h5py、
PyYAML、Pillow、pyarrow、matplotlib、tifffile 或 Scirpy 等；R 工具需要各自
的 Seurat、edgeR、spacexr 等包。先使用 runtime-preflight 核实精确解释器和版本。

DotPlot 的 Rscript 路径及期望 R/Seurat/SeuratObject 版本由 run YAML 指定；
共享 R 渲染器通过安装树相对路径解析，不依赖符号链接。DotPlot 仍要求注释
registry/guide 一致、精确键连接、来源值和适用门槛校验。队列限制使用 typed
`cohort_gate`，旧的文字范围须提供明确配置映射；不再推断疾病代码。

只提供方法契约而未附实现的项目 QC、DAPI、Visium HD 或结果绘图入口，已明确
说明应使用经过验证的项目 runner。可选 `xenium_downstream` recluster guard
属于消费项目，不是本仓库依赖；其启用和豁免须遵守项目政策。

图形默认 PNG，不在底部添加解释、方法或 provenance。保留轴、图例和程序条。
单组 DotPlot 保留 native log1p mean、缺失 z-score 与 `Expression (unscaled)`
图例，并把未缩放说明写入结果记录和文档。

## 校验

无需分析依赖的发布检查：

```bash
python -B scripts/build_catalog.py --check
python -B scripts/validate_repo.py
python -B -m unittest discover -s scripts -p 'test_*.py' -v
```

覆盖 44 个 skill 的元数据、taxonomy、资源链接、Python 语法、确定性 catalog、
依赖闭包、可移植性和数据/缓存排除；安装测试在临时目录验证全量、单个 DotPlot、
分类、重复、显式覆盖和 dry-run。GitHub Actions 执行这些检查。

在具备对应分析依赖的环境中运行：

```bash
export PYTHONDONTWRITEBYTECODE=1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
python -m unittest discover -s skills/_shared/dotplot/tests -v
python -m unittest discover -s skills/_shared/spatial/tests -v
python -m unittest discover -s skills/_meta/skills-taxonomy/scripts -v
python -m unittest discover -s skills/single-cell/workflow/workflow-sc-analysis/scripts -v
python -m unittest discover -s skills/single-cell/python/sc-multiomics-tcr-vdj-analysis/scripts -v
SC_BRIDGE_RSCRIPT=/path/to/Rscript python -m unittest discover -s skills/single-cell/shared/h5ad-to-r-bridge/scripts -v
/path/to/Rscript --vanilla skills/single-cell/python/sc-sample-level-comparison/scripts/edger_runtime_smoke.R
python skills/_shared/dotplot/tests/verify_style.py --rscript /path/to/Rscript --work-dir /path/to/new-test-output
```

最后一项合成测试当前使用 R 4.4.3、Seurat 5.5.1、SeuratObject 5.4.0；生产 run
使用自己的显式版本契约。可选 `--baseline /path/to/previous/render_scrna_dotplot.py`
增加新旧统计回归对照。不指定 baseline 时仍执行 Python–Seurat 数值一致性检查，
但不会声称完成历史版本对照。测试输出写在发布树之外。

本次发布的检查结果与未执行范围见 [VALIDATION.md](VALIDATION.md)。
本仓库目前未指定许可证；公开可见不等同于授予特定开源许可。
