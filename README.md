# Codex Omics Skills

一套从真实单细胞 RNA-seq 与空间组学分析流程中整理出的可复用 Codex
skills。首个版本包含 34 个 canonical skills，覆盖 scRNA-seq、Visium、
Visium HD、Xenium、CODEX/IMC，以及跨平台的环境、运行时和发表级质量控制。

This repository packages reusable Codex skills for single-cell RNA-seq and
spatial omics. The skills are organized by modality, spatial resolution,
workflow stage, platform, and algorithm.

## 发布范围

| 分组 | 数量 | 内容 |
|---|---:|---|
| `_meta` | 3 | taxonomy、运行时预检、conda/pip 环境策略 |
| `single-cell` | 7 | workflow、子聚类、轨迹、regulon、VDJ、Monocle2、h5ad/R 桥接 |
| `spatial` | 23 | shared、Visium spot、Visium HD、Xenium、CODEX/IMC |
| cross-cutting | 1 | 发表级 omics 分析与审稿可辩护性 |
| **总计** | **34** | 详见 [`catalog.json`](catalog.json) |

本仓库不包含患者数据、临床表格、H5AD/CSV/Parquet 结果、项目 annotation
registry、冻结的项目默认值、服务器绝对路径、环境锁文件或运行日志。项目特有
的路径和标签应保留在对应项目的 repository-local adapter 中。

## 目录

```text
skills/
├── _meta/
├── publication-grade-omics-analysis/
├── single-cell/
│   ├── workflow/
│   ├── python/
│   ├── r/
│   └── shared/
└── spatial/
    ├── shared/
    ├── low-resolution/visium-spot/
    └── high-resolution/
        ├── visium-hd/
        ├── xenium/
        └── codex-imc/
```

路由入口是
[`skills-taxonomy`](skills/_meta/skills-taxonomy/SKILL.md)。不确定分析阶段
时先用 workflow skill；明确指定方法时再组合对应 algorithm 或平台 skill。

## 安装

克隆仓库后先预览：

```bash
python scripts/install_skills.py --dry-run
```

安装全部 skills 到 `${CODEX_HOME}/skills`；未设置 `CODEX_HOME` 时使用
`~/.codex/skills`：

```bash
python scripts/install_skills.py
```

如果目标 skill 已存在，安装器默认跳过并返回非零状态。确认要用仓库版本更新
现有文件时显式使用：

```bash
python scripts/install_skills.py --force
```

只安装一个或一组 skill：

```bash
python scripts/install_skills.py --name sc-subcluster-analysis
python scripts/install_skills.py --category spatial
python scripts/install_skills.py --list
```

安装器只复制 catalog 中声明的 skill 目录，不会删除目标目录中的其他内容。

## 校验

```bash
python scripts/build_catalog.py --check
python scripts/validate_repo.py
```

校验覆盖：

- 34 个 skill 的 frontmatter、名称、长度和唯一性；
- `agents/openai.yaml` 与 `$skill-name` default prompt；
- taxonomy 完整性、退休名称和资源引用；
- Markdown 本地链接和 Python 语法；
- catalog 可重复生成；
- 本机绝对路径、项目标识、凭证特征、数据文件、缓存和大文件。

## 维护原则

- 通用 workflow 与算法逻辑进入本仓库。
- 组织、疾病、队列、样本和 repository-specific 默认值留在项目 adapter。
- 原始 counts、注释证据、投射标签和空间邻域标签保持为不同语义轨道。
- 新 skill 必须更新 taxonomy、catalog，并通过全仓校验。

当前整理副本暂未附加开源许可证。若 GitHub 仓库设为 public，应先明确选择
许可证；private 仓库可以先用于版本管理和后续审阅。
