# spatial/high-resolution/xenium

10x Xenium cell-level 空间原位分析。

| Skill | 用途 |
|-------|------|
| `xenium-segmentation-new` | XOA 2.0+ cells/cell_boundaries/nucleus_boundaries parquet 加载、QC 和 polygon 渲染。 |
| `xenium-roi-h5ad-qc-segmentation` | ROI h5ad QC、anucleate/cell_area 过滤和 segmentation method 保留。 |
| `xenium-roi-dapi-overlays` | ROI DAPI crop 上绘制 label centroid overlay。 |
| `xenium-scrna-label-transfer` | scRNA 参考标签投射到 Xenium。 |
| `xenium-scimap-neighborhood-analysis` | scimap RCN、目标细胞 NMF 邻域状态和尺度敏感性分析。 |
