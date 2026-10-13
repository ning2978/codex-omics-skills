# Environment acceptance — 2026-10-13

| Platform | Lock solve | Native installation and smoke tests |
|---|---|---|
| Linux x86_64 | PASS | PASS, fresh isolated Miniforge and Python/R environments |
| Apple Silicon macOS 14+ | PASS, cross-solved for target 14 | PENDING, not executed on macOS |

Linux passed counts → HVG/PCA/neighbors/Leiden/UMAP → H5AD round trip,
Harmony, Seurat → clusters/UMAP → RDS round trip, the existing Python–R bridge,
edgeR paired QL fit, both native DotPlot renderers and 17 PNG decodes. All 15
native DotPlot outputs were visually reviewed, with representative full-size
checks; the two basic smoke PNGs were also inspected. No bottom captions,
overlap or cropped labels/legends were observed. Explicit dimension override
fixtures intentionally keep their requested larger dimensions.

Tested core versions: Python 3.11.17, Scanpy 1.11.5, AnnData 0.12.7,
Harmony 0.0.10, Leidenalg 0.10.2, UMAP 0.5.12; R 4.4.3, Seurat 5.3.1,
SeuratObject 5.2.0, Matrix 1.7.6, edgeR 4.4.0. The conda scikit-misc 0.5.3
build reports 0.0.0 through Python metadata; reports retain both imported
metadata and the exact conda inventory. An actual LOESS fit passed.

[Machine-readable acceptance](../environments/validation.json) records config
and lock hashes, versions and platform status. Raw reports, synthetic artifacts,
failed attempts and installation logs remain outside the release tree.

The first combined cache attempt exposed a download collision; each default
creation now uses an independent cache. Explicit cache reuse is protected by
an advisory lock. Bridge acceptance exposed missing R.utils; both R profiles
now require it and Linux acceptance used a new R prefix, preserving the old one.

No macOS native run, remote CI, GPU model, specialized-method or research-data
analysis was executed. A successful platform solve alone is not runtime PASS.
