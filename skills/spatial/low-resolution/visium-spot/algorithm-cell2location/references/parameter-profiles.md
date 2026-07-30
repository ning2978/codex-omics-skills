# Cell2location Parameter Profiles

These profiles preserve two previously used configurations. They are provenance records and starting points, not universal defaults.

## Selection Rule

1. Use `recorded-author-workflow` only when reproducing the documented author-style workflow.
2. Use `historical-project-adapted` only when reproducing or extending the project that established its denser-spot assumption.
3. For a new dataset, estimate cells per spot, check chemistry and tissue density, then create an explicitly named project profile.
4. Run a small pilot before a full model and record all overrides.

## Recorded Author Workflow

Reference signature model:

```python
REF_EPOCHS = 250
BATCH_SIZE = 2500
LR = 0.002
TRAIN_SIZE = 1.0
```

Spatial model:

```python
ST_EPOCHS = 30000
N_CELLS_PER_LOCATION = 5
DETECTION_ALPHA = 20
TRAIN_SIZE = 1.0
BATCH_SIZE = None
USE_GPU = False
```

Use for author-style reproduction or a deliberate sensitivity comparison. Expect long runtime and verify that five cells per location fits the target assay.

## Historical Project-Adapted Profile

Reference signature model:

```python
REF_EPOCHS = 600
BATCH_SIZE = 2048
LR = 0.002
TRAIN_SIZE = 1.0
```

Spatial model:

```python
ST_EPOCHS = 5000
N_CELLS_PER_LOCATION = 30
DETECTION_ALPHA = 20
TRAIN_SIZE = 1.0
BATCH_SIZE = None
USE_GPU = True
```

This profile trades strict reproduction for runtime and preserves a project-specific dense-spot assumption. Do not select it merely because the task is ambiguous.

## Smoke Tests

Short runs such as `ST_EPOCHS=800` or altered `DETECTION_ALPHA` may be useful for plumbing checks. Label them as smoke tests and do not present them as final high-confidence profiles.

## Required Metadata

Record:

- profile name and source;
- package versions and accelerator;
- reference and spatial epochs;
- learning rate and batch size;
- `N_cells_per_location` and its justification;
- `detection_alpha`;
- train size, seed, and every override.
