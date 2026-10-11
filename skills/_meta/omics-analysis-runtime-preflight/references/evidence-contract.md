# Runtime Preflight Evidence Contract

## Status meanings

| Status | Meaning |
|---|---|
| PASS | Direct evidence satisfies every required gate for the scoped claim. |
| WARN | The check ran, but evidence is partial, sampled, or non-blocking. |
| FAIL | A required gate failed, timed out, or contradicted the claim. |
| SKIP | The check was intentionally not requested; it provides no evidence. |

Never convert SKIP to PASS. Keep a timed-out probe as FAIL and preserve its timeout and captured stderr.

## GPU claim ladder

To say "this model can train on this server GPU now," require all gates:

1. **Physical**: PCI discovery identifies an NVIDIA compute-capable device.
2. **Driver**: `nvidia-smi` returns successfully and reports the device.
3. **Device exposure**: the expected NVIDIA character devices are present and accessible.
4. **Framework**: the exact environment reports `torch.cuda.is_available() == True`, at least one device, and completes a CUDA tensor operation.
5. **Named model**: the model trainer reports a CUDA root device, training completes, numerical outputs are finite, and peak CUDA memory allocation is nonzero.

Record a separately tested row for scVI and cell2location. Different environments, Torch builds, Lightning versions, or model APIs can produce different outcomes.

## Data-contract gates

For all AnnData inputs require:

- successful direct HDF5 reopen without loading all layers;
- nonzero shape and unique observation/variable names;
- an explicit count source;
- a bounded sample of finite, nonnegative, integer-like values for count-model inputs;
- required sample, batch, label, or condition keys.

For Xenium or Visium additionally require:

- `obsm["spatial"]` with at least two columns;
- finite sampled coordinates;
- a stable sample/library boundary that downstream spatial graphs can respect.

The checker streams dense tiles or sparse data/index chunks directly with h5py, never through AnnData backed layer loading. Matrix chunk memory is bounded by `--chunk-values`; key-uniqueness metadata remains O(keys). Default sampling is a partial integrity gate; `--full-count-scan` streams all stored count values. Neither proves provenance or biological validity. Auto source selection warns and training still requires a declared source.

## Smoke-test boundary

Synthetic smoke tests prove package compatibility, accelerator routing, and a minimal forward/backward optimization loop. They do not prove:

- full-dataset GPU memory fit;
- acceptable throughput;
- convergence under final epochs and batch size;
- posterior-export correctness unless explicitly tested;
- biological validity or batch-correction quality.

Before a full run, estimate memory on a representative pilot, record peak allocation, and keep a margin for posterior export, neighbors, UMAP, and checkpoint writes.

## Cell2location-specific boundary

The bundled cell2location smoke test trains both `RegressionModel` and spatial `Cell2location` using synthetic integer counts and known signatures. It tests both trainer-device paths but does not replace the algorithm contract for signature export, `N_cells_per_location`, `detection_alpha`, posterior abundance keys, or convergence review.

## Minimum report table

| Scope | Interpreter | Packages | Data contract | Framework CUDA | Actual trainer device | Result | Limitation/next action |
|---|---|---|---|---|---|---|---|

Link each row to its JSON evidence file. Do not summarize only from terminal text when a machine-readable report exists.
