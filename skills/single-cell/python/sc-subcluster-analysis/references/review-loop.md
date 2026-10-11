# Reviewed round protocol and portable tools

Tools use the existing Python environment with AnnData, Scanpy, NumPy, SciPy and pandas. They import the relative shared `workflow-sc-analysis/scripts/sc_contracts.py`; retain that skills-tree structure when copying. They do not import a project package or implement a new clustering engine. Use a validated project entrypoint or the integration skill between rounds.

## Evidence review

Run `sc-marker-annotation/scripts/collect_evidence.py --input full.h5ad --config evidence.json --out new-evidence-directory` from the selected leaf's location. Read [the evidence interface](../../sc-marker-annotation/references/evidence-interface.md) for its schema. Feed its tables plus complete markers into agent biological review (and human review when requested); the script neither labels nor deletes cells.

## Round manager

Use this skill's `scripts/round_manager.py`. All paths in JSON are explicit paths; relative evidence/approval paths resolve against the invocation directory. Source files and old rounds are immutable. Do not edit decision/approval JSON to change scope. For a revised proposal before application, initialize a fresh run referencing the same source; for subsequent filtering rounds, name the applied predecessor. Independent sessions must not write one run concurrently.

Context JSON:

```json
{"platform":"xenium_segmented", "species":"declared species", "tissue":"declared tissue",
 "target_lineage":"declared target", "counts_origin":{"kind":"raw_counts","source":"identified artifact"},
 "integration_method":"harmony", "batch_key":"explicit_technical_batch",
 "cleanup_mode":"auto_supported", "cleanup_request_record":"initial-cleanup-request.txt",
 "sample_key":"sample_id", "unit_key":"donor_id", "gene_scope":"all measured panel genes"}
```

New runs use subcluster_loop_v2 and freeze subcluster_round_v2. Context must declare integration_method="harmony" and batch_key; both are bound into proposal/approval hashes. v1 runs support read-only summary/validation; continuing requires a fresh v2 run and new decision. Missing cleanup_mode in existing runs defaults to manual_per_round; init also defaults to manual when no explicit mode is supplied. For an explicit new cleanup-loop request set auto_supported and save its exact initial request/reference; do not infer authorization from a plotting or pipeline-development request. When requested set manual_per_round. The request file is hash-bound to the run.

Allowed cell models: `scrna`, `snrna`, `xenium_segmented`, `segmented_cells`. A missing biological-unit field is `unit_key: null` with the limitation recorded. `synthetic: true` is reserved for generated test data and permits explicitly simulated confirmations; never use it on a real analysis.

1. `init --run new-run --source full.h5ad --context context.json`
2. Review or cluster. `freeze --run new-run --round r01 --source reviewed.h5ad --cluster-key leiden --proposals actions.json` records existing-cluster review by default. Add `--recomputed --execution-manifest engine/run_params.json` only after a real mandatory-Harmony expression run passes its artifact checks (the default manifest location is alongside the source H5AD). Freeze verifies actual X_harmony, its fingerprints, main graph use_rep, stages and hash-verified execution output.
3. Present/archive `proposals.json`, `cells.csv`, `proposed_impact.csv` and `proposal_hash` in `frozen.json`, including target loss, alternatives, sample and cumulative attrition. New freezes record `evidence_rule_version=foreign_program_v1`.
4. In auto mode create a separate biological review JSON after freezing, then `auto-review --run new-run --round r01 --expected-hash HASH --record automatic-review.json`. Record every exclusion proposal as supported or review; all supported proposals are applied together. In manual mode save the actual user message/reference and use `approve --run new-run --round r01 --expected-hash HASH --select p1 p3 --approver USER --record approval-message.txt`. Wait for actual approval only in this requested manual mode; an identical already approved list needs no repeat question. Automatic decisions never call approve or its simulated interface.
5. `apply --run new-run --round r01` verifies source, frozen files, evidence and the mode-specific decision, then writes a new `retained_counts.h5ad`, exact excluded list, sample/cumulative impact and full initial-cell ledger. Missing automatic review or manual approval prevents filtering. Unselected/review cells remain. It drops stale representations/graphs and target cluster, retaining old labels as history. It does not claim the derivative is reclustered. If retained cells lack two valid batches, preserve counts/ledger and stop as technical_failure; no fallback. If the engine later fails, record technical_failure against the applied terminal round and retain counts.
6. Recluster the retained counts using the validated engine, review all genes again, then `freeze ... --round r02 --source next-reviewed.h5ad --cluster-key next_key --previous r01 --recomputed`. New rounds must conserve the previous retained counts/IDs, including shuffled input order. Old cluster numbers are not identities across rounds.
7. `validate --run new-run --round r02 [--artifacts artifacts.json]` distinguishes frozen-integrity, approved-filter and full-round artifact checks. `summary --run new-run` reports progression. `stop --run new-run --round r02 --reason evidence_insufficient --note '...'` records why the loop ended while retaining unresolved cells. `awaiting_approval` creates a pause snapshot and does not lock future approval; other stop reasons close the run.

Proposals JSON is a list with exactly one item per observed cluster. Cluster values are strings. Minimum retained/review item: `{"id":"p1","cluster":"0","action":"review","reason":"Insufficient competing-lineage evidence"}`. Actions: `retain`, `exclude_foreign`, `exclude_mixed`, `review`.

Every exclusion also requires `foreign_program` (schema below), `target_program`, `foreign_lineage`, `foreign_markers` (at least two distinct coordinated markers from that lineage), `unit_support`, `qc_review`, `interpretation`, and `evidence_files` (paths). `exclude_mixed` additionally requires `target_cell_loss`, describing target-positive counts/fractions and possible lost target biology. These fields retain expert reasoning; presence validation is not scientific endorsement. Record marker availability, recurrence/sparsity and alternatives explicitly, never invent support to satisfy a field.

Hashes bind the proposal to its run/context, round, source artifact, exact clusters/cells, actions, impact and evidence. Changing any of those requires a fresh review/decision. Decisions cannot migrate across rounds. Do not convert `review` to an exclusion by filling missing labels.

## Artifact validation interface

`artifacts.json` names `annotated_h5ad`, `hvg_h5ad`, `all_markers`, `top20`, `evidence` (portable collector manifest; retain expert review tables separately), `sidecar`, `palette`, `pngs` (list), `displayed_labels` (union), `categorical_fields` (plotted field to label-prefix mapping, e.g. `{"leiden":"r02:leiden:","display_name":""}`) and `execution_manifest`. Paths are explicit. Full annotation requires separate `cell_identity`, `cell_subtype`, `cell_state`, `display_name`, `annotation_status`; unknown values remain explicit. The exact-cell sidecar uses `obs_name` and the same fields. Earlier source-sidecar joins still follow the project registry's declared key.

The validator reopens counts/HVG/sidecar, aligns by keys, validates current expression artifacts when recomputed, checks full marker gene/cluster coverage against the evidence manifest, unique colors and PNG decoding, hashes outputs and requires a successful execution manifest with actual input/stages and a hash-verified output artifact for recomputed rounds. It does not infer scientific identity or prove provenance flags alone. Review stage execution/input hashes under the selected engine; manually inspect PNGs and verify points/legends share the saved mapping. Technical failures remain failures, not convergence.

Use one explicit palette for the union of semantic labels across a figure collection; extend/reuse it in subsequent rounds. Prefix cluster-only labels with round/cluster-key namespace. Never cycle a short palette. PNG only; a verified same-basename PNG permits removal of its paired PDF, never of standalone reports. Keep output manifests current.

Two simulated-confirmation switches are needed for tests: `approve --simulated` and `apply --allow-simulated`, with `context.synthetic=true`. This prevents accidentally treating a fixture's confirmation as real authorization.

## Program and automatic decision schemas

Two same-lineage foreign markers are screening evidence only. Removal requires a sourced, lineage-specific foreign program with multiple complementary markers, measured coverage and audited feature mapping; there is no universal gene-count or score cutoff. Review target programs, cell-level distributions and coordinated coexpression, biological-unit support/limitations and platform QC. Neither averages, module scores, UMAP nor seed instability alone justify deletion. Do not require every cell to express every program gene or treat panel/detection gaps as negative. Dual program positivity alone is not contamination. If genuine target subtype, biological mixed state or potential target-cell loss remains unresolved, retain as `review`. Whole-cluster removal needs supported technical mixture or an explicit outside-target-scope interpretation and quantified target-positive count/fraction and potential loss; partial-cell evidence cannot be expanded to a whole cluster. Retain or use the applicable at-most-one mixed-parent split and review a supported foreign child.

`foreign_program` requires nonempty `id`, `source` (reference/version/context), `genes` (unique audited var_names IDs), `feature_mapping` (identity mapping from these IDs to source var_names), `measured_genes` (exact intersection with source var_names), and substantive strings `specificity`, `complementarity`, `cell_coordination`, `target_comparison`, `coverage_limitations`. Resolve symbols to unique feature IDs before writing this contract; never guess aliases. `coordinated_program_supported: true` records the expert conclusion, not a two-marker count threshold. `foreign_markers` must be measured members of that program. Source, full marker/panel distributions, cell-level coordination, target comparison, unit and QC evidence must be linked in `evidence_files`; manifest tables are recursively bound by hashes. Blank/missing programs, mismapped markers, missing evidence or wrong coverage cannot freeze. Schema checks cannot authenticate sources or establish scientific truth.

Automatic review JSON:

```json
{"round":"r01", "proposal_hash":"HASH", "reviewer":"operating agent identifier",
 "decisions":[{"proposal_id":"p1", "decision":"supported",
   "reason":"evidence-based removal rationale", "alternatives":"assessed alternative explanations",
   "target_loss_assessment":"potential loss and justification",
   "target_loss_concern":false, "biological_mixed_state_unresolved":false,
   "target_subtype_unresolved":false, "partial_cluster_only":false,
   "removal_basis":"outside_target_scope",
   "target_positive":{"count":12,"fraction":0.12,"definition":"declared target-positive rule"}}]}
```

The example count/fraction must be replaced with measurements and match frozen cluster size. `removal_basis` is `technical_mixture` or `outside_target_scope`, supported by evidence. Every exclusion ID must appear exactly once. `decision: review` needs reason, alternatives and target_loss_assessment but cannot delete cells. Never set concerns false merely to pass the interface. If none are supported, retain all cells and stop with evidence_insufficient when doubts remain. After any applied removal, scientific stopping requires the next fully recomputed and reviewed round. Legacy freezes lacking this evidence rule cannot be automatically reviewed; historical records remain unchanged. The immutable automatic record and exact-cell digest are checked again on apply. It is an agent attestation distinct from user approval, not an authentication system.

## Final handoff selection

All per-round artifacts and the cumulative ledger above remain internal. For the
final user handoff apply [the shared delivery contract](../../../workflow/workflow-sc-analysis/references/delivery-contract.md):
select only final-round UMAP/DotPlot PNGs, per-cluster review CSV and README.
Summarize attrition and stopping in README with internal evidence paths. Unresolved
retained groups stay explicit; delivery selection never implies convergence.
