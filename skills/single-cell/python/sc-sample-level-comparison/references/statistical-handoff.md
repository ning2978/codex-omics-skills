# R statistical handoff

Hand off genes x unit-condition raw pseudobulk counts, sample table, eligible-unit/cell ledger, feature map and explicit contrast. Validate exact keys/order. Reassess eligibility/design per type. Fractional corrections are not raw counts.

## edgeR skeleton

Assume counts and corresponding sample rows already pass independent-unit/paired gates. Set reference condition and model formula explicitly: independent `~ condition + technical_batch` if justified/estimable; paired `~ donor + condition`. Contrasts derive from actual design column names.

```r
samples$condition <- relevel(factor(samples$condition), ref = reference_condition)
design <- model.matrix(model_formula, data = samples)
stopifnot(identical(colnames(counts), rownames(samples)))
if (qr(design)$rank < ncol(design) || nrow(design) <= qr(design)$rank)
  stop("Nonestimable design or no residual df; descriptive-only")
stopifnot(length(contrast) == ncol(design), any(contrast != 0))
y <- edgeR::DGEList(counts = counts)
keep <- edgeR::filterByExpr(y, design = design)
if (!any(keep)) stop("No eligible expressed genes")
y <- y[keep, , keep.lib.sizes = FALSE]
y <- edgeR::calcNormFactors(y, method = "TMM")
y <- edgeR::estimateDisp(y, design)
fit <- edgeR::glmQLFit(y, design, robust = TRUE)
test <- edgeR::glmQLFTest(fit, contrast = contrast)
result <- edgeR::topTags(test, n = Inf, sort.by = "none")$table
```

Under exact `Rscript --vanilla`, log `sessionInfo()`, `.libPaths()`, `find.package()` and executable. A small real fit verifies dispersion/robust dependencies (statmod where needed), not just imports. Technical gates do not prove power; inspect dispersion/library sizes/influential units.

## Paired/covariate propeller

Check installed `args(speckle::getTransformedProps)`/`args(speckle::propeller)` and return fields. Do not assume convenience `propeller()` supports paired designs. For versions with `getTransformedProps(clusters, sample, transform=...)`, `TransformedProps` is types x samples and `Proportions` holds raw proportions. Align columns to design explicitly; use limma `lmFit`, `contrasts.fit`, `eBayes`, `topTable` on transformed proportions with `~ donor + condition` for pairs. Export transformed-scale inference separately from raw-proportion effects. Fail on missing/duplicate units/nonestimable designs; adapt and log API differences, never substitute cell-level testing.

Absent types in a sampled denominator get zero proportions; they do not create zero-expression pseudobulk DEG columns. Record multiple testing/sample support separately for DEG, proportions and enrichment families.
