#!/usr/bin/env Rscript
# A real synthetic QL fit; verifies package plumbing, not biological validity/power.
suppressPackageStartupMessages(library(edgeR))
set.seed(917)
samples <- data.frame(condition=factor(rep(c("control","case"),4),levels=c("control","case")), donor=factor(rep(seq_len(4),each=2)))
rownames(samples) <- paste0("s",seq_len(nrow(samples)))
counts <- matrix(rnbinom(200*nrow(samples),mu=100,size=10),nrow=200,
                 dimnames=list(paste0("g",seq_len(200)),rownames(samples)))
counts[1:20,samples$condition=="case"] <- counts[1:20,samples$condition=="case"]+100
design <- model.matrix(~donor+condition,samples)
stopifnot(qr(design)$rank==ncol(design),nrow(design)>ncol(design))
y <- DGEList(counts=counts)
y <- y[filterByExpr(y,design=design),,keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
y <- estimateDisp(y,design)
fit <- glmQLFit(y,design,robust=TRUE)
contrast <- rep(0,ncol(design)); contrast[match("conditioncase",colnames(design))] <- 1
if (anyNA(contrast) || all(contrast==0)) stop("Expected design contrast absent")
test <- glmQLFTest(fit,contrast=contrast)
result <- topTags(test,n=Inf,sort.by="none")$table
stopifnot(all(is.finite(result$logFC)),all(is.finite(result$PValue)),all(is.finite(result$FDR)),
          median(result[1:20,"logFC"])>0)
cat("PASS edgeR paired QL fit; tested genes:",nrow(result)," residual df:",nrow(design)-qr(design)$rank,"\n")
cat("R:",R.version.string,"\nedgeR:",as.character(packageVersion("edgeR"))," path:",find.package("edgeR"),"\n")
print(.libPaths())
print(sessionInfo())
