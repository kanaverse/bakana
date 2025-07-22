# This generates alabaster-formatted object files for testing the alabaster readers.

library(alabaster.sce)
unlink("alabaster", recursive=TRUE)
dir.create("alabaster", showWarnings=FALSE)

sce <- readRDS("datasets/zeisel-brain.rds")
assay(sce) <- as(assay(sce), "dgCMatrix")
saveObject(sce, "alabaster/zeisel-brain")

mainExpName(sce) <- NULL 
altExps(sce) <- list()
saveObject(sce, "alabaster/zeisel-brain-stripped")

sce <- readRDS("datasets/zeisel-brain-dense.rds")
saveObject(sce, "alabaster/zeisel-brain-dense")

sce <- readRDS("datasets/zeisel-brain-with-results.rds")
saveObject(sce, "alabaster/zeisel-brain-with-results")
