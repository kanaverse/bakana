library(alabaster.base)
library(SingleCellExperiment)
library(testthat)

test_that("simple runs can reload correctly", {
    validateObject(file.path("from-tests", "10X"))
    X <- readObject(file.path("from-tests", "10X"))
    expect_identical(as.character(sort(unique(X$clusters))), names(rowData(X)$marker_detection))

    res.dir <- file.path("from-tests", "10X_genes")
    validateObject(file.path(res.dir, "feature_selection"))
    marker.dir <- file.path(res.dir, "marker_detection")
    for (mod in list.files(marker.dir)) {
        moddir <- file.path(marker.dir, mod)
        for (clust in list.files(moddir)) {
            validateObject(file.path(moddir, clust))
        }
    }
    validateObject(file.path("from-tests", "10X"))

    info <- acquireMetadata("from-tests", "MatrixMarket")
    X <- readObject(info, "from-tests")

    info <- acquireMetadata("from-tests", "H5AD")
    X <- readObject(info, "from-tests")

    info <- acquireMetadata("from-tests", "se")
    X <- readObject(info, "from-tests")
})

test_that("block information is correctly saved", {
    info <- acquireMetadata("from-tests", "combined")
    X <- readObject(info, "from-tests")
    expect_type(X$block, "integer")
    expect_identical(metadata(X)$block_levels[X$block+1], X$named_block)
})

test_that("custom selections are correctly saved", {
    info <- acquireMetadata("from-tests", "custom")
    X <- readObject(info, "from-tests")
    expect_type(metadata(X)$custom_selections$evens, "integer")
    expect_identical(colnames(rowData(X)$custom_selections), names(metadata(X)$custom_selections))
})

test_that("ADTs are correctly saved", {
    info <- acquireMetadata("from-tests", "adt")
    X <- readObject(info, "from-tests")
    expect_true(nrow(altExp(X, "ADT")) < nrow(X))
})

test_that("CRISPR guides are correctly saved", {
    info <- acquireMetadata("from-tests", "crispr")
    X <- readObject(info, "from-tests")
    expect_true(nrow(altExp(X, "CRISPR")) < nrow(X))
})
