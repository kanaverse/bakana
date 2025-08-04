library(alabaster.base)
library(SingleCellExperiment)
library(testthat)

check_genes <- function(dir) {
    res.dir <- file.path("from-tests", dir)
    expect_error(validateObject(file.path(res.dir, "feature_selection")), NA)

    marker.dir <- file.path(res.dir, "marker_detection")
    for (mod in list.files(marker.dir)) {
        moddir <- file.path(marker.dir, mod)
        for (clust in list.files(moddir)) {
            validateObject(file.path(moddir, clust))
        }
    }

    custom.dir <- file.path(res.dir, "custom_selections")
    if (file.exists(custom.dir)) {
        for (mod in list.files(marker.dir)) {
            moddir <- file.path(marker.dir, mod)
            for (sel in list.files(moddir)) {
                validateObject(file.path(moddir, sel))
            }
        }
    }
}

test_that("basic PBMC output is correct", {
    expect_error(validateObject(file.path("from-tests", "PBMC")), NA)
    check_genes("PBMC_genes")
})

test_that("basic zeisel output is correct", {
    expect_error(validateObject(file.path("from-tests", "zeisel")), NA)
    check_genes("zeisel_genes")
})

test_that("block information is correctly saved", {
    expect_error(validateObject(file.path("from-tests", "combined")), NA)
    check_genes("combined_genes")

    X <- readObject(file.path("from-tests", "combined"))
    expect_type(colData(X)[["kana::block"]], "character")
})

test_that("custom selections are correctly saved", {
    expect_error(validateObject(file.path("from-tests", "custom")), NA)
    check_genes("custom_genes")

    X <- readObject(file.path("from-tests", "custom"))
    expect_type(colData(X)[["kana::custom_selections::evens"]], "logical")
})

test_that("ADTs are correctly saved", {
    expect_error(validateObject(file.path("from-tests", "adt")), NA)
    check_genes("adt_genes")

    X <- readObject(file.path("from-tests", "adt"))
    expect_true(nrow(altExp(X, "ADT")) < nrow(X))
    expect_true(any(grepl("kana::ADT::quality_control", colnames(colData(X)))))
    expect_false(any(grepl("kana::quality_control", colnames(colData(altExp(X, "ADT"))))))

    X <- readObject(file.path("from-tests", "adt_split"))
    expect_false(any(grepl("kana::ADT::quality_control", colnames(colData(X)))))
    expect_true(any(grepl("kana::quality_control", colnames(colData(altExp(X, "ADT"))))))
})

test_that("CRISPR guides are correctly saved", {
    expect_error(validateObject(file.path("from-tests", "crispr")), NA)
    check_genes("crispr_genes")

    X <- readObject(file.path("from-tests", "crispr"))
    expect_true(nrow(altExp(X, "CRISPR")) < nrow(X))
})
