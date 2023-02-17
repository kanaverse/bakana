library(alabaster.base)
library(SingleCellExperiment)
library(testthat)

test_that("simple runs can reload correctly", {
    info <- acquireMetadata("from-tests", "10X")
    X <- loadObject(info, "from-tests")
    expect_identical(as.character(sort(unique(X$clusters))), names(rowData(X)$marker_detection))

    info <- acquireMetadata("from-tests", "MatrixMarket")
    X <- loadObject(info, "from-tests")

    info <- acquireMetadata("from-tests", "H5AD")
    X <- loadObject(info, "from-tests")

    info <- acquireMetadata("from-tests", "se")
    X <- loadObject(info, "from-tests")
})

test_that("block information is correctly saved", {
    info <- acquireMetadata("from-tests", "combined")
    X <- loadObject(info, "from-tests")
    expect_type(X$block, "integer")
    expect_identical(metadata(X)$block_levels[X$block+1], X$named_block)
})

test_that("custom selections are correctly saved", {
    info <- acquireMetadata("from-tests", "custom")
    X <- loadObject(info, "from-tests")
    expect_type(metadata(X)$custom_selections$evens, "integer")
    expect_identical(colnames(rowData(X)$custom_selections), names(metadata(X)$custom_selections))
})

test_that("ADTs are correctly saved", {
    info <- acquireMetadata("from-tests", "adt")
    X <- loadObject(info, "from-tests")
    expect_true(nrow(altExp(X, "ADT")) < nrow(X))
})

test_that("CRISPR guides are correctly saved", {
    info <- acquireMetadata("from-tests", "crispr")
    X <- loadObject(info, "from-tests")
    expect_true(nrow(altExp(X, "CRISPR")) < nrow(X))
})

test_that("all metadata files and MD5sums are valid", {
    all.json <- list.files("from-tests", pattern="\\.json$", recursive=TRUE)
    for (x in all.json) {
        xpath <- file.path("from-tests", x)
        info <- jsonlite::fromJSON(xpath)
        target <- system.file("schemas", info[["$schema"]], package="alabaster.schemas")

        expect_error(jsonvalidate::json_validate(xpath, target, engine="ajv", error=TRUE), NA, label=paste0("failed on '", x, "'"))
        if ("md5sum" %in% names(info)) {
            expect_identical(info$md5sum, digest::digest(file=file.path("from-tests", info$path)))
        }
    }
})
