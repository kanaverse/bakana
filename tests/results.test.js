import * as bakana from "../src/index.js";
import * as utils from "./utils.js";
import * as bioc from "bioconductor";

beforeAll(utils.initializeAll);
afterAll(async () => await bakana.terminate());

function has_noninteger(x) {
    for (const y of x) {
        if (y !== Math.round(y)) {
            return true;
        }
    }
    return false;
}

test("H5AD result readers work correctly with log-count loading", async () => {
    let info = new bakana.H5adResult("files/datasets/zeisel-brain-with-results.h5ad");

    let details = info.summary();
    expect(details.all_features.numberOfRows()).toBeGreaterThan(0);
    expect(details.cells.numberOfRows()).toBeGreaterThan(0);
    expect(details.all_assay_names.length).toEqual(2);
    expect(details.reduced_dimension_names.length).toEqual(3);

    info.setPrimaryMatrixName("layers/logcounts");
    info.setIsPrimaryNormalized(true);
    info.setReducedDimensionNames(["TSNE", "UMAP"]);

    let payload = info.load();
    expect(payload.matrix.numberOfColumns()).toEqual(details.cells.numberOfRows());
    expect(payload.matrix.get("").numberOfRows()).toEqual(payload.features[""].numberOfRows());
    expect(Object.keys(payload.reduced_dimensions)).toEqual(["TSNE", "UMAP"]);
    expect(payload.reduced_dimensions["TSNE"].length).toEqual(2);
    expect(payload.reduced_dimensions["TSNE"][0].length).toEqual(details.cells.numberOfRows());

    let col0 = payload.matrix.get("").column(0);
    expect(has_noninteger(col0)).toBe(true);

    info.clear();
    payload.matrix.free();
})

test("H5AD result readers work correctly with manual count normalization", async () => {
    let info = new bakana.H5adResult("files/datasets/zeisel-brain-with-results.h5ad");
    info.setPrimaryMatrixName("X");

    let payload = info.load();
    expect(payload.matrix.numberOfColumns()).toEqual(payload.cells.numberOfRows());
    expect(payload.matrix.get("").numberOfRows()).toEqual(payload.features[""].numberOfRows());
    expect(Object.keys(payload.reduced_dimensions)).toEqual(["PCA", "TSNE", "UMAP"]);

    let col0 = payload.matrix.get("").column(0);
    expect(has_noninteger(col0)).toBe(false);

    // Trying again after flagging it as normalizable.
    info.setIsPrimaryNormalized(false);
    let payload2 = info.load();
    let col0_2 = payload2.matrix.get("").column(0);
    expect(has_noninteger(col0_2)).toBe(true);

    info.clear();
    payload.matrix.free();
    payload2.matrix.free();
})

test("SummarizedExperiment result readers work correctly with log-count loading", async () => {
    let info = new bakana.SummarizedExperimentResult("files/datasets/zeisel-brain-with-results.rds");

    let details = info.summary();
    expect(details.modality_features[""].numberOfRows()).toBeGreaterThan(0);
    expect(details.modality_features["ERCC"].numberOfRows()).toBeGreaterThan(0);

    expect(details.cells.numberOfRows()).toBeGreaterThan(0);
    expect(details.cells.column("label") instanceof Array).toBe(true); // factors are correctly expanded.

    expect(details.modality_assay_names[""].length).toEqual(2);
    expect(details.reduced_dimension_names.length).toEqual(3);

    info.setPrimaryAssay("logcounts");
    info.setIsPrimaryNormalized(true);
    info.setReducedDimensionNames(["TSNE", "UMAP"]);

    let payload = info.load();
    expect(payload.matrix.numberOfColumns()).toEqual(details.cells.numberOfRows());
    expect(payload.matrix.get("").numberOfRows()).toEqual(payload.features[""].numberOfRows());
//    expect(payload.matrix.has("ERCC")).toBe(false); // ignored if it doesn't have log-counts.

    expect(Object.keys(payload.reduced_dimensions)).toEqual(["TSNE", "UMAP"]);
    expect(payload.reduced_dimensions["TSNE"].length).toEqual(2);
    expect(payload.reduced_dimensions["TSNE"][0].length).toEqual(details.cells.numberOfRows());

    let col0 = payload.matrix.get("").column(0);
    expect(has_noninteger(col0)).toBe(true);

    info.clear();
    payload.matrix.free();
})

test("SummarizedExperiment result readers work correctly with count normalization", async () => {
    let info = new bakana.SummarizedExperimentResult("files/datasets/zeisel-brain-with-results.rds");

    info.setPrimaryAssay({ "": "counts", "ERCC": "counts" });
    info.setIsPrimaryNormalized({ "ERCC": false });

    let payload = info.load();
    expect(payload.matrix.numberOfColumns()).toEqual(payload.cells.numberOfRows());
    expect(payload.matrix.get("").numberOfRows()).toEqual(payload.features[""].numberOfRows());
    expect(payload.matrix.get("ERCC").numberOfRows()).toEqual(payload.features["ERCC"].numberOfRows());

    // Checking that naked gene counts are loaded, as a control.
    {
        let col0 = payload.matrix.get("").column(0);
        expect(has_noninteger(col0)).toBe(false);
    }

    // Checking that the ERCC count matrix is correctly normalized.
    {
        let col0 = payload.matrix.get("ERCC").column(0);
        expect(has_noninteger(col0)).toBe(true);
    }

    info.clear();
    payload.matrix.free();
})
