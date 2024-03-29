import * as bakana from "../src/index.js";
import * as utils from "./utils.js";
import * as bioc from "bioconductor";
import * as nav from "./navigator.js";
import * as scran from "scran.js";

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

    info.setOptions({
        primaryMatrixName: "layers/logcounts",
        isPrimaryNormalized: true,
        reducedDimensionNames: ["TSNE", "UMAP"]
    });

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
    info.setOptions({ primaryMatrixName: "X" });

    let payload = info.load();
    expect(payload.matrix.numberOfColumns()).toEqual(payload.cells.numberOfRows());
    expect(payload.matrix.get("").numberOfRows()).toEqual(payload.features[""].numberOfRows());
    expect(Object.keys(payload.reduced_dimensions)).toEqual(["PCA", "TSNE", "UMAP"]);

    let col0 = payload.matrix.get("").column(0);
    expect(has_noninteger(col0)).toBe(false);

    // Trying again after flagging it as normalizable.
    info.setOptions({ isPrimaryNormalized: false });
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

    info.setOptions({ 
        primaryAssay: { "": "logcounts" },
        isPrimaryNormalized: true,
        reducedDimensionNames: ["TSNE", "UMAP"]
    });

    let payload = info.load();
    expect(payload.matrix.numberOfColumns()).toEqual(details.cells.numberOfRows());
    expect(payload.matrix.get("").numberOfRows()).toEqual(payload.features[""].numberOfRows());
    expect(payload.matrix.has("ERCC")).toBe(false); // ignored if it doesn't have log-counts.

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

    info.setOptions({
        primaryAssay: "counts",
        isPrimaryNormalized: { "ERCC": false }
    });

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

class LocalArtifactdbResult extends bakana.AbstractArtifactdbResult {
    constructor(path, dir, options={}) {
        super(path, new nav.LocalProjectDirectoryNavigator(dir), options);
    }
}

/***********************************************/

let target_simple = nav.pathExists("H5AD");
let action_simple = (target_simple == null ? test.skip : test);

action_simple("local ArtifactDB result readers work correctly with log-count loading", async () => {
    let info = new LocalArtifactdbResult(target_simple, nav.baseDirectory);

    let details = await info.summary();
    expect(details.modality_features[""].numberOfRows()).toBeGreaterThan(0);
    expect(details.modality_features[""].rowNames().length).toBeGreaterThan(0);

    expect(details.cells.numberOfRows()).toBeGreaterThan(0);
    expect(details.cells.hasColumn("kana::RNA::quality_control::sums")).toBe(true);
    expect(details.cells.hasColumn("kana::RNA::size_factors")).toBe(true);
    expect(details.cells.hasColumn("kana::clusters")).toBe(true);

    expect(details.modality_assay_names[""].length).toEqual(2);
    expect(details.reduced_dimension_names.length).toEqual(3);

    expect("custom_selections" in details.other_metadata).toBe(true);

    // Actually loading the log-counts.
    info.setOptions({
        primaryAssay: "logcounts",
        isPrimaryNormalized: true,
        reducedDimensionNames: ["TSNE", "UMAP"]
    });

    let payload = await info.load();
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

action_simple("local ArtifactDB result readers work correctly with count normalization", async () => {
    let info = new LocalArtifactdbResult(target_simple, nav.baseDirectory);

    info.setOptions({
        primaryAssay: "counts",
        isPrimaryNormalized: false
    });

    let payload = await info.load();
    expect(payload.matrix.numberOfColumns()).toEqual(payload.cells.numberOfRows());
    expect(payload.matrix.get("").numberOfRows()).toEqual(payload.features[""].numberOfRows());

    let col0 = payload.matrix.get("").column(0);
    expect(has_noninteger(col0)).toBe(true);

    info.clear();
    payload.matrix.free();
})

let target_adt = nav.pathExists("adt");
let action_adt = (target_adt == null ? test.skip : test);

action_adt("local ArtifactDB result readers work correctly with multiple modalities", async () => {
    let info = new LocalArtifactdbResult(target_adt, nav.baseDirectory);

    let details = await info.summary();
    expect(details.modality_features[""].numberOfRows()).toBeGreaterThan(0);
    expect(details.modality_features["ADT"].numberOfRows()).toBeGreaterThan(0);
    expect(details.cells.hasColumn("kana::ADT::quality_control::sums")).toBe(true);

    let payload = await info.load();
    let nRna = payload.features[""].numberOfRows();
    expect(nRna).toBeGreaterThan(0);
    expect(payload.matrix.get("").numberOfRows()).toEqual(nRna);

    let nAdt = payload.features["ADT"].numberOfRows();
    expect(nAdt).toBeGreaterThan(0);
    expect(nAdt).toBeLessThan(nRna);
    expect(payload.matrix.get("ADT").numberOfRows()).toEqual(nAdt);

    expect(payload.matrix.numberOfColumns()).toEqual(payload.cells.numberOfRows());
    payload.matrix.free();

    // What if it's split by modality?
    {
        let info = new LocalArtifactdbResult(target_adt + "_split", nav.baseDirectory);
        let details = await info.summary();
        expect(details.cells.hasColumn("kana::quality_control::sums")).toBe(true);
    }
})

/***********************************************/

action_simple("Zipped ArtifactDB result summary and loading works correctly", async () => {
    // First, zipping the contents of the target directory.
    let zipped = await nav.zipDirectory(nav.baseDirectory, [ target_simple, target_simple + ".json" ]);
    scran.writeFile("miscellaneous/bundle.zip", zipped); // for manual inspection.
    let zipfile = new bakana.SimpleFile(zipped, { name: "bundle.zip" });

    // Checking for consistency between zipped and native.
    let info = new bakana.ZippedArtifactdbResult(target_simple, zipfile);
    let ref = new LocalArtifactdbResult(target_simple, nav.baseDirectory);

    let details = await info.summary();
    let rdetails = await ref.summary();
    expect(Object.keys(details.modality_features)).toEqual(Object.keys(rdetails.modality_features));
    expect(details.modality_features[""].numberOfRows()).toEqual(rdetails.modality_features[""].numberOfRows());
    expect(details.modality_features[""].columnNames()).toEqual(rdetails.modality_features[""].columnNames());

    expect(details.cells.numberOfRows()).toEqual(rdetails.cells.numberOfRows());
    expect(details.cells.columnNames()).toEqual(rdetails.cells.columnNames());
    expect(details.modality_assay_names).toEqual(rdetails.modality_assay_names);
    expect(details.reduced_dimension_names).toEqual(rdetails.reduced_dimension_names);

    expect(Object.keys(details.other_metadata)).toEqual(Object.keys(rdetails.other_metadata));

    // Actually loading the log-counts.
    info.setOptions({
        primaryAssay: "logcounts",
        isPrimaryNormalized: true,
        reducedDimensionNames: ["TSNE", "UMAP"]
    });

    let payload = await info.load();
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
