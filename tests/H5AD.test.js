import * as bakana from "../src/index.js";
import * as scran from "scran.js";
import * as utils from "./utils.js";
import * as bioc from "bioconductor";

beforeAll(utils.initializeAll);
afterAll(async () => await bakana.terminate());

test("H5AD dataset loader works correctly", async () => {
    const fpath = "files/datasets/zeisel-brain.h5ad";
    let ds = new bakana.H5adDataset(fpath);
    await utils.checkDatasetGeneral(ds);

    let summ = await utils.checkDatasetSummary(ds);
    expect("all_features" in summ).toBe(true);

    let loaded = await utils.checkDatasetLoad(ds);
    expect(loaded.matrix.available()).toEqual(["RNA"]);

    let copy = await utils.checkDatasetSerialize(ds);
    utils.sameDatasetSummary(summ, await copy.summary());
    utils.sameDatasetLoad(loaded, await copy.load());

    // Input reorganization is done correctly.
    {
        let mat = loaded.matrix.get("RNA");
        let names = loaded.primary_ids["RNA"];
        let simple = scran.initializeSparseMatrixFromHdf5(fpath, "X", { layered: false });
        let simple_names = (new scran.H5File(fpath)).open("var").open("_index", { load: true }).values;
        utils.checkMatrixContents(simple, simple_names, mat, names);
        simple.free();
    }

    ds.clear();
})

test("H5AD result readers work correctly with log-count loading", async () => {
    let res = new bakana.H5adResult("files/datasets/zeisel-brain-with-results.h5ad");

    let summ = await utils.checkResultSummary(res);
    expect(summ.all_assay_names).toEqual(["X", "layers/logcounts"]);
    expect(summ.reduced_dimension_names).toEqual(["PCA", "TSNE", "UMAP"]);

    res.setOptions({
        primaryMatrixName: "layers/logcounts",
        isPrimaryNormalized: true,
        reducedDimensionNames: ["TSNE", "UMAP"]
    });

    let loaded = await utils.checkResultLoad(res);
    expect(Object.keys(loaded.reduced_dimensions)).toEqual(["TSNE", "UMAP"]);
    expect(loaded.reduced_dimensions["TSNE"].length).toEqual(2);
    expect(loaded.reduced_dimensions["UMAP"].length).toEqual(2);

    let col0 = loaded.matrix.get("").column(0);
    expect(utils.hasNonInteger(col0)).toBe(true);

    res.clear();
})

test("H5AD result readers work correctly with manual count normalization", async () => {
    let res = new bakana.H5adResult("files/datasets/zeisel-brain-with-results.h5ad");
    res.setOptions({ primaryMatrixName: "X" });

    let loaded = res.load();
    let col0 = loaded.matrix.get("").column(0);
    expect(utils.hasNonInteger(col0)).toBe(false);

    // Trying again after flagging it as normalizable.
    res.setOptions({ isPrimaryNormalized: false });
    let loaded2 = res.load();
    let col0_2 = loaded2.matrix.get("").column(0);
    expect(utils.hasNonInteger(col0_2)).toBe(true);

    res.clear();
})
