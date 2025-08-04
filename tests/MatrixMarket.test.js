import * as bakana from "../src/index.js";
import * as scran from "scran.js";
import * as utils from "./utils.js";
import * as bioc from "bioconductor";

beforeAll(utils.initializeAll);
afterAll(async () => await bakana.terminate());

test("MatrixMarket reader works with all files", async () => {
    let mtx_file = "files/datasets/pbmc3k-matrix.mtx.gz";
    let feat_file = "files/datasets/pbmc3k-features.tsv.gz";
    let ds = new bakana.TenxMatrixMarketDataset(mtx_file, feat_file, "files/datasets/pbmc3k-barcodes.tsv.gz")
    await utils.checkDatasetGeneral(ds);

    let summ = await utils.checkDatasetSummary(ds);
    expect(Object.keys(summ.modality_features)).toEqual(["Gene Expression"]);

    let loaded = await utils.checkDatasetLoad(ds);
    expect(loaded.matrix.available()).toEqual(["RNA"]);

    let copy = await utils.checkDatasetSerialize(ds);
    utils.sameDatasetSummary(summ, await copy.summary());
    utils.sameDatasetLoad(loaded, await copy.load());

    // Input reorganization is done correctly.
    {
        let mat = loaded.matrix.get("RNA");
        let names = loaded.primary_ids["RNA"];
        let simple = scran.initializeSparseMatrixFromMatrixMarket(mtx_file, { layered: false });
        let parsed = bakana.readTable((new bakana.SimpleFile(feat_file)).buffer(), { compression: "gz" });
        let simple_names = parsed.map(x => x[0]);
        utils.checkMatrixContents(simple, simple_names, mat, names);
        simple.free();
    }

    ds.clear();
})

test("MatrixMarket reader works with minimal files", async () => {
    let mtx_file = "files/datasets/pbmc3k-matrix.mtx.gz";
    let ds = new bakana.TenxMatrixMarketDataset(mtx_file, null, null);
    await utils.checkDatasetGeneral(ds);

    let summ = await utils.checkDatasetSummary(ds);
    expect(Object.keys(summ.modality_features)).toEqual([""]);
    expect(summ.modality_features[""] instanceof bioc.DataFrame).toBe(true);
    expect(summ.modality_features[""].numberOfColumns()).toBe(0);
    expect(summ.cells instanceof bioc.DataFrame).toBe(true);
    expect(summ.cells.numberOfColumns()).toBe(0);

    let loaded = await utils.checkDatasetLoad(ds);
    expect(loaded.matrix.available()).toEqual(["RNA"]);
    expect(loaded.features["RNA"] instanceof bioc.DataFrame).toBe(true);
    expect(loaded.features["RNA"].numberOfColumns()).toBe(0);
    expect(loaded.cells instanceof bioc.DataFrame).toBe(true);
    expect(loaded.cells.numberOfColumns()).toBe(0);

    ds.clear();
})

test("MatrixMarket reader works with multiple modalities", async () => {
    let mtx_file = "files/datasets/immune_3.0.0-matrix.mtx.gz";
    let feat_file = "files/datasets/immune_3.0.0-features.tsv.gz";
    let ds = new bakana.TenxMatrixMarketDataset(mtx_file, feat_file, "files/datasets/immune_3.0.0-barcodes.tsv.gz")
    await utils.checkDatasetGeneral(ds);

    let summ = await utils.checkDatasetSummary(ds);
    expect(Object.keys(summ.modality_features)).toEqual(["Gene Expression", "Antibody Capture"]);

    let loaded = await utils.checkDatasetLoad(ds);
    expect(loaded.matrix.available()).toEqual(["RNA", "ADT"]);

    let copy = await utils.checkDatasetSerialize(ds);
    utils.sameDatasetSummary(summ, await copy.summary());
    utils.sameDatasetLoad(loaded, await copy.load());

    ds.clear();
})
