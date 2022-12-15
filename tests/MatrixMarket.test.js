import * as bakana from "../src/index.js";
import * as butils from "../src/steps/utils/general.js";
import * as scran from "scran.js";
import * as utils from "./utils.js";
import * as bioc from "bioconductor";

beforeAll(utils.initializeAll);
afterAll(async () => await bakana.terminate());

let mtx_file = "files/datasets/pbmc3k-matrix.mtx.gz";
let feat_file = "files/datasets/pbmc3k-features.tsv.gz";
let files = { 
    default: new bakana.TenxMatrixMarketDataset(mtx_file, feat_file, "files/datasets/pbmc3k-barcodes.tsv.gz")
};

test("MatrixMarket summary works correctly", async () => {
    let summ = await files.default.summary();
    expect(summ.all_features instanceof bioc.DataFrame).toBe(true);
    expect(summ.all_features.numberOfColumns()).toBeGreaterThan(0);
    expect(summ.cells instanceof bioc.DataFrame).toBe(true);
    expect(summ.cells.numberOfColumns()).toBeGreaterThan(0);
})

test("runAnalysis works correctly (MatrixMarket)", async () => {
    let attempts = new Set;
    let started = step => {
        attempts.add(step);
    };

    let completed = new Set;
    let finished = (step) => {
        completed.add(step);
    };

    let state = await bakana.createAnalysis();
    let params = utils.baseParams();
    let res = await bakana.runAnalysis(state, files, params, { startFun: started, finishFun: finished });

    // Check that the callbacks are actually called.
    expect(attempts.has("quality_control")).toBe(true);
    expect(attempts.has("pca")).toBe(true);
    expect(completed.has("pca")).toBe(true);
    expect(completed.has("feature_selection")).toBe(true);
    expect(completed.has("cell_labelling")).toBe(true);

    // Input reorganization is done correctly.
    {
        let loaded = state.inputs.fetchCountMatrix().get("RNA");
        let loaded_names = state.inputs.fetchFeatureAnnotations()["RNA"].column("id");
        let loaded_ids = state.inputs.fetchRowIds()["RNA"];

        let simple = scran.initializeSparseMatrixFromMatrixMarket(mtx_file, { layered: false });
        let parsed = bakana.readTable((new bakana.SimpleFile(feat_file)).buffer(), { compression: "gz" });
        let simple_names = parsed.map(x => x[0]);

        utils.checkReorganization(simple.matrix, simple.row_ids, simple_names, loaded, loaded_ids, loaded_names);
        simple.matrix.free();
    }

    // Basic consistency checks.
    await utils.checkStateResultsSimple(state);
    utils.checkClusterVersusMode(state);
    await utils.triggerAnimation(state);

    // Saving and loading.
    const path = "TEST_state_MatrixMarket.h5";
    let collected = await bakana.saveAnalysis(state, path);
    utils.validateState(path);
    expect(collected.collected.length).toBe(3);
    expect(typeof(collected.collected[0])).toBe("string");

    let offsets = utils.mockOffsets(collected.collected);
    let reloaded = await bakana.loadAnalysis(
        path, 
        (offset, size) => offsets[offset]
    );

    let new_params = bakana.retrieveParameters(reloaded);
    expect(butils.changedParameters(new_params, params)).toBe(false); // should be the same after a roundtrip.

    // Checking that the results are the same.
    await utils.compareStates(state, reloaded);

    // While the ADT itself is a no-op, the parameters should still be okay.
    expect(new_params.adt_normalization.num_pcs).toBeGreaterThan(0);
    expect(new_params.adt_normalization.num_clusters).toBeGreaterThan(0);
    expect(new_params.combine_embeddings.weights).toBeNull();
    expect(new_params.batch_correction.num_neighbors).toBeGreaterThan(0);

    // Release me!
    await bakana.freeAnalysis(state);
    await bakana.freeAnalysis(reloaded);
})

let minimal_files = { 
    default: new bakana.TenxMatrixMarketDataset("files/datasets/pbmc3k-matrix.mtx.gz", null, null)
};

test("MatrixMarket summary works correctly with the bare minimum", async () => {
    let summ = await minimal_files.default.summary();
    expect(summ.all_features instanceof bioc.DataFrame).toBe(true);
    expect(summ.all_features.numberOfColumns()).toBe(0);
    expect(summ.cells instanceof bioc.DataFrame).toBe(true);
    expect(summ.cells.numberOfColumns()).toBe(0);
})

test("runAnalysis works correctly with the bare minimum (MatrixMarket)", async () => {
    let state = await bakana.createAnalysis();
    let params = utils.baseParams();
    let res = await bakana.runAnalysis(state, minimal_files, params);

    // Basic consistency checks.
    await utils.checkStateResultsSimple(state);

    // No annotations, so no mitochondrial proportions.
    expect(state.inputs.fetchFeatureAnnotations()["RNA"].numberOfColumns()).toBe(0);
    expect(state.quality_control.fetchFilters().thresholdsSubsetProportions()[0]).toBe(0);

    // Saving and loading.
    const path = "TEST_state_MatrixMarket.h5";
    let collected = await bakana.saveAnalysis(state, path);
    utils.validateState(path);
    expect(collected.collected.length).toBe(1);
    expect(typeof(collected.collected[0])).toBe("string");

    let offsets = utils.mockOffsets(collected.collected);
    let reloaded = await bakana.loadAnalysis(
        path, 
        (offset, size) => offsets[offset]
    );

    await utils.compareStates(state, reloaded);

    let new_params = bakana.retrieveParameters(reloaded);
    expect(new_params).toEqual(params);

    // Release me!
    await bakana.freeAnalysis(state);
    await bakana.freeAnalysis(reloaded);
})
