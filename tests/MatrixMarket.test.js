import * as bakana from "../src/index.js";
import * as butils from "../src/steps/utils/general.js";
import * as scran from "scran.js";
import * as utils from "./utils.js";

beforeAll(utils.initializeAll);
afterAll(async () => await bakana.terminate());

test("runAnalysis works correctly (MatrixMarket)", async () => {
    let attempts = new Set;
    let started = step => {
        attempts.add(step);
    };

    let completed = new Set;
    let finished = (step) => {
        completed.add(step);
    };

    let mtx_file = "files/datasets/pbmc3k-matrix.mtx.gz";
    let feat_file = "files/datasets/pbmc3k-features.tsv.gz";
    let files = { 
        default: new bakana.TenxMatrixMarketDataset(mtx_file, feat_file, "files/datasets/pbmc3k-barcodes.tsv.gz")
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
//    for (const [k, v] of Object.entries(params)) {
//        if (butils.changedParameters(new_params[k], v)) {
//            console.log([k, new_params[k], v]);
//        }
//    }
    expect(butils.changedParameters(new_params, params)).toBe(false); // should be the same after a roundtrip.

    {
        // Check that steps unserialize correctly.
        let old_keys = Object.keys(state);
        old_keys.sort();
        let new_keys = Object.keys(reloaded);
        new_keys.sort();
        expect(old_keys).toEqual(new_keys);

        for (const step of old_keys) {
            let qc_deets = reloaded[step].fetchParameters();
            let ref = state[step].fetchParameters();
            expect(ref).toEqual(ref);
        }

        // Check that we still get some markers.
        let reloaded_markers = reloaded.marker_detection.fetchResults();
        let ref_markers = state.marker_detection.fetchResults();
        expect(reloaded_markers["RNA"].cohen(0)).toEqual(ref_markers["RNA"].cohen(0));
    }

    {
        // While the ADT itself is a no-op, the parameters should still be okay.
        expect(new_params.adt_normalization.num_pcs).toBeGreaterThan(0);
        expect(new_params.adt_normalization.num_clusters).toBeGreaterThan(0);
        expect(new_params.combine_embeddings.weights).toBeNull();
        expect(new_params.batch_correction.num_neighbors).toBeGreaterThan(0);
    }

    // Checking that the permutation is unchanged on reload.
    let old_ids = state.inputs.fetchRowIds()["RNA"];
    let new_ids = reloaded.inputs.fetchRowIds()["RNA"];
    expect(old_ids.length).toBeGreaterThan(0);
    expect(old_ids).toEqual(new_ids);

    let old_ids2 = state.inputs.fetchFeatureAnnotations()["RNA"].column("id");
    let new_ids2 = reloaded.inputs.fetchFeatureAnnotations()["RNA"].column("id");
    expect(old_ids2.length).toBeGreaterThan(0);
    expect(old_ids2).toEqual(new_ids2);

    let old_fres = state.feature_selection.fetchResults();
    let new_fres = reloaded.feature_selection.fetchResults();
    expect(old_fres.means()).toEqual(new_fres.means());

    // Release me!
    await bakana.freeAnalysis(state);
    await bakana.freeAnalysis(reloaded);
})

test("runAnalysis works correctly with the bare minimum (MatrixMarket)", async () => {
    let state = await bakana.createAnalysis();
    let params = utils.baseParams();
    let res = await bakana.runAnalysis(state, 
        { 
            default: new bakana.TenxMatrixMarketDataset("files/datasets/pbmc3k-matrix.mtx.gz", null, null)
        },
        params
    );

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

    let new_params = bakana.retrieveParameters(reloaded);
    expect(new_params.quality_control instanceof Object).toBe(true);
    expect(new_params.pca instanceof Object).toBe(true);

    // Release me!
    await bakana.freeAnalysis(state);
    await bakana.freeAnalysis(reloaded);
})
