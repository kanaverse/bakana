import * as bakana from "../src/index.js";
import * as butils from "../src/steps/utils/general.js";
import * as scran from "scran.js";
import * as utils from "./utils.js";

beforeAll(utils.initializeAll);
afterAll(async () => await bakana.terminate());

test("runAnalysis works correctly (MatrixMarket)", async () => {
    let attempts = new Set();
    let started = step => {
        attempts.add(step);
    };

    let contents = {};
    let finished = (step, res) => {
        contents[step] = res;
    };

    let mtx_file = "files/datasets/pbmc3k-matrix.mtx.gz";
    let feat_file = "files/datasets/pbmc3k-features.tsv.gz";
    let files = { 
        default: new bakana.TenxMatrixMarketDataset(mtx_file, feat_file, "files/datasets/pbmc3k-barcodes.tsv.gz")
    };

    let state = await bakana.createAnalysis();
    let params = utils.baseParams();
    let res = await bakana.runAnalysis(state, files, params, { startFun: started, finishFun: finished });

    expect(attempts.has("pca")).toBe(true);
    expect(contents.pca instanceof Object).toBe(true);
    expect(contents.feature_selection instanceof Object).toBe(true);
    expect(contents.cell_labelling instanceof Object).toBe(true);

    // Input reorganization is done correctly.
    {
        let loaded = state.inputs.fetchCountMatrix();
        let loaded_names = state.inputs.fetchGenes().column("id");
        let loaded_ids = state.inputs.fetchRowIds();

        let simple = scran.initializeSparseMatrixFromMatrixMarket(mtx_file, { layered: false });
        let parsed = bakana.readTable((new bakana.SimpleFile(feat_file)).buffer(), { compression: "gz" });
        let simple_names = parsed.map(x => x[0]);

        utils.checkReorganization(simple.matrix, simple.row_ids, simple_names, loaded, loaded_ids, loaded_names);
        simple.matrix.free();
    }

    // Quality control.
    {
        expect(attempts.has("quality_control")).toBe(true);
        expect(contents.quality_control instanceof Object).toBe(true);
        expect(contents.quality_control.thresholds.default.proportion).toBeGreaterThan(0);

        // Undoing the filtering works as expected.
        let last_filtered = contents.cell_filtering.retained - 1;
        let idx = [0, last_filtered];
        state.cell_filtering.undoFiltering(idx);
        expect(idx[1]).toBeGreaterThan(last_filtered);
        expect(state.inputs.fetchCountMatrix().column(idx[0])).toEqual(state.cell_filtering.fetchFilteredMatrix().column(0));
        expect(state.inputs.fetchCountMatrix().column(idx[1])).toEqual(state.cell_filtering.fetchFilteredMatrix().column(last_filtered));
    }

    // Markers.
    {
        expect(contents.marker_detection instanceof Object).toBe(true);
        expect(state.marker_detection.numberOfGroups()).toBeGreaterThan(0);

        let res = state.marker_detection.fetchGroupResults(0, "cohen-mean", "RNA");
        expect("ordering" in res).toBe(true);
        expect("means" in res).toBe(true);
        expect("lfc" in res).toBe(true);
    }

    // In versus mode.
    {
        let vres = state.marker_detection.computeVersus(3, 0, "auc", "RNA");
        expect("ordering" in vres).toBe(true);
        expect("lfc" in vres).toBe(true);

        let vres2 = state.marker_detection.computeVersus(0, 3, "auc", "RNA");
        expect("ordering" in vres2).toBe(true);
        expect("lfc" in vres2).toBe(true);

        let lfcs = new Array(vres.lfc.length);
        vres.ordering.forEach((x, i) => {
            lfcs[x] = vres.lfc[i];
        });

        let lfcs2 = new Array(vres2.lfc.length);
        vres2.ordering.forEach((x, i) => {
            lfcs2[x] = -vres2.lfc[i];
        });

        expect(lfcs).toEqual(lfcs2);
    }

    // Normalized expression.
    {
        let exprs = state.normalization.fetchExpression(0);
        let nfiltered = state.cell_filtering.fetchFilteredMatrix().numberOfColumns();
        expect(exprs.length).toBe(nfiltered);
    }

    // Clustering.
    {
        let filtered_cells = contents.cell_filtering.retained;
        expect(contents.choose_clustering.clusters.length).toBe(filtered_cells);

        let nclusters = (new Set(contents.choose_clustering.clusters)).size;
        let expected = new Int32Array(filtered_cells);
        for (var c = 0; c < nclusters; c++) {
            let idx = state.choose_clustering.fetchClusterIndices(c);
            idx.forEach(x => { expected[x] = c; });
        }
        expect(contents.choose_clustering.clusters).toEqual(expected);
    }

    // ADTs are no-ops.
    {
        expect(state.adt_quality_control.summary()).toBeNull();
        expect(state.adt_normalization.summary()).toBeNull();
        expect(state.adt_pca.summary()).toBeNull();
    }

    // Animation catcher workers correctly.
    {
        let collected = { tsne: [], umap: [] }
        let fun = bakana.setVisualizationAnimate((type, x, y, iterations) => {
            collected[type].push(iterations);
        });

        let p = [state.tsne.animate(), state.umap.animate()];
        await Promise.all(p);
        bakana.setVisualizationAnimate(null)

        expect(collected.tsne.length).toBeGreaterThan(0);
        expect(collected.umap.length).toBeGreaterThan(0);
    }

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
            let qc_deets = reloaded[step].summary();
            let ref = state[step].summary();
            expect(ref).toEqual(ref);
        }

        // Check that we still get some markers.
        let reloaded_markers = reloaded.marker_detection.fetchGroupResults(0, "cohen-mean", "RNA");
        let ref_markers = reloaded.marker_detection.fetchGroupResults(0, "cohen-mean", "RNA");
        expect(reloaded_markers).toEqual(ref_markers);
    }

    {
        // While the ADT itself is a no-op, the parameters should still be okay.
        expect(new_params.adt_normalization.num_pcs).toBeGreaterThan(0);
        expect(new_params.adt_normalization.num_clusters).toBeGreaterThan(0);
        expect(new_params.combine_embeddings.weights).toBeNull();
        expect(new_params.batch_correction.num_neighbors).toBeGreaterThan(0);
    }

    // Checking that the permutation is unchanged on reload.
    let old_ids = state.inputs.summary()["genes"]["RNA"]["id"];
    let new_ids = reloaded.inputs.summary()["genes"]["RNA"]["id"];
    expect(old_ids.length).toBeGreaterThan(0);
    expect(old_ids).toEqual(new_ids);

    let old_res = state.feature_selection.summary();
    let new_res = reloaded.feature_selection.summary();
    expect("means" in old_res).toBe(true);
    expect(old_res["means"]).toEqual(new_res["means"]);

    // Release me!
    await bakana.freeAnalysis(state);
    await bakana.freeAnalysis(reloaded);
})

test("runAnalysis works correctly with the bare minimum (MatrixMarket)", async () => {
    let contents = {};
    let finished = (step, res) => {
        contents[step] = res;
    };

    let state = await bakana.createAnalysis();
    let params = utils.baseParams();
    let res = await bakana.runAnalysis(state, 
        { 
            default: new bakana.TenxMatrixMarketDataset("files/datasets/pbmc3k-matrix.mtx.gz", null, null)
        },
        params,
        {
            finishFun: finished,
        }
    );

    expect(contents.quality_control.thresholds.default.proportion).toBe(0);

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
