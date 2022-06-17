import * as bakana from "../src/index.js";
import * as scran from "scran.js";
import * as utils from "./utils.js";

beforeAll(async () => await bakana.initialize({ localFile: true }));
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

    let state = await bakana.createAnalysis();
    let params = utils.baseParams();
    let res = await bakana.runAnalysis(state, 
        { 
            default: {
                format: "MatrixMarket",
                mtx: "files/datasets/pbmc3k-matrix.mtx.gz",
                genes: "files/datasets/pbmc3k-features.tsv.gz",
                annotations: "files/datasets/pbmc3k-barcodes.tsv.gz"
            }
        },
        params,
        {
            startFun: started,
            finishFun: finished
        }
    );

    expect(attempts.has("pca")).toBe(true);
    expect(contents.pca instanceof Object).toBe(true);
    expect(contents.feature_selection instanceof Object).toBe(true);
    expect(contents.cell_labelling instanceof Object).toBe(true);

    // Quality control.
    {
        expect(attempts.has("quality_control")).toBe(true);
        expect(contents.quality_control instanceof Object).toBe(true);
        expect(contents.quality_control.thresholds.default.proportion).toBeGreaterThan(0);
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

    // Normalized expression.
    {
        let exprs = state.normalization.fetchExpression(0);
        let nfiltered = state.cell_filtering.fetchFilteredMatrix().numberOfColumns();
        expect(exprs.length).toBe(nfiltered);
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
    expect(collected.collected.length).toBe(3);
    expect(typeof(collected.collected[0])).toBe("string");

    let offsets = utils.mockOffsets(collected.collected);
    let reloaded = await bakana.loadAnalysis(
        path, 
        (offset, size) => offsets[offset]
    );

    let new_params = reloaded.parameters;
    expect(new_params.quality_control instanceof Object).toBe(true);
    expect(new_params.pca instanceof Object).toBe(true);

    {
        // Check that steps unserialize correctly.
        let old_keys = Object.keys(state);
        old_keys.sort();
        let new_keys = Object.keys(reloaded.state);
        new_keys.sort();
        expect(old_keys).toEqual(new_keys);

        test.each(old_keys, step => {
            let qc_deets = reloaded[step].summary();
            let ref = state[step].summary();
            expect(ref).toEqual(ref);
        });
    }

    // Checking that the permutation is unchanged on reload.
    let old_ids = state.inputs.summary()["genes"]["RNA"]["id"];
    let new_ids = reloaded.state.inputs.summary()["genes"]["RNA"]["id"];
    expect(old_ids.length).toBeGreaterThan(0);
    expect(old_ids).toEqual(new_ids);

    let old_res = state.feature_selection.summary();
    let new_res = reloaded.state.feature_selection.summary();
    expect("means" in old_res).toBe(true);
    expect(old_res["means"]).toEqual(new_res["means"]);

    // Release me!
    await bakana.freeAnalysis(state);
    await bakana.freeAnalysis(reloaded.state);
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
            default: {
                format: "MatrixMarket",
                mtx: "files/datasets/pbmc3k-matrix.mtx.gz"
            }
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
    expect(collected.collected.length).toBe(1);
    expect(typeof(collected.collected[0])).toBe("string");

    let offsets = utils.mockOffsets(collected.collected);
    let reloaded = await bakana.loadAnalysis(
        path, 
        (offset, size) => offsets[offset]
    );

    let new_params = reloaded.parameters;
    expect(new_params.quality_control instanceof Object).toBe(true);
    expect(new_params.pca instanceof Object).toBe(true);

    // Release me!
    await bakana.freeAnalysis(state);
    await bakana.freeAnalysis(reloaded.state);
})
