import * as bakana from "../src/index.js";
import * as scran from "scran.js";
import * as utils from "./utils.js"

beforeAll(utils.initializeAll);
afterAll(async () => await bakana.terminate());

test("multi-matrix analyses work correctly", async () => {
    let contents = {};
    let finished = (step, res) => {
        contents[step] = res;
    };

    let files = { 
        "4K": {
            format: "10X",
            h5: "files/datasets/pbmc4k-tenx.h5"
        },
        "3K": {
            format: "MatrixMarket",
            mtx: "files/datasets/pbmc3k-matrix.mtx.gz",
            genes: "files/datasets/pbmc3k-features.tsv.gz",
            annotations: "files/datasets/pbmc3k-barcodes.tsv.gz"
        }
    };

    let paramcopy = utils.baseParams();
    let state = await bakana.createAnalysis();
    let res = await bakana.runAnalysis(
        state,
        files,
        paramcopy,
        {
            finishFun: finished,
        }
    );

    expect(contents.quality_control instanceof Object).toBe(true);
    expect("3K" in contents.quality_control.thresholds).toBe(true);
    expect("4K" in contents.quality_control.thresholds).toBe(true);

    {
        // Checking the blocking factors.
        expect(state.inputs.fetchBlockLevels()).toEqual(["3K", "4K"]); 
        let ids = state.inputs.fetchBlock();
        let counts = {};
        ids.forEach((x, i) => {
            if (x in counts) {
                counts[x] = 1;
            } else {
                counts[x]++;
            }
        });
        expect(Object.keys(counts).length).toBe(2);
        expect(counts[0]).toBeGreaterThan(0);
        expect(counts[1]).toBeGreaterThan(0);
    }

    {
        // Check that some correction actually occured.
        let corrected_pcs = state.batch_correction.fetchPCs();
        let original_pcs = state.pca.fetchPCs();
        expect(corrected_pcs.length).toEqual(original_pcs.length);
        expect(corrected_pcs.pcs.slice()).not.toEqual(original_pcs.pcs.slice());
    }

    // Saving and loading.
    const path = "TEST_state_multi-matrix.h5";
    let collected = await bakana.saveAnalysis(state, path);
    utils.validateState(path);
    expect(collected.collected.length).toBe(4);
    expect(typeof(collected.collected[0])).toBe("string");

    {
        let handle = new scran.H5File(path);
        let ihandle = handle.open("inputs");
        expect(ihandle.open("results").open("num_samples", { load: true }).values[0]).toEqual(2);

        let phandle = ihandle.open("parameters");
        let sample_names = phandle.open("sample_names", { load: true }).values;
        expect(sample_names).toEqual(["3K", "4K"]); // again, should be sorted.

        let sample_groups = phandle.open("sample_groups", { load: true }).values;
        expect(Array.from(sample_groups)).toEqual([3, 1]); // 3 MatrixMarket files, one 10X file.

        let formats = phandle.open("format", { load: true }).values;
        expect(formats).toEqual(["MatrixMarket", "10X"]); 
    }

    let offsets = utils.mockOffsets(collected.collected);
    let reloaded = await bakana.loadAnalysis(
        path, 
        (offset, size) => offsets[offset]
    );

    let new_params = reloaded.parameters;
    expect(new_params.inputs.sample_factor).toBeNull();
    expect(new_params.quality_control instanceof Object).toBe(true);
    expect(new_params.pca instanceof Object).toBe(true);

    // Checking that the permutation is unchanged on reload.
    let old_ids = state.inputs.summary()["genes"]["RNA"]["id"];
    let new_ids = reloaded.state.inputs.summary()["genes"]["RNA"]["id"];
    expect(old_ids.length).toBeGreaterThan(0);
    expect(old_ids).toEqual(new_ids);

    let old_res = state.feature_selection.summary();
    let new_res = reloaded.state.feature_selection.summary();
    expect("means" in old_res).toBe(true);
    expect(old_res["means"]).toEqual(new_res["means"]);

    {
        // Check that the PCs are correctly recovered.
        let corrected_pcs = reloaded.state.batch_correction.fetchPCs();
        expect(corrected_pcs.pcs.owner).toBeNull();
        let original_pcs = reloaded.state.pca.fetchPCs();
        expect(corrected_pcs.length).toEqual(original_pcs.length);
        expect(corrected_pcs.pcs.slice()).not.toEqual(original_pcs.pcs.slice());
    }

    // Freeing.
    await bakana.freeAnalysis(state);
    await bakana.freeAnalysis(reloaded.state);
})

test("single-matrix multi-sample analyses work correctly", async () => {
    let contents = {};
    let finished = (step, res) => {
        contents[step] = res;
    };
    
    let paramcopy = utils.baseParams();
    paramcopy.inputs = {
        sample_factor: "3k"
    };
    let state = await bakana.createAnalysis();
    let res = await bakana.runAnalysis(
        state,
        { 
            "combined": {
                format: "MatrixMarket",
                mtx: "files/datasets/pbmc-combined-matrix.mtx.gz",
                genes: "files/datasets/pbmc-combined-features.tsv.gz",
                annotations: "files/datasets/pbmc-combined-barcodes.tsv.gz"
            }
        },
        paramcopy,
        {
            finishFun: finished,
        }
    );

    expect(contents.quality_control instanceof Object).toBe(true);
    expect("3k" in contents.quality_control.thresholds).toBe(true);
    expect("4k" in contents.quality_control.thresholds).toBe(true);

    // Saving and loading.
    const path = "TEST_state_multi-sample.h5";
    let collected = await bakana.saveAnalysis(state, path);
    utils.validateState(path);
    expect(collected.collected.length).toBe(3);
    expect(typeof(collected.collected[0])).toBe("string");
    
    let offsets = utils.mockOffsets(collected.collected);
    let reloaded = await bakana.loadAnalysis(
        path, 
        (offset, size) => offsets[offset]
    );

    let new_params = reloaded.parameters;
    expect(new_params.inputs.sample_factor).toBe("3k");
    expect(new_params.quality_control instanceof Object).toBe(true);
    expect(new_params.pca instanceof Object).toBe(true);

    // Freeing.
    await bakana.freeAnalysis(state);
    await bakana.freeAnalysis(reloaded.state);
})
