import * as bakana from "../src/index.js";
import * as scran from "scran.js";
import * as utils from "./utils.js";

beforeAll(async () => await bakana.initialize({ localFile: true }));
afterAll(async () => await bakana.terminate());

test("reanalysis from a reloaded analysis works correctly", async () => {
    let files = { 
        default: {
            format: "MatrixMarket",
            mtx: "files/datasets/pbmc3k-matrix.mtx.gz",
            genes: "files/datasets/pbmc3k-features.tsv.gz",
            annotations: "files/datasets/pbmc3k-barcodes.tsv.gz"
        }
    }

    let params = utils.baseParams();
    let state = await bakana.createAnalysis();
    await bakana.runAnalysis(state, files, params);

    // Saving and reloading.
    const path = "TEST_state_reloaded.h5";
    let collected = await bakana.saveAnalysis(state, path);
    await bakana.freeAnalysis(state);

    let offsets = utils.mockOffsets(collected.collected);
    let reloaded = await bakana.loadAnalysis(
        path, 
        (offset, size) => offsets[offset]
    );

    // Re-analyzing without change; should be a no-op.
    let new_state = reloaded.state;
    let new_params = reloaded.parameters;

    let contents = {};
    let finished = (step, res) => {
        if (typeof res !== "undefined") {
            contents[step] = res;
        }
    };
    await bakana.runAnalysis(new_state, null, new_params, { finishFun: finished });
    expect(Object.entries(contents).length).toBe(0);

    await bakana.freeAnalysis(new_state);

    // Re-analyzing with a change.
    let reloaded2 = await bakana.loadAnalysis(
        path, 
        (offset, size) => offsets[offset]
    );
    let new_state2 = reloaded2.state;
    let new_params2 = reloaded2.parameters;
    new_params2.quality_control.nmads = 2.5;

    contents = {};
    await bakana.runAnalysis(new_state2, null, new_params2, { finishFun: finished });

    expect(contents.inputs).toBeUndefined();
    expect(contents.quality_control instanceof Object).toBe(true);

    await bakana.freeAnalysis(new_state2);
})
