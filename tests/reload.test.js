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
    bakana.freeAnalysis(state);

    let offsets = utils.mockOffsets(collected.collected);
    let reloaded = await bakana.loadAnalysis(
        path, 
        (offset, size) => offsets[offset]
    );

    // Further re-analyzing from the first analysis.
    let new_state = reloaded.state;
    let new_params = reloaded.parameters;
    new_params.quality_control.nmads = 2.5;

    let contents = {};
    let finished = (step, res) => {
        contents[step] = res;
    };
    await bakana.runAnalysis(new_state, null, new_params, { finishFun: finished });

    expect(contents.inputs).toBeUndefined();
    expect(contents.quality_control instanceof Object).toBe(true);

    // Liberation.
    bakana.freeAnalysis(new_state);
})
