import * as bakana from "../src/index.js";
import * as scran from "scran.js";
import * as utils from "./utils.js";

beforeAll(async () => await bakana.initialize({ localFile: true }));
afterAll(async () => await bakana.terminate());

test("runAnalysis works correctly when QC method is 'none'", async () => {
    let state = await bakana.createAnalysis();
    let params = utils.baseParams();
    params.quality_control.method = "none";

    let res = await bakana.runAnalysis(state, 
        { 
            default: {
                format: "MatrixMarket",
                mtx: "files/datasets/pbmc3k-matrix.mtx.gz",
                genes: "files/datasets/pbmc3k-features.tsv.gz",
                annotations: "files/datasets/pbmc3k-barcodes.tsv.gz"
            }
        },
        params
    );

    // Checking that nothing was discarded.
    let disc = state.quality_control.fetchDiscards();
    let lost = 0;
    disc.forEach(x => lost += x);
    expect(lost).toBe(0);
    
    // Saving and loading.
    const path = "TEST_state_no_qc.h5";
    let collected = await bakana.saveAnalysis(state, path);
    expect(collected.collected.length).toBe(3);
    expect(typeof(collected.collected[0])).toBe("string");

    let offsets = utils.mockOffsets(collected.collected);
    let reloaded = await bakana.loadAnalysis(
        path, 
        (offset, size) => offsets[offset]
    );

    let new_params = reloaded.parameters;
    expect(new_params.quality_control.method).toBe("none");

    // Release me!
    await bakana.freeAnalysis(state);
})
