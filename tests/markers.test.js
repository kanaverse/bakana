import * as bakana from "../src/index.js";
import * as scran from "scran.js";
import * as utils from "./utils.js"

beforeAll(async () => await bakana.initialize({ localFile: true }));
afterAll(async () => await bakana.terminate());

let files = {
    default: {
        format: "10X",
        h5: "files/datasets/pbmc4k-tenx.h5"
    }
};

test("fetching of clustered markers works correctly", async () => {
    let state = await bakana.createAnalysis();
    let params = utils.baseParams();
    await bakana.runAnalysis(state, files, params);

    // Markers.
    expect(state.marker_detection.numberOfGroups()).toBeGreaterThan(0);

    let res = state.marker_detection.fetchGroupResults(0, "cohen-mean");
    expect("ordering" in res).toBe(true);
    expect("means" in res).toBe(true);
    expect("lfc" in res).toBe(true);

    bakana.freeAnalysis(state);
})
