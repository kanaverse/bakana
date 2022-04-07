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

test("addition, fetching and removal of custom selections works correctly", async () => {
    let state = await bakana.createAnalysis();
    let params = utils.baseParams();
    await bakana.runAnalysis(state, files, params);

    state.custom_selections.addSelection("evens", [0,2,4,6,8]);
    let res = state.custom_selections.fetchResults("evens", "cohen");
    expect("ordering" in res).toBe(true);
    expect("means" in res).toBe(true);
    expect("lfc" in res).toBe(true);

    state.custom_selections.addSelection("odds", [1,3,5,7,9]);
    let res2 = state.custom_selections.fetchResults("odds", "cohen");
    expect("ordering" in res2).toBe(true);
    expect("means" in res2).toBe(true);
    expect("lfc" in res2).toBe(true);

    state.custom_selections.removeSelection("odds");

    // Saving and loading works correctly.
    const path = "TEST_state_custom.h5";
    let collected = await bakana.saveAnalysis(state, path);

    let offsets = utils.mockOffsets(collected.collected);
    let reloaded = await bakana.loadAnalysis(
        path, 
        (offset, size) => offsets[offset]
    );

    expect(reloaded.parameters.custom_selections.selections.evens.length).toBe(5);
    let reres = reloaded.state.custom_selections.fetchResults("evens", "cohen");
    expect(reres.ordering).toEqual(res.ordering);

    bakana.freeAnalysis(state);
    bakana.freeAnalysis(reloaded.state);
})
