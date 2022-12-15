import * as bakana from "../src/index.js";
import * as scran from "scran.js";
import * as utils from "./utils.js";
import * as bioc from "bioconductor";

beforeAll(utils.initializeAll);
afterAll(async () => await bakana.terminate());

let h5path = "files/datasets/pbmc4k-tenx.h5";
let files = {
    default: new bakana.TenxHdf5Dataset(h5path)
};

test("10X summary works correctly", async () => {
    let summ = files.default.summary();
    expect(summ.all_features instanceof bioc.DataFrame).toBe(true);
    expect(summ.cells instanceof bioc.DataFrame).toBe(true);
})

test("runAnalysis works correctly (10X)", async () => {
    let state = await bakana.createAnalysis();
    let params = utils.baseParams();
    let res = await bakana.runAnalysis(state, files, params);

    // Input reorganization is done correctly.
    {
        let loaded = state.inputs.fetchCountMatrix().get("RNA");
        let loaded_ids = state.inputs.fetchRowIds()["RNA"];
        let loaded_names = state.inputs.fetchFeatureAnnotations()["RNA"].column("id");

        let simple = scran.initializeSparseMatrixFromHDF5(h5path, "matrix", { layered: false });
        let simple_names = (new scran.H5File(h5path)).open("matrix").open("features").open("id", { load: true }).values;

        utils.checkReorganization(simple.matrix, simple.row_ids, simple_names, loaded, loaded_ids, loaded_names);
        simple.matrix.free();
    }

    // Basic consistency checks.
    await utils.checkStateResultsSimple(state);
    utils.checkClusterVersusMode(state);
    await utils.triggerAnimation(state);

    // Saving and loading.
    const path = "TEST_state_10X.h5";
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

    // Freeing.
    await bakana.freeAnalysis(state);
    await bakana.freeAnalysis(reloaded);
})
