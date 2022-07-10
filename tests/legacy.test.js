import * as bakana from "../src/index.js";
import * as scran from "scran.js";
import * as utils from "./utils.js";
import * as valkana from "valkana";

beforeAll(utils.initializeAll);
afterAll(async () => await bakana.terminate());

function createFinisher(contents) {
    return (step, res) => {
        if (typeof res !== "undefined") {
            contents[step] = res;
        }
    };
}

test("reanalysis from a v0 analysis works correctly (10X)", async () => {
    const h5path = "TEST_v0_state.h5";
    let loader = await bakana.parseKanaFile("files/legacy/zeisel_tenx_20220307.kana", h5path);
    // valkana.validateState(h5path, true, 0); // Gave up, doesn't look like a valid v0 file, actually.

    let contents = {};
    let reloaded = await bakana.loadAnalysis(h5path, loader, { finishFun: createFinisher(contents) });
    let new_params = bakana.retrieveParameters(reloaded);

    // Checking that some of the summaries are correctly reported.
    expect(contents.cell_filtering.retained).toBeGreaterThan(0);
    expect(contents.combine_embeddings).toEqual({});

    // Missing steps are filled in.
    expect("cell_labelling" in new_params).toBe(true);
    expect("kmeans_cluster" in new_params).toBe(true);

    // Cleaning up.
    await bakana.freeAnalysis(reloaded);
})

test("reanalysis from a v0 analysis works correctly (MatrixMarket)", async () => {
    const h5path = "TEST_v0_state.h5";
    let loader = await bakana.parseKanaFile("files/legacy/zeisel_mtx_20220306.kana", h5path);
    // valkana.validateState(h5path, true, 0); // Gave up, doesn't look like a valid v0 file, actually.

    let contents = {};
    let reloaded = await bakana.loadAnalysis(h5path, loader, { finishFun: createFinisher(contents) });
    let new_params = bakana.retrieveParameters(reloaded);

    // Checking that some of the summaries are correctly reported.
    expect(contents.cell_filtering.retained).toBeGreaterThan(0);
    expect(contents.combine_embeddings).toEqual({});

    // Missing steps are filled in.
    expect("cell_labelling" in new_params).toBe(true);
    expect("kmeans_cluster" in new_params).toBe(true);

    // Cleaning up.
    await bakana.freeAnalysis(reloaded);
})

test("reanalysis from a v1.0 analysis works correctly (10X)", async () => {
    const h5path = "TEST_v1.0_state.h5";
    let loader = await bakana.parseKanaFile("files/legacy/zeisel_tenx_20220318.kana", h5path);
    // valkana.validateState(h5path, true, 1000000); // Gave up, doesn't look like a valid v1.0 file, actually.

    let contents = {};
    let reloaded = await bakana.loadAnalysis(h5path, loader, { finishFun: createFinisher(contents) });
    let new_params = bakana.retrieveParameters(reloaded);

    // Checking that some of the summaries are correctly reported.
    expect(contents.cell_filtering.retained).toBeGreaterThan(0);
    expect(contents.combine_embeddings).toEqual({});

    // Missing steps are filled in.
    expect("cell_labelling" in new_params).toBe(true);
    expect("kmeans_cluster" in new_params).toBe(true);

    expect(new_params.snn_graph_cluster.scheme).toBe("rank"); // recover something from mis-formatted files.

    // Cleaning up.
    await bakana.freeAnalysis(reloaded);
})

test("reanalysis from a v1.1 analysis works correctly (10X combined)", async () => {
    const h5path = "TEST_v1.1_state.h5";
    let loader = await bakana.parseKanaFile("files/legacy/pbmc-combined_tenx_20220401.kana", h5path);
    valkana.validateState(h5path, true, 1001000);

    let contents = {};
    let reloaded = await bakana.loadAnalysis(h5path, loader, { finishFun: createFinisher(contents) });
    let new_params = bakana.retrieveParameters(reloaded);

    // Checking that some of the summaries are correctly reported.
    expect(contents.cell_filtering.retained).toBeGreaterThan(0);
    expect(contents.combine_embeddings).toEqual({});

    // Missing steps are filled in.
    expect("adt_normalization" in new_params).toBe(true);

    // Cleaning up.
    await bakana.freeAnalysis(reloaded);
})

test("reanalysis from a v1.1 analysis works correctly (MatrixMarket)", async () => {
    const h5path = "TEST_v1.1_state.h5";
    let loader = await bakana.parseKanaFile("files/legacy/pbmc4k-with-custom_mtx_20220408.kana", h5path);
    valkana.validateState(h5path, true, 1001000);

    let contents = {};
    let reloaded = await bakana.loadAnalysis(h5path, loader, { finishFun: createFinisher(contents) });
    let new_params = bakana.retrieveParameters(reloaded);

    // Checking that some of the summaries are correctly reported.
    expect(contents.cell_filtering.retained).toBeGreaterThan(0);
    expect(contents.combine_embeddings).toEqual({});

    // Missing steps are filled in.
    expect("adt_normalization" in new_params).toBe(true);

    // Cleaning up.
    await bakana.freeAnalysis(reloaded);
})

test("reanalysis from a v1.2 analysis works correctly (10X combined)", async () => {
    const h5path = "TEST_v1.2_state.h5";
    let loader = await bakana.parseKanaFile("files/legacy/pbmc-combined-with-kmeans-and-custom_tenx_20220525.kana", h5path);
    // valkana.validateState(h5path, true, 1002000); // Gave up, doesn't look like a valid 1.2 file, actually.

    let contents = {};
    let reloaded = await bakana.loadAnalysis(h5path, loader, { finishFun: createFinisher(contents) });
    let new_params = bakana.retrieveParameters(reloaded);

    // Checking that some of the summaries are correctly reported.
    expect(contents.cell_filtering.retained).toBeGreaterThan(0);
    expect(contents.combine_embeddings).toEqual({});

    // Missing steps are filled in.
    expect("adt_normalization" in new_params).toBe(true);

    // Cleaning up.
    await bakana.freeAnalysis(reloaded);
})

test("reanalysis from a v1.2 analysis works correctly (MatrixMarket)", async () => {
    const h5path = "TEST_v1.2_state.h5";
    let loader = await bakana.parseKanaFile("files/legacy/pbmc4k-with-kmeans-and-custom_mtx_20220525.kana", h5path);
    // valkana.validateState(h5path, true, 1002000); // Gave up, doesn't look like a valid 1.2 file, actually.

    let contents = {};
    let reloaded = await bakana.loadAnalysis(h5path, loader, { finishFun: createFinisher(contents) });
    let new_params = bakana.retrieveParameters(reloaded);

    // Checking that some of the summaries are correctly reported.
    expect(contents.cell_filtering.retained).toBeGreaterThan(0);
    expect(contents.combine_embeddings).toEqual({});

    // Missing steps are filled in.
    expect("adt_normalization" in new_params).toBe(true);

    // Cleaning up.
    await bakana.freeAnalysis(reloaded);
})
