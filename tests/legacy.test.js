import * as bakana from "../src/index.js";
import * as scran from "scran.js";
import * as utils from "./utils.js";

beforeAll(async () => await bakana.initialize({ localFile: true }));
afterAll(async () => await bakana.terminate());

test("reanalysis from a v0 analysis works correctly (10X)", async () => {
    const h5path = "TEST_v0_state.h5";
    let loader = await bakana.parseKanaFile("files/legacy/zeisel_tenx_20220307.kana", h5path);
    let reloaded = await bakana.loadAnalysis(h5path, loader);

    let new_state = reloaded.state;
    let new_params = reloaded.parameters;

    // Missing steps are filled in.
    expect("cell_labelling" in new_params).toBe(true);
    expect("kmeans_cluster" in new_params).toBe(true);
})

test("reanalysis from a v0 analysis works correctly (MatrixMarket)", async () => {
    const h5path = "TEST_v0_state.h5";
    let loader = await bakana.parseKanaFile("files/legacy/zeisel_mtx_20220306.kana", h5path);
    let reloaded = await bakana.loadAnalysis(h5path, loader);

    let new_state = reloaded.state;
    let new_params = reloaded.parameters;

    // Missing steps are filled in.
    expect("cell_labelling" in new_params).toBe(true);
    expect("kmeans_cluster" in new_params).toBe(true);
})
