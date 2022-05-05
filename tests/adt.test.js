import * as bakana from "../src/index.js";
import * as utils from "./utils.js";
import * as scran from "scran.js";
import * as fs from "fs";

beforeAll(async () => await bakana.initialize({ localFile: true }));
afterAll(async () => await bakana.terminate());

test("runAnalysis works correctly (MatrixMarket)", async () => {
    let state = await bakana.createAnalysis();
    let params = utils.baseParams();
    let mtx = "files/datasets/immune_3.0.0-matrix.mtx.gz";
    let res = await bakana.runAnalysis(state, 
        { 
            default: {
                format: "MatrixMarket",
                mtx: mtx,
                genes: "files/datasets/immune_3.0.0-features.tsv.gz",
                annotations: "files/datasets/immune_3.0.0-barcodes.tsv.gz"
            }
        },
        params
    );

    // Checking that the ADTs were split out.
    let input = state.inputs.summary();
    let cd3 = state.inputs.fetchAnnotations("Antibody Capture: CD3");
    expect(cd3.values.length).toBe(input.dimensions.num_cells);

    let all_types = Array.from(new Set(input.genes.type));
    expect(all_types.length).toBe(1);
    expect(all_types[0]).toBe("Gene Expression");

    // Check that the subsetting was done correctly.
    let f = fs.readFileSync(mtx);
    let buff = f.buffer.slice(f.byteOffset, f.byteOffset + f.byteLength);
    let mat = scran.initializeSparseMatrixFromMatrixMarketBuffer(new Uint8Array(buff), { compressed: true });

    let loaded = state.inputs.fetchCountMatrix();
    expect(mat.numberOfRows()).toBeGreaterThan(loaded.numberOfRows());

    let loaded_ids = loaded.identities();
    let full_ids = mat.identities();
    expect(loaded_ids).toEqual(full_ids.slice(0, loaded_ids.length));

    let last = loaded_ids.length - 1;
    expect(loaded.row(last)).toEqual(mat.row(last));
    mat.free();

    // Saving and loading.
    const path = "TEST_state_adt.h5";
    let collected = await bakana.saveAnalysis(state, path);
    expect(collected.collected.length).toBe(3);
    expect(typeof(collected.collected[0])).toBe("string");

    let offsets = utils.mockOffsets(collected.collected);
    let reloaded = await bakana.loadAnalysis(
        path, 
        (offset, size) => offsets[offset]
    );

    let new_params = reloaded.parameters;

    // Checking that the permutation is unchanged on reload, 
    // even when identities are subsetted.
    let old_ids = state.inputs.summary()["genes"]["id"];
    let new_ids = reloaded.state.inputs.summary()["genes"]["id"];
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

test("runAnalysis works correctly (10X)", async () => {
    let state = await bakana.createAnalysis();
    let params = utils.baseParams();
    let res = await bakana.runAnalysis(state, 
        { 
            default: {
                format: "10X",
                h5: "files/datasets/immune_3.0.0-tenx.h5"
            }
        },
        params
    );

    // Checking that the ADTs were split out.
    let input = state.inputs.summary();
    let cd3 = state.inputs.fetchAnnotations("Antibody Capture: CD3");
    expect(cd3.values.length).toBe(input.dimensions.num_cells);

    let all_types = Array.from(new Set(input.genes.type));
    expect(all_types.length).toBe(1);
    expect(all_types[0]).toBe("Gene Expression");

    // Release me!
    await bakana.freeAnalysis(state);
})
