import * as bakana from "../src/index.js";
import * as scran from "scran.js";
import * as utils from "./utils.js";

beforeAll(async () => await bakana.initialize({ localFile: true }));
afterAll(async () => await bakana.terminate());

test("reanalysis from a reloaded analysis works correctly", async () => {
    let files = { 
        default: new bakana.TenxMatrixMarketDataset(
                "files/datasets/pbmc3k-matrix.mtx.gz",
                "files/datasets/pbmc3k-features.tsv.gz",
                "files/datasets/pbmc3k-barcodes.tsv.gz"
            )
    }

    let params = utils.baseParams();
    let state = await bakana.createAnalysis();
    await bakana.runAnalysis(state, files, params);

    // Saving and reloading.
    const path = "TEST_state_reloaded.h5";
    let collected = await bakana.saveAnalysis(state, path);

    let offsets = utils.mockOffsets(collected.collected);
    let reloaded = await bakana.loadAnalysis(
        path, 
        (offset, size) => offsets[offset]
    );

    await utils.compareStates(state, reloaded);

    // Any runAnalysis should trigger a complete re-analysis.
    let new_params = bakana.retrieveParameters(reloaded);
    await bakana.runAnalysis(reloaded, null, new_params);
    expect(reloaded.inputs.changed).toBe(true);
    expect(reloaded.rna_quality_control.changed).toBe(true);
    expect(reloaded.snn_graph_cluster.changed).toBe(true);
    await utils.compareStates(state, reloaded);

    // Any further runs should still be a no-op, though.
    await bakana.runAnalysis(reloaded, null, new_params);
    for (const x of Object.values(reloaded)) {
        expect(x.changed).toBe(false);
    }

    // Freeing everything.
    await bakana.freeAnalysis(reloaded);
    await bakana.freeAnalysis(state);
})
