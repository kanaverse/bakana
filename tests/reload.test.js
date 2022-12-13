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
    await bakana.freeAnalysis(state);

    let offsets = utils.mockOffsets(collected.collected);
    let reloaded = await bakana.loadAnalysis(
        path, 
        (offset, size) => offsets[offset]
    );

    // Re-analyzing without change; should be a no-op.
    let new_params = bakana.retrieveParameters(reloaded);
    await bakana.runAnalysis(reloaded, null, new_params);

    for (const x of Object.values(reloaded)) {
        expect(x.changed).toBe(false);
    }

    // Animations still work after reloading.
    await utils.triggerAnimation(reloaded);

    await bakana.freeAnalysis(reloaded);

    // Re-analyzing with a change.
    {
        let reloaded2 = await bakana.loadAnalysis(
            path, 
            (offset, size) => offsets[offset]
        );
        let new_params2 = bakana.retrieveParameters(reloaded2);
        new_params2.quality_control.nmads = 2.5;

        await bakana.runAnalysis(reloaded2, null, new_params2);

        expect(reloaded2.inputs.changed).toBe(false);
        expect(reloaded2.quality_control.changed).toBe(true);

        await bakana.freeAnalysis(reloaded2);
    }
})
