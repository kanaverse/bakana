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
    await bakana.runAnalysis(files, utils.baseParams);

    // Saving.
    const path = "TEST_state_reloaded.h5";
    let collected = await bakana.saveAnalysis(path);

    // Doing a completely different analysis.
    await bakana.runAnalysis(
        { 
            default: {
                format: "10X",
                h5: "files/datasets/pbmc4k-tenx.h5"
            }
        },
        utils.baseParams
    );

    // Now reloading from the first analysis. 
    let offsets = utils.mockOffsets(collected.collected);
    let new_params = await bakana.loadAnalysis(
        path, 
        (offset, size) => offsets[offset]
    );

    // Further re-analyzing from the first analysis.
    let paramcopy = bakana.analysisDefaults();
    paramcopy.quality_control.nmads = 2.5;

    let contents = {};
    let finished = (step, res) => {
        contents[step] = res;
    };
    await bakana.runAnalysis({ default: { format: "kana" } }, paramcopy, { finishFun: finished });

    expect(contents.inputs).toBeUndefined();
    expect(contents.quality_control instanceof Object).toBe(true);
})
