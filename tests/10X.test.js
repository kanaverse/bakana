import * as bakana from "../src/index.js";
import * as scran from "scran.js";
import * as utils from "./utils.js"

beforeAll(async () => await bakana.initialize({ localFile: true }));
afterAll(async () => await bakana.terminate());

test("runAnalysis works correctly (10X)", async () => {
    let contents = {};
    let finished = (step, res) => {
        contents[step] = res;
    };

    let res = await bakana.runAnalysis(
        {
            default: {
                format: "10X",
                h5: "files/datasets/pbmc4k-tenx.h5"
            }
        },
        utils.baseParams,
        {
            finishFun: finished,
        }
    );

    expect(contents.quality_control instanceof Object).toBe(true);
    expect(contents.pca instanceof Object).toBe(true);
    expect(contents.feature_selection instanceof Object).toBe(true);
    expect(contents.cell_labelling instanceof Object).toBe(true);
    expect(contents.marker_detection instanceof Object).toBe(true);

    // Saving and loading.
    const path = "TEST_state_10X.h5";
    let collected = await bakana.saveAnalysis(path);
    expect(collected.collected.length).toBe(1);
    expect(typeof(collected.collected[0])).toBe("string");
    
    let offsets = utils.mockOffsets(collected.collected);
    let new_params = await bakana.loadAnalysis(
        path, 
        (offset, size) => offsets[offset]
    );

    expect(new_params.quality_control instanceof Object).toBe(true);
    expect(new_params.pca instanceof Object).toBe(true);
})
