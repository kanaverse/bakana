import * as bakana from "../src/index.js";
import * as scran from "scran.js";
import * as utils from "./utils.js"

beforeAll(async () => { await scran.initialize({ localFile: true }) });
afterAll(async () => { await scran.terminate() });

test("runAnalysis works correctly (MatrixMarket)", async () => {
    let contents = {};
    let finished = (res, step, msg) => {
        contents[step] = res;
    };

    let res = await bakana.runAnalysis(
        [
            {
                format: "mtx",
                mtx: [ "files/datasets/pbmc3k-matrix.mtx.gz" ],
                gene: [ "files/datasets/pbmc3k-features.tsv.gz" ],
                barcode: [ "files/datasets/pbmc3k-barcodes.tsv.gz" ]
            }
        ],
        utils.baseParams,
        {
            finished: finished,
            download: utils.downloadReference
        }
    );

    expect(contents.quality_control instanceof Object).toBe(true);
    expect(contents.pca instanceof Object).toBe(true);
    expect(contents.feature_selection instanceof Object).toBe(true);
    expect(contents.cell_labelling instanceof Object).toBe(true);
    expect(contents.marker_detection instanceof Object).toBe(true);

    // Saving and loading.
    const path = "TEST_state_MatrixMarket.h5";
    let collected = await bakana.saveAnalysis(path);
    expect(collected.collected.length).toBe(3);
    expect(typeof(collected.collected[0])).toBe("string");
    
    let offsets = utils.mockOffsets(collected.collected);
    let new_params = await bakana.loadAnalysis(
        path, 
        true,
        {
            finished: (x, y, z) => null,
            loader: (offset, size) => offsets[offset]
        }
    );

    expect(new_params.quality_control instanceof Object).toBe(true);
    expect(new_params.pca instanceof Object).toBe(true);
})
