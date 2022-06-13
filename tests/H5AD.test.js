import * as bakana from "../src/index.js";
import * as scran from "scran.js";
import * as utils from "./utils.js";

beforeAll(async () => await bakana.initialize({ localFile: true }));
afterAll(async () => await bakana.terminate());

test("runAnalysis works correctly (H5AD)", async () => {
    let contents = {};
    let finished = (step, res) => {
        contents[step] = res;
    };

    let state = await bakana.createAnalysis();
    let params = utils.baseParams();
    let res = await bakana.runAnalysis(
        state,
        { 
            default: {
                format: "H5AD",
                h5: "files/datasets/zeisel-brain.h5ad"
            }
        },
        params,
        {
            finishFun: finished,
        }
    );

    expect(contents.quality_control instanceof Object).toBe(true);
    expect(contents.pca instanceof Object).toBe(true);
    expect(contents.feature_selection instanceof Object).toBe(true);
    expect(contents.cell_labelling instanceof Object).toBe(true);
    expect(contents.marker_detection instanceof Object).toBe(true);

    // Computations are done correctly.
    {
        // Markers.
        expect(state.marker_detection.numberOfGroups()).toBeGreaterThan(0);

        let res = state.marker_detection.fetchGroupResults(0, "cohen-mean", "RNA");
        expect("ordering" in res).toBe(true);
        expect("means" in res).toBe(true);
        expect("lfc" in res).toBe(true);

        // Normalized expression.
        let exprs = state.normalization.fetchExpression(0);
        let nfiltered = state.cell_filtering.fetchFilteredMatrix().numberOfColumns();
        expect(exprs.length).toBe(nfiltered);

        // Factor annotations, with and without filtering.
        let cell_anno = state.inputs.fetchAnnotations("level1class");
        expect(cell_anno.length).toBe(state.inputs.fetchCountMatrix().numberOfColumns());
        let filtered_anno = state.cell_filtering.fetchFilteredAnnotations("level1class");
        expect(filtered_anno.length).toBe(nfiltered);

        // Non-factor annotations, with and without filtering.
        let sex_anno = state.inputs.fetchAnnotations("sex");
        expect(sex_anno.length).toBe(state.inputs.fetchCountMatrix().numberOfColumns());
        let filtered_sex = state.cell_filtering.fetchFilteredAnnotations("sex");
        expect(filtered_sex.length).toBe(nfiltered);
    }

    // Saving and loading.
    const path = "TEST_state_H5AD.h5";
    let collected = await bakana.saveAnalysis(state, path);
    expect(collected.collected.length).toBe(1);
    expect(typeof(collected.collected[0])).toBe("string");
    
    let offsets = utils.mockOffsets(collected.collected);
    let reloaded = await bakana.loadAnalysis(
        path, 
        (offset, size) => offsets[offset]
    );

    let new_params = reloaded.parameters;
    expect(new_params.quality_control instanceof Object).toBe(true);
    expect(new_params.pca instanceof Object).toBe(true);

    // Freeing.
    await bakana.freeAnalysis(state);
})

