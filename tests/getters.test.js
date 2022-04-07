import * as bakana from "../src/index.js";
import * as scran from "scran.js";
import * as utils from "./utils.js"

beforeAll(async () => await bakana.initialize({ localFile: true }));
afterAll(async () => await bakana.terminate());

let files = {
    default: {
        format: "H5AD",
        h5: "files/datasets/zeisel-brain.h5ad"
    }
};

test("fetching of various bits and pieces is done correctly", async () => {
    let state = await bakana.createAnalysis();
    let params = utils.baseParams();
    await bakana.runAnalysis(state, files, params);

    // Markers.
    expect(state.marker_detection.numberOfGroups()).toBeGreaterThan(0);

    let res = state.marker_detection.fetchGroupResults(0, "cohen-mean");
    expect("ordering" in res).toBe(true);
    expect("means" in res).toBe(true);
    expect("lfc" in res).toBe(true);

    // Normalized expression.
    let exprs = state.normalization.fetchExpression(0);
    let nfiltered = state.quality_control.fetchFilteredMatrix().numberOfColumns();
    expect(exprs.length).toBe(nfiltered);

    // Annotations, with and without filtering.
    let cell_anno = state.inputs.fetchAnnotations("level1class");
    expect(cell_anno.index.length).toBe(state.inputs.fetchCountMatrix().numberOfColumns());
    let filtered_anno = state.quality_control.fetchFilteredAnnotations("level1class");
    expect(filtered_anno.index.length).toBe(nfiltered);

    await bakana.freeAnalysis(state);
})
