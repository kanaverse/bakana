import * as bakana from "../src/index.js";
import * as scran from "scran.js";
import * as utils from "./utils.js"

beforeAll(utils.initializeAll);
afterAll(async () => await bakana.terminate());

let files = {
    default: new bakana.TenxHdf5Dataset("files/datasets/pbmc4k-tenx.h5")
};

test("addition, fetching and removal of custom selections works correctly", async () => {
    let state = await bakana.createAnalysis();
    let params = utils.baseParams();
    await bakana.runAnalysis(state, files, params);

    state.custom_selections.addSelection("evens", [0,2,4,6,8]);
    let res = state.custom_selections.fetchResults("evens");
    expect(res.RNA instanceof scran.ScoreMarkersResults).toBe(true);

    state.custom_selections.addSelection("odds", [1,3,5,7,9]);
    let res2 = state.custom_selections.fetchResults("odds");
    expect(res2.RNA instanceof scran.ScoreMarkersResults).toBe(true);

    expect(state.custom_selections.fetchSelectionIndices("evens")).toEqual([0,2,4,6,8]);
    expect(state.custom_selections.fetchSelectionIndices("odds")).toEqual([1,3,5,7,9]);
    {
        let all = state.custom_selections.fetchSelections();
        expect(all["evens"]).toEqual([0,2,4,6,8]);
        expect(all["odds"]).toEqual([1,3,5,7,9]);
    }

    // Versus mode works correctly.
    {
        let vres = state.custom_selections.computeVersus("odds", "evens");
        expect(vres.results.RNA instanceof scran.ScoreMarkersResults).toBe(true);

        let vres2 = state.custom_selections.computeVersus("evens", "odds");
        expect(vres2.results.RNA instanceof scran.ScoreMarkersResults).toBe(true);

        let lfcs = vres.results.RNA.lfc(vres.left);
        let lfcs2 = vres2.results.RNA.lfc(vres2.left);
        lfcs2.forEach((x, i) => {
            lfcs2[i] *= -1;
        });

        expect(lfcs).toEqual(lfcs2);
    }

    // Removal of selections works correctly.
    state.custom_selections.removeSelection("odds");
    {
        let all = state.custom_selections.fetchSelections();
        expect("odds" in all).toBe(false)
        expect("evens" in all).toBe(true)
    }

    // Saving and loading works correctly.
    const path = "TEST_state_custom.h5";
    let collected = await bakana.saveAnalysis(state, path);
    utils.validateState(path);

    let offsets = utils.mockOffsets(collected.collected);
    let reloaded = await bakana.loadAnalysis(
        path, 
        (offset, size) => offsets[offset]
    );

    let new_params = bakana.retrieveParameters(reloaded);
    expect(new_params.custom_selections).toEqual({});
    expect(state.custom_selections.fetchSelectionIndices("evens")).toEqual([0,2,4,6,8]);

    let reres = reloaded.custom_selections.fetchResults("evens");
    expect(reres.RNA.cohen(1)).toEqual(res.RNA.cohen(1));
    expect(reres.RNA.auc(1)).toEqual(res.RNA.auc(1));
    expect(reres.RNA.means(1)).toEqual(res.RNA.means(1));

    await bakana.freeAnalysis(state);
    await bakana.freeAnalysis(reloaded);
})
