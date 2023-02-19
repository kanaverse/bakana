import * as bakana from "../src/index.js";
import * as butils from "../src/steps/utils/general.js";
import * as scran from "scran.js";
import * as utils from "./utils.js";
import * as bioc from "bioconductor";

beforeAll(utils.initializeAll);
afterAll(async () => await bakana.terminate());

let hs_files = { 
    default: new bakana.TenxMatrixMarketDataset(
        "files/datasets/pbmc3k-matrix.mtx.gz",
        "files/datasets/pbmc3k-features.tsv.gz",
        "files/datasets/pbmc3k-barcodes.tsv.gz")
};

test("feature set enrichment works correctly for humans", async () => {
    let params = utils.baseParams();
    params.feature_set_enrichment.dataset_id_column = "id";
    params.feature_set_enrichment.collections = [ "mouse-GO_3.16.0", "human-GO_3.16.0" ];

    let state = await bakana.createAnalysis();
    await bakana.runAnalysis(state, hs_files, params);

    let deets = await state.feature_set_enrichment.fetchCollectionDetails();
    let nmouse = deets["mouse-GO_3.16.0"].names.length;
    expect(nmouse).toBeGreaterThan(0);
    expect(deets["mouse-GO_3.16.0"].universe).toBe(0);

    let nhuman = deets["human-GO_3.16.0"].names.length;
    expect(nhuman).toBeGreaterThan(0);
    expect(deets["human-GO_3.16.0"].universe).toBeGreaterThan(0);

    let ngroups = state.marker_detection.fetchResults().RNA.numberOfGroups();
    let clust = ngroups - 1;

    {
        let stats = await state.feature_set_enrichment.fetchGroupResults(clust, "cohen", "mean");
        expect(stats["human-GO_3.16.0"].counts.length).toBeGreaterThan(0);
        expect(stats["human-GO_3.16.0"].num_markers).toBeLessThan(100 * 1.5); // a bit of leeway for ties.
        expect(stats["human-GO_3.16.0"].counts).not.toEqual(new Int32Array(nhuman));
        expect(stats["mouse-GO_3.16.0"].counts).toEqual(new Int32Array(nmouse));
    }

    {
        let stats = await state.feature_set_enrichment.fetchGroupResults(clust, "cohen", "min_rank");
        expect(stats["human-GO_3.16.0"].counts.length).toBeGreaterThan(0);
        expect(stats["human-GO_3.16.0"].num_markers).toBeLessThan(100 * 1.5); // a bit of leeway for ties.
        expect(stats["human-GO_3.16.0"].counts).not.toEqual(new Int32Array(nhuman));
        expect(stats["mouse-GO_3.16.0"].counts).toEqual(new Int32Array(nmouse));
    }

    {
        let stats = await state.feature_set_enrichment.fetchGroupResults(clust, "cohen", "min_rank");
        expect(stats["human-GO_3.16.0"].counts.length).toBeGreaterThan(0);
        expect(stats["human-GO_3.16.0"].num_markers).toBeLessThan(100 * 1.5); // a bit of leeway for ties.
        expect(stats["human-GO_3.16.0"].counts).not.toEqual(new Int32Array(nhuman));
        expect(stats["mouse-GO_3.16.0"].counts).toEqual(new Int32Array(nmouse));
    }

    // Check that the scores work correctly.
    {
        let chosen = -1;
        let human_sizes = deets["human-GO_3.16.0"].sizes;
        for (var i = 0; i < human_sizes.length; i++) {
            if (human_sizes[i] > 20) {
                chosen = i;
                break;
            }
        }

        expect(chosen).toBeGreaterThanOrEqual(0);
        let output = await state.feature_set_enrichment.fetchPerCellScores("human-GO_3.16.0", chosen);
        expect(output.weights.length).toEqual(human_sizes[chosen]);
        expect(output.scores.length).toEqual(state.cell_filtering.fetchFilteredMatrix().numberOfColumns());
    }


    // Release me!
    await bakana.freeAnalysis(state);
})

let mm_files = { 
    default: new bakana.H5adDataset("files/datasets/zeisel-brain.h5ad")
};

test("feature set enrichment works correctly for mice", async () => {
    let params = utils.baseParams();
    params.feature_set_enrichment.collections = [ "mouse-GO_3.16.0", "human-GO_3.16.0" ];
    params.feature_set_enrichment.dataset_id_column = null;
    params.feature_set_enrichment.reference_id_column = "SYMBOL";

    let state = await bakana.createAnalysis();
    await bakana.runAnalysis(state, mm_files, params);

    let deets = await state.feature_set_enrichment.fetchCollectionDetails();
    let nmouse = deets["mouse-GO_3.16.0"].names.length;
    expect(nmouse).toBeGreaterThan(0);
    expect(deets["mouse-GO_3.16.0"].universe).toBeGreaterThan(0);

    {
        let stats = await state.feature_set_enrichment.fetchGroupResults(0, "cohen", "mean");
        expect(stats["mouse-GO_3.16.0"].counts.length).toBeGreaterThan(0);
        expect(stats["mouse-GO_3.16.0"].num_markers).toBeLessThan(100 * 1.5); // a bit of leeway for ties.
        expect(stats["mouse-GO_3.16.0"].counts).not.toEqual(new Int32Array(nmouse));
    }

    // Check that the scores work correctly.
    {
        let chosen = -1;
        let mouse_sizes = deets["mouse-GO_3.16.0"].sizes;
        for (var i = 0; i < mouse_sizes.length; i++) {
            if (mouse_sizes[i] > 5) {
                chosen = i;
                break;
            }
        }

        expect(chosen).toBeGreaterThanOrEqual(0);
        let output = await state.feature_set_enrichment.fetchPerCellScores("mouse-GO_3.16.0", chosen);
        expect(output.weights.length).toEqual(mouse_sizes[chosen]);
        expect(output.scores.length).toEqual(state.cell_filtering.fetchFilteredMatrix().numberOfColumns());
    }

    // Release me!
    await bakana.freeAnalysis(state);
})
