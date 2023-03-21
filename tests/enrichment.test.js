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
    params.feature_set_enrichment.skip = false;
    let state = await bakana.createAnalysis();
    await bakana.runAnalysis(state, hs_files, params);

    let collections = state.feature_set_enrichment.fetchCollectionDetails();
    expect(collections.names.length).toBeGreaterThan(0);
    expect(Array.from(new Set(collections.species))).toEqual(["9606"]);

    let sets = state.feature_set_enrichment.fetchSetDetails();
    expect(sets.names.length).toBeGreaterThan(0);
    expect(sets.sizes instanceof Int32Array).toBe(true);

    expect(state.feature_set_enrichment.fetchUniverseSize()).toBeGreaterThan(10000);

    let ngroups = state.marker_detection.fetchResults()["RNA"].numberOfGroups();
    let clust = ngroups - 1;
    let ntop = params.feature_set_enrichment.top_markers;

    {
        let stats = state.feature_set_enrichment.computeEnrichment(state.marker_detection.fetchResults()["RNA"], clust, "cohen", "mean");
        expect(stats.counts.length).toBeGreaterThan(0);
        expect(stats.set_ids.every(x => x >= 0 && x < sets.names.length)).toBe(true);
        expect(stats.num_markers).toBeGreaterThan(0); 
        expect(stats.num_markers).toBeLessThan(ntop * 1.5); // a bit of leeway for ties.
        expect(stats.counts.every(x => x > 0)).toBe(true); // all counts should be positive.
        expect(stats.pvalues.some(Number.isNaN)).toBe(false); // make sure we computed a valid p-value
    }

    {
        // Goes in the right direction with min-rank.
        let stats = state.feature_set_enrichment.computeEnrichment(state.marker_detection.fetchResults()["RNA"], clust, "cohen", "min_rank");
        expect(stats.counts.length).toBeGreaterThan(0);
        expect(stats.num_markers).toBeGreaterThan(0); 
        expect(stats.num_markers).toBeLessThan(ntop * 1.5); 
        expect(stats.counts.every(x => x > 0)).toBe(true);
        expect(stats.pvalues.some(Number.isNaN)).toBe(false);
    }

    {
        // Also works correctly for AUC.
        let stats = state.feature_set_enrichment.computeEnrichment(state.marker_detection.fetchResults()["RNA"], clust, "auc", "mean");
        expect(stats.counts.length).toBeGreaterThan(0);
        expect(stats.num_markers).toBeGreaterThan(0); 
        expect(stats.num_markers).toBeLessThan(ntop * 1.5);
        expect(stats.counts.every(x => x > 0)).toBe(true);
        expect(stats.pvalues.some(Number.isNaN)).toBe(false);
    }

    // Check that the scores work correctly.
    {
        let chosen = -1;
        for (var i = 0; i < sets.sizes.length; i++) {
            if (sets.sizes[i] > 20) {
                chosen = i;
                break;
            }
        }

        expect(chosen).toBeGreaterThanOrEqual(0);
        let output = state.feature_set_enrichment.computePerCellScores(chosen);
        expect(output.weights.length).toEqual(sets.sizes[chosen]);
        expect(output.scores.length).toEqual(state.cell_filtering.fetchFilteredMatrix().numberOfColumns());
    }

    // Forcing us to load the mice.
    params.feature_set_enrichment.automatic = false;
    params.feature_set_enrichment.species = ["10090"];
    params.feature_set_enrichment.gene_id_column = "id";

    {
        await bakana.runAnalysis(state, hs_files, params);
        expect(state.marker_detection.changed).toBe(false);
        expect(state.feature_set_enrichment.changed).toBe(true);

        let collections = state.feature_set_enrichment.fetchCollectionDetails();
        expect(collections.names.length).toBeGreaterThan(0);
        expect(Array.from(new Set(collections.species))).toEqual(["10090"]);
        expect(state.feature_set_enrichment.fetchUniverseSize()).toBe(0);

        let stats = state.feature_set_enrichment.computeEnrichment(state.marker_detection.fetchResults()["RNA"], 0, "cohen", "mean");
        expect(stats.counts.length).toBe(0);
    }

    // Checking that reloading works.
    params.feature_set_enrichment.automatic = true;
    {
        await bakana.runAnalysis(state, hs_files, params);
        expect(state.marker_detection.changed).toBe(false);
        expect(state.feature_set_enrichment.changed).toBe(true);

        let collections = state.feature_set_enrichment.fetchCollectionDetails();
        expect(Array.from(new Set(collections.species))).toEqual(["9606"]);
    }

    params.feature_set_enrichment.species = ["9606"];
    {
        await bakana.runAnalysis(state, hs_files, params);
        expect(state.marker_detection.changed).toBe(false);
        expect(state.feature_set_enrichment.changed).toBe(false); // no effect when automatic=true.
    }

    // Release me!
    await bakana.freeAnalysis(state);
})

let mm_files = { 
    default: new bakana.H5adDataset("files/datasets/zeisel-brain.h5ad")
};

test("feature set enrichment works correctly for mice", async () => {
    let params = utils.baseParams();
    params.feature_set_enrichment.skip = false;
    let state = await bakana.createAnalysis();
    await bakana.runAnalysis(state, mm_files, params);

    let collections = state.feature_set_enrichment.fetchCollectionDetails();
    expect(collections.names.length).toBeGreaterThan(0);
    expect(Array.from(new Set(collections.species))).toEqual(["10090"]);
    expect(state.feature_set_enrichment.fetchUniverseSize()).toBeGreaterThan(10000);
    
    let sets = state.feature_set_enrichment.fetchSetDetails();
    let ntop = params.feature_set_enrichment.top_markers;

    {
        let stats = state.feature_set_enrichment.computeEnrichment(state.marker_detection.fetchResults()["RNA"], 0, "cohen", "mean");
        expect(stats.counts.length).toBeGreaterThan(0);
        expect(stats.num_markers).toBeGreaterThan(0); 
        expect(stats.num_markers).toBeLessThan(ntop * 1.5); // a bit of leeway for ties.
        expect(stats.counts.every(x => x > 0)).toBe(true); // all counts should be positive.
        expect(stats.pvalues.some(Number.isNaN)).toBe(false); // make sure we computed a valid p-value
    }

    // Check that the scores work correctly.
    {
        let chosen = -1;
        for (var i = 0; i < sets.sizes.length; i++) {
            if (sets.sizes[i] > 5) {
                chosen = i;
                break;
            }
        }

        expect(chosen).toBeGreaterThanOrEqual(0);
        let output = state.feature_set_enrichment.computePerCellScores(chosen);
        expect(output.weights.length).toEqual(sets.sizes[chosen]);
        expect(output.scores.length).toEqual(state.cell_filtering.fetchFilteredMatrix().numberOfColumns());
    }

    // Checking that automatic discovery actually has an effect.
    params.feature_set_enrichment.automatic = false;
    params.feature_set_enrichment.gene_id_type = "ENSEMBL";
    params.feature_set_enrichment.species = ["10090"];
    {
        await bakana.runAnalysis(state, mm_files, params);
        let collections = state.feature_set_enrichment.fetchCollectionDetails();
        expect(Array.from(new Set(collections.species))).toEqual(["10090"]);
        expect(state.feature_set_enrichment.fetchUniverseSize()).toEqual(0);
    }

    params.feature_set_enrichment.gene_id_type = "SYMBOL";
    {
        await bakana.runAnalysis(state, mm_files, params);
        let stats = state.feature_set_enrichment.computeEnrichment(state.marker_detection.fetchResults()["RNA"], 0, "cohen", "mean");
        expect(Array.from(new Set(collections.species))).toEqual(["10090"]);
        expect(state.feature_set_enrichment.fetchUniverseSize()).toBeGreaterThan(10000);
    }

    // Checking that versus mode works correctly.
    {
        let v01 = state.marker_detection.computeVersus(0, 1);
        let e01 = state.feature_set_enrichment.computeEnrichment(v01.results.RNA, v01.left, "cohen", "mean");
        expect(e01.counts instanceof Int32Array).toBe(true);
        expect(e01.counts.length).toBeGreaterThan(0);

        let v10 = state.marker_detection.computeVersus(1, 0);
        let e10 = state.feature_set_enrichment.computeEnrichment(v10.results.RNA, v10.left, "cohen", "mean");
        expect(e10.counts.length).toBeGreaterThan(0);
    }

    // Release me!
    await bakana.freeAnalysis(state);
})
