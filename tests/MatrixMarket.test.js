import * as bakana from "../src/index.js";
import * as butils from "../src/steps/utils/general.js";
import * as scran from "scran.js";
import * as utils from "./utils.js";
import * as bioc from "bioconductor";

beforeAll(utils.initializeAll);
afterAll(async () => await bakana.terminate());

let mtx_file = "files/datasets/pbmc3k-matrix.mtx.gz";
let feat_file = "files/datasets/pbmc3k-features.tsv.gz";
let files = { 
    default: new bakana.TenxMatrixMarketDataset(mtx_file, feat_file, "files/datasets/pbmc3k-barcodes.tsv.gz")
};

test("MatrixMarket summary works correctly", async () => {
    let summ = await files.default.summary();
    expect(summ.modality_features["Gene Expression"] instanceof bioc.DataFrame).toBe(true);
    expect(summ.modality_features["Gene Expression"].numberOfColumns()).toBeGreaterThan(0);
    expect(summ.cells instanceof bioc.DataFrame).toBe(true);
    expect(summ.cells.numberOfColumns()).toBeGreaterThan(0);

    let preview = await files.default.previewPrimaryIds();
    expect("RNA" in preview).toBe(true);
    expect(preview.RNA.length).toBeGreaterThan(0);
})

test("runAnalysis works correctly (MatrixMarket)", async () => {
    let attempts = new Set;
    let started = step => {
        attempts.add(step);
    };

    let completed = new Set;
    let finished = (step) => {
        completed.add(step);
    };

    let state = await bakana.createAnalysis();
    let params = utils.baseParams();
    let res = await bakana.runAnalysis(state, files, params, { startFun: started, finishFun: finished });

    // Check that the callbacks are actually called.
    expect(attempts.has("rna_quality_control")).toBe(true);
    expect(attempts.has("rna_pca")).toBe(true);
    expect(completed.has("rna_pca")).toBe(true);
    expect(completed.has("feature_selection")).toBe(true);
    expect(completed.has("cell_labelling")).toBe(true);

    // Input reorganization is done correctly.
    {
        let loaded = state.inputs.fetchCountMatrix().get("RNA");
        let loaded_names = state.inputs.fetchFeatureAnnotations()["RNA"].column("id");
        let loaded_ids = state.inputs.fetchRowIds()["RNA"];

        let simple = scran.initializeSparseMatrixFromMatrixMarket(mtx_file, { layered: false });
        let parsed = bakana.readTable((new bakana.SimpleFile(feat_file)).buffer(), { compression: "gz" });
        let simple_names = parsed.map(x => x[0]);

        utils.checkReorganization(simple.matrix, simple.row_ids, simple_names, loaded, loaded_ids, loaded_names);
        simple.matrix.free();
    }

    // Basic consistency checks.
    await utils.overlordCheckStandard(state);
    utils.checkClusterVersusMode(state);
    await utils.triggerAnimation(state);

    // Check saving of results.
    await bakana.saveSingleCellExperiment(state, "MatrixMarket", { directory: "miscellaneous/from-tests" });

    // Check reloading of the parameters/datasets.
    {
        let saved = [];
        let saver = (n, k, f) => {
            saved.push(f.content());
            return String(saved.length);
        };

        let serialized = await bakana.serializeConfiguration(state, saver);
        let reloaded = bakana.unserializeDatasets(serialized.datasets, x => saved[Number(x) - 1]); 
        expect(reloaded.default instanceof bakana.TenxMatrixMarketDataset);
        expect(serialized.parameters).toEqual(bakana.retrieveParameters(state));
    }

    // Release me!
    await bakana.freeAnalysis(state);
})

let minimal_files = { 
    default: new bakana.TenxMatrixMarketDataset("files/datasets/pbmc3k-matrix.mtx.gz", null, null)
};

test("MatrixMarket summary works correctly with the bare minimum", async () => {
    let summ = await minimal_files.default.summary();
    expect(summ.modality_features[""] instanceof bioc.DataFrame).toBe(true);
    expect(summ.modality_features[""].numberOfColumns()).toBe(0);
    expect(summ.cells instanceof bioc.DataFrame).toBe(true);
    expect(summ.cells.numberOfColumns()).toBe(0);
})

test("runAnalysis works correctly with the bare minimum (MatrixMarket)", async () => {
    let state = await bakana.createAnalysis();
    let params = utils.baseParams();
    let res = await bakana.runAnalysis(state, minimal_files, params);

    // Basic consistency checks.
    await utils.checkStateResultsMinimal(state);
    await utils.checkStateResultsRna(state, { exclusive: true, hasMito: false });
    await utils.checkStateResultsUnblocked(state);

    // No annotations, so no mitochondrial proportions.
    expect(state.inputs.fetchFeatureAnnotations()["RNA"].numberOfColumns()).toBe(0);
    expect(state.rna_quality_control.fetchFilters().thresholdsSubsetProportions(0)[0]).toBe(0);

    // Feature set enrichment behaves.
    params.feature_set_enrichment.skip = false;
    await bakana.runAnalysis(state, minimal_files, params);
    expect(state.marker_detection.changed).toBe(false);
    expect(state.feature_set_enrichment.changed).toBe(true);
    expect(state.feature_set_enrichment.fetchCollectionDetails().names.length).toBe(0);

    // Release me!
    await bakana.freeAnalysis(state);
})
