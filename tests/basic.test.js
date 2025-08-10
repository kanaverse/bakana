import * as bakana from "../src/index.js";
import * as scran from "scran.js";
import * as utils from "./utils.js";
import * as bioc from "bioconductor";

beforeAll(utils.initializeAll);
afterAll(async () => await bakana.terminate());

test("basic runAnalysis works correctly (PBMC)", async () => {
    let mtx_file = "files/datasets/pbmc3k-matrix.mtx.gz";
    let feat_file = "files/datasets/pbmc3k-features.tsv.gz";
    let files = { 
        default: new bakana.TenxMatrixMarketDataset(mtx_file, feat_file, "files/datasets/pbmc3k-barcodes.tsv.gz")
    };

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

    // Basic consistency checks.
    await utils.overlordCheckStandard(state);
    utils.checkClusterVersusMode(state);
    await utils.triggerAnimation(state);

    // Check saving of results.
    utils.purgeDirectory("miscellaneous/from-tests/pbmc");
    await bakana.saveSingleCellExperiment(state, "pbmc", { directory: "miscellaneous/from-tests" });
    await utils.checkSavedExperiment("miscellaneous/from-tests/pbmc", state);
    utils.purgeDirectory("miscellaneous/from-tests/pbmc_genes");
    await bakana.saveGenewiseResults(state, "pbmc_genes", { directory: "miscellaneous/from-tests" });

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

test("basic runAnalysis works correctly (Zeisel)", async () => {
    let fpath = "files/datasets/zeisel-brain.rds";
    let files = { 
        default: new bakana.SummarizedExperimentDataset(fpath)
    };
    files.default.setOptions({ rnaExperiment: "endogenous" });

    let state = await bakana.createAnalysis();
    let params = utils.baseParams();
    let res = await bakana.runAnalysis(state, files, params);

    // Basic checks.
    await utils.overlordCheckStandard(state);
    utils.checkClusterVersusMode(state);
    await utils.triggerAnimation(state);

    // Annotations, with and without filtering.
    {
        let nfull = state.inputs.fetchCountMatrix().numberOfColumns();
        let nfilt = state.cell_filtering.fetchFilteredMatrix().numberOfColumns();

        let cell_anno = state.inputs.fetchCellAnnotations().column("level1class");
        expect(cell_anno.length).toBe(nfull);
        let filtered_cell_anno = state.cell_filtering.applyFilter(cell_anno);
        expect(filtered_cell_anno.length).toBe(nfilt);

        let sex_anno = state.inputs.fetchCellAnnotations().column("sex");
        expect(sex_anno.length).toBe(nfull);
        let filtered_sex_anno = state.cell_filtering.applyFilter(sex_anno);
        expect(filtered_sex_anno.length).toBe(nfilt);
    }

    utils.purgeDirectory("miscellaneous/from-tests/zeisel");
    await bakana.saveSingleCellExperiment(state, "zeisel", { directory: "miscellaneous/from-tests" });
    await utils.checkSavedExperiment("miscellaneous/from-tests/zeisel", state);
    utils.purgeDirectory("miscellaneous/from-tests/zeisel_results");
    await bakana.saveGenewiseResults(state, "zeisel_results", { directory: "miscellaneous/from-tests" });

    // Release me!
    await bakana.freeAnalysis(state);
})
