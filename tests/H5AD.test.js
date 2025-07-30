import * as bakana from "../src/index.js";
import * as scran from "scran.js";
import * as utils from "./utils.js";
import * as bioc from "bioconductor";

beforeAll(utils.initializeAll);
afterAll(async () => await bakana.terminate());

let fpath = "files/datasets/zeisel-brain.h5ad";
let files = { 
    default: new bakana.H5adDataset(fpath)
};

test("H5AD summary works correctly", async () => {
    let summ = files.default.summary();
    expect(summ.all_features instanceof bioc.DataFrame).toBe(true);
    expect(summ.cells instanceof bioc.DataFrame).toBe(true);
    expect(summ.all_assay_names.length).toBeGreaterThan(0);

    let preview = files.default.previewPrimaryIds();
    expect("RNA" in preview).toBe(true);
    expect(preview.RNA.length).toBeGreaterThan(0);
})

test("runAnalysis works correctly (H5AD)", async () => {
    let state = await bakana.createAnalysis();
    let params = utils.baseParams();
    let res = await bakana.runAnalysis(state, files, params);

    // Input reorganization is done correctly.
    {
        let loaded = state.inputs.fetchCountMatrix().get("RNA");
        let loaded_names = state.inputs.fetchFeatureAnnotations()["RNA"].rowNames();

        let simple = scran.initializeSparseMatrixFromHdf5(fpath, "X", { layered: false });
        let simple_names = (new scran.H5File(fpath)).open("var").open("_index", { load: true }).values;

        utils.checkMatrixContents(simple, simple_names, loaded, loaded_names);
        simple.free();
    }

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

    // Check saving of results.
    utils.purgeDirectory("miscellaneous/from-tests/H5AD");
    await bakana.saveSingleCellExperiment(state, "H5AD", { directory: "miscellaneous/from-tests" });
    utils.purgeDirectory("miscellaneous/from-tests/H5AD_genes");
    await bakana.saveGenewiseResults(state, "H5AD_genes", { directory: "miscellaneous/from-tests" });

    // Check reloading of the parameters/datasets.
    {
        let saved = [];
        let saver = (n, k, f) => {
            saved.push(f.content());
            return String(saved.length);
        };

        let serialized = await bakana.serializeConfiguration(state, saver);
        let reloaded = bakana.unserializeDatasets(serialized.datasets, x => saved[Number(x) - 1]); 
        expect(reloaded.default instanceof bakana.H5adDataset);
        expect(serialized.parameters).toEqual(bakana.retrieveParameters(state));
    }

    // Freeing.
    await bakana.freeAnalysis(state);
})

