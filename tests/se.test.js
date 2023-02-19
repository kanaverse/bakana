import * as bakana from "../src/index.js";
import * as inputs from "../src/steps/inputs.js";
import * as scran from "scran.js";
import * as utils from "./utils.js";
import * as bioc from "bioconductor";

beforeAll(utils.initializeAll);
afterAll(async () => await bakana.terminate());

test("runAnalysis works correctly (RDS containing SingleCellExperiment)", async () => {
    let fpath = "files/datasets/zeisel-brain.rds";
    let files = { 
        default: new bakana.SummarizedExperimentDataset(fpath)
    };

    let state = await bakana.createAnalysis();
    let params = utils.baseParams();
    let res = await bakana.runAnalysis(state, files, params);

    // Input reorganization is done correctly. 
    {
        let loaded = state.inputs.fetchCountMatrix().get("RNA");
        let loaded_ids = state.inputs.fetchRowIds()["RNA"];
        let loaded_names = state.inputs.fetchFeatureAnnotations()["RNA"].rowNames();
        expect(loaded_names.length).toBeGreaterThan(0);

        let rdshandle = scran.readRds(fpath);
        let rhandle = rdshandle.value();
        expect(rhandle.className()).toBe("SingleCellExperiment");

        let rrhandle = rhandle.attribute("rowRanges");
        expect(rrhandle.className()).toMatch("GRangesList");
        let nhandle = rrhandle.attribute("partitioning").attribute("NAMES"); // technically a memory leak, but finalizers should handle it.
        let simple_names = nhandle.values();

        let ahandle = rhandle.attribute("assays").attribute("data").attribute("listData").load(0); // again, technically memory leaks here.
        let simple = scran.initializeSparseMatrixFromRds(ahandle, { layered: false, consume: true });

        ahandle.free();
        nhandle.free();
        rrhandle.free();
        rhandle.free();
        rdshandle.free();

        utils.checkReorganization(simple.matrix, simple.row_ids, simple_names, loaded, loaded_ids, loaded_names);
        simple.matrix.free();
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
    await bakana.saveSingleCellExperiment(state, "se", { directory: "miscellaneous/from-tests" });

    // Check reloading of the parameters/datasets.
    {
        let saved = [];
        let saver = (n, k, f) => {
            saved.push(f.content());
            return String(saved.length);
        };

        let serialized = await bakana.serializeConfiguration(state, saver);
        let reloaded = bakana.unserializeDatasets(serialized.datasets, x => saved[Number(x) - 1]); 
        expect(reloaded.default instanceof bakana.SummarizedExperimentDataset);
        expect(serialized.parameters).toEqual(bakana.retrieveParameters(state));
    }

    // Freeing.
    await bakana.freeAnalysis(state);
})

test("RDS loaders work correctly for a base SummarizedExperiment", async () => {
    let contents = {};
    let finished = (step, res) => {
        contents[step] = res;
    };

    let fpath = "files/datasets/paul-hsc.rds";
    let files = { 
        default: new bakana.SummarizedExperimentDataset(fpath)
    };

    // Check we're dealing with a pure SE.
    {
        let rdshandle = scran.readRds(fpath);
        let rhandle = rdshandle.value();
        expect(rhandle.className()).toBe("SummarizedExperiment");
        rhandle.free();
        rdshandle.free();
    }

    let fullstate = new inputs.InputsState;
    await fullstate.compute(files, inputs.InputsState.defaults());

    let loaded = fullstate.fetchCountMatrix();
    expect(loaded.numberOfColumns()).toBeGreaterThan(0);

    let loaded_names = fullstate.fetchFeatureAnnotations()["RNA"].rowNames();
    expect(loaded_names.length).toBeGreaterThan(0);

    fullstate.free();
})

test("RDS loaders work correctly for GRanges SummarizedExperiment", async () => {
    let contents = {};
    let finished = (step, res) => {
        contents[step] = res;
    };

    let fpath = "files/datasets/zeisel-brain-dense.rds";
    let files = { 
        default: new bakana.SummarizedExperimentDataset(fpath)
    };

    // Check we're dealing with a GRanges-containing SCE.
    {
        let rdshandle = scran.readRds(fpath);
        let rhandle = rdshandle.value();
        expect(rhandle.className()).toBe("SingleCellExperiment");

        let rrhandle = rhandle.attribute("rowRanges");
        expect(rrhandle.className()).toBe("GRanges");

        rrhandle.free();
        rhandle.free();
        rdshandle.free();
    }

    let fullstate = new inputs.InputsState;
    await fullstate.compute(files, inputs.InputsState.defaults());

    let loaded = fullstate.fetchCountMatrix();
    expect(loaded.numberOfColumns()).toBeGreaterThan(0);

    let loaded_names = fullstate.fetchFeatureAnnotations()["RNA"].rowNames();
    expect(loaded_names.length).toBeGreaterThan(0);

    fullstate.free();
})

test("RDS loaders work correctly for a SingleCellExperiment with altExps", async () => {
    let contents = {};
    let finished = (step, res) => {
        contents[step] = res;
    };

    let fpath = "files/datasets/immune_3.0.0-tenx.rds";
    let files = { 
        default: new bakana.SummarizedExperimentDataset(fpath)
    };

    // Summary works correctly.
    let summ = files.default.summary();
    expect(summ.modality_features[""] instanceof bioc.DataFrame).toBe(true);
    expect(summ.modality_features["Antibody Capture"] instanceof bioc.DataFrame).toBe(true);
    expect(summ.cells instanceof bioc.DataFrame).toBe(true);
    expect(summ.modality_assay_names[""].length).toBeGreaterThan(0);
    expect(summ.modality_assay_names["Antibody Capture"].length).toBeGreaterThan(0);

    let fullstate = new inputs.InputsState;
    await fullstate.compute(files, inputs.InputsState.defaults());

    let loaded = fullstate.fetchCountMatrix();
    let loaded_rna = loaded.get("RNA");
    expect(loaded_rna.numberOfColumns()).toBeGreaterThan(0);
    let loaded_adt = loaded.get("ADT");
    expect(loaded_adt.numberOfColumns()).toBeGreaterThan(0);

    let feat_anno = fullstate.fetchFeatureAnnotations();

    let loaded_names = feat_anno["RNA"].rowNames();
    expect(loaded_names.length).toBeGreaterThan(0);
    expect(loaded_names.length).toEqual(loaded_rna.numberOfRows());
    expect(loaded_names.some(x => x.match(/^ENSG/))).toBe(true);

    let loaded_names_adt = feat_anno["ADT"].rowNames();
    expect(loaded_names_adt.length).toBeGreaterThan(0);
    expect(loaded_names_adt.length).toEqual(loaded_adt.numberOfRows());
    expect(loaded_names_adt.some(x => x.match(/^ENSG/))).toBe(false);

    fullstate.free();
})
