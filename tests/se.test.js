import * as bakana from "../src/index.js";
import * as inputs from "../src/steps/inputs.js";
import * as scran from "scran.js";
import * as utils from "./utils.js";

beforeAll(utils.initializeAll);
afterAll(async () => await bakana.terminate());

test("runAnalysis works correctly (RDS containing SingleCellExperiment)", async () => {
    let contents = {};
    let finished = (step, res) => {
        contents[step] = res;
    };

    let files = { 
        default: {
            format: "SummarizedExperiment",
            rds: "files/datasets/zeisel-brain.rds"
        }
    };

    let state = await bakana.createAnalysis();
    let params = utils.baseParams();
    let res = await bakana.runAnalysis(state, files, params, { finishFun: finished });

    expect(contents.quality_control instanceof Object).toBe(true);
    expect(contents.pca instanceof Object).toBe(true);
    expect(contents.feature_selection instanceof Object).toBe(true);
    expect(contents.cell_labelling instanceof Object).toBe(true);
    expect(contents.marker_detection instanceof Object).toBe(true);

    // Input reorganization is done correctly. We compare it against the 
    // H5AD file, to avoid having to regurgitate the RDS loading code.
    {
        let loaded = state.inputs.fetchCountMatrix();
        let loaded_ids = state.inputs.fetchRowIds();
        let loaded_names = state.inputs.fetchGenes().id;
        expect(loaded_names.length).toBeGreaterThan(0);

        let rdshandle = scran.readRds(files.default.rds);
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
    const path = "TEST_state_SummarizedExperiment.h5";
    let collected = await bakana.saveAnalysis(state, path);
    utils.validateState(path);
    expect(collected.collected.length).toBe(1);
    expect(typeof(collected.collected[0])).toBe("string");

    let offsets = utils.mockOffsets(collected.collected);
    let reloaded = await bakana.loadAnalysis(
        path, 
        (offset, size) => offsets[offset]
    );

    let new_params = bakana.retrieveParameters(reloaded);
    expect(new_params.quality_control instanceof Object).toBe(true);
    expect(new_params.pca instanceof Object).toBe(true);

    // Freeing.
    await bakana.freeAnalysis(state);
    await bakana.freeAnalysis(reloaded);
})

test("RDS loaders work correctly for a base SummarizedExperiment", async () => {
    let contents = {};
    let finished = (step, res) => {
        contents[step] = res;
    };

    let files = { 
        default: {
            format: "SummarizedExperiment",
            rds: "files/datasets/paul-hsc.rds"
        }
    };

    // Check we're dealing with a pure SE.
    {
        let rdshandle = scran.readRds(files.default.rds);
        let rhandle = rdshandle.value();
        expect(rhandle.className()).toBe("SummarizedExperiment");
        rhandle.free();
        rdshandle.free();
    }

    let fullstate = new inputs.InputsState;
    await fullstate.compute(files, null, null);

    let loaded = fullstate.fetchCountMatrix();
    expect(loaded.numberOfColumns()).toBeGreaterThan(0);

    let loaded_names = fullstate.fetchGenes().id;
    expect(loaded_names.length).toBeGreaterThan(0);

    fullstate.free();
})

test("RDS loaders work correctly for GRanges SummarizedExperiment", async () => {
    let contents = {};
    let finished = (step, res) => {
        contents[step] = res;
    };

    let files = { 
        default: {
            format: "SummarizedExperiment",
            rds: "files/datasets/zeisel-brain-dense.rds"
        }
    };

    // Check we're dealing with a GRanges-containing SCE.
    {
        let rdshandle = scran.readRds(files.default.rds);
        let rhandle = rdshandle.value();
        expect(rhandle.className()).toBe("SingleCellExperiment");

        let rrhandle = rhandle.attribute("rowRanges");
        expect(rrhandle.className()).toBe("GRanges");

        rrhandle.free();
        rhandle.free();
        rdshandle.free();
    }

    let fullstate = new inputs.InputsState;
    await fullstate.compute(files, null, null);

    let loaded = fullstate.fetchCountMatrix();
    expect(loaded.numberOfColumns()).toBeGreaterThan(0);

    let loaded_names = fullstate.fetchGenes().id;
    expect(loaded_names.length).toBeGreaterThan(0);

    fullstate.free();
})

test("RDS loaders work correctly for a SingleCellExperiment with altExps", async () => {
    let contents = {};
    let finished = (step, res) => {
        contents[step] = res;
    };

    let files = { 
        default: {
            format: "SummarizedExperiment",
            rds: "files/datasets/immune_3.0.0-tenx.rds"
        }
    };

    let fullstate = new inputs.InputsState;
    await fullstate.compute(files, null, null);

    let loaded = fullstate.fetchCountMatrix({ type: "RNA" });
    expect(loaded.numberOfColumns()).toBeGreaterThan(0);

    let loaded_adt = fullstate.fetchCountMatrix({ type: "ADT" });
    expect(loaded_adt.numberOfColumns()).toBeGreaterThan(0);

    let loaded_names = fullstate.fetchGenes({ type: "RNA" }).id;
    expect(loaded_names.length).toBeGreaterThan(0);
    expect(loaded_names.length).toEqual(loaded.numberOfRows());
    expect(loaded_names.some(x => x.match(/^ENSG/))).toBe(true);

    let loaded_names_adt = fullstate.fetchGenes({ type: "ADT" }).id;
    expect(loaded_names_adt.length).toBeGreaterThan(0);
    expect(loaded_names_adt.length).toEqual(loaded_adt.numberOfRows());
    expect(loaded_names_adt.some(x => x.match(/^ENSG/))).toBe(false);

    fullstate.free();
})
