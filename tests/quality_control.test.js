import * as bakana from "../src/index.js";
import * as qc from "../src/steps/rna_quality_control.js";
import * as aqc from "../src/steps/adt_quality_control.js";
import * as scran from "scran.js";
import * as utils from "./utils.js"
import * as wa from "wasmarrays.js";

beforeAll(utils.initializeAll);
afterAll(async () => await bakana.terminate());

let files = {
    default: new bakana.TenxHdf5Dataset("files/datasets/immune_3.0.0-tenx.h5")
};

test("analysis works when we skip the QC steps", async () => {
    // Running the analysis with skipped QC.
    let state = await bakana.createAnalysis();
    let paramcopy = utils.baseParams();
    paramcopy.cell_filtering.use_rna = false;
    paramcopy.cell_filtering.use_adt = false;

    await bakana.runAnalysis(state, files, paramcopy);
    expect(state.rna_quality_control.fetchKeep().length).toBeGreaterThan(0);
    expect(state.adt_quality_control.fetchKeep().length).toBeGreaterThan(0);

    expect(state.cell_filtering.fetchKeep()).toBeNull();
    let ncells = state.inputs.fetchCountMatrix().numberOfColumns();
    expect(state.cell_filtering.fetchFilteredMatrix().numberOfColumns()).toBe(ncells);

    // Check saving of results.
    utils.purgeDirectory("miscellaneous/from-tests/no-qc");
    await bakana.saveSingleCellExperiment(state, "no-qc", { directory: "miscellaneous/from-tests" });
    await utils.checkSavedExperiment("miscellaneous/from-tests/no-qc", state);

    // Just applying the RNA filtering.
    {
        paramcopy.cell_filtering.use_rna = true;
        paramcopy.cell_filtering.use_adt = false;

        await bakana.runAnalysis(state, null, paramcopy);
        expect(state.rna_quality_control.changed).toBe(false);
        expect(state.adt_quality_control.changed).toBe(false);
        expect(state.cell_filtering.changed).toBe(true);

        let ncells = state.inputs.fetchCountMatrix().numberOfColumns();
        expect(state.cell_filtering.fetchKeep().owner).toBe(state.rna_quality_control.fetchKeep()); // i.e. a view
        expect(state.cell_filtering.fetchKeep().length).toBe(ncells);
        expect(state.cell_filtering.fetchFilteredMatrix().numberOfColumns()).toBeLessThan(ncells);
    }

    await bakana.freeAnalysis(state);
})

test("different options for choosing mitochondrial genes in RNA QC works as expected", async () => {
    let state = await bakana.createAnalysis();
    let paramcopy = utils.baseParams();
    await state.inputs.compute(files, paramcopy.inputs);
    await state.rna_quality_control.compute(paramcopy.rna_quality_control);
    expect(state.rna_quality_control.fetchMetrics().subsetProportion(0).reduce((a, b) => a + b)).toBeGreaterThan(0); // ID guessing works.

    await state.inputs.compute(files, paramcopy.inputs); // running it again to set changed = false.

    paramcopy.rna_quality_control.guess_ids = false;
    paramcopy.rna_quality_control.gene_id_column = "id";
    paramcopy.rna_quality_control.species = [ "9606" ];
    paramcopy.rna_quality_control.gene_id_type = "SYMBOL";
    await state.rna_quality_control.compute(paramcopy.rna_quality_control);
    expect(state.rna_quality_control.fetchMetrics().subsetProportion(0).reduce((a, b) => a + b)).toEqual(0);

    paramcopy.rna_quality_control.gene_id_type = "ENSEMBL";
    await state.rna_quality_control.compute(paramcopy.rna_quality_control);
    expect(state.rna_quality_control.fetchMetrics().subsetProportion(0).reduce((a, b) => a + b)).toBeGreaterThan(0);

    // Using the prefix.
    paramcopy.rna_quality_control.use_reference_mito = false;
    await state.rna_quality_control.compute(paramcopy.rna_quality_control);
    expect(state.rna_quality_control.fetchMetrics().subsetProportion(0).reduce((a, b) => a + b)).toEqual(0);

    paramcopy.rna_quality_control.gene_id_column = "name";
    await state.rna_quality_control.compute(paramcopy.rna_quality_control);
    expect(state.rna_quality_control.fetchMetrics().subsetProportion(0).reduce((a, b) => a + b)).toBeGreaterThan(0);

    paramcopy.rna_quality_control.mito_prefix = null;
    await state.rna_quality_control.compute(paramcopy.rna_quality_control);
    expect(state.rna_quality_control.fetchMetrics().subsetProportion(0).reduce((a, b) => a + b)).toBe(0);
    paramcopy.rna_quality_control.mito_prefix = "mt-";  // restoring.

    // When ID guessing is enabled, changes to the other parameters have no effect.
    paramcopy.rna_quality_control.guess_ids = true;
    await state.rna_quality_control.compute(paramcopy.rna_quality_control);
    expect(state.rna_quality_control.fetchMetrics().subsetProportion(0).reduce((a, b) => a + b)).toBeGreaterThan(0);
    expect(state.rna_quality_control.changed).toBe(true);
    expect(state.rna_quality_control.fetchParameters().gene_id_column).toBe("name");

    paramcopy.rna_quality_control.gene_id_column = "foobar"; // introducing a change.
    await state.rna_quality_control.compute(paramcopy.rna_quality_control);
    expect(state.rna_quality_control.fetchMetrics().subsetProportion(0).reduce((a, b) => a + b)).toBeGreaterThan(0);
    expect(state.rna_quality_control.changed).toBe(false);
    expect(state.rna_quality_control.fetchParameters().gene_id_column).toBe("foobar");

    await bakana.freeAnalysis(state);
})

test("different options for choosing IgG tags in ADT QC works as expected", async () => {
    let state = await bakana.createAnalysis();
    let paramcopy = utils.baseParams();
    await state.inputs.compute(files, paramcopy.inputs);
    await state.adt_quality_control.compute(paramcopy.adt_quality_control);
    expect(state.adt_quality_control.fetchMetrics().subsetSum(0).reduce((a, b) => a + b)).toBeGreaterThan(0); // ID guessing works.

    await state.inputs.compute(files, paramcopy.inputs); // running it again to set changed = false.

    // Disabling everything.
    paramcopy.adt_quality_control.guess_ids = false;
    paramcopy.adt_quality_control.tag_id_column = "id";
    paramcopy.adt_quality_control.igg_prefix = "igg1_"; // underscore is deliberate, otherwise we'd just match all the IDs.
    await state.adt_quality_control.compute(paramcopy.adt_quality_control);
    expect(state.adt_quality_control.fetchMetrics().subsetSum(0).reduce((a, b) => a + b)).toEqual(0); 

    paramcopy.adt_quality_control.tag_id_column = "name";
    await state.adt_quality_control.compute(paramcopy.adt_quality_control);
    expect(state.adt_quality_control.fetchMetrics().subsetSum(0).reduce((a, b) => a + b)).toBeGreaterThan(0); 

    paramcopy.adt_quality_control.igg_prefix = null;
    await state.adt_quality_control.compute(paramcopy.adt_quality_control);
    expect(state.adt_quality_control.fetchMetrics().subsetSum(0).reduce((a, b) => a + b)).toEqual(0); 
    paramcopy.adt_quality_control.igg_prefix = "igg";

    // When ID guessing is enabled, changes to the other parameters have no effect.
    paramcopy.adt_quality_control.guess_ids = true;
    await state.adt_quality_control.compute(paramcopy.adt_quality_control);
    expect(state.adt_quality_control.fetchMetrics().subsetSum(0).reduce((a, b) => a + b)).toBeGreaterThan(0); 
    expect(state.adt_quality_control.changed).toBe(true);
    expect(state.adt_quality_control.fetchParameters().tag_id_column).toBe("name");

    paramcopy.adt_quality_control.tag_id_column = "foobar"; // introducing a change.
    await state.adt_quality_control.compute(paramcopy.adt_quality_control);
    expect(state.adt_quality_control.fetchMetrics().subsetSum(0).reduce((a, b) => a + b)).toBeGreaterThan(0);
    expect(state.adt_quality_control.changed).toBe(false);
    expect(state.adt_quality_control.fetchParameters().tag_id_column).toBe("foobar");
})

test("we can switch to manual QC for RNA and ADTs", async () => {
    let state = await bakana.createAnalysis();
    let paramcopy = utils.baseParams();
    await state.inputs.compute(files, paramcopy.inputs);

    // For RNA.
    {
        let rna_param = paramcopy.rna_quality_control;
        rna_param.filter_strategy = "manual";
        rna_param.detected_threshold = 1234;
        await state.rna_quality_control.compute(rna_param);

        let filt = state.rna_quality_control.fetchFilters();
        let old_sums = filt.sum();
        expect(Array.from(old_sums)).toEqual([100]);
        let old_detected = filt.detected();
        expect(Array.from(old_detected)).toEqual([1234]);
        let old_mito = filt.subsetProportion(0);
        expect(Array.from(old_mito)).toEqual([0.2]);

        let params = state.rna_quality_control.fetchParameters();
        expect(params.filter_strategy).toEqual("manual");
        expect(params.detected_threshold).toEqual(1234);

        // Switching back.
        rna_param.filter_strategy = "automatic";
        await state.rna_quality_control.compute(rna_param);
        expect(state.rna_quality_control.changed).toBe(true);

        params = state.rna_quality_control.fetchParameters();
        filt = state.rna_quality_control.fetchFilters();
        expect(params.filter_strategy).toEqual("automatic");
        expect(filt.sum()).not.toEqual(old_sums);
        expect(filt.detected()).not.toEqual(old_detected);
        expect(filt.subsetProportion(0)).not.toEqual(old_mito);
    }

    // For ADT.
    {
        let adt_param = paramcopy.adt_quality_control;
        adt_param.filter_strategy = "manual";
        adt_param.detected_threshold = 1234;
        await state.adt_quality_control.compute(adt_param);

        let filt = state.adt_quality_control.fetchFilters();
        let old_detected = filt.detected();
        expect(Array.from(old_detected)).toEqual([1234]);
        let old_igg = filt.subsetSum(0);
        expect(Array.from(old_igg)).toEqual([1]);

        let params = state.adt_quality_control.fetchParameters();
        expect(params.filter_strategy).toEqual("manual");
        expect(params.detected_threshold).toEqual(1234);

        // Switching back.
        adt_param.filter_strategy = "automatic";
        await state.adt_quality_control.compute(adt_param);
        expect(state.adt_quality_control.changed).toBe(true);

        params = state.adt_quality_control.fetchParameters();
        filt = state.adt_quality_control.fetchFilters();
        expect(params.filter_strategy).toEqual("automatic");
        expect(filt.detected()).not.toEqual(old_detected);
        expect(filt.subsetSum(0)).not.toEqual(old_igg);
    }
})

test("we can switch to manual QC for CRISPR", async () => {
    const h5file = "files/datasets/crispr_6.0.0-tenx.h5";
    let files = { default: new bakana.TenxHdf5Dataset(h5file) };

    let state = await bakana.createAnalysis();
    let paramcopy = utils.baseParams();
    await state.inputs.compute(files, paramcopy.inputs);

    let crispr_param = paramcopy.crispr_quality_control;
    crispr_param.filter_strategy = "manual";
    crispr_param.max_threshold = 1234;
    await state.crispr_quality_control.compute(crispr_param);

    let filt = state.crispr_quality_control.fetchFilters();
    let old_max = filt.maxValue();
    expect(Array.from(old_max)).toEqual([1234]);

    let params = state.crispr_quality_control.fetchParameters();
    expect(params.filter_strategy).toEqual("manual");
    expect(params.max_threshold).toEqual(1234);

    // Switching back.
    crispr_param.filter_strategy = "automatic";
    await state.crispr_quality_control.compute(crispr_param);
    expect(state.crispr_quality_control.changed).toBe(true);

    params = state.crispr_quality_control.fetchParameters();
    filt = state.crispr_quality_control.fetchFilters();
    expect(params.filter_strategy).toEqual("automatic");
    expect(filt.maxValue()).not.toEqual(old_max);
})
