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
    expect(state.rna_quality_control.fetchDiscards().length).toBeGreaterThan(0);
    expect(state.adt_quality_control.fetchDiscards().length).toBeGreaterThan(0);

    expect(state.cell_filtering.fetchDiscards()).toBeNull();
    let ncells = state.inputs.fetchCountMatrix().numberOfColumns();
    expect(state.cell_filtering.fetchFilteredMatrix().numberOfColumns()).toBe(ncells);

    // Just applying the RNA filtering.
    {
        paramcopy.cell_filtering.use_rna = true;
        paramcopy.cell_filtering.use_adt = false;

        await bakana.runAnalysis(state, null, paramcopy);
        expect(state.rna_quality_control.changed).toBe(false);
        expect(state.adt_quality_control.changed).toBe(false);
        expect(state.cell_filtering.changed).toBe(true);

        let ncells = state.inputs.fetchCountMatrix().numberOfColumns();
        expect(state.cell_filtering.fetchDiscards().owner).toBe(state.rna_quality_control.fetchDiscards()); // i.e. a view
        expect(state.cell_filtering.fetchDiscards().length).toBe(ncells);
        expect(state.cell_filtering.fetchFilteredMatrix().numberOfColumns()).toBeLessThan(ncells);
    }

    await bakana.freeAnalysis(state);
})

test("different options for choosing mitochondrial genes in RNA QC works as expected", async () => {
    let state = await bakana.createAnalysis();
    let paramcopy = utils.baseParams();
    await state.inputs.compute(files, paramcopy.inputs);
    await state.rna_quality_control.compute(paramcopy.rna_quality_control);
    expect(state.rna_quality_control.fetchMetrics().subsetProportions(0).reduce((a, b) => a + b)).toBeGreaterThan(0); // automatic choice works.

    await state.inputs.compute(files, paramcopy.inputs); // running it again to set changed = false.

    paramcopy.rna_quality_control.automatic = false;
    paramcopy.rna_quality_control.gene_id_column = "id";
    paramcopy.rna_quality_control.species = [ "9606" ];
    paramcopy.rna_quality_control.gene_id_type = "SYMBOL";
    await state.rna_quality_control.compute(paramcopy.rna_quality_control);
    expect(state.rna_quality_control.fetchMetrics().subsetProportions(0).reduce((a, b) => a + b)).toEqual(0);

    paramcopy.rna_quality_control.gene_id_type = "ENSEMBL";
    await state.rna_quality_control.compute(paramcopy.rna_quality_control);
    expect(state.rna_quality_control.fetchMetrics().subsetProportions(0).reduce((a, b) => a + b)).toBeGreaterThan(0);

    // Using the prefix.
    paramcopy.rna_quality_control.use_reference_mito = false;
    await state.rna_quality_control.compute(paramcopy.rna_quality_control);
    expect(state.rna_quality_control.fetchMetrics().subsetProportions(0).reduce((a, b) => a + b)).toEqual(0);

    paramcopy.rna_quality_control.gene_id_column = "name";
    await state.rna_quality_control.compute(paramcopy.rna_quality_control);
    expect(state.rna_quality_control.fetchMetrics().subsetProportions(0).reduce((a, b) => a + b)).toBeGreaterThan(0);

    paramcopy.rna_quality_control.mito_prefix = null;
    await state.rna_quality_control.compute(paramcopy.rna_quality_control);
    expect(state.rna_quality_control.fetchMetrics().subsetProportions(0).reduce((a, b) => a + b)).toBe(0);
    paramcopy.rna_quality_control.mito_prefix = "mt-";  // restoring.

    // When automatic discovery is enabled, changes to the other parameters have no effect.
    paramcopy.rna_quality_control.automatic = true;
    await state.rna_quality_control.compute(paramcopy.rna_quality_control);
    expect(state.rna_quality_control.fetchMetrics().subsetProportions(0).reduce((a, b) => a + b)).toBeGreaterThan(0);
    expect(state.rna_quality_control.changed).toBe(true);
    expect(state.rna_quality_control.fetchParameters().gene_id_column).toBe("name");

    paramcopy.rna_quality_control.gene_id_column = "foobar"; // introducing a change.
    await state.rna_quality_control.compute(paramcopy.rna_quality_control);
    expect(state.rna_quality_control.fetchMetrics().subsetProportions(0).reduce((a, b) => a + b)).toBeGreaterThan(0);
    expect(state.rna_quality_control.changed).toBe(false);
    expect(state.rna_quality_control.fetchParameters().gene_id_column).toBe("foobar");

    await bakana.freeAnalysis(state);
})

test("different options for choosing IgG tags in ADT QC works as expected", async () => {
    let state = await bakana.createAnalysis();
    let paramcopy = utils.baseParams();
    await state.inputs.compute(files, paramcopy.inputs);
    await state.adt_quality_control.compute(paramcopy.adt_quality_control);
    expect(state.adt_quality_control.fetchMetrics().subsetTotals(0).reduce((a, b) => a + b)).toBeGreaterThan(0); // automatic choice works.

    await state.inputs.compute(files, paramcopy.inputs); // running it again to set changed = false.

    // Disabling everything.
    paramcopy.adt_quality_control.automatic = false;
    paramcopy.adt_quality_control.tag_id_column = "id";
    paramcopy.adt_quality_control.igg_prefix = "igg1_"; // underscore is deliberate, otherwise we'd just match all the IDs.
    await state.adt_quality_control.compute(paramcopy.adt_quality_control);
    expect(state.adt_quality_control.fetchMetrics().subsetTotals(0).reduce((a, b) => a + b)).toEqual(0); 

    paramcopy.adt_quality_control.tag_id_column = "name";
    await state.adt_quality_control.compute(paramcopy.adt_quality_control);
    expect(state.adt_quality_control.fetchMetrics().subsetTotals(0).reduce((a, b) => a + b)).toBeGreaterThan(0); 

    paramcopy.adt_quality_control.igg_prefix = null;
    await state.adt_quality_control.compute(paramcopy.adt_quality_control);
    expect(state.adt_quality_control.fetchMetrics().subsetTotals(0).reduce((a, b) => a + b)).toEqual(0); 
    paramcopy.adt_quality_control.igg_prefix = "igg";

    // When automatic discovery is enabled, changes to the other parameters have no effect.
    paramcopy.adt_quality_control.automatic = true;
    await state.adt_quality_control.compute(paramcopy.adt_quality_control);
    expect(state.adt_quality_control.fetchMetrics().subsetTotals(0).reduce((a, b) => a + b)).toBeGreaterThan(0); 
    expect(state.adt_quality_control.changed).toBe(true);
    expect(state.adt_quality_control.fetchParameters().tag_id_column).toBe("name");

    paramcopy.adt_quality_control.tag_id_column = "foobar"; // introducing a change.
    await state.adt_quality_control.compute(paramcopy.adt_quality_control);
    expect(state.adt_quality_control.fetchMetrics().subsetTotals(0).reduce((a, b) => a + b)).toBeGreaterThan(0);
    expect(state.adt_quality_control.changed).toBe(false);
    expect(state.adt_quality_control.fetchParameters().tag_id_column).toBe("foobar");
})
