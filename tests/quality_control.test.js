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

test("guessing of the mitochondrial genes for RNA QC works as expected", async () => {
    let state = await bakana.createAnalysis();
    let paramcopy = utils.baseParams();
    state.inputs.compute(files, paramcopy.inputs);

    paramcopy.rna_quality_control.automatic = false;
    paramcopy.rna_quality_control.gene_id_column = "id";
    paramcopy.rna_quality_control.reference_mito_species = [ "human" ];
    paramcopy.rna_quality_control.reference_mito_type = "SYMBOL";
    await state.rna_quality_control.compute(paramcopy.rna_quality_control);
    expect(state.rna_quality_control.fetchMetrics().subsetProportions(0).reduce((a, b) => a + b)).toEqual(0);

    paramcopy.rna_quality_control.reference_mito_type = "ENSEMBL";
    await state.rna_quality_control.compute(paramcopy.rna_quality_control);
    expect(state.rna_quality_control.fetchMetrics().subsetProportions(0).reduce((a, b) => a + b)).toBeGreaterThan(0);

    paramcopy.rna_quality_control.use_reference_mito = false;
    await state.rna_quality_control.compute(paramcopy.rna_quality_control);
    expect(state.rna_quality_control.fetchMetrics().subsetProportions(0).reduce((a, b) => a + b)).toEqual(0);

    paramcopy.rna_quality_control.gene_id_column = "symbol";
    await state.rna_quality_control.compute(paramcopy.rna_quality_control);
    expect(state.rna_quality_control.fetchMetrics().subsetProportions(0).reduce((a, b) => a + b)).toBeGreaterThan(0);

    await bakana.freeAnalysis(state);
})
