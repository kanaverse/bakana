import * as bakana from "../src/index.js";
import * as qc from "../src/steps/rna_quality_control.js";
import * as aqc from "../src/steps/adt_quality_control.js";
import * as scran from "scran.js";
import * as utils from "./utils.js"
import * as wa from "wasmarrays.js";

beforeAll(utils.initializeAll);
afterAll(async () => await bakana.terminate());

test("analysis works when we skip the QC steps", async () => {
    let files = {
        default: new bakana.TenxHdf5Dataset("files/datasets/immune_3.0.0-tenx.h5")
    };

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
