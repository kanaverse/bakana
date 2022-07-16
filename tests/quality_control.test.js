import * as bakana from "../src/index.js";
import * as qc from "../src/steps/quality_control.js";
import * as aqc from "../src/steps/adt_quality_control.js";
import * as scran from "scran.js";
import * as valkana from "valkana";
import * as utils from "./utils.js"

beforeAll(utils.initializeAll);
afterAll(async () => await bakana.terminate());

test("analysis works when we skip the QC steps", async () => {
    let files = {
        default: {
            format: "10X",
            h5: "files/datasets/immune_3.0.0-tenx.h5"
        }
    };

    let contents = {};
    let finished = (step, res) => {
        contents[step] = res;
    };

    // Running the analysis with skipped QC.
    let state = await bakana.createAnalysis();
    let paramcopy = utils.baseParams();
    paramcopy.quality_control.skip = true;
    paramcopy.adt_quality_control.skip = true;

    await bakana.runAnalysis(state, files, paramcopy, { finishFun: finished });
    expect(state.quality_control.valid()).toBe(true);
    expect(state.quality_control.skipped()).toBe(true);
    expect(contents.quality_control).toBeNull();
    expect(state.adt_quality_control.valid()).toBe(true);
    expect(state.adt_quality_control.skipped()).toBe(true);
    expect(contents.adt_quality_control).toBeNull();

    expect(state.cell_filtering.fetchDiscards()).toBeNull();
    let ncells = state.inputs.fetchCountMatrix().numberOfColumns();
    expect(state.cell_filtering.fetchFilteredMatrix().numberOfColumns()).toBe(ncells);

    const path = "TEST_state_qc-skip.h5";
    let collected = await bakana.saveAnalysis(state, path);
    utils.validateState(path);
    {
        let handle = new scran.H5File(path);
        let qhandle = handle.open("quality_control");
        let qrhandle = qhandle.open("results");
        expect("metrics" in qrhandle.children).toBe(false);
        expect("discards" in qrhandle.children).toBe(false);

        let aqhandle = handle.open("adt_quality_control");
        let aqrhandle = aqhandle.open("results");
        expect("metrics" in aqrhandle.children).toBe(false);
        expect("discards" in aqrhandle.children).toBe(false);
    }

    {
        let offsets = utils.mockOffsets(collected.collected);
        let reloaded = await bakana.loadAnalysis(
            path, 
            (offset, size) => offsets[offset]
        );
        let new_params = bakana.retrieveParameters(reloaded);

        expect(reloaded.quality_control.skipped()).toBe(true);
        expect(new_params.quality_control.skip).toBe(true);

        expect(reloaded.adt_quality_control.skipped()).toBe(true);
        expect(new_params.adt_quality_control.skip).toBe(true);

        await bakana.freeAnalysis(reloaded);
    }

    // Forcing the upstream to be non-changed.
    await state.inputs.compute(files, null, null);

    for (const stepname of [ "quality_control", "adt_quality_control" ]) {
        let curstep = state[stepname];

        // Setting up some functions to run things automatically for RNA and ADT.
        let rerun;
        let rerun_change;
        let defaults = state[stepname].constructor.defaults();
        let validate;

        if (stepname == "adt_quality_control") {
            rerun = skip => {
                curstep.compute(skip, defaults.igg_prefix, defaults.nmads, defaults.min_detected_drop);
            };

            rerun_change = skip => {
                curstep.compute(skip, defaults.igg_prefix, defaults.nmads + 0.1, defaults.min_detected_drop);
            };

            validate = path => {
                valkana.validateAdtQualityControlState(path, ncells, 1, true, bakana.kanaFormatVersion);
            };

        } else {
            rerun = skip => {
                curstep.compute(skip, defaults.use_mito_default, defaults.mito_prefix, defaults.nmads);
            };

            rerun_change = skip => {
                curstep.compute(skip, defaults.use_mito_default, defaults.mito_prefix, defaults.nmads + 0.1);
            };

            validate = path => {
                valkana.validateQualityControlState(path, ncells, 1, bakana.kanaFormatVersion);
            };
        }
        
        // Changing parameters have no effect when we're still skipping.
        rerun_change(true);
        expect(curstep.skipped()).toBe(true);
        expect(curstep.changed).toBe(false);

        // Until we reactivate everything... 
        rerun(false);
        expect(curstep.changed).toBe(true);
        expect(curstep.skipped()).toBe(false);

        {
            let has_discards = curstep.fetchDiscards();
            expect(has_discards.length).toBe(ncells);

            state.cell_filtering.compute();
            has_discards = state.cell_filtering.fetchDiscards();
            expect(has_discards.length).toBe(ncells);
            let remaining = state.cell_filtering.summary().retained;
            expect(remaining).toBeLessThan(ncells);
            expect(state.cell_filtering.fetchFilteredMatrix().numberOfColumns()).toEqual(remaining);

            let lost = 0;
            has_discards.forEach(x => { lost += x; });
            expect(lost).toBeGreaterThan(0);
        }

        // Skipping again, but the results are still valid, so we can save and load them.
        rerun(true);
        expect(curstep.changed).toBe(true);

        {
            let has_discards = curstep.fetchDiscards();
            expect(has_discards.length).toBe(ncells);

            state.cell_filtering.compute();
            expect(state.cell_filtering.fetchDiscards()).toBeNull();
            expect(state.cell_filtering.summary().retained).toEqual(ncells);
            expect(state.cell_filtering.fetchFilteredMatrix().numberOfColumns()).toEqual(ncells);

            let handle = scran.createNewHDF5File(path);
            curstep.serialize(handle);
            validate(path);
        }

        {
            let handle = new scran.H5File(path);
            let qhandle = handle.open(stepname);
            let qrhandle = qhandle.open("results");
            expect("metrics" in qrhandle.children).toBe(true);
            expect("discards" in qrhandle.children).toBe(true);

            let reloaded;
            if (stepname == "adt_quality_control") {
                reloaded = new aqc.unserialize(handle, state.inputs);
            } else {
                reloaded = new qc.unserialize(handle, state.inputs);
            }
            expect(reloaded.skipped()).toBe(true);

            let has_discards = reloaded.fetchDiscards();
            expect(has_discards.length).toBe(ncells);
            reloaded.free();
        }

        // Unskipping and just using the cached results.
        {
            let has_discards = curstep.fetchDiscards();
            has_discards.array()[0] = 2;
        }

        rerun(false);
        expect(curstep.changed).toBe(true);

        {
            let has_discards = curstep.fetchDiscards();
            expect(has_discards.array()[0]).toBe(2); // i.e., using the cached result.
        }

        // Reskipping, but with altered parameters, so we wipe the existing values. 
        rerun_change(true);
        expect(curstep.changed).toBe(true);

        {
            let handle = scran.createNewHDF5File(path);
            curstep.serialize(handle);
            validate(path);
        }

        {
            let handle = new scran.H5File(path);
            let qhandle = handle.open(stepname);
            let qrhandle = qhandle.open("results");
            expect("metrics" in qrhandle.children).toBe(true);
            expect("thresholds" in qrhandle.children).toBe(false);
            expect("discards" in qrhandle.children).toBe(false);
        }
    }

    await bakana.freeAnalysis(state);
})
