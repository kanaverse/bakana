import * as bakana from "../src/index.js";
import * as qc from "../src/steps/quality_control.js";
import * as aqc from "../src/steps/adt_quality_control.js";
import * as scran from "scran.js";
import * as utils from "./utils.js"

beforeAll(async () => await bakana.initialize({ localFile: true }));
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
    expect(state.quality_control.skipped()).toBe(true);
    expect(contents.quality_control).toEqual({});
    expect(state.adt_quality_control.skipped()).toBe(true);
    expect(contents.adt_quality_control).toEqual({});

    expect(state.cell_filtering.fetchDiscards()).toBeNull();
    let ncells = state.inputs.fetchCountMatrix().numberOfColumns();
    expect(state.cell_filtering.fetchFilteredMatrix().numberOfColumns()).toBe(ncells);

    const path = "TEST_state_qc-skip.h5";
    let collected = await bakana.saveAnalysis(state, path);
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

        expect(reloaded.state.quality_control.skipped()).toBe(true);
        expect(reloaded.parameters.quality_control.skip).toBe(true);

        expect(reloaded.state.adt_quality_control.skipped()).toBe(true);
        expect(reloaded.parameters.adt_quality_control.skip).toBe(true);

        await bakana.freeAnalysis(reloaded.state);    
    }

    // Forcing the upstream to be non-changed.
    state.inputs.compute(files, null);

    for (const stepname of [ "quality_control", "adt_quality_control" ]) {
        let curstep = state[stepname];

        // Setting up some functions to run things automatically for RNA and ADT.
        let rerun;
        let rerun_change;
        let defaults = state[stepname].constructor.defaults();

        if (stepname == "adt_quality_control") {
            rerun = skip => {
                curstep.compute(skip, defaults.igg_prefix, defaults.nmads, defaults.min_detected_drop);
            };

            rerun_change = skip => {
                curstep.compute(skip, defaults.igg_prefix, defaults.nmads + 0.1, defaults.min_detected_drop);
            };
        } else {
            rerun = skip => {
                curstep.compute(skip, defaults.use_mito_default, defaults.mito_prefix, defaults.nmads);
            };

            rerun_change = skip => {
                curstep.compute(skip, defaults.use_mito_default, defaults.mito_prefix, defaults.nmads + 0.1);
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

            let handle = scran.createNewHDF5File(path);
            curstep.serialize(handle);
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
            let state2 = reloaded.state;
            expect(state2.skipped()).toBe(true);

            let has_discards = state2.fetchDiscards();
            expect(has_discards.length).toBe(ncells);
            state2.free();
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
