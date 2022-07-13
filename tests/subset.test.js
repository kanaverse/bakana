import * as bakana from "../src/index.js";
import * as inputs from "../src/steps/inputs.js";
import * as scran from "scran.js";
import * as valkana from "valkana";
import * as utils from "./utils.js"
import * as fs from "fs";

beforeAll(utils.initializeAll);
afterAll(async () => await bakana.terminate());

function createSaver(saved) {
    return (serialized, size) => {
        saved.collected.push(serialized);
        let current = saved.total;
        saved.total += size;
        return {
            "offset": current,
            "size": size
        };
    };
}

test("subsetting behaves correctly with indices", async () => {
    let files = {
        default: {
            format: "H5AD",
            h5: "files/datasets/zeisel-brain.h5ad"
        }
    };

    let fullstate = new inputs.InputsState;
    await fullstate.compute(files, null, null);
    
    let ncells = fullstate.fetchCountMatrix().numberOfColumns();
    let subset = [];
    for (var i = 0; i < ncells; i += 2) {
        subset.push(i); // all evens.
    }

    let istate = new inputs.InputsState;
    istate.setDirectSubset(new Int32Array(subset));
    await istate.compute(files, null, null);

    expect(istate.fetchCountMatrix().numberOfColumns()).toBe(subset.length);
    expect(istate.summary().num_cells).toBe(subset.length);
    expect(istate.fetchCountMatrix().column(2)).toEqual(fullstate.fetchCountMatrix().column(4));
    expect(istate.fetchCountMatrix().column(5)).toEqual(fullstate.fetchCountMatrix().column(10));

    {
        let indices = [0,2,4,8];
        fullstate.undoSubset(indices);
        expect(indices).toEqual([0,2,4,8]);
        istate.undoSubset(indices);
        expect(indices).toEqual([0,4,8,16]);
    }

    {
        let subanno = istate.fetchAnnotations("level1class");
        expect(subanno.length).toBe(subset.length);

        let fullanno = fullstate.fetchAnnotations("level1class");
        let expected = subset.map(i => fullanno[i]);
        expect(subanno).toEqual(expected);
    }

    // Changes to the subset is a no-op when indices are around.
    {
        expect(istate.changed).toBe(true);
        await istate.compute(files, null, { field: "level1class", values: ["FOO", "BLAH"] });
        expect(istate.changed).toBe(false);
    }

    // Checking the serialization and unserialization.
    const path = "TEST_subset-inputs.h5";

    let saved = { collected: [], total: 0 };
    {
        let handle = scran.createNewHDF5File(path);
        await istate.serialize(handle, createSaver(saved));
        valkana.validateInputsState(path, true, bakana.kanaFormatVersion);
    }

    let offsets = utils.mockOffsets(saved.collected);
    {
        let handle = new scran.H5File(path);
        let restate = await inputs.unserialize(handle, (offset, size) => offsets[offset]);
        expect(restate.state.fetchCountMatrix().numberOfColumns()).toBe(subset.length);

        let indices = restate.state.fetchDirectSubset();
        expect(Array.from(indices)).toEqual(subset);
        expect(restate.state.fetchParameters().subset).toBeNull();

        restate.state.free();
    }

    // Unsetting works as expected.
    istate.setDirectSubset(null);
    expect(istate.changed).toBe(true);
    expect(istate.summary().num_cells).toBe(fullstate.summary().num_cells);
    expect(istate.fetchCountMatrix().column(1)).toEqual(fullstate.fetchCountMatrix().column(1));

    // Resetting works as expected.
    istate.setDirectSubset(subset);
    expect(istate.changed).toBe(true);
    expect(istate.summary().num_cells).toBe(subset.length);
    expect(istate.fetchCountMatrix().column(3)).toEqual(fullstate.fetchCountMatrix().column(6));

    istate.free();
    fullstate.free();
})

test("subsetting behaves correctly with a factor", async () => {
    let files = {
        default: {
            format: "H5AD",
            h5: "files/datasets/zeisel-brain.h5ad"
        }
    };

    let istate = new inputs.InputsState;
    await istate.compute(files, null, {
        field: "level1class",
        values: [ "microglia", "interneurons" ]
    });

    let expected_num = 388;
    expect(istate.fetchCountMatrix().numberOfColumns()).toEqual(expected_num);
    expect(istate.summary().num_cells).toBe(expected_num);

    {
        let subanno = istate.fetchAnnotations("level1class");
        expect(subanno.length).toBe(expected_num);

        let uniq = Array.from(new Set(subanno));
        uniq.sort();
        expect(uniq).toEqual(["interneurons", "microglia"]);
    }

    // Checking the serialization and unserialization.
    const path = "TEST_subset-inputs.h5";

    let saved = { collected: [], total: 0 };
    {
        let handle = scran.createNewHDF5File(path);
        await istate.serialize(handle, createSaver(saved));
        valkana.validateInputsState(path, true, bakana.kanaFormatVersion);
    }

    let offsets = utils.mockOffsets(saved.collected);
    {
        let handle = new scran.H5File(path);
        let restate = await inputs.unserialize(handle, (offset, size) => offsets[offset]);
        expect(restate.state.fetchCountMatrix().numberOfColumns()).toBe(expected_num);

        let new_params = restate.state.fetchParameters();
        expect(new_params.subset.field).toEqual("level1class");
        expect(new_params.subset.values).toEqual(["microglia", "interneurons"]);

        restate.state.free();
    }

    // Masked and unmasked when the direct subset changes.
    let subset = [1,10,100];
    let expected = {};
    for (const i of subset) {
        expected[i] = istate.fetchCountMatrix().column(i); 
    }

    istate.setDirectSubset(subset);
    expect(istate.summary().num_cells).toBe(subset.length);
    for (const [i, j] of Object.entries(subset)) {
        expect(istate.fetchCountMatrix().column(i)).toEqual(expected[j]);
    }

    istate.setDirectSubset(null);
    expect(istate.summary().num_cells).toBe(expected_num);
    for (const i of subset) {
        expect(istate.fetchCountMatrix().column(i)).toEqual(expected[i]);
    }

    istate.free();
})

test("subsetting processes the blocking factors", async () => {
    let files = { 
        "4K": {
            format: "10X",
            h5: "files/datasets/pbmc4k-tenx.h5"
        },
        "3K": {
            format: "MatrixMarket",
            mtx: "files/datasets/pbmc3k-matrix.mtx.gz",
            genes: "files/datasets/pbmc3k-features.tsv.gz",
            annotations: "files/datasets/pbmc3k-barcodes.tsv.gz"
        }
    };

    let fullstate = new inputs.InputsState;
    await fullstate.compute(files, null, null);

    let ncells = fullstate.fetchCountMatrix().numberOfColumns();
    let subset = [];
    for (var i = 1; i < ncells; i += 2) {
        subset.push(i); // all odds.
    }

    let istate = new inputs.InputsState;
    istate.setDirectSubset(new Int32Array(subset));
    await istate.compute(files, null, null);

    expect(istate.fetchCountMatrix().numberOfColumns()).toBe(subset.length);
    expect(istate.fetchCountMatrix().column(2)).toEqual(fullstate.fetchCountMatrix().column(5));
    expect(istate.fetchCountMatrix().column(5)).toEqual(fullstate.fetchCountMatrix().column(11));

    {
        let subanno = istate.fetchAnnotations("__batch__");
        expect(subanno.length).toBe(subset.length);

        let fullanno = fullstate.fetchAnnotations("__batch__");
        let expected = subset.map(i => fullanno[i]);
        expect(subanno).toEqual(expected);
    }

    {
        let fullblock = fullstate.fetchBlock();
        let fullarr = fullblock.array();
        let expected = subset.map(i => fullarr[i]);

        let subblock = istate.fetchBlock();
        expect(Array.from(subblock.array())).toEqual(expected);
        expect(istate.fetchBlockLevels()).toEqual(fullstate.fetchBlockLevels());
    }

    istate.free();
    fullstate.free();
})

test("end-to-end run works with subsetting", async () => {
    let files = { 
        "default": {
            format: "MatrixMarket",
            mtx: "files/datasets/pbmc-combined-matrix.mtx.gz",
            genes: "files/datasets/pbmc-combined-features.tsv.gz",
            annotations: "files/datasets/pbmc-combined-barcodes.tsv.gz"
        }
    };

    let anno_buffer = fs.readFileSync(files.default.annotations);
    let offset = anno_buffer.byteOffset;
    let buf = anno_buffer.buffer.slice(offset, offset + anno_buffer.byteLength);
    let ncells = bakana.readTable(new Uint8Array(buf), { compression: "gz" }).length - 1;

    let subset = [];
    for (var i = 1; i < ncells; i += 10) { 
        subset.push(i); // every 10th cell.
    }
    
    let state = await bakana.createAnalysis();
    state.inputs.setDirectSubset(subset);
    let params = utils.baseParams();
    params.inputs.sample_factor = "3k";

    let res = await bakana.runAnalysis(state, files, params);
    expect(state.inputs.summary().num_cells).toBe(subset.length);

    let qc_sum = state.quality_control.summary();
    expect(qc_sum.data["3k"].sums.length + qc_sum.data["4k"].sums.length).toBe(subset.length);

    // Saving and loading.
    const path = "TEST_subset-inputs.h5";
    let collected = await bakana.saveAnalysis(state, path);
    utils.validateState(path);

    let offsets = utils.mockOffsets(collected.collected);
    let reloaded = await bakana.loadAnalysis(
        path, 
        (offset, size) => offsets[offset]
    );

    expect(Array.from(state.inputs.fetchDirectSubset())).toEqual(subset);
    expect(reloaded.inputs.summary().num_cells).toEqual(subset.length);

    // subsetInputs works as expected. 
    let refcol = {};
    for (const i of [ 0, 2, 10 ]) {
        refcol[i] = state.inputs.fetchCountMatrix().column(i);
    }

    let subset2 = [];
    for (var i = 0; i < subset.length; i += 2) { 
        subset2.push(i); // every 2nd cell.
    }

    let expected_idx = [];
    for (var i = 0; i < subset.length; i += 2) { 
        expected_idx.push(subset[i]);
    }

    for (var i = 0; i < 2; i++) {
        // We do this twice to check that the subsetting didn't mutate anything
        // in 'state.inputs'.
        let subx = await bakana.subsetInputs(state, subset2);

        let ref = state.inputs.summary();
        let idetails = subx.inputs.summary();
        expect(idetails.num_cells).toEqual(subset2.length);
        expect(idetails.num_genes).toEqual(ref.num_genes);
        expect(idetails.genes).toStrictEqual(ref.genes); // exact same genes.

        for (const [k, v] of Object.entries(refcol)) {
            expect(subx.inputs.fetchCountMatrix().column(k/2)).toEqual(v);
        }

        // Indices were correctly reindexed to be relative to the original cells.
        expect(subx.inputs.fetchDirectSubset()).toEqual(expected_idx);

        // Passing of '#abbreviated' avoids extra work on compute().
        expect(subx.inputs.changed).toBe(false);
        subx.inputs.compute(files, params.inputs.sample_factor, params.inputs.subset);
        expect(subx.inputs.changed).toBe(false);

        await bakana.freeAnalysis(subx);
    }

    await bakana.freeAnalysis(state);
    await bakana.freeAnalysis(reloaded);
})
