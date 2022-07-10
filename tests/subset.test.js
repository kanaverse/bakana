import * as bakana from "../src/index.js";
import * as inputs from "../src/steps/inputs.js";
import * as scran from "scran.js";
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
    await istate.compute(files, null, {
        indices: new Int32Array(subset)
    });

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

    // Checking the serialization and unserialization.
    const path = "TEST_subset-inputs.h5";

    let saved = { collected: [], total: 0 };
    {
        let handle = scran.createNewHDF5File(path);
        await istate.serialize(handle, createSaver(saved));
    }

    let offsets = utils.mockOffsets(saved.collected);
    {
        let handle = new scran.H5File(path);
        let restate = await inputs.unserialize(handle, (offset, size) => offsets[offset]);
        expect(restate.state.fetchCountMatrix().numberOfColumns()).toBe(subset.length);

        let new_params = restate.state.fetchParameters();
        expect(Array.from(new_params.subset.indices)).toEqual(subset);

        restate.state.free();
    }

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
    await istate.compute(files, null, {
        indices: new Int32Array(subset)
    });

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
    let params = utils.baseParams();
    params.inputs.subset = { "indices": subset };
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

    let new_params = bakana.retrieveParameters(reloaded);
    expect(Array.from(new_params.inputs.subset.indices)).toEqual(subset);
    expect(reloaded.inputs.summary().num_cells).toEqual(subset.length);

    // Retrieval can also skip the indices.
    {
        let new_params = bakana.retrieveParameters(reloaded, { wipeIndices: true });
        expect(new_params.inputs.subset).toBeNull();
    }

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
        let params = subx.inputs.fetchParameters();
        expect(params.subset.indices).toEqual(expected_idx);

        // Passing of '#abbreviated' avoids extra work on compute().
        expect(subx.inputs.changed).toBe(false);
        subx.inputs.compute(files, params.sample_factor, params.subset);
        expect(subx.inputs.changed).toBe(false);

        await bakana.freeAnalysis(subx);
    }

    await bakana.freeAnalysis(state);
    await bakana.freeAnalysis(reloaded);
})

test("subset equality checks behave correctly", () => {
    expect(inputs.InputsState.areSubsetsEqual(undefined, null)).toBe(false);
    expect(inputs.InputsState.areSubsetsEqual(null, null)).toBe(true);
    expect(inputs.InputsState.areSubsetsEqual(null, {})).toBe(false);
    expect(inputs.InputsState.areSubsetsEqual(undefined, {})).toBe(false);

    expect(inputs.InputsState.areSubsetsEqual({}, {})).toBe(true);
    expect(inputs.InputsState.areSubsetsEqual({fields:"A", values:["B"]}, {fields:"A", values:["B"]})).toBe(true);
    expect(inputs.InputsState.areSubsetsEqual({fields:"A", values:["B"]}, {fields:"C", values:["B"]})).toBe(false);

    expect(inputs.InputsState.areSubsetsEqual({fields:"A", values:["B"]}, {indices:[1,2,3]})).toBe(false);
    expect(inputs.InputsState.areSubsetsEqual({indices:[1,2,3]}, {indices:[1,2,3]})).toBe(true);
    expect(inputs.InputsState.areSubsetsEqual({indices:[1,2,3]}, {indices:[1,2,3,4]})).toBe(false);

    expect(inputs.InputsState.areSubsetsEqual({indices:[1,2,3]}, {indices:[1,2,3], foo:1})).toBe(false);
    expect(inputs.InputsState.areSubsetsEqual({indices:[1,2,3], foo:1}, {indices:[1,2,3], foo:1})).toBe(true);
})

