import * as bakana from "../src/index.js";
import * as inputs from "../src/steps/inputs.js";
import * as scran from "scran.js";
import * as utils from "./utils.js"

beforeAll(async () => await bakana.initialize({ localFile: true }));
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
        let reloaded = handle.open("inputs").open("parameters").open("subset").open("indices", { load: true }).values;
        expect(Array.from(reloaded)).toEqual(subset);

        let restate = await inputs.unserialize(handle, (offset, size) => offsets[offset]);
        expect(restate.state.fetchCountMatrix().numberOfColumns()).toBe(subset.length);
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

        let subhandle = handle.open("inputs").open("parameters").open("subset");
        expect(subhandle.open("field", { load: true }).values[0]).toEqual("level1class");
        expect(subhandle.open("values", { load: true }).values).toEqual(["microglia", "interneurons"]);

        let restate = await inputs.unserialize(handle, (offset, size) => offsets[offset]);
        expect(restate.state.fetchCountMatrix().numberOfColumns()).toBe(expected_num);
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
