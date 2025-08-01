import * as bakana from "../src/index.js";
import * as inputs from "../src/steps/inputs.js";
import * as scran from "scran.js";
import * as bioc from "bioconductor";
import * as utils from "./utils.js"
import * as fs from "fs";

beforeAll(utils.initializeAll);
afterAll(async () => await bakana.terminate());

let h5path = "files/datasets/zeisel-brain.h5ad";

test("subsetting behaves correctly with indices", async () => {
    let files = { default: new bakana.H5adDataset(h5path) };

    let fullstate = new inputs.InputsState;
    let idefaults = inputs.InputsState.defaults();
    await fullstate.compute(files, idefaults);
    
    let ncells = fullstate.fetchCountMatrix().numberOfColumns();
    let subset = [];
    for (var i = 0; i < ncells; i += 2) {
        subset.push(i); // all evens.
    }

    let istate = new inputs.InputsState;
    istate.setDirectSubset(new Int32Array(subset));
    await istate.compute(files, idefaults);

    expect(istate.fetchCountMatrix().numberOfColumns()).toBe(subset.length);
    expect(istate.fetchCountMatrix().get("RNA").column(2)).toEqual(fullstate.fetchCountMatrix().get("RNA").column(4));
    expect(istate.fetchCountMatrix().get("RNA").column(5)).toEqual(fullstate.fetchCountMatrix().get("RNA").column(10));

    {
        let indices = [0,2,4,8];
        fullstate.undoSubset(indices);
        expect(indices).toEqual([0,2,4,8]);
        istate.undoSubset(indices);
        expect(indices).toEqual([0,4,8,16]);
    }

    {
        let subanno = istate.fetchCellAnnotations().column("level1class");
        expect(subanno.length).toBe(subset.length);

        let fullanno = fullstate.fetchCellAnnotations().column("level1class");
        let expected = subset.map(i => fullanno[i]);
        expect(subanno).toEqual(expected);
    }

    // Changes to the subset is a no-op when indices are around.
    {
        expect(istate.changed).toBe(true);
        await istate.compute(files, { block_factor: null, subset: { field: "level1class", values: ["FOO", "BLAH"] } });
        expect(istate.changed).toBe(false);
    }

    // Unsetting works as expected.
    istate.setDirectSubset(null);
    istate.compute(null, idefaults);
    expect(istate.changed).toBe(true);
    expect(istate.fetchCountMatrix().numberOfColumns()).toBeGreaterThan(3000);
    expect(istate.fetchCountMatrix().get("RNA").column(1)).toEqual(fullstate.fetchCountMatrix().get("RNA").column(1));

    // Resetting works as expected.
    istate.setDirectSubset(subset);
    istate.compute(null, idefaults);
    expect(istate.changed).toBe(true);
    expect(istate.fetchCountMatrix().numberOfColumns()).toBe(subset.length);
    expect(istate.fetchCountMatrix().get("RNA").column(3)).toEqual(fullstate.fetchCountMatrix().get("RNA").column(6));

    istate.free();
    fullstate.free();
})

test("subsetting behaves correctly with a factor", async () => {
    let files = { default: new bakana.H5adDataset(h5path) };

    let istate = new inputs.InputsState;
    let sub =  {
        field: "level1class",
        values: [ "microglia", "interneurons" ]
    };
    let iparams = { block_factor: null, subset: sub };
    await istate.compute(files, iparams);

    let expected_num = 0;
    {
        let handle = new scran.H5File(h5path);
        let anno = handle.open("obs").open("level1class", { load: true }).values;
        let levels = handle.open("obs").open("__categories").open("level1class", { load: true }).values;
        for (const x of anno) {
            if (levels[x] == "microglia" || levels[x] == "interneurons") {
                expected_num++;
            }
        }
    }

    expect(expected_num).toBeGreaterThan(0);
    expect(istate.fetchCountMatrix().numberOfColumns()).toEqual(expected_num);

    {
        let subanno = istate.fetchCellAnnotations().column("level1class");
        expect(subanno.length).toBe(expected_num);

        let uniq = Array.from(new Set(subanno));
        uniq.sort();
        expect(uniq).toEqual(["interneurons", "microglia"]);
    }

    // Masked and unmasked when the direct subset changes.
    let direct_subset = [1,10,100];
    let expected = {};
    for (const i of direct_subset) {
        expected[i] = istate.fetchCountMatrix().get("RNA").column(i); 
    }

    istate.setDirectSubset(direct_subset);
    istate.compute(files, iparams);
    expect(istate.fetchCountMatrix().numberOfColumns()).toBe(direct_subset.length);
    for (const [i, j] of Object.entries(direct_subset)) {
        expect(istate.fetchCountMatrix().get("RNA").column(i)).toEqual(expected[j]);
    }

    istate.setDirectSubset(null);
    istate.compute(files, iparams);
    expect(istate.fetchCountMatrix().numberOfColumns()).toBe(expected_num);
    for (const i of direct_subset) {
        expect(istate.fetchCountMatrix().get("RNA").column(i)).toEqual(expected[i]);
    }

    istate.free();
})

test("subsetting behaves correctly with ranges", async () => {
    let files = { default: new bakana.H5adDataset(h5path) };

    let istate = new inputs.InputsState;
    await istate.compute(files, { 
        block_factor: null, 
        subset: {
            field: "diameter",
            ranges: [ [5, 10], [20, 25], [100, 200] ]
        }
    });

    let expected_num = 0;
    {
        let handle = new scran.H5File(h5path);
        let diameter = handle.open("obs").open("diameter", { load: true }).values;
        diameter.forEach((x, i) => {
            if ((x >= 5 && x <= 10) || (x >= 20 && x <= 25) || (x >= 100 && x <= 200)) {
                expected_num++;
            }
        });
    }

    expect(expected_num).toBeGreaterThan(0);
    expect(istate.fetchCountMatrix().numberOfColumns()).toEqual(expected_num);

    {
        let subanno = istate.fetchCellAnnotations().column("diameter");
        expect(subanno.length).toBe(expected_num);

        let failed = 0;
        for (const x of subanno) {
            if (!(x >= 5 && x <= 10) && !(x >= 20 && x <= 25) && !(x >= 100 && x <= 200)) {
                failed++;
            }
        }
        expect(failed).toBe(0);
    }
})

test("subsetting behaves correctly with infinite ranges", async () => {
    let files = { default: new bakana.H5adDataset(h5path) };

    let istate = new inputs.InputsState;
    await istate.compute(files, {
        block_factor: null, 
        subset: {
            field: "diameter",
            ranges: [ [-Infinity, 10], [100, Infinity] ]
        }
    });
    expect(istate.fetchParameters().subset.ranges[0][0]).toEqual(-Infinity);
    expect(istate.fetchParameters().subset.ranges[1][1]).toEqual(Infinity);

    let expected_num = 0;
    {
        let handle = new scran.H5File(h5path);
        let diameter = handle.open("obs").open("diameter", { load: true }).values;
        diameter.forEach((x, i) => {
            if (x >= 100 || x <= 10) {
                expected_num++;
            }
        });
    }

    expect(expected_num).toBeGreaterThan(0);
    expect(istate.fetchCountMatrix().numberOfColumns()).toEqual(expected_num);

    istate.free();
})

test("subsetting processes the blocking factors", async () => {
    let files = { 
        "4K": new bakana.TenxHdf5Dataset("files/datasets/pbmc4k-tenx.h5"),
        "3K": new bakana.TenxMatrixMarketDataset(
                "files/datasets/pbmc3k-matrix.mtx.gz",
                "files/datasets/pbmc3k-features.tsv.gz",
                "files/datasets/pbmc3k-barcodes.tsv.gz"
            )
    };

    let fullstate = new inputs.InputsState;
    let idefaults = inputs.InputsState.defaults();
    await fullstate.compute(files, idefaults);

    let ncells = fullstate.fetchCountMatrix().numberOfColumns();
    let subset = [];
    for (var i = 1; i < ncells; i += 2) {
        subset.push(i); // all odds.
    }

    let istate = new inputs.InputsState;
    istate.setDirectSubset(new Int32Array(subset));
    await istate.compute(files, idefaults);

    expect(istate.fetchCountMatrix().numberOfColumns()).toBe(subset.length);
    expect(istate.fetchCountMatrix().get("RNA").column(2)).toEqual(fullstate.fetchCountMatrix().get("RNA").column(5));
    expect(istate.fetchCountMatrix().get("RNA").column(5)).toEqual(fullstate.fetchCountMatrix().get("RNA").column(11));

    {
        let subanno = istate.fetchCellAnnotations().column("__batch__");
        expect(subanno.length).toBe(subset.length);

        let fullanno = fullstate.fetchCellAnnotations().column("__batch__");
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

test("subsetting behaves with invalid values in the blocking factors", async () => {
    let annopath = "files/datasets/pbmc-combined-barcodes.tsv.gz";
    let files = { 
        "default": new bakana.TenxMatrixMarketDataset(
                "files/datasets/pbmc-combined-matrix.mtx.gz",
                "files/datasets/pbmc-combined-features.tsv.gz",
                annopath
            )
    };

    let deets = await files.default.load({ cache: true });
    let copy = deets.cells.column("3k").slice();
    copy[0] = null;
    deets.cells.$setColumn("dummy", copy);

    let istate = new inputs.InputsState;
    let params = inputs.InputsState.defaults();
    params.block_factor = "dummy";
    await istate.compute(files, params);

    // first column has been removed.
    expect(istate.fetchCountMatrix().numberOfColumns()).toBe(deets.cells.numberOfRows() - 1);
    expect(istate.fetchCellAnnotations().column("3k")).toEqual(deets.cells.column("3k").slice(1));
    expect(istate.fetchCountMatrix().get("RNA").column(0)).toEqual(deets.matrix.get("RNA").column(1)); 

    let sub = [ 0, 1, 2 ];
    istate.undoSubset(sub);
    expect(sub).toEqual([1,2,3]);
    expect(istate.fetchBlock().length).toEqual(deets.cells.numberOfRows() - 1);

    expect(istate.fetchFeatureAnnotations()["RNA"].numberOfRows()).toEqual(deets.features["RNA"].numberOfRows()); // RNA is unaffected.

    // Interacts correctly with additional subsetting; the first column is still dropped.
    istate.setDirectSubset(new Int32Array([0,2,4,6,8,10]), { onOriginal: true });
    istate.compute(files, params);

    expect(istate.fetchCountMatrix().numberOfColumns()).toBe(5);
    expect(istate.fetchCellAnnotations().column("3k")).toEqual(bioc.SLICE(deets.cells.column("3k"), [2, 4, 6, 8, 10]));
    expect(istate.fetchCountMatrix().get("RNA").column(0)).toEqual(deets.matrix.get("RNA").column(2));
    expect(istate.fetchCountMatrix().get("RNA").column(1)).toEqual(deets.matrix.get("RNA").column(4));

    sub = [ 0, 1, 2 ];
    istate.undoSubset(sub);
    expect(sub).toEqual([2, 4, 6]);

    istate.free();
})

test("end-to-end run works with (direct) subsetting", async () => {
    let annopath = "files/datasets/pbmc-combined-barcodes.tsv.gz";
    let files = { 
        "default": new bakana.TenxMatrixMarketDataset(
                "files/datasets/pbmc-combined-matrix.mtx.gz",
                "files/datasets/pbmc-combined-features.tsv.gz",
                annopath
            )
    };

    let anno_buffer = fs.readFileSync(annopath);
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
    params.inputs.block_factor = "3k";

    let res = await bakana.runAnalysis(state, files, params);
    expect(state.inputs.fetchCountMatrix().numberOfColumns()).toBe(subset.length);
    expect(state.rna_quality_control.fetchKeep().length).toBe(subset.length);

    // Basic checks.
    await utils.checkStateResultsMinimal(state);
    await utils.checkStateResultsRna(state, { exclusive: true });
    await utils.checkStateResultsBlocked(state);

    // Serialization and reloading works as expected.
    {
        let saved = [];
        let saver = (n, k, f) => {
            saved.push(f.content());
            return String(saved.length);
        };

        let serialized = await bakana.serializeConfiguration(state, saver);
        expect(serialized.parameters).toEqual(bakana.retrieveParameters(state));
        expect(serialized.other.inputs.direct_subset).toEqual(subset);

        let reloaded = await bakana.unserializeConfiguration(serialized, x => saved[Number(x) - 1]);
        await utils.compareStates(reloaded, state);
        await bakana.freeAnalysis(reloaded);
    }

    // subsetInputs works as expected. 
    let refcol = {};
    for (const i of [ 0, 2, 10 ]) {
        refcol[i] = state.inputs.fetchCountMatrix().get("RNA").column(i);
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

        // Exact same genes.
        let mat = state.inputs.fetchCountMatrix();
        let submat = subx.inputs.fetchCountMatrix();
        expect(mat.available()).toEqual(submat.available());
        expect(mat.get("RNA").numberOfRows()).toEqual(submat.get("RNA").numberOfRows());

        let anno = state.inputs.fetchFeatureAnnotations()["RNA"];
        let subanno = subx.inputs.fetchFeatureAnnotations()["RNA"];
        expect(anno.numberOfRows()).toEqual(subanno.numberOfRows());
        expect(anno.column(0)).toEqual(subanno.column(0));

        for (const [k, v] of Object.entries(refcol)) {
            expect(subx.inputs.fetchCountMatrix().get("RNA").column(k/2)).toEqual(v);
        }

        // Indices were correctly reindexed to be relative to the original cells.
        expect(subx.inputs.fetchDirectSubset()).toEqual(expected_idx);

        // Passing of '#abbreviated' avoids extra work on compute().
        expect(subx.inputs.changed).toBe(false);
        subx.inputs.compute(files, params.inputs);
        expect(subx.inputs.changed).toBe(false);

        await bakana.freeAnalysis(subx);
    }

    await bakana.freeAnalysis(state);
})
