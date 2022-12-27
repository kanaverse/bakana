import * as bakana from "../src/index.js";
import * as utils from "./utils.js";
import * as scran from "scran.js";
import * as valkana from "valkana";
import * as bioc from "bioconductor";
import * as fs from "fs";
import * as combine from "../src/steps/combine_embeddings.js";

beforeAll(utils.initializeAll);
afterAll(async () => await bakana.terminate());

const h5file = "files/datasets/crispr_6.0.0-tenx.h5";
let files = { default: new bakana.TenxHdf5Dataset(h5file) };

test("CRISPR MatrixMarket summary works correctly", async () => {
    let summ = await files.default.summary();
    expect(summ.modality_features["Gene Expression"] instanceof bioc.DataFrame).toBe(true);
    expect(summ.modality_features["Gene Expression"].numberOfColumns()).toBeGreaterThan(0);
    expect(summ.modality_features["CRISPR Guide Capture"] instanceof bioc.DataFrame).toBe(true);
    expect(summ.modality_features["CRISPR Guide Capture"].numberOfColumns()).toBeGreaterThan(0);
    expect(summ.cells instanceof bioc.DataFrame).toBe(true);
    expect(summ.cells.numberOfRows()).toBeGreaterThan(0);
})

test("runAnalysis works correctly (10X)", async () => {
    let state = await bakana.createAnalysis();
    let params = utils.baseParams();
    params.combine_embeddings.crispr_weight = 1; // turning it on for PCA.
    await bakana.runAnalysis(state, files, params);

    // Input reorganization is done correctly. 
    {
        let simple = scran.initializeSparseMatrixFromHDF5(h5file, "matrix", { layered: false });
        let fhandle = (new scran.H5File(h5file)).open("matrix").open("features");
        let simple_names_all = fhandle.open("id", { load: true }).values;
        let simple_types = fhandle.open("feature_type", { load: true }).values;

        const regexes = { "RNA": /Gene Expression/i, "CRISPR": /CRISPR Guide Capture/i };
        for (const [k, v] of Object.entries(regexes)) {
            let keep = [];
            simple_types.forEach((x, i) => {
                if (x.match(v)) {
                    keep.push(i);
                }
            });

            let simple_names = keep.map(x => simple_names_all[x]);
            let simple_mat = scran.subsetRows(simple.matrix, keep);
            let simple_ids = keep.slice();

            let loaded = state.inputs.fetchCountMatrix().get(k);
            let loaded_ids = state.inputs.fetchRowIds()[k];
            let loaded_names = state.inputs.fetchFeatureAnnotations()[k].column("id");

            expect(simple.matrix.numberOfRows()).toBeGreaterThan(loaded.numberOfRows());
            utils.checkReorganization(simple_mat, simple_ids, simple_names, loaded, loaded_ids, loaded_names, { referenceSubset: true }); 
            simple_mat.free();
        }

        simple.matrix.free();
    }

    await utils.checkStateResultsMinimal(state);
    await utils.checkStateResultsUnblocked(state);
    await utils.checkStateResultsRna(state);
    await utils.checkStateResultsCrispr(state);
    await utils.checkStateResultsRnaPlusCrispr(state);

    // Saving and loading.
    const path = "TEST_state_crispr.h5";
    let collected = await bakana.saveAnalysis(state, path);
    // utils.validateState(path); // TODO: kanaval doesn't know anything about CRISPR yet.
    expect(collected.collected.length).toBe(1);
    expect(typeof(collected.collected[0])).toBe("string");

    let offsets = utils.mockOffsets(collected.collected);
    let reloaded = await bakana.loadAnalysis(
        path, 
        (offset, size) => offsets[offset]
    );

    let new_params = bakana.retrieveParameters(reloaded);
    expect(new_params).toEqual(params);

    await utils.compareStates(state, reloaded, { checkCrispr: true });

    // Release me!
    await bakana.freeAnalysis(state);
    await bakana.freeAnalysis(reloaded);
})

const mtxfile = "files/datasets/crispr_6.0.0-matrix.mtx.gz";
const featfile = "files/datasets/crispr_6.0.0-features.tsv.gz";
const barfile = "files/datasets/crispr_6.0.0-barcodes.tsv.gz";

test("runAnalysis works for CRISPR (MatrixMarket) with blocking", async () => {
    let state = await bakana.createAnalysis();
    let params = utils.baseParams();

    // Mocking up a blocking file with pretend batches.
    let exfile = "TEST_crispr_block.tsv";
    let nblocks = 3;
    {
        let f = fs.readFileSync(barfile);
        let buff = f.buffer.slice(f.byteOffset, f.byteOffset + f.byteLength);
        let stuff = bakana.readTable(new Uint8Array(buff));

        let ncells = stuff.length;
        let per_block = Math.ceil(ncells / nblocks);
        let blocks = new Array(ncells);
        for (var c = 0; c < ncells; c++) {
            blocks[c] = 'A' + String(Math.floor(c / per_block));
        }

        fs.writeFileSync(exfile, blocks.join("\n"));
    }
    params.inputs.sample_factor = "A0";

    let res = await bakana.runAnalysis(state, 
        { "combined": new bakana.TenxMatrixMarketDataset(mtxfile, featfile, exfile) },
        params
    );

    await utils.checkStateResultsMinimal(state);
    await utils.checkStateResultsBlocked(state);
    await utils.checkStateResultsRna(state);
    await utils.checkStateResultsCrispr(state, { use_embeddings: false }); // CRISPR is not used in the embeddings by default.
    await utils.checkStateResultsRnaPlusCrispr(state, { use_embeddings: false });

    // However, CRISPR should still show up in the marker analyses.
    let vres = utils.checkClusterVersusMode(state);
    expect(vres.results.RNA.numberOfBlocks()).toEqual(nblocks);
    expect(vres.results.CRISPR.numberOfBlocks()).toEqual(nblocks);

    let custom = utils.launchCustomSelections(state);
    expect(custom.first.CRISPR.numberOfBlocks()).toEqual(nblocks);
    expect(custom.last.CRISPR.numberOfBlocks()).toEqual(nblocks);
    expect(custom.versus.results.CRISPR.numberOfBlocks()).toBeGreaterThan(1); // as subset might not actually have all 3 blocks.

    // Freeing everyone.
    await bakana.freeAnalysis(state);
})

test("CRISPR-only runAnalysis works correctly", async () => {
    let state = await bakana.createAnalysis();
    let params = utils.baseParams();
    params.combine_embeddings.crispr_weight = 1;
    let files = { default: new bakana.TenxHdf5Dataset(h5file, { featureTypeRnaName: null }) };
    await bakana.runAnalysis(state, files, params);

    // Check the results.
    await utils.checkStateResultsMinimal(state);
    await utils.checkStateResultsUnblocked(state);
    await utils.checkStateResultsCrispr(state, { exclusive: true });

    // Can save and reload.
    const path = "TEST_state_crispr_only.h5";
    let collected = await bakana.saveAnalysis(state, path);
    //utils.validateState(path);
    expect(collected.collected.length).toBe(1);
    expect(typeof(collected.collected[0])).toBe("string");

    let offsets = utils.mockOffsets(collected.collected);
    let reloaded = await bakana.loadAnalysis(
        path, 
        (offset, size) => offsets[offset]
    );

    await utils.compareStates(state, reloaded, { checkRna: false, checkAdt: false, checkCrispr: true });

    await bakana.freeAnalysis(state);
    await bakana.freeAnalysis(reloaded);
})

