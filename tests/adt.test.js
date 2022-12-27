import * as bakana from "../src/index.js";
import * as utils from "./utils.js";
import * as scran from "scran.js";
import * as valkana from "valkana";
import * as bioc from "bioconductor";
import * as fs from "fs";
import * as combine from "../src/steps/combine_embeddings.js";

beforeAll(utils.initializeAll);
afterAll(async () => await bakana.terminate());

let mtx = "files/datasets/immune_3.0.0-matrix.mtx.gz";
let feats = "files/datasets/immune_3.0.0-features.tsv.gz";
let bars = "files/datasets/immune_3.0.0-barcodes.tsv.gz";
let files = { default: new bakana.TenxMatrixMarketDataset(mtx, feats, bars) };

async function overlord(state) {
    await utils.checkStateResultsMinimal(state);
    await utils.checkStateResultsRna(state);
    await utils.checkStateResultsAdt(state);
    await utils.checkStateResultsRnaPlusAdt(state);
}

test("ADT MatrixMarket summary works correctly", async () => {
    let summ = await files.default.summary();
    expect(summ.modality_features["Gene Expression"] instanceof bioc.DataFrame).toBe(true);
    expect(summ.modality_features["Gene Expression"].numberOfColumns()).toBeGreaterThan(0);
    expect(summ.modality_features["Antibody Capture"] instanceof bioc.DataFrame).toBe(true);
    expect(summ.modality_features["Antibody Capture"].numberOfColumns()).toBeGreaterThan(0);
    expect(summ.cells instanceof bioc.DataFrame).toBe(true);
    expect(summ.cells.numberOfColumns()).toBeGreaterThan(0);
})

test("runAnalysis works correctly (MatrixMarket)", async () => {
    let state = await bakana.createAnalysis();
    let params = utils.baseParams();
    let res = await bakana.runAnalysis(state, files, params);

    // Input reorganization is done correctly. 
    {
        let simple = scran.initializeSparseMatrixFromMatrixMarket(mtx, { layered: false });
        let parsed = bakana.readTable((new bakana.SimpleFile(feats)).buffer(), { compression: "gz" });
        let simple_names_all = parsed.map(x => x[0]);
        let simple_types = parsed.map(x => x[2]);

        const regexes = { "RNA": /Gene Expression/i, "ADT": /Antibody Capture/i };
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

    // Running through some basic checks.
    await overlord(state);
    await utils.checkStateResultsUnblocked(state);

    let vres = utils.checkClusterVersusMode(state);
    expect(vres.results.ADT.numberOfGroups()).toEqual(vres.results.RNA.numberOfGroups());

    let custom = utils.launchCustomSelections(state);
    expect(custom.first.ADT.numberOfGroups()).toEqual(custom.first.RNA.numberOfGroups());
    expect(custom.last.ADT.numberOfGroups()).toEqual(custom.last.RNA.numberOfGroups());
    expect(custom.versus.results.ADT.numberOfGroups()).toEqual(custom.versus.results.RNA.numberOfGroups());

    // Saving and loading.
    const path = "TEST_state_adt.h5";
    let collected = await bakana.saveAnalysis(state, path);
    utils.validateState(path);
    expect(collected.collected.length).toBe(3);
    expect(typeof(collected.collected[0])).toBe("string");

    let offsets = utils.mockOffsets(collected.collected);
    let reloaded = await bakana.loadAnalysis(
        path, 
        (offset, size) => offsets[offset]
    );

    let new_params = bakana.retrieveParameters(reloaded);
    expect(new_params).toEqual(params);

    await utils.compareStates(state, reloaded, { checkAdt: true });

    // Release me!
    await bakana.freeAnalysis(state);
    await bakana.freeAnalysis(reloaded);
})

test("runAnalysis works correctly (10X)", async () => {
    let state = await bakana.createAnalysis();
    let params = utils.baseParams();
    let res = await bakana.runAnalysis(state, 
        { default: new bakana.TenxHdf5Dataset("files/datasets/immune_3.0.0-tenx.h5") },
        params
    );

    await overlord(state);
    await utils.checkStateResultsUnblocked(state);

    // What happens when only one of the modalities has non-zero weight?
    {
        state.combine_embeddings.compute(/* rna_weight = */ 1, /* adt_weight = */ 0, /* crispr_weight = */ 0, true);
        let pcs = state.combine_embeddings.fetchCombined();
        expect(pcs.owner !== null).toBe(true);
        expect(state.combine_embeddings.fetchNumberOfDimensions()).toBe(state.rna_pca.fetchPCs().numberOfPCs());

        const path = "TEST_state_combine-embed.h5";
        let fhandle = scran.createNewHDF5File(path);
        state.combine_embeddings.serialize(fhandle);

        {
            let ncells = state.cell_filtering.fetchFilteredMatrix().numberOfColumns();
            let npcs_rna = state.rna_pca.fetchPCs().numberOfPCs();
            let npcs_adt = state.adt_pca.fetchPCs().numberOfPCs();
            valkana.validateCombineEmbeddingsState(path, ncells, ["RNA"], npcs_rna + npcs_adt, bakana.kanaFormatVersion);
        }

        let reloaded = combine.unserialize(fhandle, {"RNA": state.rna_pca, "ADT": state.adt_pca, "CRISPR": state.crispr_pca });
        let repcs = reloaded.fetchCombined();
        expect(repcs.owner !== null).toBe(true);
        expect(reloaded.fetchNumberOfDimensions()).toBe(state.rna_pca.fetchPCs().numberOfPCs());

        reloaded.free();
    }

    // Release me!
    await bakana.freeAnalysis(state);
})

test("runAnalysis works for ADTs with blocking", async () => {
    let state = await bakana.createAnalysis();
    let params = utils.baseParams();

    // Mocking up a blocking file with pretend batches.
    let exfile = "TEST_adt_block.tsv";
    let nblocks = 3;
    {
        let f = fs.readFileSync(bars);
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
        { "combined": new bakana.TenxMatrixMarketDataset(mtx, feats, exfile) },
        params
    );

    await overlord(state);
    await utils.checkStateResultsBlocked(state);

    let vres = utils.checkClusterVersusMode(state);
    expect(vres.results.RNA.numberOfBlocks()).toEqual(nblocks);
    expect(vres.results.ADT.numberOfBlocks()).toEqual(nblocks);

    let custom = utils.launchCustomSelections(state);
    expect(custom.first.ADT.numberOfBlocks()).toEqual(nblocks);
    expect(custom.last.ADT.numberOfBlocks()).toEqual(nblocks);
    expect(custom.versus.results.ADT.numberOfBlocks()).toBeGreaterThan(1); // as subset might not actually have all 3 blocks.

    // Freeing everyone.
    await bakana.freeAnalysis(state);
})

test("ADT-only runAnalysis works correctly", async () => {
    let state = await bakana.createAnalysis();
    let params = utils.baseParams();
    let files = { default: new bakana.TenxMatrixMarketDataset(mtx, feats, bars, { featureTypeRnaName: null }) };
    await bakana.runAnalysis(state, files, params);

    // Check the results.
    await utils.checkStateResultsAdt(state, { exclusive: true });
    await utils.checkStateResultsUnblocked(state);

    // Can save and reload.
    const path = "TEST_state_adt_only.h5";
    let collected = await bakana.saveAnalysis(state, path);
    //utils.validateState(path);
    expect(collected.collected.length).toBe(3);
    expect(typeof(collected.collected[0])).toBe("string");

    let offsets = utils.mockOffsets(collected.collected);
    let reloaded = await bakana.loadAnalysis(
        path, 
        (offset, size) => offsets[offset]
    );

    await utils.compareStates(state, reloaded, { checkRna: false, checkAdt: true });

    await bakana.freeAnalysis(state);
    await bakana.freeAnalysis(reloaded);
})

