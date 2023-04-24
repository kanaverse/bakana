import * as bakana from "../src/index.js";
import * as utils from "./utils.js";
import * as scran from "scran.js";
import * as bioc from "bioconductor";
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

    let preview = await files.default.previewPrimaryIds();
    expect("RNA" in preview).toBe(true);
    expect("ADT" in preview).toBe(true);
    expect(preview.RNA.length).toBeGreaterThan(0);
    expect(preview.ADT.length).toBeGreaterThan(0);
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

    // Check saving of results.
    await bakana.saveSingleCellExperiment(state, "adt", { directory: "miscellaneous/from-tests" });
    await bakana.saveSingleCellExperiment(state, "adt_split", { directory: "miscellaneous/from-tests", storeModalityColumnData: true });
    await bakana.saveGenewiseResults(state, "adt_genes", { directory: "miscellaneous/from-tests" });

    // Release me!
    await bakana.freeAnalysis(state);
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
        params.combine_embeddings.rna_weight = 1;
        params.combine_embeddings.adt_weight = 0;
        params.combine_embeddings.crispr_weight = 0;

        await bakana.runAnalysis(state, null, params);
        expect(state.rna_pca.changed).toBe(false);
        expect(state.adt_pca.changed).toBe(false);
        expect(state.crispr_pca.changed).toBe(false);
        expect(state.combine_embeddings.changed).toBe(true);
        expect(state.snn_graph_cluster.changed).toBe(true); // downstream is changed.

        // We should end up with a view.
        let pcs = state.combine_embeddings.fetchCombined();
        expect(pcs.owner !== null).toBe(true);
        expect(state.combine_embeddings.fetchNumberOfDimensions()).toBe(state.rna_pca.fetchPCs().numberOfPCs());
    }

    // Checking that we can successfully turn off the bias removal.
    {
        let old_sf = state.adt_normalization.fetchSizeFactors().slice();
        params.adt_normalization.remove_bias = false;
        state.adt_normalization.compute(params.adt_normalization);
        expect(state.adt_normalization.changed).toBe(true);

        let new_sf = state.adt_normalization.fetchSizeFactors().slice();
        expect(new_sf.length).toEqual(old_sf.length);
        expect(new_sf).not.toEqual(old_sf);

        // Parameter changes have no effect when bias removal is disabled.
        params.adt_normalization.num_pcs = 50;
        state.adt_normalization.compute(params.adt_normalization);
        expect(state.adt_normalization.changed).toBe(false);
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
    params.inputs.block_factor = utils.mockBlocks(bars, exfile, nblocks);

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
    let files = { default: new bakana.TenxMatrixMarketDataset(mtx, feats, bars) };
    files.default.setOptions({ featureTypeRnaName: null }); 
    await bakana.runAnalysis(state, files, params);

    // Check the results.
    await utils.checkStateResultsAdt(state, { exclusive: true });
    await utils.checkStateResultsUnblocked(state);

    // Freeing.
    await bakana.freeAnalysis(state);
})

