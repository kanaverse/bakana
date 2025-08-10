import * as bakana from "../src/index.js";
import * as utils from "./utils.js";
import * as scran from "scran.js";
import * as bioc from "bioconductor";
import * as combine from "../src/steps/combine_embeddings.js";

beforeAll(utils.initializeAll);
afterAll(async () => await bakana.terminate());

const h5file = "files/datasets/crispr_6.0.0-tenx.h5";
let files = { default: new bakana.TenxHdf5Dataset(h5file) };

async function overlord(state, use_embeddings) {
    await utils.checkStateResultsMinimal(state);
    await utils.checkStateResultsRna(state);
    await utils.checkStateResultsCrispr(state, { use_embeddings });
    await utils.checkStateResultsRnaPlusCrispr(state, { use_embeddings });
}

test("CRISPR MatrixMarket summary works correctly", async () => {
    let summ = await files.default.summary();
    expect(summ.modality_features["Gene Expression"] instanceof bioc.DataFrame).toBe(true);
    expect(summ.modality_features["Gene Expression"].numberOfColumns()).toBeGreaterThan(0);
    expect(summ.modality_features["CRISPR Guide Capture"] instanceof bioc.DataFrame).toBe(true);
    expect(summ.modality_features["CRISPR Guide Capture"].numberOfColumns()).toBeGreaterThan(0);
    expect(summ.cells instanceof bioc.DataFrame).toBe(true);
    expect(summ.cells.numberOfRows()).toBeGreaterThan(0);

    let preview = files.default.previewPrimaryIds({ cache: true });
    expect("RNA" in preview).toBe(true);
    expect("CRISPR" in preview).toBe(true);
    expect(preview.RNA.length).toBeGreaterThan(0);
    expect(preview.CRISPR.length).toBeGreaterThan(0);
})

test("runAnalysis works correctly (10X)", async () => {
    let state = await bakana.createAnalysis();
    let params = utils.baseParams();
    params.combine_embeddings.crispr_weight = 1; // turning it on for PCA.
    await bakana.runAnalysis(state, files, params);

    // Input reorganization is done correctly. 
    {
        let simple = scran.initializeSparseMatrixFromHdf5(h5file, "matrix", { layered: false });
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
            let simple_mat = scran.subsetRows(simple, keep);

            let loaded = state.inputs.fetchCountMatrix().get(k);
            let loaded_names = state.inputs.fetchFeatureAnnotations()[k].column("id");

            expect(simple.numberOfRows()).toBeGreaterThan(loaded.numberOfRows());
            utils.checkMatrixContents(simple_mat, simple_names, loaded, loaded_names);
            simple_mat.free();
        }

        simple.free();
    }

    // Simple checks, where CRISPR is used in the combined embeddings.
    await overlord(state, true);
    await utils.checkStateResultsUnblocked(state);

    // Check saving of results.
    utils.purgeDirectory("miscellaneous/from-tests/crispr");
    await bakana.saveSingleCellExperiment(state, "crispr", { directory: "miscellaneous/from-tests" });
    await utils.checkSavedExperiment("miscellaneous/from-tests/crispr", state);
    utils.purgeDirectory("miscellaneous/from-tests/crispr_genes");
    await bakana.saveGenewiseResults(state, "crispr_genes", { directory: "miscellaneous/from-tests" });

    // Release me!
    await bakana.freeAnalysis(state);
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
    params.inputs.block_factor = utils.mockBlocks(barfile, exfile, nblocks);

    let res = await bakana.runAnalysis(state, 
        { "combined": new bakana.TenxMatrixMarketDataset(mtxfile, featfile, exfile) },
        params
    );

    // Simple checks. Remember, CRISPR is not used in the embeddings by default.
    await overlord(state, false);
    await utils.checkStateResultsBlocked(state);

    // However, CRISPR should still show up in the marker analyses.
    let vres = utils.checkClusterVersusMode(state);
    expect(vres.results.CRISPR instanceof scran.ScoreMarkersResults).toBe(true);

    let custom = utils.launchCustomSelections(state);
    expect(custom.first.CRISPR instanceof scran.ScoreMarkersResults).toBe(true);
    expect(custom.last.CRISPR instanceof scran.ScoreMarkersResults).toBe(true);
    expect(custom.versus.results.CRISPR instanceof scran.ScoreMarkersResults).toBe(true);

    // Freeing everyone.
    await bakana.freeAnalysis(state);
})

test("CRISPR-only runAnalysis works correctly", async () => {
    let state = await bakana.createAnalysis();
    let params = utils.baseParams();
    params.combine_embeddings.crispr_weight = 1;
    let files = { default: new bakana.TenxHdf5Dataset(h5file) };
    files.default.setOptions({ featureTypeRnaName: null });
    await bakana.runAnalysis(state, files, params);

    // Check the results.
    await utils.checkStateResultsMinimal(state);
    await utils.checkStateResultsUnblocked(state);
    await utils.checkStateResultsCrispr(state, { exclusive: true });

    // Freeing.
    await bakana.freeAnalysis(state);
})

