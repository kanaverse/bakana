import * as bakana from "../src/index.js";
import * as scran from "scran.js";
import * as utils from "./utils.js"

beforeAll(utils.initializeAll);
afterAll(async () => await bakana.terminate());

test("multi-matrix analyses work correctly", async () => {
    let fpath4k = "files/datasets/pbmc4k-tenx.h5";
    let fpath3k_mat = "files/datasets/pbmc3k-matrix.mtx.gz";
    let fpath3k_feat = "files/datasets/pbmc3k-features.tsv.gz";
    let files = { 
        "4K": new bakana.TenxHdf5Dataset(fpath4k),
        "3K": new bakana.TenxMatrixMarketDataset(fpath3k_mat, fpath3k_feat, "files/datasets/pbmc3k-barcodes.tsv.gz")
    };

    let paramcopy = utils.baseParams();
    let state = await bakana.createAnalysis();
    let res = await bakana.runAnalysis(state, files, paramcopy);

    {
        // Checking that the matrix was correctly loaded.
        let loaded3k = scran.initializeSparseMatrixFromMatrixMarket(fpath3k_mat, { layered: false });
        let parsed3k = bakana.readTable((new bakana.SimpleFile(fpath3k_feat)).buffer(), { compression: "gz" });
        let names3k = {};
        parsed3k.forEach((x, i) => { names3k[x[0]] = i; })

        let loaded4k = scran.initializeSparseMatrixFromHDF5(fpath4k, "matrix", { layered: false });
        let info4k = (new scran.H5File(fpath4k)).open("matrix").open("features").open("id", { load: true }).values;
        let names4k = {};
        info4k.forEach((x, i) => { names4k[x] = i; }); 

        // Randomly picking every 100th gene and checking we get the same results.
        let combined_names = state.inputs.fetchFeatureAnnotations()["RNA"].column("id");
        let commat = state.inputs.fetchCountMatrix().get("RNA");
        for (var i = 0; i < combined_names.length; i+=100) {
            let i3 = names3k[combined_names[i]];
            let x3 = loaded3k.matrix.row(i3);
            let i4 = names4k[combined_names[i]];
            let x4 = loaded4k.matrix.row(i4);

            let expected = new Float64Array(x3.length + x4.length);
            expected.set(x3, 0);
            expected.set(x4, x3.length);
            expect(commat.row(i)).toEqual(expected);
        }

        // Checking that the intersetions is as expected.
        let common = [];
        for (const x of parsed3k) {
            if (x[0] in names4k) {
                common.push(x[0]);
            }
        }
        expect(common.sort()).toEqual(combined_names.slice().sort());
    }

    await utils.overlordCheckBlocked(state);
    expect(state.inputs.fetchBlockLevels()).toEqual(["3K", "4K"]); 

    let vres = utils.checkClusterVersusMode(state);
    expect(vres.results.RNA.numberOfBlocks()).toEqual(2);

    let custom = utils.launchCustomSelections(state);
    expect(custom.first.RNA.numberOfBlocks()).toEqual(2);
    expect(custom.last.RNA.numberOfBlocks()).toEqual(2);
    expect(custom.versus.results.RNA.numberOfBlocks()).toEqual(2);

    // Check saving of results.
    await bakana.saveSingleCellExperiment(state, "combined", { directory: "miscellaneous/from-tests" });
    await bakana.saveGenewiseResults(state, "combined_genes", { directory: "miscellaneous/from-tests" });

    // Check reloading of the parameters/datasets.
    {
        let saved = [];
        let saver = (n, k, f) => {
            saved.push(f.content());
            return String(saved.length);
        };

        let serialized = await bakana.serializeConfiguration(state, saver);
        let reloaded = bakana.unserializeDatasets(serialized.datasets, x => saved[Number(x) - 1]); 
        expect(reloaded["4K"] instanceof bakana.TenxHdf5Dataset);
        expect(reloaded["3K"] instanceof bakana.TenxMatrixMarketDataset);
        expect(serialized.parameters).toEqual(bakana.retrieveParameters(state));
    }

    // Freeing.
    await bakana.freeAnalysis(state);
})

test("single-matrix multi-sample analyses work correctly", async () => {
    let paramcopy = utils.baseParams();
    paramcopy.inputs.block_factor = "3k";
    let state = await bakana.createAnalysis();
    let res = await bakana.runAnalysis(
        state,
        { 
            "combined": new bakana.TenxMatrixMarketDataset(
                    "files/datasets/pbmc-combined-matrix.mtx.gz", 
                    "files/datasets/pbmc-combined-features.tsv.gz", 
                    "files/datasets/pbmc-combined-barcodes.tsv.gz"
                )
        },
        paramcopy
    );

    await utils.overlordCheckBlocked(state);

    // Freeing.
    await bakana.freeAnalysis(state);
})
