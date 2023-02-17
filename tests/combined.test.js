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

        // IDs should be the same as the first matrix.
        let expected_ids = new Int32Array(combined_names.length);
        combined_names.forEach((x, i) => { expected_ids[i] = names3k[x]; });
        expect(state.inputs.fetchRowIds()["RNA"]).toEqual(expected_ids);
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

    // Saving and loading.
    const path = "TEST_state_multi-matrix.h5";
    let collected = await bakana.saveAnalysis(state, path);
    utils.validateState(path);
    expect(collected.collected.length).toBe(4);
    expect(typeof(collected.collected[0])).toBe("string");

    {
        let handle = new scran.H5File(path);
        let ihandle = handle.open("inputs");
        expect(ihandle.open("results").open("num_blocks", { load: true }).values[0]).toEqual(2);

        let phandle = ihandle.open("parameters");
        let dhandle = phandle.open("datasets");
        expect(dhandle.open("0").open("name", { load: true }).values[0]).toBe("3K"); // should be sorted: 3K, then 4K.
        expect(dhandle.open("1").open("name", { load: true }).values[0]).toBe("4K");

        expect(dhandle.open("0").open("format", { load: true }).values[0]).toBe("MatrixMarket"); 
        expect(dhandle.open("1").open("format", { load: true }).values[0]).toBe("10X");

        expect(Object.keys(dhandle.open("0").open("files").children).length).toBe(3);
        expect(Object.keys(dhandle.open("1").open("files").children).length).toBe(1);
    }

    let offsets = utils.mockOffsets(collected.collected);
    let reloaded = await bakana.loadAnalysis(
        path, 
        (offset, size) => offsets[offset]
    );

    let new_params = bakana.retrieveParameters(reloaded);
    expect(new_params).toEqual(paramcopy);

    await utils.compareStates(state, reloaded);

    {
        // Check that the filtered blocks are correctly restored. This checks
        // for bugs when the blocks are requested before the filtered matrix is
        // restored, given that the former depends on the latter.
        let fblocks = reloaded.cell_filtering.fetchFilteredBlock();
        let fmat = reloaded.cell_filtering.fetchFilteredMatrix();
        expect(fblocks.length).toBe(fmat.numberOfColumns());
    }

    // Freeing.
    await bakana.freeAnalysis(state);
    await bakana.freeAnalysis(reloaded);
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

    // Saving and loading.
    const path = "TEST_state_multi-sample.h5";
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
    expect(new_params).toEqual(paramcopy);

    await utils.compareStates(state, reloaded);

    // Freeing.
    await bakana.freeAnalysis(state);
    await bakana.freeAnalysis(reloaded);
})
