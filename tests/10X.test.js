import * as bakana from "../src/index.js";
import * as scran from "scran.js";
import * as utils from "./utils.js";
import * as bioc from "bioconductor";

beforeAll(utils.initializeAll);
afterAll(async () => await bakana.terminate());

let h5path = "files/datasets/pbmc4k-tenx.h5";
let files = {
    default: new bakana.TenxHdf5Dataset(h5path)
};

test("10X summary works correctly", async () => {
    let summ = files.default.summary();
    expect(summ.modality_features[""] instanceof bioc.DataFrame).toBe(true);
    expect(summ.modality_features[""].numberOfRows()).toBeGreaterThan(0);
    expect(summ.cells instanceof bioc.DataFrame).toBe(true);

    let preview = files.default.previewPrimaryIds();
    expect("RNA" in preview).toBe(true);
    expect(preview.RNA.length).toBeGreaterThan(0);

    // Trying what happens if we force it to be another modality.
    files.default.setOptions({ featureTypeAdtName: "", featureTypeCrisprName: null, featureTypeRnaName: null });
    let preview2 = files.default.previewPrimaryIds();
    expect("ADT" in preview2).toBe(true);
    expect("RNA" in preview2).toBe(false);
    expect(preview2.ADT.length).toBeGreaterThan(0);
})

test("runAnalysis works correctly (10X)", async () => {
    let state = await bakana.createAnalysis();
    let params = utils.baseParams();
    let res = await bakana.runAnalysis(state, files, params);

    // Input reorganization is done correctly.
    {
        let loaded = state.inputs.fetchCountMatrix().get("RNA");
        let loaded_names = state.inputs.fetchFeatureAnnotations()["RNA"].column("id");

        let simple = scran.initializeSparseMatrixFromHdf5(h5path, "matrix", { layered: false });
        let simple_names = (new scran.H5File(h5path)).open("matrix").open("features").open("id", { load: true }).values;

        utils.checkMatrixContents(simple, simple_names, loaded, loaded_names);
        simple.free();
    }

    // Basic consistency checks.
    await utils.overlordCheckStandard(state);
    utils.checkClusterVersusMode(state);
    await utils.triggerAnimation(state);

    // Check saving of results.
    utils.purgeDirectory("miscellaneous/from-tests/10X");
    await bakana.saveSingleCellExperiment(state, "10X", { directory: "miscellaneous/from-tests" });
    utils.purgeDirectory("miscellaneous/from-tests/10X_genes");
    await bakana.saveGenewiseResults(state, "10X_genes", { directory: "miscellaneous/from-tests" });

    // Check reloading of the parameters/datasets.
    {
        let saved = [];
        let saver = (n, k, f) => {
            saved.push(f.content());
            return String(saved.length);
        };

        let serialized = await bakana.serializeConfiguration(state, saver);
        let reloaded = bakana.unserializeDatasets(serialized.datasets, x => saved[Number(x) - 1]); 
        expect(reloaded.default instanceof bakana.TenxHdf5Dataset);
        expect(serialized.parameters).toEqual(bakana.retrieveParameters(state));
    }

    // Freeing.
    await bakana.freeAnalysis(state);
})
