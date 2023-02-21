import * as bakana from "../src/index.js";
import * as scran from "scran.js";
import * as utils from "./utils.js";
import * as bioc from "bioconductor";
import * as fs from "fs";
import * as nav from "./navigator.js";

beforeAll(utils.initializeAll);
afterAll(async () => await bakana.terminate());

class ArtifactDbLocalDirectoryDataset extends bakana.ArtifactDbSummarizedExperimentDatasetBase {
    constructor(path, dir, options={}) {
        super(path, new nav.LocalProjectDirectoryNavigator(dir), options);
    }

    static format() {
        return "ArtifactDB-local-directory";
    }

    static unserialize(files, options) {
        let dec = new TextDecoder;
        let path = dec.decode(files[0].file.buffer());
        let dir = dec.decode(files[1].file.buffer());
        return new ArtifactDbLocalDirectoryDataset(path, dir, options);
    }
}

test("ArtifactDB summary works correctly", async () => {
    let target = await nav.obtainLocalTestPath();
    let files = { default: new ArtifactDbLocalDirectoryDataset(target, nav.baseDirectory) };
    let summ = await files.default.summary({ cache: true });

    expect(summ.modality_features[""] instanceof bioc.DataFrame).toBe(true);
    expect(summ.modality_features[""].numberOfRows()).toBeGreaterThan(0);
    expect(summ.cells instanceof bioc.DataFrame).toBe(true);

    // Clear the cache.
    files.default.clear();
})

test("runAnalysis works correctly (ArtifactDB)", async () => {
    let target = await nav.obtainLocalTestPath();
    let files = { default: new ArtifactDbLocalDirectoryDataset(target, nav.baseDirectory) };
    let state = await bakana.createAnalysis();
    let params = utils.baseParams();
    await bakana.runAnalysis(state, files, params);

    // Input reorganization is done correctly.
    {
        let loaded = state.inputs.fetchCountMatrix().get("RNA");
        let loaded_ids = state.inputs.fetchRowIds()["RNA"];
        let loaded_names = state.inputs.fetchFeatureAnnotations()["RNA"].rowNames();

        let fpath = nav.baseDirectory + "/" + target + "/assay-counts/matrix.h5";
        let simple = scran.initializeSparseMatrixFromHDF5(fpath, "matrix", { layered: false });
        let rpath = nav.baseDirectory + "/" + target + "/rowdata/simple.h5";
        let simple_names = (new scran.H5DataSet(rpath, "data/row_names", { load: true })).values;

        utils.checkReorganization(simple.matrix, simple.row_ids, simple_names, loaded, loaded_ids, loaded_names, { mustDiffer: false });
        simple.matrix.free();
    }

    // Basic checks.
    await utils.overlordCheckStandard(state);

    // Check reloading of the parameters/datasets.
    {
        let saved = [];
        let saver = (n, k, f) => {
            saved.push(f.content());
            return String(saved.length);
        };

        let serialized = await bakana.serializeConfiguration(state, saver);
        bakana.availableReaders["ArtifactDB-local-directory"] = ArtifactDbLocalDirectoryDataset;
        let reloaded = bakana.unserializeDatasets(serialized.datasets, x => saved[Number(x) - 1]); 
        expect(reloaded.default instanceof ArtifactDbLocalDirectoryDataset);
        expect(serialized.parameters).toEqual(bakana.retrieveParameters(state));
    }

    await bakana.freeAnalysis(state);
})

