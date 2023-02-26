import * as bakana from "../src/index.js";
import * as scran from "scran.js";
import * as utils from "./utils.js";
import * as bioc from "bioconductor";
import * as fs from "fs";
import * as nav from "./navigator.js";

beforeAll(utils.initializeAll);
afterAll(async () => await bakana.terminate());

class ArtifactDbLocalDirectoryDataset extends bakana.ArtifactDbSummarizedExperimentDatasetBase {
    #path;
    #dir;

    constructor(path, dir, options={}) {
        super(path, new nav.LocalProjectDirectoryNavigator(dir), options);
        this.#path = path;
        this.#dir = dir;
    }

    static format() {
        return "ArtifactDB-local-directory";
    }

    static unserialize(files, options) {
        let dec = new TextDecoder;
        let payload = files[0].file.buffer();
        let info = JSON.stringify(dec.decode(payload));
        return new ArtifactDbLocalDirectoryDataset(info.path, info.directory, options);
    }

    abbreviate() {
        return {
            path: this.#path,
            directory: this.#dir,
            options: this.options()
        };
    }

    serialize() {
        let enc = new TextEncoder;
        let config = { path: this.#path, directory: this.#dir };
        let buffer = enc.encode(JSON.stringify(config));
        return {
            files: [ { file: new bakana.SimpleFile(buffer, { name: "config.json" }), type: "config" } ],
            options: this.options()
        }
    }
}

let target_simple = nav.pathExists("H5AD");
let action_simple = (target_simple == null ? test.skip : test);

action_simple("ArtifactDB summary works correctly", async () => {
    let files = { default: new ArtifactDbLocalDirectoryDataset(target_simple, nav.baseDirectory) };
    let summ = await files.default.summary({ cache: true });

    expect(summ.modality_features[""] instanceof bioc.DataFrame).toBe(true);
    expect(summ.modality_features[""].numberOfRows()).toBeGreaterThan(0);
    expect(summ.cells instanceof bioc.DataFrame).toBe(true);

    // Clear the cache.
    files.default.clear();
})

action_simple("runAnalysis works correctly (ArtifactDB)", async () => {
    let files = { default: new ArtifactDbLocalDirectoryDataset(target_simple, nav.baseDirectory) };
    let state = await bakana.createAnalysis();
    let params = utils.baseParams();
    await bakana.runAnalysis(state, files, params);

    // Input reorganization is done correctly.
    {
        let loaded = state.inputs.fetchCountMatrix().get("RNA");
        let loaded_ids = state.inputs.fetchRowIds()["RNA"];
        let loaded_names = state.inputs.fetchFeatureAnnotations()["RNA"].rowNames();

        let fpath = nav.baseDirectory + "/" + target_simple + "/assay-counts/matrix.h5";
        let simple = scran.initializeSparseMatrixFromHDF5(fpath, "matrix", { layered: false });
        let rpath = nav.baseDirectory + "/" + target_simple + "/rowdata/simple.h5";
        let simple_names = (new scran.H5DataSet(rpath, "data/row_names", { load: true })).values;

        utils.checkReorganization(simple.matrix, simple.row_ids, simple_names, loaded, loaded_ids, loaded_names, { mustDiffer: false });
        simple.matrix.free();
    }

    // Basic checks.
    await utils.checkStateResultsMinimal(state, { mustFilter: false });
    await utils.checkStateResultsRna(state, { exclusive: true, mustFilter: false });
    await utils.checkStateResultsUnblocked(state);

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

let target_adt = nav.pathExists("adt");
let action_adt = (target_adt == null ? test.skip : test);

action_adt("ArtifactDB summary and loading works with multiple modalities", async () => {
    let files = { super: new ArtifactDbLocalDirectoryDataset(target_adt, nav.baseDirectory) };

    let summ = await files.super.summary();
    expect(Object.keys(summ.modality_features).sort()).toEqual(["", "ADT"]);
    expect(summ.modality_assay_names).toEqual({ "": ["counts", "logcounts"], "ADT": [ "counts", "logcounts" ] });

    // Trying with a name for the experiment.
    files.super.setAdtExperiment("ADT");
    {
        let everything = await files.super.load();

        let nRna = everything.features["RNA"].numberOfRows();
        expect(nRna).toBeGreaterThan(0);
        expect(everything.matrix.get("RNA").numberOfRows()).toEqual(nRna);

        let nAdt = everything.features["ADT"].numberOfRows();
        expect(nAdt).toBeGreaterThan(0);
        expect(nAdt).toBeLessThan(nRna);
        expect(everything.matrix.get("ADT").numberOfRows()).toEqual(nAdt);

        expect(everything.matrix.numberOfColumns()).toEqual(everything.cells.numberOfRows());
        everything.matrix.free();
    }

    // Trying with a numeric index for the ADT experiment.
    files.super.setAdtExperiment(0);
    {
        let everything = await files.super.load();

        let nRna = everything.features["RNA"].numberOfRows();
        expect(nRna).toBeGreaterThan(0);
        expect(everything.matrix.get("RNA").numberOfRows()).toEqual(nRna);

        let nAdt = everything.features["ADT"].numberOfRows();
        expect(nAdt).toBeGreaterThan(0);
        expect(nAdt).toBeLessThan(nRna);
        expect(everything.matrix.get("ADT").numberOfRows()).toEqual(nAdt);

        expect(everything.matrix.numberOfColumns()).toEqual(everything.cells.numberOfRows());
        everything.matrix.free();
    }
})

