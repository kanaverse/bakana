import * as bakana from "../src/index.js";
import * as scran from "scran.js";
import * as utils from "./utils.js";
import * as bioc from "bioconductor";
import * as fs from "fs";
import * as nav from "./navigator.js";
import jszip from "jszip";

beforeAll(utils.initializeAll);
afterAll(async () => await bakana.terminate());

class LocalArtifactdbDataset extends bakana.AbstractArtifactdbDataset {
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
        return new LocalArtifactdbDataset(info.path, info.directory, options);
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

/***********************************************/

let target_simple = nav.pathExists("H5AD");
let action_simple = (target_simple == null ? test.skip : test);

test("ArtifactDB summary works correctly", async () => {
    let files = { default: new LocalArtifactdbDataset("experiment.json", nav.baseDirectory + "/zeisel-brain-stripped") };
    let summ = await files.default.summary({ cache: true });

    expect(summ.modality_features[""] instanceof bioc.DataFrame).toBe(true);
    expect(summ.modality_features[""].numberOfRows()).toBeGreaterThan(0);
    expect(summ.cells instanceof bioc.DataFrame).toBe(true);

    let preview = await files.default.previewPrimaryIds({ cache: true });
    expect("RNA" in preview).toBe(true);
    expect(preview.RNA.length).toBeGreaterThan(0);

    // Clear the cache.
    files.default.clear();
})

test("runAnalysis works correctly (ArtifactDB)", async () => {
    let files = { default: new LocalArtifactdbDataset("experiment.json", nav.baseDirectory + "/zeisel-brain-stripped") };
    let state = await bakana.createAnalysis();
    let params = utils.baseParams();
    await bakana.runAnalysis(state, files, params);

    // Input reorganization is done correctly.
    {
        let loaded = state.inputs.fetchCountMatrix().get("RNA");
        let loaded_names = state.inputs.fetchFeatureAnnotations()["RNA"].rowNames();

        let fpath = nav.baseDirectory + "/" + target_simple + "/assay-counts/matrix.h5";
        let simple = scran.initializeSparseMatrixFromHdf5(fpath, "matrix", { layered: false });
        let rpath = nav.baseDirectory + "/" + target_simple + "/rowdata/simple.h5";
        let simple_names = (new scran.H5DataSet(rpath, "data/row_names", { load: true })).values;

        utils.checkMatrixContents(simple, simple_names, loaded, loaded_names);
        simple.free();
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
        bakana.availableReaders["ArtifactDB-local-directory"] = LocalArtifactdbDataset;
        let reloaded = bakana.unserializeDatasets(serialized.datasets, x => saved[Number(x) - 1]); 
        expect(reloaded.default instanceof LocalArtifactdbDataset);
        expect(serialized.parameters).toEqual(bakana.retrieveParameters(state));
    }

    await bakana.freeAnalysis(state);
})

/***********************************************/

test("ArtifactDB summary and loading works with multiple modalities", async () => {
    let files = { super: new LocalArtifactdbDataset("experiment.json", nav.baseDirectory + "/zeisel-brain-sparse") };

    let summ = await files.super.summary();
    expect(Object.keys(summ.modality_features).sort()).toEqual(["", "ERCC","repeat"]);
    expect(summ.modality_assay_names).toEqual({ "": ["counts"], "ERCC": ["counts"], "repeat": ["counts"] });

    // Trying with a name for the experiment.
    files.super.setOptions({ adtExperiment: "ERCC" });
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
    files.super.setOptions({ adtExperiment: 0 });
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

/***********************************************/

test("Zipped ArtifactDB dataset summary and loading works correctly", async () => {
    // First, zipping the contents of the target directory.
    let zipped = await nav.zipDirectory(nav.baseDirectory + "/zeisel-brain-stripped", [ "experiment.json" ]);
    let zipfile = new bakana.SimpleFile(zipped, { name: "bundle.zip" });

    let files = { zipped: new bakana.ZippedArtifactdbDataset(target_simple, zipfile) };
    expect(files.zipped.constructor.format()).toEqual("ArtifactDB-zipped");
    let abbr = files.zipped.abbreviate();
    expect(abbr.files[0].type).toEqual("zip");
    expect(abbr.options.datasetName).toEqual(target_simple);

    let summ = await files.zipped.summary();
    expect(summ.modality_features[""].numberOfRows()).toBeGreaterThan(0);
    expect(summ.modality_features[""].rowNames()).not.toBeNull();
    expect(summ.cells.numberOfColumns()).toBeGreaterThan(0);
    expect(summ.modality_assay_names[""]).toEqual(["counts", "logcounts"]);

    // Loading everything.
    {
        let loaded = await files.zipped.load();
        expect(loaded.matrix.get("RNA").numberOfRows()).toEqual(loaded.features["RNA"].numberOfRows());
        expect(loaded.matrix.numberOfColumns()).toEqual(loaded.cells.numberOfRows());
        loaded.matrix.free();
    }

    // Running through a serialization cycle.
    {
        let dump = await files.zipped.serialize();
        expect(dump.options.datasetName).toEqual(target_simple);
        expect(dump.files.length).toEqual(1);

        let reloaded = files.zipped.constructor.unserialize(dump.files, dump.options);
        expect(reloaded instanceof bakana.ZippedArtifactdbDataset).toBe(true);
        expect(reloaded.abbreviate()).toEqual(abbr);
    }

    // Checking other input modes.
    let handle = await jszip.loadAsync(await zipfile.buffer());
    {
        let files2 = { zipped: new bakana.ZippedArtifactdbDataset(target_simple, zipfile, { existingHandle: handle }) };
        let summ2 = await files2.zipped.summary();
        expect(summ2.modality_features[""].rowNames()).toEqual(summ2.modality_features[""].rowNames());
        files2.zipped.clear();
    }

    // Inspection works as expected.
    {
        let results = await bakana.searchZippedArtifactdb(handle);
        expect(results.get(target_simple).length).toEqual(2);
    }
})
