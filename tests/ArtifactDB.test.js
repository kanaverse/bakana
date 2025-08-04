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

    constructor(path, dir) {
        super(path, new nav.LocalProjectDirectoryNavigator(dir));
        this.#path = path;
        this.#dir = dir;
    }

    static format() {
        return "ArtifactDB-local-directory";
    }

    static unserialize(files, options) {
        let dec = new TextDecoder;
        let payload = files[0].file.buffer();
        let info = JSON.parse(dec.decode(payload));
        let output = new LocalArtifactdbDataset(info.path, info.directory);
        output.setOptions(options);
        return output;
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

test("ArtifactDB readers work correctly", async () => {
    let ds = new LocalArtifactdbDataset("experiment.json", nav.baseDirectory + "/zeisel-brain-stripped"); 
    await utils.checkDatasetGeneral(ds);

    let summ = await utils.checkDatasetSummary(ds);
    expect(Object.keys(summ.modality_features)).toEqual([""]);

    let loaded = await utils.checkDatasetLoad(ds);
    expect(loaded.matrix.available()).toEqual(["RNA"]);

    let copy = await utils.checkDatasetSerialize(ds);
    utils.sameDatasetSummary(summ, await copy.summary());
    utils.sameDatasetLoad(loaded, await copy.load());

    // Input reorganization is done correctly.
    {
        let fpath = nav.baseDirectory + "/zeisel-brain-stripped/assay-1/matrix.h5";
        let simple = scran.initializeSparseMatrixFromHdf5(fpath, "sparse", { layered: false });
        let rpath = nav.baseDirectory + "/zeisel-brain-stripped/rowdata/simple.h5";
        let simple_names = (new scran.H5DataSet(rpath, "contents/row_names", { load: true })).values;
        utils.checkMatrixContents(simple, simple_names, loaded.matrix.get("RNA"), loaded.features["RNA"].rowNames());
        simple.free();
    }

    ds.clear();
})

test("ArtifactDB readers work with multiple modalities", async () => {
    let ds = new LocalArtifactdbDataset("experiment.json", nav.baseDirectory + "/zeisel-brain-sparse");
    await utils.checkDatasetGeneral(ds);

    let summ = await utils.checkDatasetSummary(ds);
    expect(Object.keys(summ.modality_features)).toEqual(["", "repeat", "ERCC"]);
    expect(summ.modality_assay_names).toEqual({ "": ["counts"], "ERCC": ["counts"], "repeat": ["counts"] });

    // Trying with a name for the experiment.
    ds.setOptions({ adtExperiment: "ERCC" });
    {
        let everything = await utils.checkDatasetLoad(ds);
        expect(Object.keys(everything.features)).toEqual(["RNA", "ADT"]);
        let nRna = everything.features["RNA"].numberOfRows();
        let nAdt = everything.features["ADT"].numberOfRows();
        expect(nAdt).toBeLessThan(nRna);
    }

    // Trying with a numeric index for the ADT experiment.
    ds.setOptions({ adtExperiment: 0 });
    {
        let everything = await utils.checkDatasetLoad(ds);
        expect(Object.keys(everything.features)).toEqual(["RNA", "ADT"]);
        let nRna = everything.features["RNA"].numberOfRows();
        let nAdt = everything.features["ADT"].numberOfRows();
        expect(nAdt).toBeLessThan(nRna);
    }
})

test("Zipped ArtifactDB readers work correctly", async () => {
    let zipped = await nav.zipDirectory(nav.baseDirectory + "/zeisel-brain-stripped", [ ".", "experiment.json" ]);
    let zipfile = new bakana.SimpleFile(zipped, { name: "bundle.zip" });
    let ds = new bakana.ZippedArtifactdbDataset("experiment.json", zipfile);

    let summ = await utils.checkDatasetSummary(ds);
    expect(Object.keys(summ.modality_features)).toEqual([""]);
    expect(summ.modality_assay_names).toEqual({ "": ["counts"] });

    let everything = await utils.checkDatasetLoad(ds);
    expect(Object.keys(everything.features)).toEqual(["RNA"]);

    let copy = await utils.checkDatasetSerialize(ds);
    utils.sameDatasetSummary(summ, await copy.summary());
    utils.sameDatasetLoad(everything, await copy.load());

    // Checking other input modes.
    let handle = await jszip.loadAsync(await zipfile.buffer());
    {
        let zipped = new bakana.ZippedArtifactdbDataset("experiment.json", zipfile, { existingHandle: handle });
        let summ2 = await zipped.summary();
        expect(summ2.modality_features[""].rowNames()).toEqual(summ2.modality_features[""].rowNames());
        zipped.clear();
    }

    // Inspection works as expected.
    {
        let results = await bakana.searchZippedArtifactdb(handle);
        expect(results.get("experiment.json").length).toEqual(2);
    }
})

/***********************************************/

class LocalArtifactdbResult extends bakana.AbstractArtifactdbResult {
    constructor(path, dir) {
        super(path, new nav.LocalProjectDirectoryNavigator(dir));
    }
}

test("local ArtifactDB result readers work correctly with log-count loading", async () => {
    let res = new LocalArtifactdbResult("experiment.json", nav.baseDirectory + "/zeisel-brain-sparse-results");

    let details = await utils.checkResultSummary(res);
    expect(Object.keys(details.modality_features)).toEqual([""]);
    expect(details.modality_assay_names[""]).toEqual(["filtered", "normalized"]);
    expect(details.reduced_dimension_names).toEqual(["pca", "tsne", "umap"]);

    // Actually loading the log-counts.
    res.setOptions({
        primaryAssay: "normalized",
        isPrimaryNormalized: true,
        reducedDimensionNames: ["tsne", "umap"]
    });
    {
        let payload = await utils.checkResultLoad(res);
        expect(Object.keys(payload.reduced_dimensions)).toEqual(["tsne", "umap"]);
        expect(payload.reduced_dimensions["tsne"].length).toEqual(2);
        expect(payload.reduced_dimensions["umap"].length).toEqual(2);

        let col0 = payload.matrix.get("").column(0);
        expect(utils.hasNonInteger(col0)).toBe(true);
    }

    // Works correctly with manual normalization.
    res.setOptions({
        primaryAssay: "filtered",
        isPrimaryNormalized: false
    });
    {
        let payload = await utils.checkResultLoad(res);
        let col0 = payload.matrix.get("").column(0);
        expect(utils.hasNonInteger(col0)).toBe(true);
    }

    res.clear();
})

test("local ArtifactDB result readers work correctly with multiple modalities", async () => {
    let res = new LocalArtifactdbResult("experiment.json", nav.baseDirectory + "/zeisel-brain-dense-multimodal-results");

    let details = await utils.checkResultSummary(res);
    expect(Object.keys(details.modality_features)).toEqual(["", "adt"]);
    expect(details.modality_features[""].numberOfRows()).toBeGreaterThan(0);
    expect(details.modality_features["adt"].numberOfRows()).toBeGreaterThan(0);

    let payload = await utils.checkResultLoad(res);
    let nRna = payload.features[""].numberOfRows();
    let nAdt = payload.features["adt"].numberOfRows();
    expect(nAdt).toBeLessThan(nRna);

    res.clear();
})

/***********************************************/

test("Zipped ArtifactDB result summary and loading works correctly", async () => {
    // First, zipping the contents of the target directory.
    let zipped = await nav.zipDirectory(nav.baseDirectory + "/zeisel-brain-sparse-results", [ ".", "experiment.json" ]);
    scran.writeFile("miscellaneous/bundle.zip", zipped); // for manual inspection.
    let zipfile = new bakana.SimpleFile(zipped, { name: "bundle.zip" });
    let res = new bakana.ZippedArtifactdbResult("experiment.json", zipfile);

    let details = await utils.checkResultSummary(res);
    expect(Object.keys(details.modality_features)).toEqual([""]);
    expect(details.modality_assay_names[""]).toEqual(["filtered", "normalized"]);
    expect(details.reduced_dimension_names).toEqual(["pca", "tsne", "umap"]);

    // Actually loading the log-counts.
    res.setOptions({
        primaryAssay: "normalized",
        isPrimaryNormalized: true,
        reducedDimensionNames: ["tsne", "umap"]
    });

    let payload = await utils.checkResultLoad(res);
    expect(Object.keys(payload.reduced_dimensions)).toEqual(["tsne", "umap"]);
    expect(payload.reduced_dimensions["tsne"].length).toEqual(2);
    expect(payload.reduced_dimensions["umap"].length).toEqual(2);

    let col0 = payload.matrix.get("").column(0);
    expect(utils.hasNonInteger(col0)).toBe(true);

    res.clear();
})
