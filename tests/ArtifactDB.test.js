import * as bakana from "../src/index.js";
import * as scran from "scran.js";
import * as utils from "./utils.js";
import * as bioc from "bioconductor";
import * as fs from "fs";

beforeAll(utils.initializeAll);
afterAll(async () => await bakana.terminate());

const base = "miscellaneous/from-tests";

async function obtain_test_path() {
    let h5path = base + "/H5AD";
    if (fs.existsSync(h5path)) {
        return "H5AD";
    }

    let fpath = "files/datasets/zeisel-brain.h5ad";
    let files = { 
        default: new bakana.H5adDataset(fpath)
    };

    let altpath = "ArtifactDB-testing";
    let fullalt = base + "/" + altpath;
    if (fs.existsSync(fullalt)) {
        return altpath;
    }

    let state = await bakana.createAnalysis();
    try {
        let params = utils.baseParams();
        await bakana.runAnalysis(state, files, params);
        await bakana.saveSingleCellExperiment(state, altpath, { directory: base });
    } finally {
        await bakana.freeAnalysis(state);
    }

    return altpath;
}

class ProjectDirectoryNavigator {
    #directory;

    constructor(d) {
        this.#directory = d;
    }

    serialize() {
        let enc = new TextEncoder;
        let contents = enc.encode(this.#directory);
        let f = new bakana.SimpleFile(contents, { name: this.#directory });
        return [{ type: "directory", file: f }];
    }

    metadata(p) {
        let fullpath = this.#directory + "/" + p;

        while (1) {
            if (!fullpath.endsWith(".json")) { 
                fullpath += ".json";
            }
            let contents = fs.readFileSync(fullpath);

            let dec = new TextDecoder;
            let json = dec.decode(contents);
            let values = JSON.parse(json);

            if (values["$schema"].startsWith("redirection/")){
                fullpath = this.#directory + "/" + values.redirection.targets[0].location;
            } else {
                return values;
            }
        }
    }

    file(p) {
        return this.#directory + "/" + p;
    }
}

class ArtifactDbProjectDataset extends bakana.ArtifactDbSummarizedExperimentDatasetBase {
    constructor(path, dir, options={}) {
        super(path, new ProjectDirectoryNavigator(dir), options);
    }

    static format() {
        return "ArtifactDB-directory";
    }

    static unserialize(files, options) {
        let dec = new TextDecoder;
        let path = dec.decode(files[0].file.buffer());
        let dir = dec.decode(files[1].file.buffer());
        return new ArtifactDbProjectDataset(path, dir, options);
    }
}

test("ArtifactDB summary works correctly", async () => {
    let target = await obtain_test_path();
    let files = { default: new ArtifactDbProjectDataset(target, base) };
    let summ = await files.default.summary({ cache: true });

    expect(summ.modality_features[""] instanceof bioc.DataFrame).toBe(true);
    expect(summ.modality_features[""].numberOfRows()).toBeGreaterThan(0);
    expect(summ.cells instanceof bioc.DataFrame).toBe(true);

    // Clear the cache.
    files.default.clear();
})

test("runAnalysis works correctly (ArtifactDB)", async () => {
    let target = await obtain_test_path();
    let files = { default: new ArtifactDbProjectDataset(target, base) };
    let state = await bakana.createAnalysis();
    let params = utils.baseParams();
    await bakana.runAnalysis(state, files, params);

    // Input reorganization is done correctly.
    {
        let loaded = state.inputs.fetchCountMatrix().get("RNA");
        let loaded_ids = state.inputs.fetchRowIds()["RNA"];
        let loaded_names = state.inputs.fetchFeatureAnnotations()["RNA"].rowNames();

        let fpath = base + "/" + target + "/assay-counts/matrix.h5";
        let simple = scran.initializeSparseMatrixFromHDF5(fpath, "matrix", { layered: false });
        let rpath = base + "/" + target + "/rowdata/simple.h5";
        let simple_names = (new scran.H5File(rpath)).open("data").open("row_names", { load: true }).values;

        utils.checkReorganization(simple.matrix, simple.row_ids, simple_names, loaded, loaded_ids, loaded_names);
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
        bakana.availableReaders["ArtifactDB-directory"] = ArtifactDbProjectDataset;
        let reloaded = bakana.unserializeDatasets(serialized.datasets, x => saved[Number(x) - 1]); 
        expect(reloaded.default instanceof ArtifactDbProjectDataset);
        expect(serialized.parameters).toEqual(bakana.retrieveParameters(state));
    }

    await bakana.freeAnalysis(state);
})

