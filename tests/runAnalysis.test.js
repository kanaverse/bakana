import * as bakana from "../src/index.js";
import * as scran from "scran.js";
import * as path from "path";
import * as fs from "fs";

beforeAll(async () => { await scran.initialize({ localFile: true }) });
afterAll(async () => { await scran.terminate() });

let params = {
    inputs: {
        sample_factor: null
    },
    quality_control: {
        use_mito_default: true,
        mito_prefix: "mt-",
        nmads: 3
    },
    feature_selection: {
        span: 0.3
    },
    pca: {
        num_hvgs: 2000,
        num_pcs: 10,
        block_method: "none"
    },
    neighbor_index: {
        approximate: true
    },
    kmeans_cluster: {
        k: 10
    },
    snn_graph_cluster: {
        k: 10,
        scheme: 0,
        resolution: 1
    },
    choose_clustering: {
        method: "snn_graph"
    },
    cell_labelling: {
        mouse_references: [],
        human_references: [ "BlueprintEncode" ],
    }
};

let ref_download = (url) => {
    let fpath = path.basename(decodeURIComponent(url));
    let obj = fs.readFileSync("files/references/" + fpath);
    return obj.buffer.slice(obj.byteOffset, obj.byteOffset + obj.byteLength);
};

test("runAnalysis works correctly (MatrixMarket)", async () => {
    let contents = {};
    let finished = (res, step, msg) => {
        contents[step] = res;
    };

    let res = await bakana.runAnalysis(
        [
            {
                format: "mtx",
                mtx: [ "files/datasets/pbmc3k-matrix.mtx.gz" ],
                gene: [ "files/datasets/pbmc3k-features.tsv.gz" ],
                annotation: [ "files/datasets/pbmc3k-barcodes.tsv.gz" ]
            }
        ],
        params,
        {
            finished: finished,
            download: ref_download
        }
    );

    expect(contents.quality_control instanceof Object).toBe(true);
    expect(contents.pca instanceof Object).toBe(true);
    expect(contents.feature_selection instanceof Object).toBe(true);
    expect(contents.cell_labelling instanceof Object).toBe(true);
    expect(contents.marker_detection instanceof Object).toBe(true);
})

test("runAnalysis works correctly (10X)", async () => {
    let contents = {};
    let finished = (res, step, msg) => {
        contents[step] = res;
    };

    let res = await bakana.runAnalysis(
        [
            {
                format: "tenx",
                file: [ "files/datasets/pbmc4k-tenx.h5" ]
            }
        ],
        params,
        {
            finished: finished,
            download: ref_download
        }
    );

    expect(contents.quality_control instanceof Object).toBe(true);
    expect(contents.pca instanceof Object).toBe(true);
    expect(contents.feature_selection instanceof Object).toBe(true);
    expect(contents.cell_labelling instanceof Object).toBe(true);
    expect(contents.marker_detection instanceof Object).toBe(true);
})

test("runAnalysis works correctly (combined)", async () => {
    let contents = {};
    let finished = (res, step, msg) => {
        contents[step] = res;
    };

    let res = await bakana.runAnalysis(
        [
            {
                format: "tenx",
                file: [ "files/datasets/pbmc4k-tenx.h5" ]
            },
            {
                format: "mtx",
                mtx: [ "files/datasets/pbmc3k-matrix.mtx.gz" ],
                gene: [ "files/datasets/pbmc3k-features.tsv.gz" ],
                annotation: [ "files/datasets/pbmc3k-barcodes.tsv.gz" ]
            }
        ],
        params,
        {
            finished: finished,
            download: ref_download
        }
    );

    expect(contents.quality_control instanceof Object).toBe(true);
    expect(contents.pca instanceof Object).toBe(true);
    expect(contents.feature_selection instanceof Object).toBe(true);
    expect(contents.cell_labelling instanceof Object).toBe(true);
    expect(contents.marker_detection instanceof Object).toBe(true);
})
