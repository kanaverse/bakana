import * as bakana from "../src/index.js";
import * as scran from "scran.js";
import * as utils from "./utils.js"

beforeAll(utils.initializeAll);
afterAll(async () => await bakana.terminate());

test("switching between clustering methods (SNN first)", async () => {
    let files = {
        default: {
            format: "10X",
            h5: "files/datasets/pbmc4k-tenx.h5"
        }
    };

    // First running the analysis with SNN graph,
    // and confirming that only SNN graph results are reported.
    let contents = {};
    let finished = (step, res) => {
        contents[step] = res;
    };

    let state = await bakana.createAnalysis();
    let paramcopy = utils.baseParams();
    await bakana.runAnalysis(state, files, paramcopy, { finishFun: finished });

    expect(contents.pca instanceof Object).toBe(true);
    expect(contents.snn_graph_cluster instanceof Object).toBe(true);
    expect(contents.kmeans_cluster instanceof Object).toBe(true);

    const path = "TEST_state_clusters.h5";
    await bakana.saveAnalysis(state, path);
    utils.validateState(path);
    {
        let handle = new scran.H5File(path);
        let khandle = handle.open("kmeans_cluster");
        let krhandle = khandle.open("results");
        expect("clusters" in krhandle.children).toBe(false);

        let shandle = handle.open("snn_graph_cluster");
        let srhandle = shandle.open("results");
        expect("clusters" in srhandle.children).toBe(true);
    }

    // Now trying with k-means. This should cause both sets of
    // results to be saved, as both clusterings are still valid.
    paramcopy.choose_clustering.method = "kmeans";

    contents = {};
    await bakana.runAnalysis(state, files, paramcopy, { finishFun: finished });
    expect(contents.inputs).toBeUndefined();
    expect(contents.pca).toBeUndefined();
    expect(contents.snn_graph_cluster).toBeUndefined();
    expect(contents.kmeans_cluster instanceof Object).toBe(true);

    await bakana.saveAnalysis(state, path);
    utils.validateState(path);
    {
        let handle = new scran.H5File(path);
        let khandle = handle.open("kmeans_cluster");
        let krhandle = khandle.open("results");
        expect("clusters" in krhandle.children).toBe(true);

        let shandle = handle.open("snn_graph_cluster");
        let srhandle = shandle.open("results");
        expect("clusters" in srhandle.children).toBe(true);
    }

    // Checking that invalidation of the results behaves correctly.
    // If we change the parameters but we're still on the old set of
    // results, the SNN graph results don't get rerun and the results get wiped.
    paramcopy.snn_graph_cluster.resolution = 0.77;

    contents = {};
    await bakana.runAnalysis(state, files, paramcopy, { finishFun: finished });
    expect(contents.snn_graph_cluster instanceof Object).toBe(true);
    expect(contents.kmeans_cluster).toBeUndefined();
    expect(contents.choose_clustering).toBeUndefined();

    await bakana.saveAnalysis(state, path);
    utils.validateState(path);
    {
        let handle = new scran.H5File(path);
        let khandle = handle.open("kmeans_cluster");
        let krhandle = khandle.open("results");
        expect("clusters" in krhandle.children).toBe(true);

        let shandle = handle.open("snn_graph_cluster");
        let sphandle = shandle.open("parameters");
        expect(sphandle.open("resolution", { load: true }).values[0]).toEqual(0.77);
        let srhandle = shandle.open("results");
        expect("clusters" in srhandle.children).toBe(false);
    }

    // Freeing all states.
    await bakana.freeAnalysis(state);
});

test("switching between clustering methods (k-means first)", async () => {
    let files = {
        default: {
            format: "MatrixMarket",
            mtx: "files/datasets/pbmc3k-matrix.mtx.gz",
            genes: "files/datasets/pbmc3k-features.tsv.gz"
        }
    };

    let contents = {};
    let finished = (step, res) => {
        contents[step] = res;
    };

    let state = await bakana.createAnalysis();
    let paramcopy = utils.baseParams();
    paramcopy.choose_clustering = { method: "kmeans" };

    await bakana.runAnalysis(state, files, paramcopy, { finishFun: finished });
    expect(contents.snn_graph_cluster instanceof Object).toBe(true);
    expect(contents.kmeans_cluster instanceof Object).toBe(true);

    const path = "TEST_state_clusters.h5";
    await bakana.saveAnalysis(state, path);
    utils.validateState(path);
    {
        let handle = new scran.H5File(path);
        let khandle = handle.open("kmeans_cluster");
        let krhandle = khandle.open("results");
        expect("clusters" in krhandle.children).toBe(true);

        let shandle = handle.open("snn_graph_cluster");
        let srhandle = shandle.open("results");
        expect("clusters" in srhandle.children).toBe(false);
    }

    // Trying with a different clustering method.
    paramcopy.choose_clustering.method = "snn_graph"; 
    contents = {};
    await bakana.runAnalysis(state, files, paramcopy, { finishFun: finished });

    expect(contents.kmeans_cluster).toBeUndefined();
    expect(contents.snn_graph_cluster instanceof Object).toBe(true);

    await bakana.saveAnalysis(state, path);
    utils.validateState(path);
    {
        let handle = new scran.H5File(path);
        let khandle = handle.open("kmeans_cluster");
        let krhandle = khandle.open("results");
        expect("clusters" in krhandle.children).toBe(true);

        let shandle = handle.open("snn_graph_cluster");
        let srhandle = shandle.open("results");
        expect("clusters" in srhandle.children).toBe(true);
    }

    // Checking that invalidation works.
    paramcopy.kmeans_cluster.k = 7;
    contents = {};
    await bakana.runAnalysis(state, files, paramcopy, { finishFun: finished });

    expect(contents.kmeans_cluster instanceof Object).toBe(true);
    expect(contents.snn_graph_cluster).toBeUndefined();
    expect(contents.choose_clustering).toBeUndefined();

    await bakana.saveAnalysis(state, path);
    utils.validateState(path);
    {
        let handle = new scran.H5File(path);
        let khandle = handle.open("kmeans_cluster");
        let kphandle = khandle.open("parameters");
        expect(kphandle.open("k", { load: true }).values[0]).toEqual(7);
        let krhandle = khandle.open("results");
        expect("clusters" in krhandle.children).toBe(false);

        let shandle = handle.open("snn_graph_cluster");
        let srhandle = shandle.open("results");
        expect("clusters" in srhandle.children).toBe(true);
    }

    // Freeing all states.
    await bakana.freeAnalysis(state);
});
