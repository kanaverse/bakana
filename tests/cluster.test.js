import * as bakana from "../src/index.js";
import * as scran from "scran.js";
import * as utils from "./utils.js"

beforeAll(utils.initializeAll);
afterAll(async () => await bakana.terminate());

test("switching between clustering methods (SNN first)", async () => {
    let files = {
        default: new bakana.TenxHdf5Dataset("files/datasets/pbmc4k-tenx.h5")
    };

    // First running the analysis with SNN graph,
    // and confirming that only SNN graph results are reported.
    let state = await bakana.createAnalysis();
    let paramcopy = utils.baseParams();
    await bakana.runAnalysis(state, files, paramcopy);
    expect(state.snn_graph_cluster.fetchClusters().length).toBeGreaterThan(0);
    expect(() => state.kmeans_cluster.fetchClusters()).toThrow("cannot fetch");

    // Now trying with k-means. This should cause both sets of
    // results to be saved, as both clusterings are still valid.
    paramcopy.choose_clustering.method = "kmeans";

    await bakana.runAnalysis(state, files, paramcopy);
    expect(state.rna_pca.changed).toBe(false);
    expect(state.snn_graph_cluster.changed).toBe(false);
    expect(state.kmeans_cluster.changed).toBe(true);
    expect(state.choose_clustering.changed).toBe(true);
    expect(state.snn_graph_cluster.fetchClusters().length).toBeGreaterThan(0);
    expect(state.kmeans_cluster.fetchClusters().length).toBeGreaterThan(0);

    // Checking that invalidation of the results behaves correctly.
    // If we change the parameters but we're still on the old set of
    // results, the SNN graph results don't get rerun and the results get wiped.
    paramcopy.snn_graph_cluster.multilevel_resolution = 0.77;

    await bakana.runAnalysis(state, files, paramcopy);
    expect(state.rna_pca.changed).toBe(false);
    expect(state.snn_graph_cluster.changed).toBe(true);
    expect(state.kmeans_cluster.changed).toBe(false);
    expect(state.choose_clustering.changed).toBe(false);
    expect(() => state.snn_graph_cluster.fetchClusters()).toThrow("cannot fetch");
    expect(state.kmeans_cluster.fetchClusters().length).toBeGreaterThan(0);

    // Freeing all states.
    await bakana.freeAnalysis(state);
});

let files_3k = {
    default: new bakana.TenxMatrixMarketDataset(
        "files/datasets/pbmc3k-matrix.mtx.gz", 
        "files/datasets/pbmc3k-features.tsv.gz",
        null
    )
};

test("switching between clustering methods (k-means first)", async () => {
    let state = await bakana.createAnalysis();
    let paramcopy = utils.baseParams();
    paramcopy.choose_clustering = { method: "kmeans" };

    await bakana.runAnalysis(state, files_3k, paramcopy);
    expect(state.snn_graph_cluster.changed).toBe(true);
    expect(state.kmeans_cluster.changed).toBe(true);
    expect(() => state.snn_graph_cluster.fetchClusters()).toThrow("cannot fetch");
    expect(state.kmeans_cluster.fetchClusters().length).toBeGreaterThan(0);

    // Trying with a different clustering method.
    paramcopy.choose_clustering.method = "snn_graph"; 

    await bakana.runAnalysis(state, files_3k, paramcopy);
    expect(state.snn_graph_cluster.changed).toBe(true);
    expect(state.kmeans_cluster.changed).toBe(false);
    expect(state.snn_graph_cluster.fetchClusters().length).toBeGreaterThan(0);
    expect(state.kmeans_cluster.fetchClusters().length).toBeGreaterThan(0);

    // Checking that invalidation works.
    paramcopy.kmeans_cluster.k = 7;

    await bakana.runAnalysis(state, files_3k, paramcopy);
    expect(state.snn_graph_cluster.changed).toBe(false);
    expect(state.kmeans_cluster.changed).toBe(true);
    expect(state.snn_graph_cluster.fetchClusters().length).toBeGreaterThan(0);
    expect(() => state.kmeans_cluster.fetchClusters().length).toThrow("cannot fetch");

    // Freeing all states.
    await bakana.freeAnalysis(state);
});

test("trying some of the other SNN algorithms", async () => {
    let state = await bakana.createAnalysis();
    let paramcopy = utils.baseParams();

    paramcopy.snn_graph_cluster.algorithm = "walktrap";
    await bakana.runAnalysis(state, files_3k, paramcopy);

    let copy = state.snn_graph_cluster.fetchClusters().slice();
    expect(copy.length).toEqual(state.cell_filtering.fetchFilteredMatrix().numberOfColumns());
    expect(state.snn_graph_cluster.fetchParameters().algorithm).toEqual("walktrap");
    expect(state.marker_detection.fetchResults().RNA.numberOfGroups()).toBeGreaterThan(1);

    // Also trying with Leiden.
    paramcopy.snn_graph_cluster.algorithm = "leiden";
    await bakana.runAnalysis(state, files_3k, paramcopy);

    expect(state.neighbor_index.changed).toBe(false);
    expect(state.tsne.changed).toBe(false);
    expect(state.marker_detection.changed).toBe(true);
    expect(state.snn_graph_cluster.changed).toBe(true);

    let copy2 = state.snn_graph_cluster.fetchClusters().slice();
    expect(copy.length).toEqual(copy2.length);
    expect(copy).not.toEqual(copy2);

    // Freeing states.
    await bakana.freeAnalysis(state);
})

