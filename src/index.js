import * as scran from "scran.js";
import * as inputs from "./inputs.js";
import * as preflight from "./preflight.js";
import * as qc from "./quality_control.js";
import * as normalization from "./normalization.js";
import * as variance from "./feature_selection.js";
import * as pca from "./pca.js";
import * as index from "./neighbor_index.js";
import * as cluster_choice from "./choose_clustering.js";
import * as kmeans_cluster from "./kmeans_cluster.js";
import * as snn_cluster from "./snn_graph_cluster.js";
import * as tsne from "./tsne.js";
import * as umap from "./umap.js";
import * as cluster_markers from "./marker_detection.js";
import * as label_cells from "./cell_labelling.js";
import * as custom_markers from "./custom_selections.js";
export { analysisDefaults } from "./defaults.js";

import * as utils from "./utils/general.js";
import * as serialize_utils from "./utils/serialize.js";
import * as rutils from "./utils/reader.js";

const step_inputs = "inputs";
const step_qc = "quality_control";
const step_norm = "normalizaton";
const step_feat = "feature_selection";
const step_pca = "pca";
const step_neighbors = "neighbor_index";
const step_tsne = "tsne";
const step_umap = "umap";
const step_kmeans = "kmeans_cluster";
const step_snn = "snn_graph_cluster";
const step_choice = "choose_clustering";
const step_markers = "marker_detection";
const step_labels = "cell_labelling";
const step_custom = "custom_selections";

/**
 * Initialize the backend for computation.
 * This is required prior to running any other **bakana** function.
 *
 * @param {object} [options] - Optional parameters.
 * @param {number} [options.numberOfThreads] - Number of threads used by **scran.js**.
 * @param {boolean} [options.localFile] - Whether to use local file paths for imported modules in **scran.js**.
 * This only needs to be `true` for old Node versions that do not support file URIs.
 * @param {function} [options.animateFun] - Function to handle animation iterations for t-SNE or UMAP.
 * If `null`, a no-op function is automatically created.
 * Otherwise, a supplied function should accept three arguments, in the following order:
 * - A `Float64Array` containing the x-coordinates for each cell.
 * - A `Float64Array` containing the y-coordinates for each cell.
 * - An integer specifying the iteration number.
 * 
 * @return A promise that resolves to `null` when initialization is complete.
 */
export function initialize({ numberOfThreads = 1, localFile = false, animateFun = null } = {}) {
    let s = scran.initialize({ 
        numberOfThreads: numberOfThreads,
        localFile: localFile
    });

    // Making a dummy iterator, if none is supplied.
    if (animateFun === null) {
        animateFun = (x, y, it) => true;
    }
    let t = tsne.initialize(animateFun, { localFile: localFile });
    let u = umap.initialize(animateFun, { localFile: localFile });

    return Promise.all([s, t, u]).then(x => null); 
}

/**
 * Terminate the backend, in particular shutting down all workers.
 * This is typically necessary for a clean shutdown in Node.js applications.
 *
 * @return A promise that resolves to `null` when all workers are terminated.
 */
export function terminate() {
    let s = scran.terminate();
    let t = tsne.terminate();
    let u = umap.terminate();
    return Promise.all([s, t, u]).then(x => null);
}

/**
 * Run a basic single-cell RNA-seq analysis with the specified files and parameters.
 * This will cache the results from each step so that, if the parameters change, only the affected steps will be rerun.
 *
 * @param {Array} matrices - An array of input matrices.
 * Each matrix is represented as an object with property `type` and additional properties referring to individual data files.
 * In the browser, each data file manifests as a `File` object; for Node.js, each data file should be represented as a string containing a file path.
 * - If `type: "MatrixMarket"`, the object should contain an `mtx` property, referring to a (possibly Gzipped) Matrix Market file containing a count matrix.
 *   The object may contain a `genes` property, referring to a (possibly Gzipped) tab-separated file with the gene ID and symbols for each row of the count matrix.
 *   The object may contain a `annotation` property, referring to a (possibly Gzipped) tab-separated file with the gene ID and symbols for each row of the count matrix.
 * - If `type: "10X"`, the object should contain an `h5` property, referring to a HDF5 file following the 10X Genomics feature-barcode matrix format.
 *   It is assumed that the matrix has already been filtered to contain only the cell barcodes.
 * - If `type: "H5AD"`, the object should contain an `h5` property, referring to a H5AD file.
 * @param {object} params - An object containing parameters for all steps.
 * See {@linkcode analysisDefaults} for more details.
 * @param {function} downloadFun - Function that accepts a single URL, downloads the requested resource, and returns an ArrayBuffer of the contents.
 * @param {object} [options] - Optional parameters.
 * @param {function} [options.finishFun] - Function that is called on successful execution of each step.
 * This should accept two arguments - the name of the step and an object containing the results of that step.
 * If `null`, a no-op function is automatically created.
 * 
 * @return A promise that resolves to `null` when all asynchronous analysis steps are complete.
 */
export function runAnalysis(matrices, params, downloadFun, { finishFun = null } = {}) {
    if (finishFun === null) {
        finishFun = (name, results) => true;
    }

    let promises = [];

    inputs.compute(
        matrices, 
        params[step_inputs]["sample_factor"]
    );
    if (inputs.changed) {
        finishFun(step_inputs, inputs.results());
    }

    qc.compute(
        params[step_qc]["use_mito_default"], 
        params[step_qc]["mito_prefix"], 
        params[step_qc]["nmads"]
    );
    if (qc.changed) {
        finishFun(step_qc, qc.results());
    }
 
    normalization.compute();
    if (normalization.changed) {
        finishFun(step_norm, normalization.results());
    }

    variance.compute(
        params[step_feat]["span"]
    );
    if (variance.changed) {
        finishFun(step_feat, variance.results());
    }

    pca.compute(
        params[step_pca]["num_hvgs"],
        params[step_pca]["num_pcs"],
        params[step_pca]["block_method"]
    );
    if (pca.changed) {
        finishFun(step_pca, pca.results());
    }

    index.compute(
        params[step_neighbors]["approximate"]
    );
    if (index.changed) {
        finishFun(step_neighbors, index.results());
    }

    tsne.compute(
        params[step_tsne]["perplexity"],
        params[step_tsne]["iterations"], 
        params[step_tsne]["animate"]
    );
    if (tsne.changed) {
        promises.push(tsne.results().then(res => finishFun(step_tsne, res)));
    }

    umap.compute(
        params[step_umap]["num_neighbors"], 
        params[step_umap]["num_epochs"], 
        params[step_umap]["min_dist"], 
        params[step_umap]["animate"]
    );
    if (umap.changed) {
        promises.push(umap.results().then(res => finishFun(step_umap, res)));
    }

    let method = params[step_choice]["method"];

    kmeans_cluster.compute(
        method == "kmeans", 
        params[step_kmeans]["k"]
    );
    if (kmeans_cluster.changed) {
        finishFun(step_kmeans, kmeans_cluster.results());
    }

    snn_cluster.compute(
        method == "snn_graph", 
        params[step_snn]["k"], 
        params[step_snn]["scheme"], 
        params[step_snn]["resolution"]
    );
    if (snn_cluster.changed) {
        finishFun(step_snn, snn_cluster.results());
    }
  
    cluster_choice.compute(
        method
    );
    if (cluster_choice.changed) {
        finishFun(step_choice, cluster_choice.results());
    }

    cluster_markers.compute();
    if (cluster_markers.changed) {
        finishFun(step_markers, cluster_markers.results());
    }

    label_cells.compute(
        params[step_labels]["human_references"],
        params[step_labels]["mouse_references"],
        downloadFun
    );
    if (label_cells.changed) {
        promises.push(label_cells.results().then(res => finishFun(step_labels, res)));
    }

    custom_markers.compute();
    if (custom_markers.changed) {
        finishFun(step_custom, custom_markers.results());
    }

    return Promise.all(promises).then(x => null);
}

/**
 * Save the current analysis state into a HDF5 file.
 * This HDF5 file can then be embedded into a `*.kana` file for distribution.
 *
 * @param {string} path - Path to the output HDF5 file.
 * On browsers, this will lie inside the virtual file system of the **scran.js** module.
 * @param {object} [options] - Optional parameters.
 * @param {?function} [options.linkFun] - Function that returns a linking idenfier to a data file.
 * Each link is stored in the state file, which is ultimately used to create a `*.kana` file where the data files are linked rather than embedded.
 * The function should accept an ArrayBuffer containing the file content (for browsers) or a string containing the file path (for Node.js),
 * and return a string containing some unique identifier to the file.
 * This is most typically used to register the file with some user-specified database system for later retrieval.
 * If `null`, we assume that embedding of the file is desired instead.
 * 
 * @return A HDF5 file is created at `path` containing the analysis parameters and results - see https://ltla.github.io/kanaval for more details on the structure.
 * An object is returned containing:
 * - `collected`: an array of length equal to the number of data files.
 *   If `linkFun: null`, each element is an ArrayBuffer containing the file contents, which can be used to assemble an embedded `*.kana` file.
 *   Otherwise, if `linkFun` is supplied, each element is a string containing the linking identifier to the corresponding file.
 * - `total`: an integer containing the total length of all ArrayBuffers in `collected`.
 *   This will only be present if `linkFun` is not supplied.
 */
export async function saveAnalysis(path, { linkFun = null } = {}) {
    let saver;
    let saved;

    if (linkFun === null) {
        saved = serialize_utils.embedFiles();
        saver = saved.saver;
        delete saved.saver;
    } else {
        saved = { collected: [] };
        saver = (obj) => {
            let id = linkFun(obj);
            saved.collected.push(id);
            return id;
        }
    }

    let handle = scran.createNewHDF5File(path);

    await inputs.serialize(handle, saver);

    qc.serialize(handle);

    normalization.serialize(handle);

    variance.serialize(handle);

    pca.serialize(handle);

    index.serialize(handle);

    await tsne.serialize(handle);

    await umap.serialize(handle);

    kmeans_cluster.serialize(handle);

    snn_cluster.serialize(handle);

    cluster_choice.serialize(handle);

    cluster_markers.serialize(handle);

    await label_cells.serialize(handle);

    custom_markers.serialize(handle);

    return saved;
}

/**
 * Save the current analysis state into a HDF5 file.
 * This HDF5 file can then be embedded into a `*.kana` file for distribution.
 *
 * @param {string} path - Path to the HDF5 file containing the analysis state.
 * On browsers, this should lie inside the virtual file system of the **scran.js** module.
 * @param {function} loadFun - Function to load the individual data file.
 * - If the state file contains links to files, the function should accept a single string.
 *   Otherwise, the function should accept two arguments - an offset to the start of the file in the embedded file buffer, and the size of the file.
 * - In the browser, the function should return an ArrayBuffer containing the contents of the file.
 *   For Node.js, the function should return a string containing a path to the file.
 * @param {object} [options] - Optional parameters.
 * @param {function} [options.finishFun] - Function that is called on after extracting results for each step.
 * This should accept two arguments - the name of the step and an object containing the results of that step.
 * If `null`, a no-op function is automatically created.
 *
 * @return An object containing the analysis parameters, similar to that created by {@linkcode analysisDefaults}.
 */
export async function loadAnalysis(path, loadFun, { finishFun = null } = {}) {
    if (finishFun === null) {
        finishFun = (name, results) => true;
    }

    let response = {};
    let handle = new scran.H5File(path);

    let permuter = await inputs.unserialize(handle, loadFun);
    finishFun(step_inputs, inputs.results());

    response[step_qc] = qc.unserialize(handle);
    finishFun(step_qc, qc.results());

    response[step_norm] = normalization.unserialize(handle);
    finishFun(step_norm, normalization.results());

    response[step_feat] = variance.unserialize(handle, permuter);
    finishFun(step_feat, variance.results());

    response[step_pca] = pca.unserialize(handle);
    finishFun(step_pca, pca.results());

    response[step_neighbors] = index.unserialize(handle);
    finishFun(step_neighbors, index.results());

    // Note that all awaits here are trivial, and just occur because results()
    // is async for general usage.  So we can chuck them in without really
    // worrying that they're blocking anything here.
    response[step_tsne] = tsne.unserialize(handle);
    finishFun(step_tsne, await tsne.results());

    response[step_umap] = umap.unserialize(handle);
    finishFun(step_umap, await umap.results());

    response[step_kmeans] = kmeans_cluster.unserialize(handle);
    finishFun(step_kmeans, kmeans_cluster.results());

    response[step_snn] = snn_cluster.unserialize(handle);
    finishFun(step_snn, snn_cluster.results());

    response[step_choice] = cluster_choice.unserialize(handle);
    finishFun(step_choice, cluster_choice.results());

    response[step_markers] = cluster_markers.unserialize(handle, permuter);
    finishFun(step_markers, cluster_markers.results());

    response[step_labels] = label_cells.unserialize(handle);
    finishFun(step_labels, await label_cells.results());

    response[step_custom] = custom_markers.unserialize(handle, permuter);
    finishFun(step_custom, custom_markers.results()); 

    return response;
}
