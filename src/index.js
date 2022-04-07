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
export { validateAnnotations } from "./preflight.js";
export { availableReaders } from "./utils/inputs.js";
export { setVisualizationAnimate } from "./utils/viz_parent.js";
export { setCreateLink, setResolveLink } from "./utils/reader.js";
export { setCellLabellingDownload } from "./cell_labelling.js";

import * as utils from "./utils/general.js";
import * as serialize_utils from "./utils/serialize.js";
import * as rutils from "./utils/reader.js";
import * as vizutils from "./utils/viz_parent.js";

import * as aserialize from "./abstract/serialize.js";

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
 * 
 * @return A promise that resolves to `null` when initialization is complete.
 */
export function initialize({ numberOfThreads = 1, localFile = false } = {}) {
    let s = scran.initialize({ 
        numberOfThreads: numberOfThreads,
        localFile: localFile
    });
    vizutils.scranOptions.localFile = localFile;
    return s.then(x => null); 
}

/**
 * Terminate the backend, in particular shutting down all workers.
 * This is typically necessary for a clean shutdown in Node.js applications.
 *
 * @return A promise that resolves to `null` when all workers are terminated.
 */
export function terminate() {
    let s = scran.terminate();
    let w = vizutils.killAllWorkers();
    return Promise.all([s, w]).then(x => null);
}

/**
 * Create a new analysis state in preparation for calling {@linkcode runAnalysis}.
 * Multiple states can be created and used interchangeably within the same Javascript runtime.
 *
 * @return A promise that resolves to an object containing states for all analysis steps.
 * This object can be used as input into {@linkcode runAnalysis}.
 */
export async function createAnalysis() {
    let output = {};
    output[step_inputs] = new inputs.State;
    output[step_qc] = new qc.State(output[step_inputs]);
    output[step_norm] = new normalization.State(output[step_qc]);
    output[step_feat] = new variance.State(output[step_qc], output[step_norm]);
    output[step_pca] = new pca.State(output[step_qc], output[step_norm], output[step_feat]);
    output[step_neighbors] = new index.State(output[step_pca]);

    output[step_tsne] = new tsne.State(output[step_neighbors]);
    output[step_umap] = new umap.State(output[step_neighbors]);

    output[step_kmeans] = new kmeans_cluster.State(output[step_pca]);
    output[step_snn] = new snn_cluster.State(output[step_neighbors]);
    output[step_choice] = new cluster_choice.State(output[step_snn], output[step_kmeans]);
    output[step_markers] = new cluster_markers.State(output[step_qc], output[step_norm], output[step_choice]);
    output[step_labels] = new label_cells.State(output[step_inputs], output[step_markers]);
    output[step_custom] = new custom_markers.State(output[step_qc], output[step_norm]);

    return Promise.all([output[step_tsne].ready(), output[step_umap].ready()]).then(val => output);
}

/**
 * Free the contents of an analysis state.
 * This releases memory on the **scran.js** Wasm heap and terminates any workers associated with this analysis.
 *
 * @param state An existing analysis state, produced by {@linkcode createAnalysis} or {@linkcode loadAnalysis}.
 *
 * @return A promise that resolves to `null` when all states are freed.
 */
export function freeAnalysis(state) {
    let promises = [];
    for (const [k, v] of Object.entries(state)) {
        let p = v.free();
        if (p !== null) {
            promises.push(p); 
        }
    }
    return Promise.all(promises).then(x => null);
}

/**
 * Run a basic single-cell RNA-seq analysis with the specified files and parameters.
 * This will cache the results from each step so that, if the parameters change, only the affected steps will be rerun.
 *
 * @param {object} state - Object containing the analysis state, produced by {@linkcode createAnalysis} or {@linkcode loadAnalysis}.
 * @param {Array} matrices - Object where each (arbitrarily named) property corresponds to an input matrix. 
 * Each matrix should be an object with `type` string property and any number of additional properties referring to individual data files.
 *
 * - If `type: "MatrixMarket"`, the object should contain an `mtx` property, referring to a (possibly Gzipped) Matrix Market file containing a count matrix.
 *   The object may contain a `genes` property, referring to a (possibly Gzipped) tab-separated file with the gene ID and symbols for each row of the count matrix.
 *   The object may contain a `annotation` property, referring to a (possibly Gzipped) tab-separated file with the gene ID and symbols for each row of the count matrix.
 * - If `type: "10X"`, the object should contain an `h5` property, referring to a HDF5 file following the 10X Genomics feature-barcode matrix format.
 *   It is assumed that the matrix has already been filtered to contain only the cell barcodes.
 * - If `type: "H5AD"`, the object should contain an `h5` property, referring to a H5AD file.
 *
 * The representation of each reference to a data file depends on the runtime.
 * In the browser, each data file manifests as a `File` object; for Node.js, each data file should be represented as a string containing a file path.
 *
 * Alternatively, `matrices` may be `null`, in which case the count matrices are extracted from `state`.
 * This assumes that the data matrices were already cached in `state`, either from a previous call to {@linkcode runAnalysis} or from @{linkcode loadAnalysis}.
 * @param {object} params - An object containing parameters for all steps.
 * See {@linkcode analysisDefaults} for more details.
 * @param {object} [options] - Optional parameters.
 * @param {function} [options.finishFun] - Function that is called on successful execution of each step.
 * This should accept two arguments - the name of the step and an object containing the results of that step.
 * If `null`, a no-op function is automatically created.
 * 
 * @return A promise that resolves to `null` when all asynchronous analysis steps are complete.
 * The contents of `state` are modified by reference to reflect the latest state of the analysis with the supplied parameters.
 */
export function runAnalysis(state, matrices, params, { finishFun = null } = {}) {
    let quickFun = step => {
        if (state[step].changed && finishFun !== null) {
            finishFun(step, state[step].results());
        }
    }

    let promises = [];

    state[step_inputs].compute(
        matrices, 
        params[step_inputs]["sample_factor"]
    );
    quickFun(step_inputs);

    state[step_qc].compute(
        params[step_qc]["use_mito_default"], 
        params[step_qc]["mito_prefix"], 
        params[step_qc]["nmads"]
    );
    quickFun(step_qc);
 
    state[step_norm].compute();
    quickFun(step_norm);

    state[step_feat].compute(
        params[step_feat]["span"]
    );
    quickFun(step_feat);

    state[step_pca].compute(
        params[step_pca]["num_hvgs"],
        params[step_pca]["num_pcs"],
        params[step_pca]["block_method"]
    );
    quickFun(step_pca);

    state[step_neighbors].compute(
        params[step_neighbors]["approximate"]
    );
    quickFun(step_neighbors);

    state[step_tsne].compute(
        params[step_tsne]["perplexity"],
        params[step_tsne]["iterations"], 
        params[step_tsne]["animate"]
    );
    if (state[step_tsne].changed && finishFun !== null) {
        promises.push(state[step_tsne].results().then(res => finishFun(step_tsne, res)));
    }

    state[step_umap].compute(
        params[step_umap]["num_neighbors"], 
        params[step_umap]["num_epochs"], 
        params[step_umap]["min_dist"], 
        params[step_umap]["animate"]
    );
    if (state[step_umap].changed && finishFun !== null) {
        promises.push(state[step_umap].results().then(res => finishFun(step_umap, res)));
    }

    let method = params[step_choice]["method"];

    state[step_kmeans].compute(
        method == "kmeans", 
        params[step_kmeans]["k"]
    );
    quickFun(step_kmeans);

    state[step_snn].compute(
        method == "snn_graph", 
        params[step_snn]["k"], 
        params[step_snn]["scheme"], 
        params[step_snn]["resolution"]
    );
    quickFun(step_snn);
  
    state[step_choice].compute(
        method
    );
    quickFun(step_choice);

    state[step_markers].compute();
    quickFun(step_markers);

    state[step_labels].compute(
        params[step_labels]["human_references"],
        params[step_labels]["mouse_references"]
    );
    if (state[step_labels].changed && finishFun !== null) {
        promises.push(state[step_labels].results().then(res => finishFun(step_labels, res)));
    }

    state[step_custom].compute();
    quickFun(step_custom);

    return Promise.all(promises).then(x => null);
}

/**
 * Save the current analysis state into a HDF5 file.
 * This HDF5 file can then be embedded into a `*.kana` file for distribution.
 *
 * @param {object} state - Object containing the analysis state, produced by {@linkcode createAnalysis} or {@linkcode loadAnalysis}.
 * If produced by {@linkcode createAnalysis}, it should have been run through {@linkcode runAnalysis} beforehand.
 * @param {string} path - Path to the output HDF5 file.
 * On browsers, this will lie inside the virtual file system of the **scran.js** module.
 * @param {object} [options] - Optional parameters.
 * @param {boolean} [options.embedded] - Whether to store information for embedded data files.
 * If `false`, links to data files are stored instead, see {@linkcode setCreateLink}.
 * 
 * @return A HDF5 file is created at `path` containing the analysis parameters and results - see https://ltla.github.io/kanaval for more details on the structure.
 * If `embedded = false`, a promise is returned that resolves to `null` when the saving is complete.
 * Otherwise, an object is returned containing:
 * - `collected`: an array of length equal to the number of data files.
 *   If `linkFun: null`, each element is an ArrayBuffer containing the file contents, which can be used to assemble an embedded `*.kana` file.
 *   Otherwise, if `linkFun` is supplied, each element is a string containing the linking identifier to the corresponding file.
 * - `total`: an integer containing the total length of all ArrayBuffers in `collected`.
 *   This will only be present if `linkFun` is not supplied.
 */
export async function saveAnalysis(state, path, { embedded = true } = {}) {
    let saver = null;
    let saved = null;

    if (embedded) {
        saved = serialize_utils.embedFiles();
        saver = saved.saver;
        delete saved.saver;
    }

    let handle = scran.createNewHDF5File(path);

    await state[step_inputs].serialize(handle, saver);

    state[step_qc].serialize(handle);

    state[step_norm].serialize(handle);

    state[step_feat].serialize(handle);

    state[step_pca].serialize(handle);

    state[step_neighbors].serialize(handle);

    await state[step_tsne].serialize(handle);

    await state[step_umap].serialize(handle);

    state[step_kmeans].serialize(handle);

    state[step_snn].serialize(handle);

    state[step_choice].serialize(handle);

    state[step_markers].serialize(handle);

    await state[step_labels].serialize(handle);

    state[step_custom].serialize(handle);

    return saved;
}

/**
 * Load an analysis state from a HDF5 state file, usually excised from a `*.kana` file.
 *
 * @param {string} path - Path to the HDF5 file containing the analysis state.
 * On browsers, this should lie inside the virtual file system of the **scran.js** module.
 * @param {function} loadFun - Function to load each embedded data file.
 * This should accept two arguments - an offset to the start of the file in the embedded file buffer, and the size of the file.
 *
 * In the browser, the function should return an ArrayBuffer containing the contents of the file.
 * For Node.js, the function should return a string containing a path to the file.
 * In both cases, the function may instead return a promise that resolves to the expected values.
 *
 * Note that this function is only used if the state file at `path` contains information for embedded files; 
 * otherwise, links are resolved using reader-specific functions (see {@linkcode setResolveLink} for the common use cases).
 * @param {object} [options] - Optional parameters.
 * @param {function} [options.finishFun] - Function that is called on after extracting results for each step.
 * This should accept two arguments - the name of the step and an object containing the results of that step.
 * If `null`, a no-op function is automatically created.
 *
 * @return An object with the following properties:
 * - `parameters`: an object containing the analysis parameters, similar to that created by {@linkcode analysisDefaults}.
 * - `state`: an object containing the loaded analysis state.
 *   This is conceptually equivalent to creating a state with {@linkcode createAnalysis} and running it through {@linkcode runAnalysis} with the associated `parameters`.
 */
export async function loadAnalysis(path, loadFun, { finishFun = null } = {}) {
    let response = {};
    let state = {};
    let handle = new scran.H5File(path);
    let quickFun = step => {
        if (finishFun !== null) {
            finishFun(step, state[step].results());
        }
    }

    let permuter;
    {
        let out = await inputs.unserialize(handle, loadFun);
        state[step_inputs] = out.state;
        response[step_inputs] = out.parameters;
        permuter = out.permuter;
        quickFun(step_inputs);
    }

    {
        let out = qc.unserialize(handle, state[step_inputs]);
        state[step_qc] = out.state;
        response[step_qc] = out.parameters;
        quickFun(step_qc);
    }

    {
        let out = normalization.unserialize(handle, state[step_qc]);
        state[step_norm] = out.state;
        response[step_norm] = out.parameters;
        quickFun(step_norm);
    }

    {
        let out = variance.unserialize(handle, permuter, state[step_qc], state[step_norm]);
        state[step_feat] = out.state;
        response[step_feat] = out.parameters;
        quickFun(step_feat);
    }

    {
        let out = pca.unserialize(handle, state[step_qc], state[step_norm], state[step_feat]);
        state[step_pca] = out.state;
        response[step_pca] = out.parameters;
        quickFun(step_pca);
    }

    {
        let out = index.unserialize(handle, state[step_pca]);
        state[step_neighbors] = out.state;
        response[step_neighbors] = out.parameters;
        quickFun(step_neighbors);
    }

    // Note that all awaits here are trivial, and just occur because results()
    // is async for general usage.  So we can chuck them in without really
    // worrying that they're blocking anything here.
    {
        let out = tsne.unserialize(handle, state[step_neighbors]);
        state[step_tsne] = out.state;
        response[step_tsne] = out.parameters;
        if (finishFun !== null) {
            finishFun(step_tsne, await state[step_tsne].results());
        }
    }

    {
        let out = umap.unserialize(handle, state[step_neighbors]);
        state[step_umap] = out.state;
        response[step_umap] = out.parameters;
        if (finishFun !== null) {
            finishFun(step_umap, await state[step_umap].results());
        }
    }

    {
        let out = kmeans_cluster.unserialize(handle, state[step_pca]);
        state[step_kmeans] = out.state;
        response[step_kmeans] = out.parameters;
        quickFun(step_kmeans);
    }

    {
        let out = snn_cluster.unserialize(handle, state[step_neighbors]);
        state[step_snn] = out.state;
        response[step_snn] = out.parameters;
        quickFun(step_snn);
    }

    {
        let out = cluster_choice.unserialize(handle, state[step_snn], state[step_kmeans]);
        state[step_choice] = out.state;
        response[step_choice] = out.parameters;
        quickFun(step_choice);
    }

    {
        let out = cluster_markers.unserialize(handle, permuter, state[step_qc], state[step_norm], state[step_choice]);
        state[step_markers] = out.state;
        response[step_markers] = out.parameters;
        quickFun(step_markers);
    }

    {
        let out = label_cells.unserialize(handle, state[step_inputs], state[step_markers]);
        state[step_labels] = out.state;
        response[step_labels] = out.parameters;
        if (finishFun !== null) {
            finishFun(step_labels, await state[step_labels].results());
        }
    }

    {
        let out = custom_markers.unserialize(handle, permuter, state[step_qc], state[step_norm]);
        state[step_custom] = out.state;
        response[step_custom] = out.parameters;
        quickFun(step_custom);
    }

    return {
        state: state,
        parameters: response
    };
}

/**
 * Create a `*.kana` file from the HDF5 state file and the various data files.
 *
 * @param {string} statePath - String containing a file path to an existing HDF5 state file.
 * This should be the same as the `path` used in {@linkcode saveAnalysis}.
 * On browsers, the path should exist inside the virtual file system of the **scran.js** module.
 * @param {Array} inputFiles - Array of files to be embedded into the `*.kana` file.
 * On Node.js, this should be an array of file paths; on browsers, this should be an array of ArrayBuffers.
 * Typically this is obtained as the resolved return value of {@linkcode saveAnalysis}.
 *
 * If `null`, it is assumed that files are linked instead of embedded.
 * @param {object} [options] - Further options. 
 * For Node.js, callers can specify `outputPath`, a string containing the output path for the newly created `*.kana` file.
 *
 * @return 
 * For Node.js, a promise is returned that resolves to a path to a new `*.kana` file. 
 * This is equal to `outputPath` if supplied, otherwise a path to a file in a temporary directory is returned.
 *
 * In browsers, an ArrayBuffer is returned containing the full contents of the new `*.kana` file.
 */
export function createKanaFile(statePath, inputFiles, options = {}) {
    return aserialize.createKanaFileInternal(statePath, inputFiles, options);
}

/**
 * Parse a `*.kana` file by extracting the HDF5 state file and returning a function to extract embeddded data files.
 *
 * @param {string|ArrayBuffer} input - The input `*.kana` file.
 * For Node.js, this should be a string containing a path to the file.
 * On browsers, this should be an ArrayBuffer containing the full file contents.
 * @param {string} statePath - String containing a file path to save the HDF5 state file.
 * This will also be the path supplied to {@linkcode loadAnalysis} to load the state into memory.
 * On browsers, this will exist inside the virtual file system of the **scran.js** module.
 * @param {object} [options] - Further options. 
 * For Node.js, callers can specify `stageDir`, a string containing a path to a staging directory for the extracted data files.
 *
 * @return The HDF5 state file is written to `statePath`.

 * In the browser, if `input` contains embedded files, a function is returned that extracts each data file given its offset and size.
 * (This should be used as `loadFun` in {@linkcode loadAnalysis}.)
 * If `input` contains linked files, `null` is returned.
 *
 * For Node.js, a promise is returned that evaluates to the aforementioned function or `null`. 
 */
export function parseKanaFile(input, statePath, options = {}) {
    return aserialize.parseKanaFileInternal(input, statePath, options);
}
