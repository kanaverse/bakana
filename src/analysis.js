import * as scran from "scran.js";

import * as inputs from "./steps/inputs.js";

import * as qc from "./steps/rna_quality_control.js";
import * as qcadt from "./steps/adt_quality_control.js";
import * as qccrispr from "./steps/crispr_quality_control.js";
import * as filters from "./steps/cell_filtering.js";

import * as normalization from "./steps/rna_normalization.js";
import * as normadt from "./steps/adt_normalization.js";
import * as normcrispr from "./steps/crispr_normalization.js";

import * as variance from "./steps/feature_selection.js";

import * as pca from "./steps/rna_pca.js";
import * as pcaadt from "./steps/adt_pca.js";
import * as pcacrispr from "./steps/crispr_pca.js";
import * as combine from "./steps/combine_embeddings.js";
import * as correct from "./steps/batch_correction.js";

import * as index from "./steps/neighbor_index.js";
import * as cluster_choice from "./steps/choose_clustering.js";
import * as kmeans_cluster from "./steps/kmeans_cluster.js";
import * as snn_cluster from "./steps/snn_graph_cluster.js";

import * as tsne from "./steps/tsne.js";
import * as umap from "./steps/umap.js";

import * as cluster_markers from "./steps/marker_detection.js";
import * as label_cells from "./steps/cell_labelling.js";
import * as custom_markers from "./steps/custom_selections.js";

import * as feature_set_enrichment from "./steps/feature_set_enrichment.js";

import { FORMAT_VERSION } from "./abstract/utils/serialize.js";
import { bakana_version } from "./version.js";

export { setCreateLink, setResolveLink } from "./steps/inputs.js";
export { MarkerDetectionState } from "./steps/marker_detection.js";
export { CustomSelectionsState } from "./steps/custom_selections.js";

const step_inputs = inputs.step_name;
const step_qc = qc.step_name;
const step_qc_adt = qcadt.step_name;
const step_qc_crispr = qccrispr.step_name;
const step_filter = filters.step_name;
const step_norm = normalization.step_name;
const step_norm_adt = normadt.step_name;
const step_norm_crispr = normcrispr.step_name;
const step_feat = "feature_selection";
const step_pca = pca.step_name;
const step_pca_adt = pcaadt.step_name;
const step_pca_crispr = pcacrispr.step_name;
const step_combine = "combine_embeddings";
const step_correct = "batch_correction";
const step_neighbors = index.step_name;
const step_tsne = "tsne";
const step_umap = "umap";
const step_kmeans = "kmeans_cluster";
const step_snn = snn_cluster.step_name;
const step_choice = "choose_clustering";
const step_markers = cluster_markers.step_name;
const step_labels = "cell_labelling";
const step_custom = custom_markers.step_name;
const step_enrichment = feature_set_enrichment.step_name;

const load_flag = "_loaded";

/**
 * Create a new analysis state in preparation for calling {@linkcode runAnalysis}.
 * Multiple states can be created and used interchangeably within the same Javascript runtime.
 *
 * @return A promise that resolves to an object containing states for all analysis steps.
 * This object can be used as input into {@linkcode runAnalysis}.
 */
export async function createAnalysis() {
    return create_analysis(new inputs.InputsState);
}

function create_analysis(input_state) {
    let output = {};
    output[step_inputs] = input_state;

    output[step_qc] = new qc.RnaQualityControlState(output[step_inputs]);
    output[step_qc_adt] = new qcadt.AdtQualityControlState(output[step_inputs]);
    output[step_qc_crispr] = new qccrispr.CrisprQualityControlState(output[step_inputs]);

    let qc_states = { "RNA": output[step_qc], "ADT": output[step_qc_adt], "CRISPR": output[step_qc_crispr] }
    output[step_filter] = new filters.CellFilteringState(output[step_inputs], qc_states);

    output[step_norm] = new normalization.RnaNormalizationState(output[step_qc], output[step_filter]);
    output[step_norm_adt] = new normadt.AdtNormalizationState(output[step_qc_adt], output[step_filter]);
    output[step_norm_crispr] = new normcrispr.CrisprNormalizationState(output[step_qc_crispr], output[step_filter]);

    output[step_feat] = new variance.FeatureSelectionState(output[step_filter], output[step_norm]);

    output[step_pca] = new pca.RnaPcaState(output[step_filter], output[step_norm], output[step_feat]);
    output[step_pca_adt] = new pcaadt.AdtPcaState(output[step_filter], output[step_norm_adt]);
    output[step_pca_crispr] = new pcacrispr.CrisprPcaState(output[step_filter], output[step_norm_crispr]);

    let pca_states = { "RNA": output[step_pca], "ADT": output[step_pca_adt], "CRISPR": output[step_pca_crispr] }
    output[step_combine] = new combine.CombineEmbeddingsState(pca_states);
    output[step_correct] = new correct.BatchCorrectionState(output[step_filter], output[step_combine]);

    output[step_neighbors] = new index.NeighborIndexState(output[step_correct]);

    output[step_tsne] = new tsne.TsneState(output[step_neighbors]);
    output[step_umap] = new umap.UmapState(output[step_neighbors]);

    output[step_kmeans] = new kmeans_cluster.KmeansClusterState(output[step_correct]);
    output[step_snn] = new snn_cluster.SnnGraphClusterState(output[step_neighbors]);
    output[step_choice] = new cluster_choice.ChooseClusteringState(output[step_snn], output[step_kmeans]);

    let norm_states = { "RNA": output[step_norm], "ADT": output[step_norm_adt], "CRISPR": output[step_norm_crispr] };
    output[step_markers] = new cluster_markers.MarkerDetectionState(output[step_filter], norm_states, output[step_choice]);
    output[step_labels] = new label_cells.CellLabellingState(output[step_inputs], output[step_markers]);
    output[step_custom] = new custom_markers.CustomSelectionsState(output[step_filter], norm_states);

    output[step_enrichment] = new feature_set_enrichment.FeatureSetEnrichmentState(output[step_inputs], output[step_filter], output[step_norm], output[step_markers]);

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
        if (k == load_flag) {
            continue;
        }
        let p = v.free();
        if (p) { // not null, not undefined.
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
 * @param {object} datasets - Object where each (arbitrarily named) property corresponds to an input dataset.
 * Each dataset should be a object that satisfies the {@linkplain Dataset} contract.
 *
 * Alternatively, `datasets` may be `null` if the input datasets were already loaded and cached in `state`.
 * This avoids the need to respecify the inputs after a previous call to {@linkcode runAnalysis} or from {@linkcode loadAnalysis}.
 * @param {object} params - An object containing parameters for all steps.
 * See {@linkcode analysisDefaults} for more details.
 * @param {object} [options] - Optional parameters.
 * @param {function} [options.startFun] - Function that is called when each step is started.
 * This should accept a single argument - the name of the step.
 * The return value is ignored, but any promises will be awaited before the analysis proceeds to the next step.
 * If `null`, nothing is executed.
 * @param {function} [options.finishFun] - Function that is called on successful execution of each step.
 * This should accept a single argument - the name of the step.
 * The return value is ignored, but any promises will be awaited before the analysis proceeds to the next step.
 * If `null`, nothing is executed.
 * 
 * @return A promise that resolves to `null` when all asynchronous analysis steps are complete.
 * The contents of `state` are modified by reference to reflect the latest state of the analysis with the supplied parameters.
 */
export async function runAnalysis(state, datasets, params, { startFun = null, finishFun = null } = {}) {
    let quickStart = async step => {
        if (startFun !== null) {
            await startFun(step);
        }
    }

    let quickFinish = async step => {
        if (finishFun !== null) {
            await finishFun(step);
        }
    }

    let promises = [];
    let deferredQuickFinish = (step, p) => {
        if (finishFun !== null) {
            if (state[step].changed) {
                p = p.then(out => finishFun(step, state[step]));
            } else {
                p = p.then(out => finishFun(step));
            }
        }
        promises.push(p);
    }

    /*** Loading ***/
    await quickStart(step_inputs);
    await state[step_inputs].compute(
        datasets,
        params[step_inputs]["block_factor"],
        params[step_inputs]["subset"]
    );
    await quickFinish(step_inputs);

    if (load_flag in state) {
        // Force recompute for all downstream steps. This avoids mixing results
        // from different versions if we're re-running off a reloaded state; if
        // some steps rerun, but others don't, we end up with a bastard state
        // from possibly different versions of this pipeline. It's also
        // difficult to guarantee that enough results were saved for use in
        // downstream steps, given that not everything is saved to file (and
        // indeed, the requirements of downstream steps may change in future
        // versions). So we just keep it simple and flush the whole state.
        state[step_inputs].changed = true;
        delete state[load_flag];
    }

    /*** Quality control ***/
    await quickStart(step_qc);
    state[step_qc].compute(
        params[step_qc]["use_mito_default"], 
        params[step_qc]["mito_prefix"], 
        params[step_qc]["nmads"]
    );
    await quickFinish(step_qc);

    await quickStart(step_qc_adt);
    state[step_qc_adt].compute(
        params[step_qc_adt]["igg_prefix"], 
        params[step_qc_adt]["nmads"],
        params[step_qc_adt]["min_detected_drop"]
    );
    await quickFinish(step_qc_adt);

    await quickStart(step_qc_crispr);
    state[step_qc_crispr].compute(
        params[step_qc_crispr]["nmads"]
    );
    await quickFinish(step_qc_crispr);

    await quickStart(step_filter);
    state[step_filter].compute(
        params[step_filter]["use_rna"], 
        params[step_filter]["use_adt"],
        params[step_filter]["use_crispr"]
    );
    await quickFinish(step_filter);

    /*** Normalization ***/
    await quickStart(step_norm);
    state[step_norm].compute();
    await quickFinish(step_norm);

    await quickStart(step_norm_adt);
    state[step_norm_adt].compute(
        params[step_norm_adt]["num_pcs"],    
        params[step_norm_adt]["num_clusters"]    
    );
    await quickFinish(step_norm_adt);

    await quickStart(step_norm_crispr);
    state[step_norm_crispr].compute();
    await quickFinish(step_norm_crispr);

    /*** Feature selection ***/
    await quickStart(step_feat);
    state[step_feat].compute(
        params[step_feat]["span"]
    );
    await quickFinish(step_feat);
  
    /*** Dimensionality reduction ***/
    await quickStart(step_pca);
    state[step_pca].compute(
        params[step_pca]["num_hvgs"],
        params[step_pca]["num_pcs"],
        params[step_pca]["block_method"]
    );
    await quickFinish(step_pca);

    await quickStart(step_pca_adt);
    state[step_pca_adt].compute(
        params[step_pca_adt]["num_pcs"],
        params[step_pca_adt]["block_method"]
    );
    await quickFinish(step_pca_adt);

    await quickStart(step_pca_crispr);
    state[step_pca_crispr].compute(
        params[step_pca_crispr]["num_pcs"],
        params[step_pca_crispr]["block_method"]
    );
    await quickFinish(step_pca_crispr);

    await quickStart(step_combine);
    state[step_combine].compute(
        params[step_combine]["rna_weight"],
        params[step_combine]["adt_weight"],
        params[step_combine]["crispr_weight"],
        params[step_combine]["approximate"]
    );
    await quickFinish(step_combine);

    await quickStart(step_correct);
    state[step_correct].compute(
        params[step_correct]["method"],
        params[step_correct]["num_neighbors"],
        params[step_correct]["approximate"]
    );
    await quickFinish(step_correct);

    /*** Nearest neighbors ***/
    await quickStart(step_neighbors);
    state[step_neighbors].compute(
        params[step_neighbors]["approximate"]
    );
    await quickFinish(step_neighbors);

    /*** Visualization ***/
    {
        await quickStart(step_tsne);
        let p = state[step_tsne].compute(
            params[step_tsne]["perplexity"],
            params[step_tsne]["iterations"], 
            params[step_tsne]["animate"]
        );
        deferredQuickFinish(step_tsne, p);
    }

    {
        await quickStart(step_umap);
        let p = state[step_umap].compute(
            params[step_umap]["num_neighbors"], 
            params[step_umap]["num_epochs"], 
            params[step_umap]["min_dist"], 
            params[step_umap]["animate"]
        );
        deferredQuickFinish(step_umap, p);
    }

    /*** Clustering ***/
    let method = params[step_choice]["method"];

    await quickStart(step_kmeans);
    state[step_kmeans].compute(
        method == "kmeans", 
        params[step_kmeans]["k"]
    );
    await quickFinish(step_kmeans);

    await quickStart(step_snn);
    state[step_snn].compute(
        method == "snn_graph", 
        params[step_snn]["k"], 
        params[step_snn]["scheme"], 
        params[step_snn]["algorithm"],
        params[step_snn]["multilevel_resolution"],
        params[step_snn]["leiden_resolution"],
        params[step_snn]["walktrap_steps"]
    );
    await quickFinish(step_snn);

    await quickStart(step_choice);
    state[step_choice].compute(
        method
    );
    await quickFinish(step_choice);

    /*** Markers and labels ***/
    await quickStart(step_markers);
    state[step_markers].compute(
        params[step_markers]["lfc_threshold"],
        params[step_markers]["compute_auc"]
    );
    await quickFinish(step_markers);

    {
        await quickStart(step_labels);
        let p = state[step_labels].compute(
            params[step_labels]["human_references"],
            params[step_labels]["mouse_references"]
        );
        deferredQuickFinish(step_labels, p);
    }

    await quickStart(step_custom);
    state[step_custom].compute(
        params[step_custom]["lfc_threshold"],
        params[step_custom]["compute_auc"]
    );
    await quickFinish(step_custom);

    await quickStart(step_enrichment);
    await state[step_enrichment].compute(
        params[step_enrichment]["collections"],
        params[step_enrichment]["dataset_id_column"],
        params[step_enrichment]["reference_id_column"],
        params[step_enrichment]["top_markers"]
    );
    await quickFinish(step_enrichment);

    await Promise.all(promises);
    return null;
}

/**
 * Retrieve analysis parameters from a state object.
 *
 * @param {object} state - Object containing the analysis state, produced by {@linkcode createAnalysis} or {@linkcode loadAnalysis}.
 *
 * @return {object} Object containing the analysis parameters for each step, similar to that created by {@linkcode analysisDefaults}.
 */
export function retrieveParameters(state) {
    let params = {};
    for (const [k, v] of Object.entries(state)) {
        if (k == load_flag) {
            continue;
        }
        params[k] = v.fetchParameters();
    }
    return params;
}

/**
 * Create a new analysis state object consisting of a subset of cells from an existing analysis state.
 * This assumes that the existing state already contains loaded matrix data in its `inputs` property,
 * which allows us to create a cheap reference without reloading the data into memory.
 *
 * @param {object} state - State object such as that produced by {@linkcode createAnalysis} or {@linkcode linkAnalysis}.
 * This should already contain loaded data, e.g., after a run of {@linkcode runAnalysis}.
 * @param {TypedArray|Array} indices - Array containing the indices for the desired subset of cells.
 * This should be sorted and non-duplicate.
 * Any existing subset in `state` will be overridden by `indices`.
 * @param {object} [options] - Optional parameters.
 * @param {boolean} [options.copy=true] - Whether to make a copy of `indices` before storing it inside the returned state object.
 * If `false`, it is assumed that the caller makes no further use of the passed `indices`.
 * @param {boolean} [options.onOriginal=false] - Whether `indices` contains indices on the original dataset or on the dataset in `state`.
 * This distinction is only relevant if `state` itself contains an analysis of a subsetted dataset.
 * If `false`, the `indices` are assumed to refer to the columns of the already-subsetted dataset that exists in `state`;
 * if `true`, the `indices` are assumed to refer to the columns of the original dataset from which the subset in `state` was created.
 *
 * @return {object} A state object containing loaded matrix data in its `inputs` property.
 * Note that the other steps do not have any results, so this object should be passed through {@linkcode runAnalysis} before it can be used.
 */
export async function subsetInputs(state, indices, { copy = true, onOriginal = false } = {}) {
    return create_analysis(state.inputs.createDirectSubset(indices, { copy: copy, onOriginal: onOriginal }));
}

/****************************
 ********* LEGACY ***********
 ****************************/

/**
 * Load an analysis state from a HDF5 state file, usually excised from a `*.kana` file.
 *
 * @param {string} path - Path to the HDF5 file containing the analysis state.
 * On browsers, this should lie inside the virtual file system of the **scran.js** module.
 * @param {function} loadFun - Function to load each embedded data file.
 * This function is typically generated by {@linkcode parseKanaFile} and has the following characteristics:
 *
 * - It should accept two arguments, `offset` and `size`.
 *   `offset` is a number containing the offset to the start of any given file in the embedded file buffer.
 *   `size` is, of course, the size of the file.
 * - It should return any argument that can be used in the {@linkplain SimpleFile} constructor.
 *   This may be a Uint8Array, a File (on browsers) or a string containing a file path (on Node.js).
 *   Alternatively, the function may return a promise that resolves to the expected values.
 *
 * Note that this function is only used if the state file at `path` contains information for embedded files; 
 * otherwise, links are resolved using reader-specific functions (see {@linkcode setResolveLink} for the common use cases).
 * @param {object} [options] - Optional parameters.
 * @param {function} [options.finishFun] - Function that is called on after extracting results for each step.
 * This should accept two arguments - the name of the step and an object containing the results of that step.
 * The function may optionally be `async`.
 * If `null`, a no-op function is automatically created.
 *
 * @return An object containing the loaded analysis state.
 * This is conceptually equivalent to creating a state with {@linkcode createAnalysis} and running it through {@linkcode runAnalysis}.
 */
export async function loadAnalysis(path, loadFun, { finishFun = null } = {}) {
    let state = {};
    let handle = new scran.H5File(path);
    let quickFun = async step => {
        if (finishFun !== null) {
            await finishFun(step);
        }
    }

    /*** Loading ***/
    let permuters;
    {
        let out = await inputs.unserialize(handle, loadFun);
        state[step_inputs] = out.state;
        permuters = out.permuters;
        await quickFun(step_inputs);
    }

    /*** Quality control ***/
    {
        state[step_qc] = qc.unserialize(handle, state[step_inputs]);
        await quickFun(step_qc);
    }

    {
        state[step_qc_adt] = qcadt.unserialize(handle, state[step_inputs]);
        await quickFun(step_qc_adt);
    }

    {
        state[step_qc_crispr] = qccrispr.unserialize(handle, state[step_inputs]);
        await quickFun(step_qc_crispr);
    }

    {
        let qc_states = { RNA: state[step_qc], ADT: state[step_qc_adt], CRISPR: state[step_qc_crispr] };
        state[step_filter] = filters.unserialize(handle, state[step_inputs], qc_states);
        await quickFun(step_filter);
    }

    /*** Normalization ***/
    {
        state[step_norm] = normalization.unserialize(handle, state[step_qc], state[step_filter]);
        await quickFun(step_norm);
    }

    {
        state[step_norm_adt] = normadt.unserialize(handle, state[step_qc_adt], state[step_filter]);
        await quickFun(step_norm_adt);
    }

    {
        state[step_norm_crispr] = normcrispr.unserialize(handle, state[step_qc_crispr], state[step_filter]);
        await quickFun(step_norm_crispr);
    }

    /*** Feature selection ***/
    {
        state[step_feat] = variance.unserialize(handle, permuters["RNA"], state[step_filter], state[step_norm]);
        await quickFun(step_feat);
    }

    /*** Dimensionality reduction ***/
    {
        state[step_pca] = pca.unserialize(handle, state[step_filter], state[step_norm], state[step_feat]);
        await quickFun(step_pca);
    }

    {
        state[step_pca_adt] = pcaadt.unserialize(handle, state[step_filter], state[step_norm_adt]);
        await quickFun(step_pca_adt);
    }

    {
        state[step_pca_crispr] = pcacrispr.unserialize(handle, state[step_filter], state[step_norm_crispr]);
        await quickFun(step_pca_crispr);
    }

    {
        let pca_states = { RNA: state[step_pca], ADT: state[step_pca_adt], CRISPR: state[step_pca_crispr] };
        state[step_combine] = combine.unserialize(handle, pca_states);
        await quickFun(step_combine);
    }

    {
        state[step_correct] = correct.unserialize(handle, state[step_filter], state[step_combine]);
        await quickFun(step_correct);
    }

    /*** Nearest neighbors ***/
    {
        state[step_neighbors] = index.unserialize(handle, state[step_correct]);
        await quickFun(step_neighbors);
    }

    /*** Visualization ***/
    {
        state[step_tsne] = await tsne.unserialize(handle, state[step_neighbors]);
        await quickFun(step_tsne);
    }

    {
        state[step_umap] = await umap.unserialize(handle, state[step_neighbors]);
        await quickFun(step_umap);
    }

    /*** Clustering ***/
    {
        state[step_kmeans] = kmeans_cluster.unserialize(handle, state[step_correct]);
        await quickFun(step_kmeans);
    }

    {
        state[step_snn] = snn_cluster.unserialize(handle, state[step_neighbors]);
        await quickFun(step_snn);
    }

    {
        state[step_choice] = cluster_choice.unserialize(handle, state[step_snn], state[step_kmeans]);
        await quickFun(step_choice);
    }

    /*** Markers and labels ***/
    let norm_states = { "RNA": state[step_norm], "ADT": state[step_norm_adt], "CRISPR": state[step_norm_crispr] };
    {
        state[step_markers] = cluster_markers.unserialize(handle, permuters, state[step_filter], norm_states, state[step_choice]);
        await quickFun(step_markers);
    }

    {
        state[step_labels] = label_cells.unserialize(handle, state[step_inputs], state[step_markers]);
        await quickFun(step_labels);
    }

    {
        state[step_custom] = custom_markers.unserialize(handle, permuters, state[step_filter], norm_states);
        await quickFun(step_custom);
    }

    /*** Faeture set enrichment ***/
    {
        state[step_enrichment] = feature_set_enrichment.unserialize(handle, state[step_inputs], state[step_filter], state[step_norm], state[step_markers]);
    }

    // Adding a tripwire for runAnalysis.
    state[load_flag] = true;

    return state;
}
