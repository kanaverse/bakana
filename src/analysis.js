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
const step_feat = variance.step_name;
const step_pca = pca.step_name;
const step_pca_adt = pcaadt.step_name;
const step_pca_crispr = pcacrispr.step_name;
const step_combine = combine.step_name;
const step_correct = correct.step_name;
const step_neighbors = index.step_name;
const step_tsne = tsne.step_name;
const step_umap = umap.step_name;
const step_kmeans = kmeans_cluster.step_name;
const step_snn = snn_cluster.step_name;
const step_choice = cluster_choice.step_name;
const step_markers = cluster_markers.step_name;
const step_labels = label_cells.step_name;
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
    output[step_labels] = new label_cells.CellLabellingState(output[step_inputs]);
    output[step_enrichment] = new feature_set_enrichment.FeatureSetEnrichmentState(output[step_inputs], output[step_filter], output[step_norm]);
    output[step_custom] = new custom_markers.CustomSelectionsState(output[step_filter], norm_states);

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
 * @param {?function} [options.startFun=null] - Function that is called when each step is started.
 * This should accept a single argument - the name of the step.
 * The return value is ignored, but any promises will be awaited before the analysis proceeds to the next step.
 * If `null`, nothing is executed.
 * @param {?function} [options.finishFun=null] - Function that is called on successful execution of each step.
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
    await state[step_inputs].compute(datasets, params[step_inputs]);
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

    /*** Preprocessing steps ***/
    let basic_steps = [
        step_qc, step_qc_adt, step_qc_crispr,
        step_filter,
        step_norm, step_norm_adt, step_norm_crispr,
        step_feat,
        step_pca, step_pca_adt, step_pca_crispr,
        step_combine,
        step_correct,
        step_neighbors
    ];

    for (const step of basic_steps) {
        await quickStart(step);
        await state[step].compute(params[step]);
        await quickFinish(step);
    }

    /*** Visualization ***/
    for (const step of [ step_tsne, step_umap ]) {
        await quickStart(step);
        let p = state[step].compute(params[step]);
        deferredQuickFinish(step, p);
    }

    /*** Clustering ***/
    let method = params[step_choice]["method"];

    await quickStart(step_kmeans);
    state[step_kmeans].compute(method == "kmeans", params[step_kmeans]);
    await quickFinish(step_kmeans);

    await quickStart(step_snn);
    state[step_snn].compute(method == "snn_graph", params[step_snn]);
    await quickFinish(step_snn);

    await quickStart(step_choice);
    state[step_choice].compute(params[step_choice]);
    await quickFinish(step_choice);

    /*** Markers and labels ***/
    let remaining = [
        step_markers,
        step_labels,
        step_custom,
        step_enrichment
    ];

    for (const step of remaining) {
        await quickStart(step);
        await state[step].compute(params[step]);
        await quickFinish(step);
    }

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
