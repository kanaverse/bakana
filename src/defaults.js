import * as inputs from "./inputs.js";
import * as qc from "./quality_control.js";
import * as qcadt from "./adt/quality_control.js";
import * as filter from "./cell_filtering.js";
import * as norm from "./normalization.js";
import * as normadt from "./adt/normalization.js";
import * as pca from "./pca.js";
import * as pcaadt from "./adt/pca.js";
import * as combine from "./combine_embeddings.js";
import * as correct from "./batch_correction.js";
import * as index from "./neighbor_index.js";
import * as markers from "./marker_detection.js";
import * as custom from "./custom_selections.js";

/**
 * Generate an object containing all of the default analysis parameters.
 *
 * @return An object where each property corresponds to an analysis step and contains the default parameters for that step.
 * See the documentation for each step's `compute` method for more details:
 * 
 * - {@linkcode InputsState#compute inputs}
 * - {@linkcode QualityControlState#compute quality_control}
 * - {@linkcode AdtQualityControlState#compute adt_quality_control}
 * - {@linkcode CellFiltering#compute cell_filtering}
 * - {@linkcode FeatureSelectionState#compute feature_selection}
 * - {@linkcode PcaState#compute pca}
 * - {@linkcode NeighborIndexState#compute neighbor_index}
 * - {@linkcode TsneState#compute tsne}
 * - {@linkcode UmapState#compute umap}
 * - {@linkcode KmeansClusterState#compute kmeans_cluster}
 * - {@linkcode SnnGraphClusterState#compute snn_graph_cluster}
 * - {@linkcode ChooseClusteringState#compute choose_clustering}
 * - {@linkcode CellLabellingState#compute cell_labelling}
 */
export function analysisDefaults() {
    var output = {
        feature_selection: {
            span: 0.3
        },
        combine_embeddings: {
            weights: null
        },
        batch_correction: {
            method: "none"
        },
        tsne: {
            perplexity: 30,
            iterations: 500,
            animate: false
        },
        umap: {
            num_neighbors: 15,
            num_epochs: 500,
            min_dist: 0.1,
            animate: false
        },
        kmeans_cluster: {
            k: 10
        },
        snn_graph_cluster: {
            k: 10,
            scheme: "rank",
            resolution: 1
        },
        choose_clustering: {
            method: "snn_graph"
        },
        cell_labelling: {
            mouse_references: [],
            human_references: []
        }
    };

    output[inputs.step_name] = inputs.InputsState.defaults();

    output[qc.step_name] = qc.QualityControlState.defaults();
    output[qcadt.step_name] = qcadt.AdtQualityControlState.defaults();
    output[filter.step_name] = filter.CellFilteringState.defaults();

    output[norm.step_name] = norm.NormalizationState.defaults();
    output[normadt.step_name] = normadt.AdtNormalizationState.defaults();

    output[pca.step_name] = pca.PcaState.defaults();
    output[pcaadt.step_name] = pcaadt.AdtPcaState.defaults();

    output[combine.step_name] = combine.CombineEmbeddingsState.defaults();
    output[correct.step_name] = correct.BatchCorrectionState.defaults();

    output[index.step_name] = index.NeighborIndexState.defaults();

    output[markers.step_name] = markers.MarkerDetectionState.defaults();
    output[custom.step_name] = custom.CustomSelectionsState.defaults();

    return output;
}

const correctible_pca_steps = [pca.step_name, pcaadt.step_name];

/**
 * Set the batch correction parameters across multiple steps.
 * This is a convenient helper as the correction process is split across the PCA and batch correction steps.
 * For MNN, we need to weight by block in PCA before performing MNN correction;
 * for linear regression, we need to regress by block in PCA without any additional correction;
 * and for no correction, we need to turn off any block handling in PCA as well as removing any additional correction.
 *
 * @param {object} parameters Object containing parameters for all steps, e.g., from {@linkcode analysisDefaults}.
 * @param {string} method Correction method to perform, one of `"mnn"`, "`regress"` or `"none"`.
 * 
 * @return `parameters` is modified with appropriate parameters in `batch_correction`, `pca` and `pcaadt`.
 */
export function configureBatchCorrection(parameters, method) {
    let correct_method;
    let pca_blocker;

    if (method == "mnn") {
        correct_method = method;
        pca_blocker = "weight";
    } else if (method == "regress") {
        correct_method = "none";
        pca_blocker = method;
    } else if (method == "none") {
        correct_method = method;
        pca_blocker = method;
    } else {
        throw new Error("unknown correction method '" + method + "'");
    }

    parameters[correct.step_name].method = correct_method;
    for (const x of correctible_pca_steps) {
        parameters[x].block_method = pca_blocker;
    }

    return parameters;
}

/**
 * Guess the `method` value from {@linkcode configureBatchCorrection} based on the parameter object.
 * This effectively consolidates the various correction parameters into a single setting.
 *
 * @param {object} parameters - Object containing parameters for all steps, typically after {@linkcode configureBatchCorrection}.
 * @param {object} [options] - Optional parameters.
 * @param {boolean} [options.strict] - Whether to only report the `method` when the set of parameters modified by {@linkcode configureBatchCorrection} are consistent.
 *
 * @return {?string} One of `"mnn"`, `"regress"` or `"none"`, based on the expected set of modifications from {@linkcode configureBatchCorrection}.
 * If `strict = false` and there is no exact match to the expected set, the most appropriate method is returned;
 * otherwise, if `strict = true`, `null` is returned.
 */
export function guessBatchCorrectionConfig(parameters, { strict = false } = {}) {
    let pca_blockers = new Set(correctible_pca_steps.map(x => parameters[x].block_method));

    let resp;
    if (parameters[correct.step_name].method == "mnn") {
        resp = "mnn";
        if (strict) {
            if (pca_blockers.size > 1 || !pca_blockers.has("weight")) {
                resp = null;
            }
        }
    } else {
        if (pca_blockers.has("regress")) {
            if (strict && pca_blockers.size > 1) {
                resp = null;
            } else {
                resp = "regress";
            }
        } else if (pca_blockers.has("none")) {
            if (strict && pca_blockers.size > 1) {
                resp = null;
            } else {
                resp = "none";
            }
        } else {
            // If pca block_methods are set to 'weight',
            // this doesn't really correspond to anything,
            // but is closest to 'none'.
            if (strict) {
                resp = null;
            } else {
                resp = "none";
            }
        }
    }

    return resp;
}

const approximatable_steps = [correct.step_name, combine.step_name, index.step_name];

/**
 * Specify whether approximate neighbor searches should be performed across all affected steps.
 * This is a convenient helper as it is generally unnecessary to switch between exact and approximate searches in different steps.
 *
 * @param {object} parameters Object containing parameters for all steps, e.g., from {@linkcode analysisDefaults}.
 * @param {boolean} approximate Whether to perform approximate nearest neighbor searces.
 * 
 * @return `parameters` is modified with appropriate parameters in relevant steps, e.g., `batch_correction`, `combine_embeddings` and `neighbor_index`.
 */
export function configureApproximateNeighbors(parameters, approximate) {
    for (const step of approximatable_steps) {
        parameters[step].approximate = approximate;
    }
    return parameters;
}

/**
 * Guess the value of `approximate` from {@linkcode configureApproximateNeighbors} based on the parameter object.
 * This effectively consolidates the various approximation parameters into a single setting.
 *
 * @param {object} parameters - Object containing parameters for all steps, typically after {@linkcode configureApproximateNeighbors}.
 * @param {object} [options] - Optional parameters.
 * @param {boolean} [options.strict] - Whether to only report `approximate` when the set of parameters modified by {@linkcode configureApproximateNeighbors} are consistent.
 *
 * @return {?boolean} Whether or not approximate neighbor search was used.
 * If `strict = false` and there is a mixture of approximate and exact matches, an approximate search is reported;
 * otherwise, if `strict = true`, `null` is returned.
 */
export function guessApproximateNeighborsConfig(parameters, { strict = false } = {}) {
    let approximates = new Set(approximatable_steps.map(x => parameters[x].approximate));
    if (strict && approximates.size > 1) {
        return null;
    } else {
        return approximates.has(true);
    }
}
