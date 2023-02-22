import * as inputs from "./steps/inputs.js";
import * as qc from "./steps/rna_quality_control.js";
import * as qcadt from "./steps/adt_quality_control.js";
import * as qccrispr from "./steps/crispr_quality_control.js";
import * as filter from "./steps/cell_filtering.js";
import * as norm from "./steps/rna_normalization.js";
import * as normadt from "./steps/adt_normalization.js";
import * as normcrispr from "./steps/crispr_normalization.js";
import * as pca from "./steps/rna_pca.js";
import * as pcaadt from "./steps/adt_pca.js";
import * as pcacrispr from "./steps/crispr_pca.js";
import * as combine from "./steps/combine_embeddings.js";
import * as correct from "./steps/batch_correction.js";
import * as index from "./steps/neighbor_index.js";
import * as snngraph from "./steps/snn_graph_cluster.js";
import * as markers from "./steps/marker_detection.js";
import * as custom from "./steps/custom_selections.js";
import * as enrichment from "./steps/feature_set_enrichment.js";
import * as labelling from "./steps/cell_labelling.js";

/**
 * Generate an object containing all of the default analysis parameters.
 *
 * @return An object where each property corresponds to an analysis step and contains the default parameters for that step.
 * See the documentation for each step's `compute` method for more details:
 * 
 * - {@linkcode InputsState#compute inputs}
 * - {@linkcode RnaQualityControlState#compute rna_quality_control}
 * - {@linkcode AdtQualityControlState#compute adt_quality_control}
 * - {@linkcode CrisprQualityControlState#compute crispr_quality_control}
 * - {@linkcode CellFiltering#compute cell_filtering}
 * - {@linkcode RnaNormalizationState#compute rna_normalization}
 * - {@linkcode AdtNormalizationState#compute adt_normalization}
 * - {@linkcode CrisprNormalizationState#compute crispr_normalization}
 * - {@linkcode FeatureSelectionState#compute feature_selection}
 * - {@linkcode RnaPcaState#compute rna_pca}
 * - {@linkcode AdtPcaState#compute adt_pca}
 * - {@linkcode CrisprPcaState#compute crispr_pca}
 * - {@linkcode NeighborIndexState#compute neighbor_index}
 * - {@linkcode TsneState#compute tsne}
 * - {@linkcode UmapState#compute umap}
 * - {@linkcode KmeansClusterState#compute kmeans_cluster}
 * - {@linkcode SnnGraphClusterState#compute snn_graph_cluster}
 * - {@linkcode ChooseClusteringState#compute choose_clustering}
 * - {@linkcode CellLabellingState#compute cell_labelling}
 * - {@linkcode FeatureSetEnrichmentState#compute feature_set_enrichment}
 *
 * See also {@linkcode configureBatchCorrection} and {@linkcode configureApproximateNeighbors} to synchronize certain parameter settings across multiple steps.
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
        choose_clustering: {
            method: "snn_graph"
        }
    };

    output[inputs.step_name] = inputs.InputsState.defaults();

    output[qc.step_name] = qc.RnaQualityControlState.defaults();
    output[qcadt.step_name] = qcadt.AdtQualityControlState.defaults();
    output[qccrispr.step_name] = qccrispr.CrisprQualityControlState.defaults();
    output[filter.step_name] = filter.CellFilteringState.defaults();

    output[norm.step_name] = norm.RnaNormalizationState.defaults();
    output[normadt.step_name] = normadt.AdtNormalizationState.defaults();
    output[normcrispr.step_name] = normcrispr.CrisprNormalizationState.defaults();

    output[pca.step_name] = pca.RnaPcaState.defaults();
    output[pcaadt.step_name] = pcaadt.AdtPcaState.defaults();
    output[pcacrispr.step_name] = pcacrispr.CrisprPcaState.defaults();

    output[combine.step_name] = combine.CombineEmbeddingsState.defaults();
    output[correct.step_name] = correct.BatchCorrectionState.defaults();

    output[index.step_name] = index.NeighborIndexState.defaults();
    output[snngraph.step_name] = snngraph.SnnGraphClusterState.defaults();

    output[markers.step_name] = markers.MarkerDetectionState.defaults();
    output[custom.step_name] = custom.CustomSelectionsState.defaults();

    output[enrichment.step_name] = enrichment.FeatureSetEnrichmentState.defaults();
    output[labelling.step_name] = labelling.CellLabellingState.defaults();

    return output;
}

const correctible_pca_steps = [pca.step_name, pcaadt.step_name, pcacrispr.step_name];

/**
 * Set the batch correction parameters across multiple steps.
 * This is a convenient helper as the correction process is split across the PCA and batch correction steps.
 * For MNN, we need to weight by block in the {@linkplain PcaState} before performing MNN correction in {@linkplain BatchCorrectionState};
 * for linear regression, we need to regress by block in {@linkplain PcaState} without any additional correctioni in {@linkplain BatchCorrectionState};
 * and for no correction, we need to turn off any block handling in {@linkplain PcaState}as well as removing any additional correction in {@linkplain BatchCorrectionState}.
 *
 * @param {object} parameters Object containing parameters for all steps, e.g., from {@linkcode analysisDefaults}.
 * @param {string} method Correction method to perform, one of `"mnn"`, "`regress"` or `"none"`.
 * 
 * @return `parameters` is modified with appropriate parameters in `batch_correction`, `pca` and `adt_pca`.
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
 * Affected steps are {@linkplain BatchCorrectionState}, {@linkplain CombineEmbeddingsState} and {@linkplain NeighborIndexState}.
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
