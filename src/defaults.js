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
    if (method == "mnn") {
        parameters[correct.step_name].method = "mnn";
        parameters[pca.step_name].block_method = "weight";
        parameters[pcaadt.step_name].block_method = "weight";
    } else if (method == "regress") {
        parameters[correct.step_name].method = "none";
        parameters[pca.step_name].block_method = "regress";
        parameters[pcaadt.step_name].block_method = "regress";
    } else if (method == "none") {
        parameters[correct.step_name].method = "none";
        parameters[pca.step_name].block_method = "none";
        parameters[pcaadt.step_name].block_method = "none";
    } else {
        throw new Error("unknown method '" + method + "'");
    }
    return parameters;
}

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
    for (const step of [correct.step_name, combine.step_name, index.step_name]) {
        parameters[step].approximate = approximate;
    }
    return parameters;
}
