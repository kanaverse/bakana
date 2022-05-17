import * as inputs from "./inputs.js";
import * as qc from "./quality_control.js";
import * as qcadt from "./adt/quality_control.js";
import * as filter from "./cell_filtering.js";
import * as norm from "./normalization.js";
import * as normadt from "./adt/normalization.js";
import * as pca from "./pca.js";
import * as pcaadt from "./adt/pca.js";

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
        adt_pca: {
            num_pcs: 20,
            block_method: "none"
        },
        combine_embeddings: {
            weights: null
        },
        batch_correction: {
            method: "none"
        },
        neighbor_index: {
            approximate: true
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

    return output;
}
