/**
 * Generate an object containing all of the default analysis parameters.
 *
 * @return An object where each property corresponds to an analysis step and contains the default parameters for that step.
 */
export function analysisDefaults() {
    return {
        inputs: {
            sample_factor: null
        },
        quality_control: {
            use_mito_default: true,
            mito_prefix: "mt-",
            nmads: 3
        },
        feature_selection: {
            span: 0.3
        },
        pca: {
            num_hvgs: 2000,
            num_pcs: 20,
            block_method: "none"
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
}
