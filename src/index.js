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

import * as utils from "./utils/general.js";
import * as serialize_utils from "./utils/serialize.js";
import * as rutils from "./utils/reader.js";

const step_inputs = "inputs";
const step_qc = "quality_control";
const step_norm = "normalizaton";
const step_feat = "feature_selecton";
const step_pca = "pca";
const step_neighbors = "neighbor_index";
const step_tsne = "tsne";
const step_umap = "umap";
const step_kmeans = "kmeans_cluster";
const step_snn = "snn_cluster_graph";
const step_choice = "choose_clustering";
const step_markers = "marker_detection";
const step_labels = "cell_labelling";
const step_custom = "custom_marker_management";

export function runAnalysis(files, params, finished, download) {
    inputs.compute(files, state.files.batch);
    if (inputs.changed) {
        finished(inputs.results(), step_inputs, "Count matrix loaded");
    }

    qc.compute(
        params[step_qc]["use_mito_default"], 
        params[step_qc]["mito_prefix"], 
        params[step_qc]["nmads"]
    );
    if (qc.changed) {
        finished(qc.results(), step_qc, "Applying quality control filters");
    }
 
    normalization.compute();
    if (normalization.changed) {
        finished(normalization.results(), step_norm, "Log-normalization completed");
    }

    variance.compute(
        params[step_feat]["span"]
    );
    if (variance.changed) {
        finished(variance.results(), step_feat, "Variance modelling completed");
    }

    pca.compute(
        params[step_pca]["num_hvgs"],
        params[step_pca]["num_pcs"],
        params[step_pca]["block_method"]
    );
    if (pca.changed) {
        finished(pca.results(), step_pca, "Principal components analysis completed");
    }

    index.compute(
        params[step_neighbors]["clus-approx"]
    );
    if (index.changed) {
        finished(index.results(), step_neighbors, "Neighbor search index constructed");
    }

    tsne.compute(
        params[step_tsne]["perplexity"],
        params[step_tsne]["iterations"], 
        params[step_tsne]["animate"]
    );
    if (tsne.changed) {
        tsne.results().then(res => finished(res, step_tsne, "t-SNE completed"));
    }

    umap.compute(
        params[step_umap]["num_neighbors"], 
        params[step_umap]["num_epochs"], 
        params[step_umap]["min_dist"], 
        params[step_umap]["animate"]
    );
    if (umap.changed) {
        umap.results().then(res => finished(umap, step_umap, "UMAP completed"));
    }

    let method = params[step_choice]["method"];

    kmeans_cluster.compute(
        method == "kmeans", 
        params[step_kmeans]["k"]
    );
    if (kmeans_cluster.changed) {
        finished(kmeans_cluster.results(), step_kmeans, "K-means clustering completed");
    }

    snn_cluster.compute(
        method == "snn_graph", 
        params[step_snn]["k"], 
        params[step_snn]["scheme"], 
        params[step_snn]["resolution"]
    );
    if (snn_cluster.changed) {
        finished(snn_cluster.results(), step_snn, "SNN graph clustering completed");
    }
  
    cluster_choice.compute(
        method
    );
    if (cluster_choice.changed) {
        finished(cluster_choice.results(), step_choice, "Clustering of interest chosen");
    }

    cluster_markers.compute();
    if (cluster_markers.changed) {
        finished(cluster_markers.results(), step_markers, "Marker detection complete");
    }

    label_cells.compute(
        params[step_labels]["human_references"],
        params[step_labels]["mouse_references"],
        download
    );
    if (label_cells.changed) {
        label_cells.results().then(res => finished(res, step_labels, "Cell type labelling complete"));
    }

    custom_markers.compute();
    if (custom_markers.changed) {
        finished(custom_markers.results(), step_custom, "Pruning of custom markers finished");
    }

    return;
}
 
export async function saveAnalysis(linker) {
    const path = generateRandomName("state", ".h5");

    let saver;
    let saved;
    let embedded = (typeof linker == "undefined");

    if (embedded) {
        saved = sutils.embedFiles();
        saver = saved.saver;
        delete saved.saver;
    } else {
        saved = { collected: [] };
        saver = (obj) => {
            let id = linker(obj);
            saved.collected.push(id);
            return id;
        }
    }

    let output;
    try {
        let handle = scran.createNewHDF5File(path);

        await inputs.serialize(handle, saver, embedded);

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

        output = scran.readFile(path);
    } finally {
        if (scran.fileExists(path)) {
            scran.removeFile(path);
        }
    }

    return {
        "state": output,
        "saved": saved
    };
}

export async function loadAnalysis(path, finished, loader, embedded) {
    let response = {};
    let handle = new scran.H5File(path);

    let permuter = await inputs.unserialize(handle, loader, embedded);
    finished(inputs.results(), step_inputs, "Reloaded count matrix");

    response[step_qc] = qc.unserialize(handle);
    finished(qc.results(), step_qc, "Reloaded QC metrics");

    response[step_norm] = normalization.unserialize(handle);
    finished(normalization.results(), step_norm, "Reloaded log-normalization");

    response[step_feat] = variance.unserialize(handle, permuter);
    finished(variance.results(), step_feat, "Reloaded variance modelling statistics");

    response[step_pca] = pca.unserialize(handle);
    finished(pca.results(), step_pca, "Reloaded principal components");

    response[step_neighbors] = index.unserialize(handle);
    finished(index.results(), step_neighbors, "Reloaded neighbor search index");

    response[step_tsne] = tsne.unserialize(handle);
    tsne.results().then(res => finished(res, step_tsne, "t-SNE reloaded"));

    response[step_umap] = umap.unserialize(handle);
    umap.results().then(res => finished(res, step_umap, "UMAP reloaded"));

    response[step_kmeans] = kmeans_cluster.unserialize(handle);
    finished(kmeans_cluster.results(), step_kmeans, "K-means clustering reloaded");

    response[step_snn] = snn_cluster.unserialize(handle);
    finished(snn_cluster.results(), step_snn, "SNN graph clustering reloaded");

    response[step_choice] = cluster_choice.unserialize(handle);
    finished(cluster_choice.results(), step_choice, "Clustering of interest chosen");

    response[step_markers] = cluster_markers.unserialize(handle, permuter);
    finished(cluster_markers.results(), step_markers, "Reloaded per-cluster markers");

    response[step_labels] = label_cells.unserialize(handle);
    labels.results().then(res => finished(res, step_labels, "Reloaded cell type labels"));

    response[step_custom] = custom_markers.unserialize(handle, permuter);
    finished(custom_markers.results(), step_custom, "Pruning of custom markers finished");

    return response;
}
