import * as path from "path";
import * as fs from "fs";

export var baseParams = {
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
        num_pcs: 10,
        block_method: "none"
    },
    neighbor_index: {
        approximate: true
    },
    tsne: {
        perplexity: 30,
        iterations: 10,
        animate: false
    },
    umap: {
        num_neighbors: 15,
        num_epochs: 10,
        min_dist: 0.1,
        animate: false
    },
    kmeans_cluster: {
        k: 10
    },
    snn_graph_cluster: {
        k: 10,
        scheme: 0,
        resolution: 1
    },
    choose_clustering: {
        method: "snn_graph"
    },
    cell_labelling: {
        mouse_references: [ "ImmGen" ],
        human_references: [ "BlueprintEncode" ]
    }
};

export function downloadReference(url) {
    let fpath = path.basename(decodeURIComponent(url));
    let obj = fs.readFileSync("files/references/" + fpath);
    return obj.buffer.slice(obj.byteOffset, obj.byteOffset + obj.byteLength);
}

export function mockOffsets(paths) {
    let offsets = {};
    let sofar = 0;
    for (const p of paths) {
        offsets[sofar] = p;
        sofar += fs.statSync(p).size;
    }
    return offsets;
}

