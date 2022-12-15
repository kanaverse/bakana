import * as workers from "worker_threads";

export function createTsneWorker() {
    return new workers.Worker(new URL("../tsne.worker.js", import.meta.url));
}

export function createUmapWorker() {
    return new workers.Worker(new URL("../umap.worker.js", import.meta.url));
}
