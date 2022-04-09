import * as workers from "worker_threads";

export function createTsneWorker() {
    return new workers.Worker(new URL("../tsne.worker.js", import.meta.url));
}

export function createUmapWorker() {
    return new workers.Worker(new URL("../umap.worker.js", import.meta.url));
}

export function registerCallback(worker, callback) {
    // Wrapping it in an extra 'data' to mimic web workers.
    worker.on("message", msg => callback({ data: msg }));
    return;
}

export function sendMessage(worker, message, transfer) {
    worker.postMessage(message, transfer);
    return;
}

export function terminateWorker(worker) {
    worker.terminate();
    return new Promise(resolve => resolve(true));
}
