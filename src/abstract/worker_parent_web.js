/**
 * Webpack (or something in kana's build system) specifically recognizes the
 * hard-coded path in this 'new Worker(new URL(...))' pattern. This is why we
 * have hard-coded creators for the workers rather than allowing callers to
 * pass in the URL as a variable, as that doesn't pack the worker's JS.
 */

export function createTsneWorker() {
    return new Worker(new URL("../tsne.worker.js", import.meta.url), { type: "module" });
}

export function createUmapWorker() {
    return new Worker(new URL("../umap.worker.js", import.meta.url), { type: "module" });
}

export function registerCallback(worker, callback) {
    worker.onmessage = callback;
    return;
}

export function sendMessage(worker, message, transfer) {
    worker.postMessage(message, transfer);
    return;
}

export function terminateWorker(worker) {
    worker.terminate();
    return;
}
