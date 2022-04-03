export function createWorker(uri) {
    return new Worker(uri, { type: "module" });
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
    return worker.terminate();
}
