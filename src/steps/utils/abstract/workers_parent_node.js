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
    return worker.terminate();
}
