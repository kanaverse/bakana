import * as workers from "worker_threads";

export function registerCallback(callback) {
    // Wrapping it in an extra 'data:' to mimic Web Workers.
    workers.parentPort.on("message", msg => callback({ data: msg }));
    return;
}

export function sendMessage(msg) {
    workers.parentPort.postMessage(msg);
    return;
}
