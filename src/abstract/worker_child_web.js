export function registerCallback(callback) {
    self.onmessage = callback;
    return;
}

export function sendMessage(message, transfer) {
    self.postMessage(message, transfer);
    return;
}
