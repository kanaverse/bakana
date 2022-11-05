export function stream(file, chunkSize, callback, resolve, reject) {
    reject("no support for file paths in the browser context");
}

export function peek(file, n) {
    throw new Error("no support for file paths in the browser context");
}
