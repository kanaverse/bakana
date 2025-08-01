// Stub file to ensure that something is loaded.

export function fsexists() {
    return false;
}

export function read(dir, path, asBuffer) {
    throw new Error("read() is not supported in a web context");
}

export function write(dir, path, x) {
    throw new Error("write() is not supported in a web context");
}

export function mkdir(dir, path) {
    throw new Error("mkdir() is not supported in a web context");
}

export function copy(dir, from, to) {
    throw new Error("copy() is not supported in a web context");
}
