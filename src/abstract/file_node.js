/**
 * 'obj' is expected to be a string containing the file path.
 */

import * as fs from "fs";
import * as path from "path";

export function rawSize(obj) {
    return fs.statSync(obj).size;
}

export function rawName(obj) {
    return path.basename(obj);
}

export class LoadedFile {
    constructor(obj) {
        this.path = obj;
    }

    buffer() {
        let f = fs.readFileSync(this.path);
        return f.buffer.slice(f.byteOffset, f.byteOffset + f.byteLength);
    }

    size() {
        return fs.statSync(this.path).size;
    }

    serialized() {
        return this.path;
    }
};

export function realizeMatrixMarket(loaded) {
    return loaded.path;
}

export function realizeH5(loaded) {
    return loaded.path;
}

// The default is allowed = false here; we are not allowed to delete stuff
// from the host filesystem, we just read from it most of the time.
export function removeH5(path, { allowed = false } = {}) {
    if (allowed && fs.existsSync(path)) {
        fs.unlinkSync(path);
    }
    return;
}
