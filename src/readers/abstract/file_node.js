import * as fs from "fs";
import * as path from "path";

/**
 * Simple wrapper to represent a file in both Node.js and the browser.
 * This enables downstream code to operate independently of the source of the file.
 */
export class SimpleFile {
    #mode;
    #path;
    #buffer;
    #name;

    /**
     * @param {Uint8Array|string|File} x - Contents of a file or a pointer to it.
     * For browsers, this may be a File object;
     * for Node.js, this may be a string containing a file path.
     * @param {object} [options={}] - Optional parameters.
     * @param {?string} [options.name=null] - Name of the file.
     * If `null`, it is determined automatically from `x` where possible (i.e., basename of the file path or extracted from the File object).
     * If `x` is a Uint8Array, an explicit name is required.
     */
    constructor(x, { name = null } = {}) {
        if (typeof x == "string") {
            this.#mode = "path";
            this.#path = x;
            if (name === null) {
                name = path.basename(x);
            }
            this.#name = name;
        } else {
            this.#mode = "buffer";
            this.#buffer = x;
            if (name === null) {
                throw new Error("'name' must be provided for Uint8Array inputs in SimpleFile constructor");
            }
            this.#name = name;
        }
    }

    #get_buffer(copy) {
        if (copy) {
            return this.#buffer.slice();
        } else {
            return this.#buffer;
        }
    }

    /**
     * @param {object} [options={}] - Optional parameters.
     * @param {boolean} [options.copy=false] - Whether to create a new copy of the buffer.
     * Callers should set this to `true` for non-read-only applications.
     * @return {Uint8Array} Contents of the file as an array.
     */
    buffer({ copy = false } = {}) {
        if (this.#mode == "path") {
            let f = fs.readFileSync(this.#path);
            let b = f.buffer.slice(f.byteOffset, f.byteOffset + f.byteLength);
            return new Uint8Array(b);
        } else {
            return this.#get_buffer(copy);
        }
    }

    /**
     * @return {string} Name of the file, for disambiguation purposes.
     */
    name() {
        return this.#name;
    }

    /**
     * @return {number} Size of the file, in bytes.
     */
    size() {
        if (this.#mode == "path") {
            return fs.statSync(this.#path).size;
        } else {
            return this.#buffer.length;
        }
    }

    /**
     * @param {object} [options={}] - Optional parameters.
     * @param {boolean} [options.copy=false] - Whether to create a new copy of the buffer.
     * Callers should set this to `true` for non-read-only applications.
     * @return {Uint8Array|string} Contents of the file, as in {@linkcode SimpleFile#buffer SimpleFile.buffer}.
     * For Node.js, a string may also be returned containing the path to the original file.
     */
    content({ copy = false } = {}) {
        if (this.#mode == "path") {
            return this.#path;
        } else {
            return this.#get_buffer(copy);
        }
    }
};
