export class SimpleFile {
    #mode;
    #buffer;
    #file;
    #name;

    constructor(x, { name = null } = {}) {
        if (x instanceof File) {
            this.#mode = "file";
            this.#file = x;
            if (name === null) {
                name = x.name;
            }
            this.#name = name;
        } else if (x instanceof Uint8Array) {
            this.#mode = "buffer";
            this.#buffer = x; 
            if (name === null) {
                throw new Error("'name' must be provided for Uint8Array inputs in SimpleFile constructor");
            }
            this.#name = name;
        } else {
            throw new Error("unknown type '" + typeof(x) + "' for SimpleFile constructor");
        }
    }

    buffer({ copy = false } = {}) {
        if (this.#mode == "file") {
            let reader = new FileReaderSync();
            let b = await reader.readAsArrayBuffer(x);
            return new Uint8Array(b);
        } else {
            if (copy) {
                return this.#buffer.slice();
            } else {
                return this.#buffer;
            }
        }
    }

    size() {
        return this.#buffer.length;
    }

    name() {
        return this.#name;
    }

    content({ copy = false } = {}) {
        return this.buffer({ copy: copy });
    }
}
