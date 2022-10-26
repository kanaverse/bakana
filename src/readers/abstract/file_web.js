export class SimpleFile {
    #buffer;
    #name;

    constructor(x, { name = null } = {}) {
        if (x instanceof File) {
            let reader = new FileReaderSync();
            let b = reader.readAsArrayBuffer(x);
            this.#buffer = new Uint8Array(b);
            this.#name = x.name;
        } else if (x instanceof Uint8Array) {
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
        if (copy) {
            return this.#buffer.slice();
        } else {
            return this.#buffer;
        }
    }

    size() {
        return this.#buffer.length;
    }

    name() {
        return this.#name;
    }

    contents({ copy = false } = {}) {
        return this.buffer({ copy: copy });
    }
}
