import * as scran from "scran.js";
import * as jsp from "jaspagate";
import * as internal from "./abstract/dump.js";

function toStringType(data, options) {
    let length = null;
    if ("maxStringLength" in options) {
        length = options.maxStringLength;
    }
    if (length == null) {
        if (data == null) {
            throw new Error("'maxStringLength' must be supplied if 'data' is null"); 
        }
        length = scran.findMaxStringLength(data, null);
    }
    return new scran.H5StringType("UTF-8", length);
}

export class AlabasterH5Group extends jsp.H5Group {
    #handle;

    constructor(handle) {
        super();
        this.#handle = handle;
    }

    writeAttribute(name, type, shape, data, options = {}) {
        if (type == "String") {
            type = toStringType(data, options);
        }
        this.#handle.writeAttribute(name, type, shape, data);
    }

    open(name) {
        let out = this.#handle.open(name);
        if (out instanceof scran.H5Group) {
            return new AlabasterH5Group(out);
        } else {
            return new AlabasterH5DataSet(out);
        }
    }

    children() {
        return this.#handle.children;
    }

    createGroup(name) {
        return new AlabasterH5Group(this.#handle.createGroup(name));
    }

    createDataSet(name, type, shape, options = {}) {
        if ("data" in options) {
            if (type == "String") {
                type = toStringType(options.data, options);
            }
            return new AlabasterH5DataSet(this.#handle.writeDataSet(name, type, shape, options.data));
        } else {
            if (type == "String") {
                type = toStringType(null, options);
            }
            return new AlabasterH5DataSet(this.#handle.createDataSet(name, type, shape));
        }
    }

    close() {}
}

export class AlabasterH5DataSet extends jsp.H5DataSet {
    #handle;

    constructor(handle) {
        super();
        this.#handle = handle;
    }

    writeAttribute(name, type, shape, data, options = {}) {
        if (type == "String") {
            type = toStringType(data, options);
        }
        this.#handle.writeAttribute(name, type, shape, data);
    }

    write(x) {
        this.#handle.write(x);
    }

    close() {}
}

export class AlabasterGlobalsInterface extends jsp.GlobalsInterface {
    #directory;
    #files;

    constructor(directory, files) {
        super();
        this.#directory = directory;
        this.#files = files;
    }

    get(path, options = {}) {
        if (this.#directory !== null) {
            const { asBuffer = true } = options;
            return internal.read(this.#directory, path, asBuffer);
        } else {
            return this.#files[path];
        }
    }

    clean(localPath) {}

    async write(path, contents) {
        if (this.#directory !== null) {
            await internal.write(this.#directory, path, contents);
        } else {
            this.#files[path] = contents;
        }
    }

    async mkdir(path) {
        if (this.#directory !== null) {
            await internal.mkdir(this.#directory, path);
        }
    }

    async copy(from, to) {
        if (this.#directory !== null) {
            await internal.copy(this.#directory, from, to);
        } else {
            this.#files[to] = this.#files[from];
        }
    }

    h5create(path) {
        if (this.#directory !== null) {
            let actual_path = jsp.joinPath(this.#directory, path);
            let latest = scran.createNewHdf5File(actual_path);
            let output = new AlabasterH5Group(latest);
            output._path = actual_path;
            return output;
        } else {
            let tmppath = scran.chooseTemporaryPath();
            let latest = scran.createNewHdf5File(tmppath);
            let output = new AlabasterH5Group(latest);
            output._path = tmppath;
            output._intended = path;
            return output;
        }
    }

    h5finish(handle, failed) {
        if (this.#directory === null) {
            if (!failed) {
                this.write(handle._intended, scran.readFile(handle._path));
            }
            scran.removeFile(handle._path);
        }
    }
}
