import * as scran from "scran";
import * as jsp from "jaspagate";

/****************************
 *** Jaspagate interfaces ***
 ****************************/

function toStringType(data, options) {
    let length = null;
    if ("maxStringLength" in options) {
        length = options.maxStringLength;
    }
    if (length == null) {
        if (data == null) {
            throw new Error("'maxStringLength' must be supplied if 'data' is null"); 
        }
        length = scran.findMaxStringLength(data);
    }
    return new scran.H5StringType("UTF-8", length);
}

class AlabasterH5Group extends jsp.H5Group {
    #handle;

    constructor(handle) {
        super();
        this.#handle = handle;
    }

    writeAttribute(name, type, shape, data, options) {
        if (type == "String") {
            type = toStringType(data, options);
        }
        this.#handle.writeAttribute(name, type, shape, data);
    }

    createGroup(name) {
        return new AlabasterH5Group(this.#handle.createGroup(name);
    }

    createDataSet(name, type, shape, options) {
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

class AlabasterH5DataSet extends jsp.H5DataSet {
    #handle;

    constructor(handle) {
        super();
        this.#handle = handle;
    }

    writeAttribute(name, type, shape, data, options) {
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

class AlabasterGlobalsInterface extends jsp.GlobalsInterface {
    #directory;
    #ziphandle;

    constructor(directory, ziphandle) {
        super();
        this.#directory = internal.initialize(directory);
        this.#ziphandle = ziphandle;
    }

    get(path, options = {}) {
        if (this.#ziphandle !== null) {
            return this.#ziphandle.file(path).async("uint8array");
        } else {
            const { asBuffer = true } = options;
            return internal.read(this.#directory, path, asBuffer);
        }
    }

    clean(localPath) {}

    async write(path, contents) {
        if (this.#ziphandle !== null) {
            handle.file(f, contents);
        } else {
            await internal.write(this.#directory, path, contents);
        }
    }

    mkdir(path) {
        if (this.#ziphandle !== null) {
            handle.folder(path);
        } else {
            await internal.mkdir(this.#directory, path);
        }
    }

    copy(from, to) {
        if (this.#ziphandle !== null) {
            // More reliable to just clone the payload, as Linux symlinks won't survive on Windows.
            this.write(to, await this.get(from, { asBuffer: true }));
        } else {
            await internal.copy(this.#directory, from, to);
        }
    }

    h5create(path) {
        if (this.#ziphandle !== null || typeof this.#directory != "string") {
            let tmppath = scran.chooseTemporaryPath();
            let latest = scran.createNewHdf5File(tmppath);
            let output = new AlabasterH5Group(latest);
            output._path = tmppath;
            output._intended = path;
            return output;
        } else {
            let actual_path = internal.join(this.#directory, path);
            let latest = scran.createNewHdf5File(actual_path);
            let output = new AlabasterH5Group(latest);
            output._path = actual_path;
            return output;
        }
    }

    h5finish(handle, failed) {
        if (this.#ziphandle !== null || typeof this.#directory != "string") {
            if (!failed) {
                this.write(handle._intended, scran.readFile(handle._path));
            }
            scran.removeFile(handle._path);
        }
    }
}
