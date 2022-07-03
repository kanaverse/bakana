import * as TENxReader from "./10x.js";
import * as H5ADReader from "./h5ad.js";
import * as MtxReader from "./mtx.js";

export function chooseReader(format) {
    if (!(format in availableReaders)) {
        throw "unknown matrix format '" + format + "'";
    }
    return availableReaders[format];
}

/**
 * List of available readers.
 * Each key specifies a matrix format; the corresponding value should be an ES6 module containing methods to interpret that format.
 * See {@tutorial custom_readers} for details on adding new readers.
 */
export var availableReaders = {
    "MatrixMarket": MtxReader,
    "10X": TENxReader,
    "H5AD": H5ADReader
};
