import * as scran from "scran.js";
import * as utils from "./../utils/general.js";
import * as rutils from "./../utils/reader.js";
import * as afile from "./../abstract/file.js";

export function format(args, signatureOnly) {
    return { 
        "format": "10X", 
        "h5": rutils.formatFile(args.h5, signatureOnly)
    };
}

function extract_features(handle) {
    let genes = null;

    if (!("matrix" in handle.children) || handle.children["matrix"] != "Group") {
        throw new Error("expected a 'matrix' group at the top level of the file");
    }

    let mhandle = handle.open("matrix");
    if ("features" in mhandle.children && mhandle.children["features"] == "Group") {
        let fhandle = mhandle.open("features");

        let ids = rutils.extractHDF5Strings(fhandle, "id");
        if (ids !== null) {
            genes = { id: ids };
            let names = rutils.extractHDF5Strings(fhandle, "name");
            if (names !== null) {
                genes.name = names;
            }
        }
    }

    return genes;
}

export function preflight(args) {
    let output = {};
    let formatted = format(args, false)

    const tmppath = afile.realizeH5(formatted.h5.content);
    try {
        let handle = new scran.H5File(tmppath);
        output.genes = extract_features(handle);
        output.annotations = null;
        // TODO: try pull out sample IDs from the 10X file, if they exist?
    } finally {
        afile.removeH5(tmppath);
    }

    return output;
}

export function load(formatted) {
    let output = {};

    const tmppath = afile.realizeH5(formatted.h5.content);
    try {
        output.matrix = scran.initializeSparseMatrixFromHDF5(tmppath, "matrix");
        let handle = new scran.H5File(tmppath);
        output.genes = extract_features(handle);
        output.annotations = null;
    } catch (e) {
        utils.freeCache(output.matrix);
        throw e;
    } finally {
        afile.removeH5(tmppath);
    }

    return output;
}

export async function serialize(formatted, embeddedSaver) {
    return [await rutils.standardSerialize(formatted.h5, "h5", embeddedSaver)];
}

export async function unserialize(values, embeddedLoader) {
    return { 
        "format": "10X",
        "h5": await rutils.standardUnserialize(values[0], embeddedLoader)
    };
}
