import * as scran from "scran.js";
import * as utils from "./../utils/general.js";
import * as rutils from "./../utils/reader.js";
import * as afile from "./../abstract/file.js";

export function abbreviate(args) {
    return { 
        "format": "10X", 
        "h5": rutils.formatFile(args.h5, true)
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
    let formatted = rutils.formatFile(args.h5, false)

    const tmppath = afile.realizeH5(formatted.content);
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

export class Reader {
    #h5;

    constructor(args, formatted = false) {
        if (formatted) {
            this.#h5 = args;
        } else {
            this.#h5 = rutils.formatFile(args.h5, false);
        }
        return;
    }

    load() {
        let output = {};

        const tmppath = afile.realizeH5(this.#h5.content);
        try {
            output.matrix = scran.initializeSparseMatrixFromHDF5(tmppath, "matrix");
            let handle = new scran.H5File(tmppath);
            output.genes = extract_features(handle);
            output.annotations = null;

            // Stop-gap solution to remove non-gene entries.
            rutils.subsetToGenes(output);

        } catch (e) {
            utils.freeCache(output.matrix);
            throw e;
        } finally {
            afile.removeH5(tmppath);
        }

        return output;
    }

    format() {
        return "10X";
    }

    async serialize(embeddedSaver) {
        return [await rutils.standardSerialize(this.#h5, "h5", embeddedSaver)];
    }
}

export async function unserialize(values, embeddedLoader) {
    return new Reader(await rutils.standardUnserialize(values[0], embeddedLoader), true);
}
