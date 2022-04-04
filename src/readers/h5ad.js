import * as scran from "scran.js";
import * as utils from "./../utils/general.js";
import * as rutils from "./../utils/reader.js";
import * as afile from "./../abstract/file.js";

export function format(args, signatureOnly) {
    return { 
        "format": "H5AD", 
        "h5": rutils.formatFile(args.h5, signatureOnly)
    };
}

function extract_features(handle) {
    let genes = null;

    if ("var" in handle.children && handle.children["var"] == "Group") {
        let vhandle = handle.open("var");
        let index = rutils.extractHDF5Strings(vhandle, "_index");
        if (index !== null) {
            genes = { "_index": index };

            for (const [key, val] of Object.entries(vhandle.children)) {
                if (val === "DataSet" && (key.match(/name/i) || key.match(/symb/i))) {
                    let dhandle2 = vhandle.open(key);
                    if (dhandle2.type == "String") {
                        genes[key] = dhandle2.load();
                    }
                }
            }
        }
    }

    return genes;
}

function extract_annotations(handle, { namesOnly = false } = {}) {
    let annotations = null;

    if ("obs" in handle.children && handle.children["obs"] == "Group") {
        let ohandle = handle.open("obs");
        annotations = {};

        // Maybe it has names, maybe not, who knows; let's just add what's there.
        let index = rutils.extractHDF5Strings(ohandle, "_index");
        if (index !== null) {
            annotations["_index"] = index;
        }

        for (const [key, val] of Object.entries(ohandle.children)) {
            if (val != "DataSet") {
                continue;
            }
            let dhandle = ohandle.open(key);
            if (dhandle.type != "Other") {
                annotations[key] = (!namesOnly ? dhandle.load() : null);
            }
        }

        if (!namesOnly && "__categories" in ohandle.children && ohandle.children["__categories"] == "Group") {
            let chandle = ohandle.open("__categories");

            for (const [key, val] of Object.entries(chandle.children)) {
                if (key in annotations) {
                    let cats = rutils.extractHDF5Strings(chandle, key);
                    if (cats.type !== null) {
                        annotations[key] = {
                            "type": "factor",
                            "index": val,
                            "factor": cats
                        }
                    }
                }
            }
        }
    }

    if (namesOnly && annotations !== null) {
        return Object.keys(annotations);
    } else {
        return annotations;
    }
}

export function preflight(args) {
    let output = {};
    let formatted = format(args, false)

    const tmppath = afile.realizeH5(formatted.h5.content);
    try {
        let handle = new scran.H5File(tmppath);
        output.genes = extract_features(handle);
        let annotations = extract_annotations(handle, { load: false });
        output.annotations = Object.keys(annotations); 
    } finally {
        afile.removeH5(tmppath);
    }

    return output;
}

export function load(formatted) {
    let output = {};

    const tmppath = afile.realizeH5(formatted.h5.content);
    try {
        output.matrix = scran.initializeSparseMatrixFromHDF5(tmppath, "X");
        let handle = new scran.H5File(tmppath);
        output.genes = extract_features(handle); 
        output.annotations = extract_annotations(handle);
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
        "format": "H5AD",
        "h5": await rutils.standardUnserialize(values[0], embeddedLoader)
    };
}
