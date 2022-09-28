import * as scran from "scran.js";
import * as rutils from "./utils/index.js";
import * as afile from "../abstract/file.js";

export function abbreviate(args) {
    return {
        "format": "SummarizedExperiment", 
        "rds": rutils.formatFile(args.rds, true)
    };
}

function load_listData_names(lhandle) {
    let ndx = lhandle.findAttribute("names");
    if (ndx < 0) {
        return null;
    }

    let nhandle = lhandle.attribute(ndx);
    let names;
    try {
        names = nhandle.values();
    } finally {
        scran.free(nhandle);
    }

    if (names.length != lhandle.length()) {
        throw new Error("expected names to have same length as listData");
    }
    return names;
}

function load_data_frame(handle) {
    check_class(handle, { "DFrame": "S4Vectors" }, "DFrame");
    let output = {};

    // Loading the atomic columns.
    let ldx = handle.findAttribute("listData");
    if (ldx < 0) {
        throw new Error("expected a listData slot in a DFrame");
    }

    let columns = {};
    let lhandle = handle.attribute(ldx);
    try {
        if (!(lhandle instanceof scran.RdsGenericVector)) {
            throw new Error("listData slot should be a generic list");
        }

        let colnames = load_listData_names(lhandle);
        if (colnames == null) {
            throw new Error("expected the listData list to be named");
        }

        for (var i = 0; i < lhandle.length(); i++) {
            let curhandle = lhandle.load(i);
            try {
                if (curhandle instanceof scran.RdsVector && !(curhandle instanceof scran.RdsGenericVector)) {
                    columns[colnames[i]] = curhandle.values();
                }
            } finally {
                scran.free(curhandle);
            }
        }
    } finally {
        scran.free(lhandle);
    }
    output.columns = columns;

    // Loading the row names.
    let rndx = handle.findAttribute("rownames");
    if (rndx < 0) {
        throw new Error("expected a rownames slot in a DFrame");
    }

    let rnhandle = handle.attribute(rndx);
    try {
        if (rnhandle.type() instanceof scran.RdsStringVector) {
            output.row_names = rhandle.values();
        }
    } finally {
        scran.free(rnhandle);
    }

    return output;
}

function check_class(handle, accepted, base) {
    if (!(handle instanceof scran.RdsS4Object)) {
        throw new Error("expected an S4 object as the data frame");
    }

    for (const [k, v] of Object.entries(accepted)) {
        if (handle.className() == k && handle.packageName() == v) {
            return;
        }
    }
    throw new Error("object is not a " + base + " or one of its recognized subclasses");
}

function extract_annotations(handle) {
    let cdx = handle.findAttribute("colData");
    if (cdx < 0) {
        throw new Error("expected a colData slot in a SummarizedExperiment");
    }

    let chandle = handle.attribute(cdx);
    let output;
    try {
        output = load_data_frame(chandle);
    } finally {
        scran.free(chandle);
    }

    return output.columns;
}

function extract_NAMES(handle) {
    let nidx = handle.findAttribute("NAMES");
    if (nidx < 0) {
        return null;
    }

    let nhandle = handle.attribute(nidx);
    let output = null;
    try {
        if (nhandle instanceof scran.RdsStringVector) {
            output = nhandle.values();
        }
    } finally {
        scran.free(nhandle);
    }

    return output;
}

function extract_features(handle) {
    let rowdata = {};
    let names = null;

    let rrdx = handle.findAttribute("rowRanges");
    if (rrdx < 0) {
        // This is a base SummarizedExperiment.
        let rdx = handle.findAttribute("elementMetadata");
        if (rdx < 0) {
            throw new Error("expected a elementMetadata slot in a SummarizedExperiment");
        }

        let rhandle = handle.attribute(rdx);
        try {
            rowdata = load_data_frame(rhandle);
        } finally {
            scran.free(rhandle);
        }

        names = extract_NAMES(handle);

    } else {
        // Everything else is assumed to be an RSE.
        if (rrdx < 0) {
            throw new Error("expected a rowRanges slot in a RangedSummarizedExperiment");
        }

        let rrhandle = handle.attribute(rrdx);
        let output;
        try {
            let edx = rrhandle.findAttribute("elementMetadata");
            if (edx < 0) {
                throw new Error("expected an elementMetadata slot in a rowRanges object");
            }

            let ehandle = rrhandle.attribute(edx);
            try {
                rowdata = load_data_frame(ehandle);
            } finally {
                scran.free(ehandle);
            }

            let pidx = rrhandle.findAttribute("partitioning");
            if (pidx < 0) { // if absent, we'll assume it's a GRanges.
                names = extract_NAMES(rrhandle);
            } else { // otherwise, it's a GRangesList.
                let phandle = rrhandle.attribute(pidx);
                try {
                    names = extract_NAMES(phandle);
                } finally {
                    scran.free(phandle);
                }
            }

        } finally {
            scran.free(rrhandle);
        }
    }

    let output = { id: names };
    for (const [k, v] of Object.entries(rowdata.columns)) {
        if (k.match(/^sym/)) {
            output[k] = v;
        }
    }

    return output;
}

function extract_counts(handle) {
    let adx = handle.findAttribute("assays");
    if (adx < 0) {
        throw new Error("expected an assays slot in a SummarizedExperiment");
    }

    let output;
    let ahandle = handle.attribute(adx);
    try {
        let ddx = ahandle.findAttribute("data");
        if (ddx < 0) {
            throw new Error("expected a data slot in the assays object");
        }

        let dhandle = ahandle.attribute(ddx);
        try {
            let ldx = dhandle.findAttribute("listData");
            if (ldx) {
                throw new Error("expected a listData slot in the SimpleList of the assays");
            }

            let lhandle = dhandle.attribute(ldx);
            try {
                let names = load_listData_names(lhandle);

                let chosen = 1;
                if (names != null) {
                    for (var n = 0; n < names.length; n++) {
                        if (names[n].match(/^count/i)) {
                            chosen = n;
                            break;
                        }
                    }
                }

                let xhandle = lhandle.load(chosen);
                try {
                    output = scran.initializeSparseMatrixFromRds(xhandle);
                } finally {
                    scran.free(xhandle);
                }

            } finally {
                scran.free(lhandle);
            }

        } finally {
            scran.free(dhandle);
        }

    } finally {
        scran.free(ahandle);
    }

    return output;
}

function preflight_raw(handle) {
    check_class(handle, { 
        "SummarizedExperiment": "SummarizedExperiment",
        "RangedSummarizedExperiment": "SummarizedExperiment",
        "SingleCellExperiment": "SingleCellExperiment",
        "SpatialExperiment": "SpatialExperiment"
    }, "SummarizedExperiment");

    let output = {};

    // Loading the gene info.
    output.genes = extract_features(handle);

    // Loading the cell annotation.
    let anno = extract_annotations(handle);
    let acols = Object.entries(anno);
    if (acols.length) { 
        output.annotations = {};
        for (const [k, v] of acols) {
            output.annotations[k] = rutils.summarizeArray(v);
        }
    } else {
        output.annotations = null;
    }

    return output;
}

export async function preflight(args) {
    let output_anno = {};
    let output_feat = {};

    let rds_data = rutils.formatFile(args.rds, false);
    let rds_stuff = afile.realizeFile(rds_data.content);

    let rdshandle = scran.readRds(rds_stuff);
    try {
        let rhandle = rdshandle.value();
        try {
            let rna = preflight_raw(rhandle);
            output_anno = rna.annotations;
            output_feat.RNA = rna.genes;
        } finally {
            scran.free(rhandle);
        }
    } finally {
        scran.free(rdshandle);
    }

    return { annotations: output_anno, genes: output_feat };
}

export class Reader {
    #rds;

    constructor(args, formatted = false) {
        if (!formatted) {
            this.#rds = rutils.formatFile(args.rds, false);
        } else {
            this.#rds = args.rds;
        }
    }

    load() {
        let rds_stuff = afile.realizeFile(this.#rds.content);
        let output = {};

        let rdshandle = scran.readRds(rds_stuff);
        try {

            let rhandle = rdshandle.value();
            try {

                output.matrix = new scran.MultiMatrix;
                try {
                    let loaded = extract_counts(rhandle);

                    let out_mat = loaded.matrix;
                    let out_ids = loaded.row_ids;
                    output.matrix.add("RNA", out_mat);

                    let raw_gene_info = extract_features(rhandle);
                    let gene_info = rutils.reorganizeGenes(out_mat.numberOfRows(), out_ids, raw_gene_info);
                    output.genes = { RNA: gene_info };
                    output.row_ids = { RNA: out_ids };
                    output.annotations = extract_annotations(rhandle) ;

                } catch (e) {
                    scran.free(output.matrix);
                    throw e;
                }

            } finally {
                scran.free(rhandle);
            }

        } finally {
            scran.free(rdshandle);
        }

        return output;
    }

    format() {
        return "SummarizedExperiment";
    }

    async serialize(embeddedSaver) {
        return [await rutils.standardSerialize(this.#rds, "rds", embeddedSaver)];
    }
}

export async function unserialize(values, embeddedLoader) {
    let args = {};

    // This should contain 'rds'.
    for (const x of values) {
        args[x.type] = await rutils.standardUnserialize(x, embeddedLoader);
    }

    return new Reader(args, true);
}
