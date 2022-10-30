import * as scran from "scran.js";
import * as bioc from "bioconductor";
import * as afile from "./abstract/file.js";
import * as eutils from "./utils/extract.js";
import * as futils from "./utils/features.js";

function load_listData_names(lhandle) {
    let ndx = lhandle.findAttribute("names");
    if (ndx < 0) {
        return null;
    }

    let nhandle;
    let names;
    try {
        nhandle = lhandle.attribute(ndx);
        names = nhandle.values();
    } catch(e) {
        throw new Error("failed to load listData names; " + e.message);
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
    let columns = {};
    let lhandle;
    try {
        lhandle = handle.attribute("listData");
        if (!(lhandle instanceof scran.RdsGenericVector)) {
            throw new Error("listData slot should be a generic list");
        }

        let colnames = load_listData_names(lhandle);
        if (colnames == null) {
            throw new Error("expected the listData list to be named");
        }

        for (var i = 0; i < lhandle.length(); i++) {
            let curhandle;
            try {
                curhandle = lhandle.load(i);
                if (curhandle instanceof scran.RdsVector && !(curhandle instanceof scran.RdsGenericVector)) {
                    let curcol = curhandle.values();
                    columns[colnames[i]] = curcol;
                    output.nrow = curcol.length;
                }
            } finally {
                scran.free(curhandle);
            }
        }
    } catch(e) {
        throw new Error("failed to retrieve data from DataFrame's listData; " + e.message);
    } finally {
        scran.free(lhandle);
    }
    output.columns = columns;

    // Loading the row names.
    let rnhandle;
    try {
        rnhandle = handle.attribute("rownames");
        if (rnhandle instanceof scran.RdsStringVector) {
            output.row_names = rnhandle.values();
            output.nrow = output.row_names.length;
        }
    } catch(e) {
        throw new Error("failed to retrieve row names from DataFrame; " + e.message);
    } finally {
        scran.free(rnhandle);
    }

    // Loading the number of rows.
    if (!("nrow" in output)) {
        let nrhandle;
        try {
            nrhandle = handle.attribute("nrows");
            if (!(nrhandle instanceof scran.RdsIntegerVector)) {
                throw new Error("expected an integer vector as the 'nrows' slot");
            }
            let NR = nrhandle.values();
            if (NR.length != 1) {
                throw new Error("expected an integer vector of length 1 as the 'nrows' slot");
            }
            output.nrow = NR[0];
        } catch (e) {
            throw new Error("failed to retrieve nrows from DataFrame; " + e.message);
        } finally {
            scran.free(nrhandle);
        }
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

function extract_NAMES(handle) {
    let nidx = handle.findAttribute("NAMES");
    if (nidx < 0) {
        return null;
    }

    let nhandle;
    let output = null;
    try {
        nhandle = handle.attribute(nidx);
        if (nhandle instanceof scran.RdsStringVector) {
            output = nhandle.values();
        }
    } catch(e) {
        throw new Error("failed to extract NAMES; " + e.message);
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
        let rhandle;
        try {
            rhandle = handle.attribute("elementMetadata");
            rowdata = load_data_frame(rhandle);
        } catch(e) {
            throw new Error("failed to extract features from the rowData; " + e.message);
        } finally {
            scran.free(rhandle);
        }
        names = extract_NAMES(handle);

    } else {
        // Everything else is assumed to be an RSE.
        let rrhandle;
        let output;
        try {
            rrhandle = handle.attribute(rrdx);
            let ehandle = rrhandle.attribute("elementMetadata");
            try {
                rowdata = load_data_frame(ehandle);
            } catch(e) {
                throw new Error("failed to extract mcols from the rowRanges; " + e.message);
            } finally {
                scran.free(ehandle);
            }

            let pidx = rrhandle.findAttribute("partitioning");
            if (pidx < 0) { // if absent, we'll assume it's a GRanges.
                let r2handle;
                try {
                    r2handle = rrhandle.attribute("ranges");
                    names = extract_NAMES(r2handle);
                } catch(e) {
                    throw new Error("failed to extract names from the rowRanges; " + e.message);
                } finally {
                    scran.free(r2handle);
                }
            } else { // otherwise, it's a GRangesList.
                let phandle;
                try {
                    phandle = rrhandle.attribute(pidx);
                    names = extract_NAMES(phandle);
                } catch(e) {
                    throw new Error("failed to extract names from the rowRanges; " + e.message);
                } finally {
                    scran.free(phandle);
                }
            }

        } catch(e) {
            throw new Error("failed to extract features from the rowRanges; " + e.message);
        } finally {
            scran.free(rrhandle);
        }
    }

    if (names == null) {
        return new bioc.DataFrame({}, { numberOfRows: rowdata.nrow });
    } else {
        let output = { id: names };
        for (const [k, v] of Object.entries(rowdata.columns)) {
            if (k.match(/^sym/)) {
                output[k] = v;
            }
        }
        return new bioc.DataFrame(output);
    }
}

function extract_counts(handle) {
    let output;
    let ahandle;
    let dhandle;
    let lhandle;

    try {
        ahandle = handle.attribute("assays");
        dhandle = ahandle.attribute("data");
        lhandle = dhandle.attribute("listData");

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

        let xhandle;
        try {
            xhandle = lhandle.load(chosen);
            output = scran.initializeSparseMatrixFromRds(xhandle);
        } catch(e) {
            throw new Error("failed to initialize sparse matrix from assay; " + e.message);
        } finally {
            scran.free(xhandle);
        }

    } catch(e) {
        throw new Error("failed to extract assay data; " + e.message);
    } finally {
        scran.free(ahandle);
        scran.free(lhandle);
        scran.free(dhandle);
    }

    return output;
}

function extract_alt_exps(handle) {
    let indx = handle.findAttribute("int_colData");
    if (indx < 0) {
        return {};
    }

    let output = {};
    let in_handle;
    let inld_handle;
    let innn_handle;
    let ae_handle;
    let aeld_handle;
    let aenn_handle;

    try {
        in_handle = handle.attribute(indx);
        let inld_dx = in_handle.findAttribute("listData");
        if (inld_dx < 0) {
            return {};
        }

        inld_handle = in_handle.attribute(inld_dx);
        let innn_dx = inld_handle.findAttribute("names");
        if (innn_dx < 0) {
            return {};
        }

        innn_handle = inld_handle.attribute(innn_dx);
        let in_names = innn_handle.values();
        let ae_dx = in_names.indexOf("altExps");
        if (ae_dx < 0) {
            return {};
        }

        ae_handle = inld_handle.load(ae_dx);
        let aeld_dx = ae_handle.findAttribute("listData");
        if (aeld_dx < 0) {
            return {};
        }

        aeld_handle = ae_handle.attribute(aeld_dx);
        let aenn_dx = aeld_handle.findAttribute("names");
        if (aenn_dx < 0) {
            return {};
        }

        aenn_handle = aeld_handle.attribute(aenn_dx);
        let ae_names = aenn_handle.values();

        for (var i = 0; i < ae_names.length; i++) {
            let curhandle;
            try {
                curhandle = aeld_handle.load(i);
                let asehandle = curhandle.attribute("se");
                output[ae_names[i]] = asehandle;
                check_for_se(asehandle);
            } catch (e) {
                throw new Error("failed to load alternative Experiment '" + ae_names[i] + "'; " + e.message);
            } finally {
                scran.free(curhandle);
            }
        }

    } catch(e) {
        for (const v of Object.values(output)) {
            scran.free(v);
        }
        throw e;

    } finally {
        scran.free(aenn_handle);
        scran.free(aeld_handle);
        scran.free(innn_handle);
        scran.free(inld_handle);
        scran.free(in_handle);
    }

    return output;
}

function check_for_se(handle) {
    check_class(handle, { 
        "SummarizedExperiment": "SummarizedExperiment",
        "RangedSummarizedExperiment": "SummarizedExperiment",
        "SingleCellExperiment": "SingleCellExperiment",
        "SpatialExperiment": "SpatialExperiment"
    }, "SummarizedExperiment");
}

/**
 * Dataset stored as a `SummarizedExperiment` object (or one of its subclasses) inside an RDS file.
 */
export class SummarizedExperimentDataset {
    #rds_file;
    #rds_handle;
    #se_handle;
    #alt_handles;

    #raw_features;
    #raw_cells;

    /**
     * @param {SimpleFile|string|Uint8Array|File} rdsFile - Contents of a RDS file.
     * On browsers, this may be a File object.
     * On Node.js, this may also be a string containing a file path.
     */
    constructor(rdsFile) {
        super();

        if (rdsFile instanceof afile.SimpleFile) {
            this.#rds_file = rdsFile;
        } else {
            this.#rds_file = new afile.SimpleFile(rdsFile);
        }

        this.#raw_features = null;
        this.#raw_cells = null;
        this.#rds_handle = null;
    }

    static format() {
        return "SummarizedExperiment";
    }

    abbreviate() {
        return {
            "rds": {
                name: this.#rds_file.name(),
                size: this.#rds_file.size()
            }
        };
    }

    #initialize() {
        if (this.#rds_handle !== null) {
            return;
        }

        this.#rds_handle = scran.readRds(this.#rds_file.content());
        this.#se_handle = this.#rds_handle.value();
        try {
            check_for_se(this.#se_handle);
            this.#alt_handles = extract_alt_exps(this.#se_handle);
        } catch (e) {
            this.#se_handle.free();
            this.#rds_handle.free();
            throw e;
        }
    }

    #features() {
        if (this.#raw_features !== null) {
            return;
        }

        this.#initialize();
        this.#raw_features = { "": extract_features(this.#se_handle) };

        for (const [k, v] of Object.entries(this.#alt_handles)) {
            try {
                this.#raw_features[k] = extract_features(v);
            } catch (e) {
                console.warn("failed to extract features for alternative Experiment '" + k + "'; " + e.message);
            }
        }

        return;
    }

    #cells() {
        if (this.#raw_cells !== null) {
            return;
        }

        this.#initialize();
        let info;
        let chandle = this.#se_handle.attribute("colData");
        try {
            info = load_data_frame(chandle);
        } catch(e) {
            throw new Error("failed to extract colData from a SummarizedExperiment; " + e.message);
        } finally {
            scran.free(chandle);
        }

        this.#raw_cells = new bioc.DataFrame(info.columns, { numberOfRows: info.nrow });
        return;
    }

    annotations() {
        this.#initialize();
        this.#features();
        this.#cells();

        let copy = {};
        for (const [k, v] of Object.entries(this.#raw_features)) {
            copy[k] = bioc.CLONE(v);
        }

        let anno = {};
        for (const k of this.#raw_cells.columnNames()) {
            anno[k] = eutils.summarizeArray(this.#raw_cells.column(k));
        }

        return {
            features: copy,
            cells: {
                number: this.#raw_cells.numberOfColumns(),
                summary: anno
            }
        };
    }

    load() {
        this.#initialize();
        this.#features();
        this.#cells();

        let output = { 
            matrix: new scran.MultiMatrix,
            row_ids: {},
            features: {},
            cells: bioc.CLONE(this.#raw_cells)
        };

        try {
            let loaded = extract_counts(this.#se_handle);
            let out_mat = loaded.matrix;
            let out_ids = loaded.row_ids;
            output.matrix.add("", out_mat);
            output.row_ids[""] = out_ids;
            output.features[""] = bioc.SLICE(this.#raw_features[""], out_ids);

            for (const [k, v] of Object.entries(this.#alt_handles)) {
                let alt = extract_counts(v);
                output.matrix.add(k, alt.matrix);
                output.row_ids[k] = alt.row_ids;
                output.features[k] = bioc.SLICE(this.#raw_features[k], alt.row_ids);
            }

        } catch (e) {
            scran.free(output.matrix);
            throw e;
        }

        return output;
    }

    free() {
        this.#se_handle.free();
        for (const v of Object.values(this.#alt_handles)) {
            v.free();
        }
        this.#rds_handle.free();
    }

    serialize() {
        return [ { type: "rds", file: this.#rds_file } ];
    }

    static async unserialize(files) {
        if (files.length != 1 || files[0].type != "rds") {
            throw new Error("expected exactly one file of type 'rds' for SummarizedExperiment unserialization");
        }
        return new SummarizedExperimentDataset(files[0].file);
    }
}
