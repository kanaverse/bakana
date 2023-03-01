import * as scran from "scran.js";
import * as bioc from "bioconductor";
import * as afile from "./abstract/file.js";
import * as eutils from "./utils/extract.js";
import * as futils from "./utils/features.js";

/**************************
 ******* Internals ********
 **************************/

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

const acceptable_df_subclasses = { "DFrame": "S4Vectors" };

function load_data_frame(handle) {
    check_class(handle, acceptable_df_subclasses, "DFrame");

    let columns = {};
    let colnames = [];
    let lhandle;
    try {
        lhandle = handle.attribute("listData");
        if (!(lhandle instanceof scran.RdsGenericVector)) {
            throw new Error("listData slot should be a generic list");
        }

        colnames = load_listData_names(lhandle);
        if (colnames == null) {
            throw new Error("expected the listData list to be named");
        }

        for (var i = 0; i < lhandle.length(); i++) {
            let curhandle;
            try {
                curhandle = lhandle.load(i);

                if (curhandle instanceof scran.RdsVector && !(curhandle instanceof scran.RdsGenericVector)) {
                    let curcol = curhandle.values();

                    // Expand factors, if we detect them.
                    if (curhandle.findAttribute("class") >= 0) {
                        let clshandle;
                        let levhandle;
                        try {
                            clshandle = curhandle.attribute("class");
                            if (clshandle.values().indexOf("factor") >= 0 && curhandle.findAttribute("levels") >= 0) {
                                levhandle = curhandle.attribute("levels");
                                let copy = curcol.slice();
                                copy.forEach((x, i) => { copy[i] = x - 1 }); // get back to 0-based indices.
                                curcol = bioc.SLICE(levhandle.values(), copy);
                            }
                        } finally {
                            scran.free(clshandle);
                            scran.free(levhandle);
                        }
                    }

                    columns[colnames[i]] = curcol;

                } else if (curhandle instanceof scran.RdsS4Object && check_acceptable_class(curhandle, acceptable_df_subclasses)) {
                    // Handle nested DataFrames.
                    columns[colnames[i]] = load_data_frame(curhandle);
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

    // Loading the row names.
    let rnhandle;
    let rownames = null;
    try {
        rnhandle = handle.attribute("rownames");
        if (rnhandle instanceof scran.RdsStringVector) {
            rownames = rnhandle.values();
        }
    } catch(e) {
        throw new Error("failed to retrieve row names from DataFrame; " + e.message);
    } finally {
        scran.free(rnhandle);
    }

    // Loading the number of rows.
    let nrows = null;
    if (colnames.length == 0 && rownames == null) {
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
            nrows = NR[0];
        } catch (e) {
            throw new Error("failed to retrieve nrows from DataFrame; " + e.message);
        } finally {
            scran.free(nrhandle);
        }
    }

    return new bioc.DataFrame(columns, { columnOrder: colnames, rowNames: rownames, numberOfRows: nrows });
}

function check_acceptable_class(handle, accepted) {
    for (const [k, v] of Object.entries(accepted)) {
        if (handle.className() == k && handle.packageName() == v) {
            return true;
        }
    }
    return false;
}

function check_class(handle, accepted, base) {
    if (!(handle instanceof scran.RdsS4Object)) {
        throw new Error("expected an S4 object as the data frame");
    }
    if (!check_acceptable_class(handle, accepted)) {
        throw new Error("object is not a " + base + " or one of its recognized subclasses");
    }
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
    let rowdata;
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

    if (names !== null) {
        rowdata.$setRowNames(names);
    }
    return rowdata;
}

function extract_assay_names(handle) {
    let output;
    let ahandle;
    let dhandle;
    let lhandle;

    try {
        ahandle = handle.attribute("assays");
        dhandle = ahandle.attribute("data");
        lhandle = dhandle.attribute("listData");

        output = load_listData_names(lhandle);
        if (output == null) {
            output = new Array(lhandle.length());
            output.fill(null);
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

function extract_assay(handle, assay, forceInteger) {
    let output;
    let ahandle;
    let dhandle;
    let lhandle;

    try {
        ahandle = handle.attribute("assays");
        dhandle = ahandle.attribute("data");
        lhandle = dhandle.attribute("listData");

        // Choosing the assay index.
        let chosen = null;
        if (typeof assay == "string") {
            let names = load_listData_names(lhandle);
            if (assay !== null && names != null) {
                for (var n = 0; n < names.length; n++) {
                    if (names[n] == assay) {
                        chosen = n;
                        break;
                    }
                }
            }
            if (chosen == null) {
                throw new Error("no assay named '" + assay + "'");
            }
        } else {
            if (assay >= lhandle.length()) {
                throw new Error("assay index " + String(assay) + " out of range");
            }
            chosen = assay;
        }

        let xhandle;
        try {
            xhandle = lhandle.load(chosen);
            output = scran.initializeSparseMatrixFromRds(xhandle, { forceInteger });
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
    let output = { handles: {}, order: [] };
    let indx = handle.findAttribute("int_colData");
    if (indx < 0) {
        return output;
    }

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
            return output;
        }

        inld_handle = in_handle.attribute(inld_dx);
        let innn_dx = inld_handle.findAttribute("names");
        if (innn_dx < 0) {
            return output;
        }

        innn_handle = inld_handle.attribute(innn_dx);
        let in_names = innn_handle.values();
        let ae_dx = in_names.indexOf("altExps");
        if (ae_dx < 0) {
            return output;
        }

        ae_handle = inld_handle.load(ae_dx);
        let aeld_dx = ae_handle.findAttribute("listData");
        if (aeld_dx < 0) {
            return output;
        }

        aeld_handle = ae_handle.attribute(aeld_dx);
        let aenn_dx = aeld_handle.findAttribute("names");
        if (aenn_dx < 0) {
            return output;
        }

        aenn_handle = aeld_handle.attribute(aenn_dx);
        let ae_names = aenn_handle.values();

        for (var i = 0; i < ae_names.length; i++) {
            let curhandle;
            try {
                curhandle = aeld_handle.load(i);
                let asehandle = curhandle.attribute("se");
                output.handles[ae_names[i]] = asehandle;
                output.order.push(ae_names[i]);
                check_for_se(asehandle);
            } catch (e) {
                throw new Error("failed to load alternative Experiment '" + ae_names[i] + "'; " + e.message);
            } finally {
                scran.free(curhandle);
            }
        }

    } catch(e) {
        for (const v of Object.values(output.handles)) {
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

function extract_reduced_dims(handle) {
    let output = { handles: {}, order: [] };
    let indx = handle.findAttribute("int_colData");
    if (indx < 0) {
        return output;
    }

    let in_handle;
    let inld_handle;
    let innn_handle;
    let rd_handle;
    let rdld_handle;
    let rdnn_handle;

    try {
        in_handle = handle.attribute(indx);
        let inld_dx = in_handle.findAttribute("listData");
        if (inld_dx < 0) {
            return output;
        }

        inld_handle = in_handle.attribute(inld_dx);
        let innn_dx = inld_handle.findAttribute("names");
        if (innn_dx < 0) {
            return output;
        }

        innn_handle = inld_handle.attribute(innn_dx);
        let in_names = innn_handle.values();
        let rd_dx = in_names.indexOf("reducedDims");
        if (rd_dx < 0) {
            return output;
        }

        rd_handle = inld_handle.load(rd_dx);
        let rdld_dx = rd_handle.findAttribute("listData");
        if (rdld_dx < 0) {
            return output;
        }

        rdld_handle = rd_handle.attribute(rdld_dx);
        let rdnn_dx = rdld_handle.findAttribute("names");
        if (rdnn_dx < 0) {
            return output;
        }

        rdnn_handle = rdld_handle.attribute(rdnn_dx);
        let rd_names = rdnn_handle.values();

        for (var i = 0; i < rd_names.length; i++) {
            let curhandle;
            try {
                curhandle = rdld_handle.load(i);
                let okay = false;

                if (curhandle.type() == "double" && curhandle.findAttribute("dim") >= 0) { // only accepting double-precision matrics.
                    let dimhandle = curhandle.attribute("dim");
                    if (dimhandle.length() == 2) {
                        output.handles[rd_names[i]] = { handle: curhandle, dimensions: dimhandle.values() };
                        output.order.push(rd_names[i]);
                        okay = true;
                    }
                }

                if (!okay) {
                    scran.free(curhandle);
                }
            } catch (e) {
                throw new Error("failed to load reduced dimension '" + rd_names[i] + "'; " + e.message);
            }
        }

    } catch(e) {
        for (const v of Object.values(output.handles)) {
            scran.free(v.handle);
        }
        throw e;

    } finally {
        scran.free(rdnn_handle);
        scran.free(rdld_handle);
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

const main_experiment_name = "";

/************************
 ******* Dataset ********
 ************************/

/**
 * Dataset stored as a SummarizedExperiment object (or one of its subclasses) inside an RDS file.
 */
export class SummarizedExperimentDataset {
    #rds_file;

    #rds_handle;
    #se_handle;
    #alt_handles;
    #alt_handle_order;

    #raw_features;
    #raw_cells;

    #rnaCountAssay;
    #adtCountAssay;
    #crisprCountAssay;

    #rnaExperiment;
    #adtExperiment;
    #crisprExperiment;

    #primaryRnaFeatureIdColumn;
    #primaryAdtFeatureIdColumn;
    #primaryCrisprFeatureIdColumn;

    #dump_summary(fun) {
        let files = [{ type: "rds", file: fun(this.#rds_file) }];
        let options = {
            rnaCountAssay: this.#rnaCountAssay,
            adtCountAssay: this.#adtCountAssay,
            crisprCountAssay: this.#crisprCountAssay,
            rnaExperiment: this.#rnaExperiment,
            adtExperiment: this.#adtExperiment,
            crisprExperiment: this.#crisprExperiment,
            primaryRnaFeatureIdColumn: this.#primaryRnaFeatureIdColumn,
            primaryAdtFeatureIdColumn: this.#primaryAdtFeatureIdColumn,
            primaryCrisprFeatureIdColumn: this.#primaryCrisprFeatureIdColumn
        };
        return { files, options };
    }

    /**
     * @param {SimpleFile|string|Uint8Array|File} rdsFile - Contents of a RDS file.
     * On browsers, this may be a File object.
     * On Node.js, this may also be a string containing a file path.
     * @param {object} [options={}] - Optional parameters.
     * @param {string|number} [options.rnaCountAssay=0] - See {@linkcode SummarizedExperimentDataset#setRnaCountAssay setRnaCountAssay}.
     * @param {string|number} [options.adtCountAssay=0] - See {@linkcode SummarizedExperimentDataset#setAdtCountAssay setAdtCountAssay}.
     * @param {string|number} [options.crisprCountAssay=0] - See {@linkcode SummarizedExperimentDataset#setCrisprCountAssay setCrisprCountAssay}.
     * @param {?(string|number)} [options.rnaExperiment=""] - See {@linkcode SummarizedExperimentDataset#setRnaExperiment setRnaExperiment}.
     * @param {?(string|number)} [options.adtExperiment="Antibody Capture"] - See {@linkcode SummarizedExperimentDataset#setAdtExperiment setAdtExperiment}.
     * @param {?(string|number)} [options.crisprExperiment="CRISPR Guide Capture"] - See {@linkcode SummarizedExperimentDataset#setCrisprExperiment setCrisprExperiment}.
     * @param {string|number} [options.primaryRnaFeatureIdColumn=null] - See {@linkcode SummarizedExperimentDataset#setPrimaryRnaFeatureIdColumn setPrimaryRnaFeatureIdColumn}.
     * @param {string|number} [options.primaryAdtFeatureIdColumn=null] - See {@linkcode SummarizedExperimentDataset#setPrimaryAdtFeatureIdColumn setPrimaryAdtFeatureIdColumn}.
     * @param {string|number} [options.primaryCrisprFeatureIdColumn=null] - See {@linkcode SummarizedExperimentDataset#setPrimaryCrisprFeatureIdColumn setPrimaryCrisprFeatureIdColumn}.
     */
    constructor(rdsFile, { 
        rnaCountAssay = 0, 
        adtCountAssay = 0, 
        crisprCountAssay = 0,
        rnaExperiment = "", 
        adtExperiment = "Antibody Capture", 
        crisprExperiment = "CRISPR Guide Capture",
        primaryRnaFeatureIdColumn = null, 
        primaryAdtFeatureIdColumn = null,
        primaryCrisprFeatureIdColumn = null 
    } = {}) {
        if (rdsFile instanceof afile.SimpleFile) {
            this.#rds_file = rdsFile;
        } else {
            this.#rds_file = new afile.SimpleFile(rdsFile);
        }

        this.#rnaCountAssay = rnaCountAssay;
        this.#adtCountAssay = adtCountAssay;
        this.#crisprCountAssay = crisprCountAssay;

        this.#rnaExperiment = rnaExperiment;
        this.#adtExperiment = adtExperiment;
        this.#crisprExperiment = crisprExperiment;

        this.#primaryRnaFeatureIdColumn = primaryRnaFeatureIdColumn;
        this.#primaryAdtFeatureIdColumn = primaryAdtFeatureIdColumn;
        this.#primaryCrisprFeatureIdColumn = primaryCrisprFeatureIdColumn;

        this.clear();
    }

    /**
     * @param {string|number} i - Name or index of the assay containing the RNA count matrix.
     */
    setRnaCountAssay(i) {
        this.#rnaCountAssay = i;
        return;
    }

    /**
     * @param {string|number} i - Name or index of the assay containing the ADT count matrix.
     */
    setAdtCountAssay(i) {
        this.#adtCountAssay = i;
        return;
    }

    /**
     * @param {string|number} i - Name or index of the assay containing the CRISPR count matrix.
     */
    setCrisprCountAssay(i) {
        this.#adtCountAssay = i;
        return;
    }

    /**
     * @param {?(string|number)} i - Name or index of the alternative experiment containing gene expression data.
     * If `i` is `null` or invalid (e.g., out of range index, unavailable name), it is ignored and no RNA data is assumed to be present.
     * If `i` is an empty string, the main experiment is assumed to contain the gene expression data.
     */
    setRnaExperiment(i) {
        this.#rnaExperiment = i;
        return;
    }

    /**
     * @param {?(string|number)} i - Name or index of the alternative experiment containing ADT data.
     * If `i` is `null` or invalid (e.g., out of range index, unavailable name), it is ignored and no ADTs are assumed to be present.
     * If `i` is an empty string, the main experiment is assumed to contain the ADT data.
     */
    setAdtExperiment(i) {
        this.#adtExperiment = i;
        return;
    }

    /**
     * @param {?(string|number)} i - Name or index of the alternative experiment containing CRISPR guide data.
     * If `i` is `null` or invalid (e.g., out of range index, unavailable name), it is ignored and no CRISPR guides are assumed to be present.
     * If `i` is an empty string, the main experiment is assumed to contain the guide data.
     */
    setCrisprExperiment(i) {
        this.#crisprExperiment = i;
        return;
    }

    /**
     * @param {string|number} i - Name or index of the column of the `features` {@linkplain external:DataFrame DataFrame} that contains the primary feature identifier for gene expression.
     *
     * If `i` is `null` or invalid (e.g., out of range index, unavailable name), it is ignored and the primary identifier is defined as the existing row names.
     * However, if no row names are present in the SummarizedExperiment, no primary identifier is defined.
     */
    setPrimaryRnaFeatureIdColumn(i) {
        this.#primaryRnaFeatureIdColumn = i;
        return;
    }

    /**
     * @param {string|number} i - Name or index of the column of the `features` {@linkplain external:DataFrame DataFrame} that contains the primary feature identifier for the ADTs.
     *
     * If `i` is `null` or invalid (e.g., out of range index, unavailable name), it is ignored and the primary identifier is defined as the existing row names.
     * However, if no row names are present in the SummarizedExperiment, no primary identifier is defined.
     */
    setPrimaryAdtFeatureIdColumn(i) {
        this.#primaryAdtFeatureIdColumn = i;
        return;
    }

    /**
     * @param {string|number} i - Name or index of the column of the `features` {@linkplain external:DataFrame DataFrame} that contains the primary feature identifier for the CRISPR guides.
     *
     * If `i` is `null` or invalid (e.g., out of range index, unavailable name), it is ignored and the existing row names (if they exist) are used as the primary identifier.
     * However, if no row names are present in the SummarizedExperiment, no primary identifier is defined.
     */
    setPrimaryCrisprFeatureIdColumn(i) {
        this.#primaryCrisprFeatureIdColumn = i;
        return;
    }

    /**
     * Destroy caches if present, releasing the associated memory.
     * This may be called at any time but only has an effect if `cache = true` in {@linkcode SummarizedExperimentDataset#load load} or {@linkcodeSummarizedExperimentDataset#summary summary}.
     */
    clear() {
        scran.free(this.#se_handle);
        if (typeof this.#alt_handles != 'undefined' && this.#alt_handles !== null) {
            for (const v of Object.values(this.#alt_handles)) {
                scran.free(v);
            }
        }
        scran.free(this.#rds_handle);

        this.#se_handle = null;
        this.#alt_handles = null;
        this.#rds_handle = null;

        this.#raw_features = null;
        this.#raw_cells = null;
    }

    /**
     * @return {string} Format of this dataset class.
     * @static
     */
    static format() {
        return "SummarizedExperiment";
    }

    /**
     * @return {object} Object containing the abbreviated details of this dataset.
     */
    abbreviate() {
        return this.#dump_summary(f => { return { name: f.name(), size: f.size() }; });
    }

    #initialize() {
        if (this.#rds_handle !== null) {
            return;
        }

        this.#rds_handle = scran.readRds(this.#rds_file.content());
        this.#se_handle = this.#rds_handle.value();
        try {
            check_for_se(this.#se_handle);
            const { handles, order } = extract_alt_exps(this.#se_handle);
            this.#alt_handles = handles;
            this.#alt_handle_order = order;
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
        this.#raw_features = {};
        this.#raw_features[main_experiment_name] = extract_features(this.#se_handle);

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
        let chandle = this.#se_handle.attribute("colData");
        try {
            this.#raw_cells = load_data_frame(chandle);
        } catch(e) {
            throw new Error("failed to extract colData from a SummarizedExperiment; " + e.message);
        } finally {
            scran.free(chandle);
        }

        return;
    }

    /**
     * @param {object} [options={}] - Optional parameters.
     * @param {boolean} [options.cache=false] - Whether to cache the results for re-use in subsequent calls to this method or {@linkcode SummarizedExperimentDataset#load load}.
     * If `true`, users should consider calling {@linkcode SummarizedExperimentDataset#clear clear} to release the memory once this dataset instance is no longer needed.
     * 
     * @return {object} Object containing the per-feature and per-cell annotations.
     * This has the following properties:
     *
     * - `modality_features`: an object where each key is a modality name and each value is a {@linkplain external:DataFrame DataFrame} of per-feature annotations for that modality.
     * - `cells`: a {@linkplain external:DataFrame DataFrame} of per-cell annotations.
     * - `modality_assay_names`: an object where each key is a modality name and each value is an Array containing the names of available assays for that modality.
     *    Unnamed assays are represented as `null` names.
     */
    summary({ cache = false } = {}) {
        this.#initialize();
        this.#features();
        this.#cells();

        let assays = {};
        assays[main_experiment_name] = extract_assay_names(this.#se_handle);
        for (const [k, v] of Object.entries(this.#alt_handles)) {
            try {
                assays[k] = extract_assay_names(v);
            } catch (e) {
                console.warn("failed to extract features for alternative Experiment '" + k + "'; " + e.message);
            }
        }

        let output = {
            modality_features: this.#raw_features,
            cells: this.#raw_cells,
            modality_assay_names: assays
        };

        if (!cache) {
            this.clear();
        }
        return output;
    }

    /**
     * @param {object} [options={}] - Optional parameters.
     * @param {boolean} [options.cache=false] - Whether to cache the results for re-use in subsequent calls to this method or {@linkcode SummarizedExperimentDataset#summary summary}.
     * If `true`, users should consider calling {@linkcode SummarizedExperimentDataset#clear clear} to release the memory once this dataset instance is no longer needed.
     *
     * @return {object} Object containing the per-feature and per-cell annotations.
     * This has the following properties:
     *
     * - `features`: an object where each key is a modality name and each value is a {@linkplain external:DataFrame DataFrame} of per-feature annotations for that modality.
     * - `cells`: a {@linkplain external:DataFrame DataFrame} containing per-cell annotations.
     * - `matrix`: a {@linkplain external:MultiMatrix MultiMatrix} containing one {@linkplain external:ScranMatrix ScranMatrix} per modality.
     * - `primary_ids`: an object where each key is a modality name and each value is an integer array containing the feature identifiers for each row in that modality.
     */
    load({ cache = false } = {}) {
        this.#initialize();
        this.#features();
        this.#cells();

        let output = { 
            matrix: new scran.MultiMatrix,
            row_ids: {},
            features: {},
            cells: this.#raw_cells
        };

        let mapping = { 
            RNA: { exp: this.#rnaExperiment, assay: this.#rnaCountAssay },
            ADT: { exp: this.#adtExperiment, assay: this.#adtCountAssay },
            CRISPR: { exp: this.#crisprExperiment, assay: this.#crisprCountAssay }
        };

        try {
            for (const [k, v] of Object.entries(mapping)) {
                if (v.exp === null) {
                    continue;
                }

                let handle;
                let name = v.exp;
                if (typeof v.exp == "string") {
                    if (v.exp === "") {
                        handle = this.#se_handle;
                    } else {
                        if (!(v.exp in this.#alt_handles)) {
                            continue;
                        }
                        handle = this.#alt_handles[v.exp];
                    }
                } else {
                    if (v.exp >= this.#alt_handle_order.length) {
                        continue;
                    }
                    name = this.#alt_handle_order[v.exp];
                    handle = this.#alt_handles[name];
                }

                let loaded = extract_assay(handle, v.assay, true);
                output.matrix.add(k, loaded.matrix);
                let out_ids = loaded.row_ids;
                output.row_ids[k] = out_ids;
                output.features[k] = bioc.SLICE(this.#raw_features[name], out_ids);
            }

            let primaries = { 
                RNA: this.#primaryRnaFeatureIdColumn, 
                ADT: this.#primaryAdtFeatureIdColumn,
                CRISPR: this.#primaryCrisprFeatureIdColumn
            };
            output.primary_ids = futils.extractPrimaryIds(output.features, primaries);

        } catch (e) {
            scran.free(output.matrix);
            throw e;
        }

        if (!cache) {
            this.clear();
        }
        return output;
    }

    /**
     * @return {object} Object describing this dataset, containing:
     *
     * - `files`: Array of objects representing the files used in this dataset.
     *   Each object corresponds to a single file and contains:
     *   - `type`: a string denoting the type.
     *   - `file`: a {@linkplain SimpleFile} object representing the file contents.
     * - `options`: An object containing additional options to saved.
     */
    serialize() {
        return this.#dump_summary(f => f);
    }

    /**
     * @param {Array} files - Array of objects like that produced by {@linkcode SummarizedExperimentDataset#serialize serialize}.
     * @param {object} options - Object containing additional options to be passed to the constructor.
     * @return {SummarizedExperimentDataset} A new instance of this class.
     * @static
     */
    static async unserialize(files, options) {
        if (files.length != 1 || files[0].type != "rds") {
            throw new Error("expected exactly one file of type 'rds' for SummarizedExperiment unserialization");
        }
        return new SummarizedExperimentDataset(files[0].file, options);
    }
}

/***********************
 ******* Result ********
 ***********************/

/**
 * Pre-computed analysis results stored as a SummarizedExperiment object (or one of its subclasses) inside an RDS file.
 */
export class SummarizedExperimentResult {
    #rds_file;

    #rds_handle;
    #se_handle;
    #alt_handles;
    #alt_handle_order;

    #raw_features;
    #raw_cells;
    #rd_handles;
    #rd_handle_order;

    #primaryAssay;
    #isPrimaryNormalized;
    #reducedDimensionNames;

    /**
     * @param {SimpleFile|string|Uint8Array|File} rdsFile - Contents of a RDS file.
     * On browsers, this may be a File object.
     * On Node.js, this may also be a string containing a file path.
     * @param {object} [options={}] - Optional parameters.
     * @param {object|string|number} [options.primaryAssay=0] - See {@linkcode SummarizedExperimentResult#setPrimaryAssay setPrimaryAssay}.
     * @param {object|boolean} [options.isPrimaryNormalized={}] - See {@linkcode SummarizedExperimentResult#setIsPrimaryNormalized setIsPrimaryNormalized}.
     * @param {?Array} [options.reducedDimensionNames=null] - See {@linkcode SummarizedExperimentResult#setReducedDimensionNames setReducedDimensionNames}.
     */
    constructor(rdsFile, { 
        primaryAssay = 0,
        isPrimaryNormalized = true,
        reducedDimensionNames = null
    } = {}) {
        if (rdsFile instanceof afile.SimpleFile) {
            this.#rds_file = rdsFile;
        } else {
            this.#rds_file = new afile.SimpleFile(rdsFile);
        }

        // Cloning to avoid pass-by-reference links.
        this.#primaryAssay = bioc.CLONE(primaryAssay);
        this.#isPrimaryNormalized = bioc.CLONE(isPrimaryNormalized);
        this.#reducedDimensionNames = bioc.CLONE(reducedDimensionNames);

        this.clear();
    }

    /**
     * @param {object|string|number} primary - Assay containing the relevant data for each modality.
     *
     * - If a string, this is used as the name of the assay across all modalities.
     * - If a number, this is used as the index of the assay across all modalities.
     * - If any object, the key should be the name of a modality and the value may be either a string or number specifying the assay to use for that modality.
     *   Modalities absent from this object will not be loaded.
     */
    setPrimaryAssay(primary) {
        this.#primaryAssay = bioc.CLONE(primary);
        return;
    }

    /**
     * @param {object|boolean} normalized - Whether or not the assay for a particular modality has already been normalized.
     *
     * - If a boolean, this is used to indicate normalization status of assays across all modalities.
     *   If `false`, that modality's assay is assumed to contain count data and is subjected to library size normalization. 
     * - If any object, the key should be the name of a modality and the value should be a boolean indicating whether that modality's assay has been normalized.
     *   Modalities absent from this object are assumed to have been normalized.
     */
    setIsPrimaryNormalized(normalized) {
        this.#isPrimaryNormalized = bioc.CLONE(normalized);
        return;
    }

    /**
     * @param {?Array} names - Array of names of the reduced dimensions to load.
     * If `null`, all reduced dimensions found in the file are loaded.
     */
    setReducedDimensionNames(names) {
        this.#reducedDimensionNames = bioc.CLONE(names);
        return;
    }

    /**
     * Destroy caches if present, releasing the associated memory.
     * This may be called at any time but only has an effect if `cache = true` in {@linkcode SummarizedExperimentDataset#load load} or {@linkcodeSummarizedExperimentDataset#summary summary}.
     */
    clear() {
        scran.free(this.#se_handle);

        if (typeof this.#alt_handles != 'undefined' && this.#alt_handles !== null) {
            for (const v of Object.values(this.#alt_handles)) {
                scran.free(v);
            }
        }

        if (typeof this.#rd_handles != 'undefined' && this.#rd_handles !== null) {
            for (const v of Object.values(this.#rd_handles)) {
                scran.free(v.handle);
            }
        }

        scran.free(this.#rds_handle);

        this.#se_handle = null;
        this.#alt_handles = null;
        this.#rds_handle = null;

        this.#raw_features = null;
        this.#raw_cells = null;
    }

    #initialize() {
        if (this.#rds_handle !== null) {
            return;
        }

        this.#rds_handle = scran.readRds(this.#rds_file.content());
        this.#se_handle = this.#rds_handle.value();
        try {
            check_for_se(this.#se_handle);

            {
                const { handles, order } = extract_alt_exps(this.#se_handle);
                this.#alt_handles = handles;
                this.#alt_handle_order = order;
            }

            {
                const { handles, order } = extract_reduced_dims(this.#se_handle);
                this.#rd_handles = handles;
                this.#rd_handle_order = order;
            }

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
        this.#raw_features = {};
        this.#raw_features[main_experiment_name] = extract_features(this.#se_handle);

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
        let chandle = this.#se_handle.attribute("colData");
        try {
            this.#raw_cells = load_data_frame(chandle);
        } catch(e) {
            throw new Error("failed to extract colData from a SummarizedExperiment; " + e.message);
        } finally {
            scran.free(chandle);
        }

        return;
    }

    /**
     * @param {object} [options={}] - Optional parameters.
     * @param {boolean} [options.cache=false] - Whether to cache the results for re-use in subsequent calls to this method or {@linkcode SummarizedExperimentDataset#load load}.
     * If `true`, users should consider calling {@linkcode SummarizedExperimentDataset#clear clear} to release the memory once this dataset instance is no longer needed.
     * 
     * @return {object} Object containing the per-feature and per-cell annotations.
     * This has the following properties:
     *
     * - `modality_features`: an object where each key is a modality name and each value is a {@linkplain external:DataFrame DataFrame} of per-feature annotations for that modality.
     * - `cells`: a {@linkplain external:DataFrame DataFrame} of per-cell annotations.
     * - `modality_assay_names`: an object where each key is a modality name and each value is an Array containing the names of available assays for that modality.
     *    Unnamed assays are represented as `null` names.
     * - `reduced_dimension_names`: an Array of strings containing names of dimensionality reduction results.
     */
    summary({ cache = false } = {}) {
        this.#initialize();
        this.#features();
        this.#cells();

        let assays = {};
        assays[main_experiment_name] = extract_assay_names(this.#se_handle);
        for (const [k, v] of Object.entries(this.#alt_handles)) {
            try {
                assays[k] = extract_assay_names(v);
            } catch (e) {
                console.warn("failed to extract features for alternative Experiment '" + k + "'; " + e.message);
            }
        }

        let output = {
            modality_features: this.#raw_features,
            cells: this.#raw_cells,
            modality_assay_names: assays,
            reduced_dimension_names: this.#rd_handle_order
        };

        if (!cache) {
            this.clear();
        }
        return output;
    }

    /**
     * @param {object} [options={}] - Optional parameters.
     * @param {boolean} [options.cache=false] - Whether to cache the results for re-use in subsequent calls to this method or {@linkcode SummarizedExperimentDataset#summary summary}.
     * If `true`, users should consider calling {@linkcode SummarizedExperimentDataset#clear clear} to release the memory once this dataset instance is no longer needed.
     *
     * @return {object} Object containing the per-feature and per-cell annotations.
     * This has the following properties:
     *
     * - `features`: an object where each key is a modality name and each value is a {@linkplain external:DataFrame DataFrame} of per-feature annotations for that modality.
     * - `cells`: a {@linkplain external:DataFrame DataFrame} containing per-cell annotations.
     * - `matrix`: a {@linkplain external:MultiMatrix MultiMatrix} containing one {@linkplain external:ScranMatrix ScranMatrix} per modality.
     * - `reduced_dimensions`: an object containing the dimensionality reduction results.
     *   Each value is an array of arrays, where each inner array contains the coordinates for one dimension.
     */
    load({ cache = false } = {}) {
        this.#initialize();
        this.#features();
        this.#cells();

        let output = { 
            matrix: new scran.MultiMatrix,
            features: {},
            cells: this.#raw_cells,
            reduced_dimensions: {}
        };

        // Fetch the reduced dimensions first.
        let reddims = this.#reducedDimensionNames;
        if (reddims == null) {
            reddims = this.#rd_handle_order;
        }

        for (const k of reddims) {
            let v = this.#rd_handles[k];
            let acquired = [];
            let dims = v.dimensions;
            let contents = v.handle.values();
            for (var d = 0; d < dims[1]; d++) {
                acquired.push(contents.slice(d * dims[0], (d + 1) * dims[0]));
            }
            output.reduced_dimensions[k] = acquired;
        }

        // Now fetching the assay matrix.
        try {
            for (const [k, v] of Object.entries(this.#raw_features)) {
                let curassay = this.#primaryAssay;
                if (typeof curassay == "object") {
                    if (k in curassay) {
                        curassay = curassay[k];
                    } else {
                        continue;
                    }
                }

                let curnormalized = this.#isPrimaryNormalized;
                if (typeof curnormalized == "object") {
                    if (k in curnormalized) {
                        curnormalized = curnormalized[k];
                    } else {
                        curnormalized = true;
                    }
                }

                let handle;
                if (k === "") {
                    handle = this.#se_handle;
                } else {
                    handle = this.#alt_handles[k];
                }

                let loaded = extract_assay(handle, curassay, !curnormalized);
                output.matrix.add(k, loaded.matrix);

                if (!curnormalized) {
                    let normed = scran.logNormCounts(loaded.matrix, { allowZeros: true });
                    output.matrix.add(k, normed);
                    output.features[k] = bioc.SLICE(this.#raw_features[k], loaded.row_ids);
                } else {
                    output.features[k] = this.#raw_features[k];
                }
            }

        } catch (e) {
            scran.free(output.matrix);
            throw e;
        }

        if (!cache) {
            this.clear();
        }
        return output;
    }
}
