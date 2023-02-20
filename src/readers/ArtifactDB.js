import * as scran from "scran.js";
import * as bioc from "bioconductor";
import * as eutils from "./utils/extract.js";
import * as futils from "./utils/features.js";
import * as afile from "./abstract/file.js";

/**************************
 ******* Internals ********
 **************************/

async function load_data_frame(info, navigator) {
    if (typeof info == "string") {
        info = await navigator.metadata(info);
    }
    let contents = await navigator.file(info.path);

    let colnames;
    let columns;
    let rownames = null;

    if (info["$schema"].startsWith("csv_data_frame/")) {
        // TODO: replace with comservatory parser.
        let parsed = await eutils.readTable2(contents, { delim: "," });
        let colnames = parsed.shift();

        columns = new Array(colnames.length);
        for (var i = 0; i < columns.length; i++) {
            let current = [];
            for (const x of parsed) {
                current.push(x[i]);
            }
            columns[i] = current;
        }

        if (info.data_frame.row_names) {
            rownames = columns.shift();
            colnames.shift();
        }

        // Mutating the type... this doesn't quite handle NAs and NaNs properly, but whatever.
        for (var i = 0; i < columns.length; i++) {
            let type = info.data_frame.columns[i].type;
            if (type == "integer") {
                columns[i] = new Int32Array(columns[i]);
            } else if (type == "number") {
                columns[i] = new Float64Array(columns[i]);
            } else if (type == "boolean") {
                columns[i] = columns[i].map(x => x == "true");
            }
        }

    } else if (info["$schema"].startsWith("hdf5_data_frame/")) {
        let out = scran.realizeFile(contents);
        try {
            let handle = new scran.H5File(out.path);
            let ghandle = handle.open(info.hdf5_data_frame.group);

            colnames = ghandle.open("column_names", { load: true }).values;
            if (info.data_frame.row_names) {
                rownames = ghandle.open("row_names", { load: true }).values;
            }

            columns = [];
            let chandle = ghandle.open("data");
            for (var i = 0; i < colnames.length; i++) {
                if (!(String(i) in chandle.children)) {
                    columns.push(null);
                    continue;
                }

                let dhandle = chandle.open(String(i), { load: true });
                let current = dhandle.values;

                let type = info.data_frame.columns[i];
                if (type == "integer") {
                    if (current instanceof Float64Array || current instanceof Float32Array) {
                        current = new Int32Array(current);
                    }

                } else if (type == "number") {
                    if (!(current instanceof Float64Array) && !(current instanceof Float32Array)) {
                        current = new Float64Array(current);
                    }

                } else if (type == "boolean") {
                    let replacement = new Array(current.length);
                    for (var i = 0; i < current.length; i++) {
                        if (current[i] == -2147483648) {
                            replacement[i] = null;
                        } else {
                            replacement[i] = current[i] != 0
                        }
                    }
                    current = replacement;

                } else if (type == "string" || type == "date") {
                    if ("missing-value-placeholder" in dhandle.attributes) {
                        let placeholder = dhandle.readAttribute("missing-value-placeholder").values[0];
                        for (var i = 0; i < current.length; i++) {
                            if (current[i] == placeholder) {
                                current[i] = null;
                            }
                        }
                    }
                }

                columns.push(current);
            }
        } finally {
            out.flush();
        }

    } else {
        throw new Error("unknown data_frame schema type '" + info["$schema"] + "'");
    }

    // Ignoring placeholder columns for the time being.
    let new_columns = {};
    let new_colnames = [];
    for (var i = 0; i < columns.length; i++) {
        if (info.data_frame.columns[i].type === "other") {
//            let nest_meta = await navigator.metadata(info.data_frame.columns[i].resource.path);
//            try {
//                new_columns[colnames[i]] = await load_data_frame(nest_meta, navigator);
//                new_colnames.push(colnames[i]);
//            } catch (e) {
//                console.warn(e);
//            }
        } else {
            new_columns[colnames[i]] = columns[i];
            new_colnames.push(colnames[i]);
        }
    }

    return new bioc.DataFrame(new_columns, { 
        columnOrder: new_colnames, 
        rowNames: rownames, 
        numberOfRows: info.data_frame.dimensions[0] 
    });
}

const main_experiment_name = "";

async function extract_all_features(path, navigator) {
    let extract_features = async se_meta => {
        if ("row_data" in se_meta.summarized_experiment) {
            let row_path = se_meta.summarized_experiment.row_data.resource.path;
            return await load_data_frame(row_path, navigator);
        } else {
            return new bioc.DataFrame({}, { numberOfRows: se_meta.summarized_experiment.dimensions[0] });
        }
    };

    let full_meta = await navigator.metadata(path);
    let output = {};
    output[main_experiment_name] = await extract_features(full_meta);

    if ("single_cell_experiment" in full_meta) {
        for (const alt of full_meta.single_cell_experiment.alternative_experiments) {
            try {
                let alt_meta = await navigator.metadata(alt.resource.path);
                output[alt.name] = await extract_features(alt_meta);
            } catch (e) {
                console.warn("failed to extract features for alternative Experiment '" + alt.name + "'; " + e.message);
            }
        }
    }

    return output;
}

async function extract_all_assay_names(path, navigator) {
    let extract_assay_names = se_meta => {
        let output = [];
        for (const ass of se_meta.summarized_experiment.assays) {
            output.push(ass.name);
        }
        return output;
    };

    let full_meta = await navigator.metadata(path);
    let assays = {};
    assays[main_experiment_name] = extract_assay_names(full_meta);

    if ("single_cell_experiment" in full_meta) {
        for (const alt of full_meta.single_cell_experiment.alternative_experiments) {
            try {
                let alt_meta = await navigator.metadata(alt.resource.path);
                assays[alt.name] = extract_assay_names(alt_meta);
            } catch (e) {
                console.warn("failed to extract features for alternative Experiment '" + alt.name + "'; " + e.message);
            }
        }
    }

    return assays;
}

// This specifically loads the log-counts created by the dumper.
// TODO: replace this with chihaya.js.
function extract_logcounts(handle, navigator) {
    if (ghandle.attribute("delayed_type").values[0] !== "operation") {
        return null;
    }
    if (ghandle.attribute("delayed_operation").values[0] !== "unary arithmetic") {
        return null;
    }
    if (Math.abs(ghandle.open("value", { load: true }).values[0] - Math.log(2)) > 0.00000001) {
        return null;
    }
    if (ghandle.open("method", { load: true }).values[0] !== "/") {
        return null;
    }
    if (ghandle.open("seed", { load: true }).values[0] !== "right") {
        return null;
    }

    let ghandle2 = ghandle.open("seed");
    if (ghandle2.attribute("delayed_type").values[0] !== "operation") {
        return null;
    }
    if (ghandle2.attribute("delayed_operation").values[0] !== "unary math") {
        return null;
    }
    if (ghandle2.open("method", { load: true }).values[0] !== "log1p") {
        return null;
    }

    let ghandle3 = ghandle2.open("seed");
    if (ghandle3.attribute("delayed_type").values[0] !== "operation") {
        return null;
    }
    if (ghandle3.attribute("delayed_operation").values[0] !== "unary arithmetic") {
        return null;
    }
    if (ghandle3.open("method", { load: true }).values[0] !== "/") {
        return null;
    }
    if (ghandle3.open("seed", { load: true }).values[0] !== "right") {
        return null;
    }
    if (ghandle3.open("along", { load: true }).values[0] !== 1) {
        return null;
    }
    let sf = ghandle3.open("value", { load: true }).values;

    let ahandle = ghandle3.open("seed");
    if (ahandle.attribute("delayed_type").values[0] !== "array") {
        return null;
    }
    if (ahandle.attribute("delayed_operation").values[0] !== "custom alabaster local array") {
        return null;
    }
    let path = ahandle.open("path", { load: true }).values[0];

    let mat;
    let output = {};
    try {
        mat = extract_assay_raw(path, navigator, false); // don't force it to be integer, but we don't mind if it is.
        output.matrix = scran.logNormCounts(mat.matrix, { sizeFactors: sf, center: false });
        output.row_ids = mat.row_ids;
    } finally {
        scran.free(mat.matrix);
        realized.flush();
    }
                
    return output;
}

async function extract_assay(meta, assay, navigator, forceInteger) {
    if (typeof assay == "string") {
        var counter = 0;
        for (const ass of meta.summarized_experiment.assays) {
            if (ass.name == assay) {
                assay = counter;
                break;
            }
            counter++;
        }
        if (counter == meta.summarized_experiment.assays.length) {
            throw new Error("assay '" + assay + "' not found");
        }
    } else {
        if (assay >= meta.summarized_experiment.assays.length) {
            throw new Error("assay " + String(assay) + " out of range");
        }
    }

    let asspath = meta.summarized_experiment.assays[assay].resource.path;
    return extract_assay_raw(asspath, navigator, forceInteger);
}

async function extract_assay_raw(asspath, navigator, forceInteger) {
    let assmeta = await navigator.metadata(asspath);
    let contents = await navigator.file(assmeta.path);
    let output;

    let schema = assmeta["$schema"];
    let is_dense = schema.startsWith("hdf5_dense_array/");
    let is_sparse = schema.startsWith("hdf5_sparse_matrix/");

    if (is_dense || is_sparse) {
        let name = (is_sparse ?  assmeta.hdf5_sparse_matrix.group : assmeta.hdf5_dense_array.dataset);
        let stuff = scran.realizeFile(contents);
        try {
            output = scran.initializeSparseMatrixFromHDF5(stuff.path, name, { forceInteger });
        } finally {
            stuff.flush();
        }

    } else if (assmeta["$schema"].startsWith("hdf5_delayed_array/")) {
        let stuff = scran.realizeFile(contents);
        try {
            let fhandle = new scran.H5File(stuff.path);
            let ghandle = fhandle.open(assmeta.hdf5_delayed_array.group);

            // TODO: replace with calls to chihaya.js.
            output = extract_logcounts(ghandle, navigator);
            if (output == null) {
                throw new Error("currently only supporting bakana-generated log-counts for delayed arrays");
            }
        } finally {
            stuff.flush();
        }

    } else {
        throw new Error("array schema '" + assmeta["$schema"] + "' is currently not supported");
    }

    return output;
}

/************************
 ******* Dataset ********
 ************************/

/**
 * Dataset stored as a SummarizedExperiment in the **ArtifactDB** format.
 * This is intended as a virtual base class; applications should define subclasses that are tied to a specific {@linkplain external:ArtifactDbProjectNavigator ArtifactDbProjectNavigator} class.
 * Subclasses should define the static `format()` and `unserialize()` methods.
 */
export class ArtifactDbSummarizedExperimentDatasetBase {
    #path;
    #navigator;

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

    #dump_options() {
        return {
            rnaCountAssay: this.#rnaCountAssay,
            adtCountAssay: this.#adtCountAssay,
            rnaExperiment: this.#rnaExperiment,
            adtExperiment: this.#adtExperiment,
            primaryRnaFeatureIdColumn: this.#primaryRnaFeatureIdColumn,
            primaryAdtFeatureIdColumn: this.#primaryAdtFeatureIdColumn,
            primaryCrisprFeatureIdColumn: this.#primaryCrisprFeatureIdColumn
        };
    }

    /**
     * @param {string} path - Path to the SummarizedExperiment in the ArtifactDB project directory.
     * @param {extern:ArtifactDbProjectNavigator} navigator - A navigator object that describes how to obtain the various assets from the project directory containing `path`.
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
    constructor(path, navigator, { 
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
        this.#path = path;
        this.#navigator = navigator;

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
    setRnaCountAssay(name) {
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
        this.#raw_features = null;
        this.#raw_cells = null;
    }

    /**
     * @return {object} Object containing the abbreviated details of this dataset.
     */
    abbreviate() {
        return {
            path: this.#path,
            options: this.#dump_options()
        };
    }

    async #features() {
        if (this.#raw_features !== null) {
            return;
        }
        this.#raw_features = await extract_all_features(this.#path, this.#navigator);
        return;
    }

    async #cells() {
        if (this.#raw_cells !== null) {
            return;
        }
        let full_meta = await this.#navigator.metadata(this.#path);
        let col_path = full_meta.summarized_experiment.column_data.resource.path;
        this.#raw_cells = await load_data_frame(col_path, this.#navigator);
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
     *
     * @async
     */
    async summary({ cache = false } = {}) {
        await this.#features();
        await this.#cells();

        let output = {
            modality_features: this.#raw_features,
            cells: this.#raw_cells,
            modality_assay_names: await extract_all_assay_names(this.#path, this.#navigator)
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
     *
     * @async
     */
    async load({ cache = false } = {}) {
        await this.#features();
        await this.#cells();

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

        let full_meta = await this.#navigator.metadata(this.#path);
        let altmap = {};
        let alts = [];
        if ("single_cell_experiment" in full_meta) {
            alts = full_meta.single_cell_experiment.alternative_experiments;
            for (const alt of alts) {
                altmap[alt.name] = alt.resource.path;
            }
        }

        try {
            for (const [k, v] of Object.entries(mapping)) {
                if (v.exp === null) {
                    continue;
                }

                let meta = null;
                let name = v.exp;
                if (typeof v.exp == "string") {
                    if (v.exp === "") {
                        meta = full_meta;
                    } else {
                        if (!(v.exp in altmap)) {
                            continue;
                        }
                        meta = await this.#navigator(altmap[v.exp]);
                    }
                } else {
                    if (v.exp >= alts.length) {
                        continue;
                    }
                    name = alts[v.exp].name;
                    meta = await this.#navigator(alts[v.exp].resource.path);
                }

                let loaded = await extract_assay(meta, v.assay, this.#navigator, true);
                output.matrix.add(k, loaded.matrix);
                let out_ids = loaded.row_ids;
                output.row_ids[k] = out_ids;
                output.features[k] = bioc.SLICE(this.#raw_features[name], out_ids);
            }
            console.log(output.features);

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
        let enc = new TextEncoder;
        let contents = enc.encode(this.#path);
        let f = new afile.SimpleFile(contents, { name: this.#path });
        let existing = this.#navigator.serialize();
        return {
            files: [ { type: "path", file: f }, ...existing ],
            options: this.#dump_options()
        };
    }
}

/***********************
 ******* Result ********
 ***********************/

/**
 * Pre-computed analysis results stored as a SummarizedExperiment object (or one of its subclasses) in the **ArtifactDB** format.
 * This is intended as a virtual base class; applications should define subclasses that are tied to a specific {@linkplain external:ArtifactDbProjectNavigator ArtifactDbProjectNavigator} class.
 */
export class ArtifactDbSummarizedExperimentResultBase {
    #path;
    #navigator;

    #raw_features;
    #raw_cells;

    #primaryAssay;
    #isPrimaryNormalized;
    #reducedDimensionNames;

    /**
     * @param {string} path - Path to the SummarizedExperiment in the ArtifactDB project directory.
     * @param {extern:ArtifactDbProjectNavigator} navigator - A navigator object that describes how to obtain the various assets from the project directory containing `path`.
     * @param {object} [options={}] - Optional parameters.
     * @param {object|string|number} [options.primaryAssay=0] - See {@linkcode SummarizedExperimentResult#setPrimaryAssay setPrimaryAssay}.
     * @param {object|boolean} [options.isPrimaryNormalized={}] - See {@linkcode SummarizedExperimentResult#setIsPrimaryNormalized setIsPrimaryNormalized}.
     * @param {?Array} [options.reducedDimensionNames=null] - See {@linkcode SummarizedExperimentResult#setReducedDimensionNames setReducedDimensionNames}.
     */
    constructor(path, navigator, { 
        primaryAssay = 0,
        isPrimaryNormalized = true,
        reducedDimensionNames = null
    } = {}) {
        this.#path = path;
        this.#navigator = navigator;

        this.#primaryAssay = primaryAssay;
        this.#isPrimaryNormalized = isPrimaryNormalized;
        this.#reducedDimensionNames = reducedDimensionNames;

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
        this.#primaryAssay = primary;
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
        this.#isPrimaryNormalized = normalized;
        return;
    }

    /**
     * @param {?Array} names - Array of names of the reduced dimensions to load.
     * If `null`, all reduced dimensions found in the file are loaded.
     */
    setReducedDimensionNames(names) {
        this.#reducedDimensionNames = names;
        return;
    }

    /**
     * Destroy caches if present, releasing the associated memory.
     * This may be called at any time but only has an effect if `cache = true` in {@linkcode SummarizedExperimentDataset#load load} or {@linkcodeSummarizedExperimentDataset#summary summary}.
     */
    clear() {
        this.#raw_features = null;
        this.#raw_cells = null;
    }

    async #features() {
        if (this.#raw_features !== null) {
            return;
        }
        this.#raw_features = extract_all_features(this.#path, this.#navigator);
        return;
    }

    async #cells() {
        if (this.#raw_cells !== null) {
            return;
        }
        let full_meta = await this.#navigator.metadata(this.#path);
        let col_path = full_meta.summarized_experiment.col_data.resource.path;
        this.#raw_cells = await load_data_frame(col_path, this.#navigator);
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
     *
     * @async 
     */
    async summary({ cache = false } = {}) {
        await this.#features();
        await this.#cells();

        let output = {
            modality_features: this.#raw_features,
            cells: this.#raw_cells,
            modality_assay_names: await extract_all_assay_names(this.#path, this.#navigator),
            reduced_dimension_names: []
        };

        let full_meta = await this.#navigator.metadata(this.#path);
        if ("single_cell_experiment" in full_meta) {
            for (const red of full_meta.single_cell_experiment.reduced_dimensions) {
                let redmeta = await this.navigator.metadata(red.resource.path);
                if (redmeta["$schema"].startsWith("hdf5_dense_array/") && redmeta.array.dimensions.length == 2) {
                    output.reduced_dimension_names.push(red.name);
                }
            }
        }

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
     *
     * @async
     */
    async load({ cache = false } = {}) {
        await this.#features();
        await this.#cells();
        let full_meta = await this.#navigator.metadata(this.#path);

        let output = { 
            matrix: new scran.MultiMatrix,
            features: {},
            cells: this.#raw_cells,
            reduced_dimensions: {}
        };

        // Fetch the reduced dimensions first.
        {
            let reddims = this.#reducedDimensionNames;
            if (reddims == null) {
                reddims = [];
                if ("single_cell_experiment" in full_meta) {
                    for (const red of full_meta.single_cell_experiment.reduced_dimensions) {
                        reddims.push(red.name);
                    }
                }
            }

            if (reddims.length > 0) {
                let redmap = {};
                for (const red of full_meta.single_cell_experiment.reduced_dimensions) {
                    redmap[red.name] = red.path;
                }

                for (const k of reddims) {
                    let redmeta = await this.#navigator.metadata(redmap[k]); // this should be only HDF5 dense matrices.
                    let dims = redmeta.array.dimensions;
                    let redcontents = await this.#navigator.file(redmeta.path); 

                    let realized = scran.realizeFile(redcontents);
                    let acquired = [];
                    try {
                        let fhandle = new scran.H5File(realized.path);
                        let dhandle = fhandle.open(redmeta.hdf5_dense_array.dataset, { load: true });
                        let contents = dhandle.values;
                        for (var d = 0; d < dims[1]; d++) {
                            acquired.push(contents.slice(d * dims[0], (d + 1) * dims[0]));
                        }
                    } finally {
                        realized.flush();
                    }

                    output.reduced_dimensions[k] = acquired;
                }
            }
        }

        // Now fetching the assay matrix.
        {
            let altmap = {};
            if ("single_cell_experiment" in full_meta) {
                for (const alt of full_meta.single_cell_experiment.alternative_experiments) {
                    altmap[alt.name] = alt.resource.path;
                }
            }

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

                    let meta;
                    if (k === "") {
                        meta = full_meta;
                    } else {
                        meta = await this.#navigator.metadata(altmap[k]);
                    }

                    let loaded = extract_assay(meta, curassay, this.#navigator, !curnormalized);
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
        }

        if (!cache) {
            this.clear();
        }
        return output;
    }
}
