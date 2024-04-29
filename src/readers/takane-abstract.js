import * as scran from "scran.js";
import * as bioc from "bioconductor";
import * as eutils from "./utils/extract.js";
import * as futils from "./utils/features.js";
import * as afile from "./abstract/file.js";

class TakaneNavigator {
    #content;
    #listing;
    #getter;
    #lister;

    constructor(getter, lister) {
        this.clear();
        this.#getter = getter;
        this.#lister = lister;
    }

    async get_json_file(path) {
        if (this.#content.has(path)) {
            return this.#content.get(path);
        } else {
            const contents = await this.#getter(path);
            const dec = new TextDecoder;
            const payload = dec.decode(contents);
            const output = JSON.parse(payload);
            this.#content.set(path, output);
            return output;
        }
    }

    async get_object_file(path) {
        return await get_json_file(path + "/OBJECT");
    }

    async list_files(path) {
        if (this.#listing.has(path)) {
            return this.#listing.get(path);
        } else {
            const listed = await this.#lister(path);
            this.#listing.set(path, listed);
            return listed;
        }
    }

    get(path) {
        return this.#getter(path);
    }

    list(path) {
        return this.#lister(path);
    }

    clear() {
        this.#content = new Map;
        this.#listing = new Map;
    }
};

/*************************
 *** Data frame loader ***
 *************************/

async function load_data_frame(df_path, navigator) {
    let colnames;
    let columns = [];
    let rownames = null;
    let nrows;

    let contents = await navigator.get(df_path + "/basic_columns.h5")
    let out = scran.realizeFile(contents);

    try {
        let handle = new scran.H5File(out.path);
        let ghandle = handle.open("data_frame");

        nrows = ghandle.readAttribute("row-count").values[0];
        colnames = ghandle.open("column_names", { load: true }).values;
        if ("row_names" in ghandle.children) {
            rownames = ghandle.open("row_names", { load: true }).values;
        }

        let chandle = ghandle.open("data");
        for (var i = 0; i < colnames.length; i++) {
            let iname = String(i);

            if (!(iname in chandle.children)) {
                // Try to load nested objects if they're DFs, but don't try too hard.
                let nested_path = df_path + "/other_columns/" + iname;
                try {
                    columns.push(await load_data_frame(nested_path, getter));
                } catch (e) {
                    console.warn("failed to extract nested DataFrame at '" + nested_path + "'; " + e.message);
                    console.warn(e)
                    columns.push(null);
                }
                continue;
            }

            let current;
            let objtype = chandle.children[iname]

            if (obtype == "DataSet") {
                let dhandle = chandle.open(iname, { load: true });
                let type = dhandle.readAttribute("type").values[0];

                if (type == "integer" || type == "string") {
                    current = dhandle.values;
                } else if (type == "number") {
                    current = dhandle.values;
                    if (!(current instanceof Float64Array) && !(current instanceof Float32Array)) {
                        current = new Float64Array(current);
                    }
                } else if (type == "boolean") {
                    current = new Array(dhandle.values.length);
                    for (const [i, x] of dhandle.values.entries()) {
                        current[i] = (x != 0);
                    }
                } else {
                    throw new Error("data frame column has unknown type '" + type + "'");
                }

                if ("missing-value-placeholder" in dhandle.attributes) {
                    let placeholder = dhandle.readAttribute("missing-value-placeholder").values[0];
                    current = Array.from(current);
                    for (var i = 0; i < current.length; i++) {
                        if (current[i] == placeholder) {
                            current[i] = null;
                        }
                    }
                }

            } else if (objtype == "Group") {
                let fhandle = chandle.open(iname);
                let type = fhandle.readAttribute("type").values[0];

                if (type == "factor") {
                    let levels = fhandle.open("levels", { load: true }).values;
                    let chandle = fhandle.open("codes", { load: true });
                    let codes = chandle.values;

                    let placeholder = -1;
                    if ("missing-value-placeholder" in chandle.attributes) {
                        placeholder = chandle.readAttribute("missing-value-placeholder").values[0];
                    }

                    current = new Array(codes.length);
                    for (const [i, x] of codes.entries()) {
                        if (x != placeholder) {
                            current[i] = levels[x];
                        } else {
                            current[i] = null;
                        }
                    }

                } else {
                    throw new Error("data frame column has unknown type '" + type + "'");
                }

            } else {
                throw new Error("data frame column is of an unknown HDF5 object type");
            }

            columns.push(current);
        }
    } finally {
        out.flush();
    }

    let new_columns = {};
    let new_colnames = [];
    for (var i = 0; i < columns.length; i++) {
        if (columns[i] != null) {
            new_columns[colnames[i]] = columns[i];
            new_colnames.push(colnames[i]);
        }
    }

    return new bioc.DataFrame(new_columns, { 
        columnOrder: new_colnames, 
        rowNames: rownames, 
        numberOfRows: nrows
    });
}

/********************
 *** Array loader ***
 ********************/

async function extract_array(arr_path, navigator, forceInteger) {
    const arrmeta = await navigator.get_object_file(arr_path);
    const arrtype = arrmeta.type;
    let output;

    if (arrtype == "dense_array") {
        const contents = await navigator.get(arr_path + "/array.h5");
        const stuff = scran.realizeFile(contents);
        try {
            output = scran.initializeSparseMatrixFromHdf5(stuff.path, name, { forceInteger });
        } finally {
            stuff.flush();
        }

    } else if (arrtype == "compressed_sparse_matrix") {
        const contents = await navigator.get(arr_path + "/matrix.h5");
        const stuff = scran.realizeFile(contents);
        try {
            output = scran.initializeSparseMatrixFromHdf5(stuff.path, name, { forceInteger });
        } finally {
            stuff.flush();
        }

    } else {
        throw new Error("assay type '" + arrtype + "' is currently not supported");
    }

    return output;
}

/*******************
 *** List loader ***
 *******************/

async function extract_list(list_path, navigator) {
    const list_meta = await navigator.get_object_file(list_path);

    let list_format = "hdf5";
    if ("format" in list_meta) {
        list_format = list_meta.format;
    }

    if (list_format == "json.gz") {
        const contents = await navigator.get(list_path + "/list.json.gz");
        const ofile = new afile.SimpleFile(contents, { name: "list.json.gz" });
        const unpacked = eutils.unpackText(ofile.buffer(), { compression: "gz" });
        const parsed = JSON.parse(unpacked);
        return extract_json_list(parsed);

    } else if (list_format == "hdf5") {
        const contents = await navigator.get(list_path + "/list.h5");
        const stuff = scran.realizeFile(contents);
        let loaded;
        try {
            const handle = new scran.H5File(stuff.path);
            const ghandle = handle.open("simple_list");
            loaded = extract_hdf5_list(ghandle);
        } finally {
            stuff.flush();
        }
        return loaded;

    } else {
        throw new Error("unknown simple_list format '" + list_format + "'");
    }
}

// I'm just going to assume all of the lists are of version >= 1.2 (JSON).
function extract_json_list(obj) {
    if (!("type" in obj)) {
        throw new Error("non-standard JSON object for JSON 'simple_list'");
    }

    if (obj.type == "list") {
        if ("names" in obj) {
            let output = {};
            for (var i = 0; i < obj.values.length; i++) {
                output[obj.names[i]] = extract_json_list(obj.values[i]);
            }
            return output;

        } else {
            let output = [];
            for (var i = 0; i < obj.values.length; i++) {
                output.push(extract_json_list(obj.values[i]));
            }
            return output;
        }

    } else if (obj.type == "nothing") {
        return null;

    } else if (obj.type == "factor") {
        return obj.values.map(i => { 
            if (i == null) {
                return null;
            } else {
                return obj.levels[i];
            }
        });

    } else if (obj.type = "number") {
        const output = obj.values.map(i => {
            if (typeof x != "string") {
                return x;
            } else if (x == "Inf") {
                return Number.POSITIVE_INFINITY;
            } else if (x == "-Inf") {
                return Number.NEGATIVE_INFINITY;
            } else {
                return Number.NaN;
            }
        });
        if (output.indexOf(null) == -1) {
            return new Float64Array(output);
        } else {
            return output;
        }

    } else if (obj.type == "integer") {
        if (obj.values.indexOf(null) == -1) {
            return new Int32Array(obj.values);
        } else {
            return obj.values;
        }

    } else if (obj.type == "string" || obj.type == "boolean") {
        return obj.values;

    } else {
        console.warn("JSON simple list containing type '" + obj.type + "' is not yet supported");
        return null;
    }
}

// I'm just going to assume all of the lists are of version >= 1.3 (HDF5).
function extract_hdf5_list(handle) {
    const objtype = handle.readAttribute("uzuki_object").values[0];

    if (objtype == "list") {
        const dhandle = handle.open("data");
        if ("names" in handle.children) {
            const list_names = handle.open("names", { load: true }).values;
            const output = {};
            for (const [i, n] of list_names.entries()) {
                output[n] = extract_hdf5_list(dhandle.open(String(i), { load: true }));
            }
            return output;

        } else {
            const children = dhandle.children.keys();
            const output = new Array(children.length);
            for (const i of children) {
                output[Number(i)] = extract_hdf5_list(dhandle.open(i, { load: true }));
            }
            return output;
        }

    } else if (objtype == "nothing") {
        return null;

    } else if (objtype == "factor") {
        const levels = handle.open("levels", { load: true }).values;
        const ihandle = handle.open("values", { load: true });
        const codes = ihandle.values;

        const output = new Array(codes.length);
        if ("missing-value-placeholder" in ihandle.attributes) {
            const placeholder = ihandle.readAttribute("missing-value-placeholder").values[0];
            for (const [i, x] of codes.entries()) {
                if (i == placeholder) {
                    output[i] = null;
                } else {
                    output[i] = levels[x];
                }
            }
        } else {
            for (const [i, x] of codes.entries()) {
                output[i] = levels[x];
            }
        }

        return output;

    } else if (objtype = "string") {
        const vhandle = handle.open("values", { load: true });
        if ("missing-value-placeholder" in vhandle.attributes) {
            const placeholder = vhandle.readAttribute("missing-value-placeholder").values[0];
            return vhandle.values.map(x => {
                if (x == placeholder) {
                    return null;
                } else {
                    return x;
                }
            });
        } else {
            return vhandle.values;
        }

    } else if (objtype == "boolean") {
        const vhandle = handle.open("values", { load: true });
        const values = vhandle.values;

        const output = new Array(values.length);
        if ("missing-value-placeholder" in vhandle.attributes) {
            const placeholder = vhandle.readAttribute("missing-value-placeholder").values[0];
            for (const [i, x] of values.entries()) {
                if (i == placeholder) {
                    output[i] = null;
                } else {
                    output[i] = (x != 0);
                }
            }
        } else {
            for (const [i, x] of values.entries()) {
                output[i] = (x != 0);
            }
        }

        return output;

    } else if (objtype = "integer" || objtype == "number") {
        const vhandle = handle.open("values", { load: true });
        let output = vhandle.values;

        if ("missing-value-placeholder" in vhandle.attributes) {
            const placeholder = vhandle.readAttribute("missing-value-placeholder").values[0];
            output = Array.from(output);
            for (const [i, x] of output.entries()) {
                if (x == placeholder) {
                    output[i] = null;
                }
            }
        } else if (objtype == "number") {
            if (!(output instanceof Float64Array) && !(output instanceof Float32Array)) {
                output = new Float64Array(output);
            }
        }

        return output;

    } else {
        console.warn("JSON simple list containing type '" + objtype + "' is not yet supported");
        return null;
    }
}

/********************
 *** S(C)E loader ***
 ********************/

async function load_sce_components(se_path, navigator) {
    const extract_features = (obj_info, listing, se_path) => {
        if (!("summarized_experiment" in obj_info)) {
            throw new Error("expected a 'summarized_experiment' object");
        }
        if (listing.indexOf("row_data") != -1) {
            return load_data_frame(se_path + "/row_data", navigator);
        } else {
            const nrows = obj_info.summarized_experiment.dimensions[0];
            return new bioc.DataFrame({}, { numberOfRows: nrows });
        }
    };

    const extract_assay_names = async (listing, path) => {
        if (listing.indexOf("assays") == -1) {
            return [];
        }

        // Only reporting the names of assays with supported types.
        const assnames = await navigator.get_json_file(path + "/assays/names.json");
        const collected = [];
        for (const [i, a] of assnames) {
            const assmeta = await navigator.get_object_file(path + "/assays/" + String(i));
            const asstype = assmeta.type;
            if (asstype != "dense_array" && asstype != "compressed_sparse_matrix" && asstype != "delayed_array") {
                collected.push(a);
            }
        }
        return collected;
    };

    const obj_info = await navigator.get_object_file(se_path);
    const listing  = await navigator.list_files(se_path);
    const main = {
        "name": "",
        "path": se_path,
        "features": await extract_features(obj_info, listing, se_path),
        "assays": await extract_assay_names(listing, se_path),
        "cells": null,
    };

    const is_sce = ("single_cell_experiment" in obj_info);
    if (is_sce) {
        const sce_meta = obj_info.single_cell_experiment;
        if ("main_experiment_name" in sce_meta) {
            main.name = sce_meta.main_experiment_name;
        }
    }

    if (listing.indexOf("column_data") != -1) {
        let col_path = this.#path + "/column_data";
        main.cells = await load_data_frame(col_path, navigator);
    } else {
        main.cells = new bioc.DataFrame({}, { numberOfRows: obj_info.summarized_experiment.dimensions[1] });
    }

    let atlernatives = [];
    if (is_sce && listing.indexOf("alternative_experiments") != -1) {
        alternatives = await navigator.get_json_file(se_path + "/alternative_experiments/names.json");

        for (const [i, alt] of alternatives.entries()) {
            let alt_path = se_path + "/alternative_experiments/" + String(i);
            try {
                const alt_obj_info = await navigator.get_object_file(alt_path);
                const alt_listing = await navigator.list_files(alt_path);
                altnames[i] = {
                    "name": alt,
                    "path": alt_path,
                    "features": await extract_features(alt_obj_info, alt_listing, alt_path),
                    "assays": await extract_assay_names(alt_listing, alt_path)
                };
            } catch (e) {
                console.warn("failed to extract features for alternative experiment at '" + alt_path + "'; " + e.message);
            }
        }
    }

    return { main: main, alternative: alternatives };
}

async function extract_assay(path, assay, navigator, forceInteger) {
    const assnames = await navigator.get_json_file(path + "/assays/names.json");

    if (typeof assay == "string") {
        const counter = assnames.indexOf(assay);
        if (counter == -1) {
            throw new Error("assay '" + assay + "' not found");
        }
        assay = counter;
    } else {
        if (assay >= assnames.length) {
            throw new Error("assay " + String(assay) + " out of range");
        }
    }

    return extract_array(path + "/assays/" + String(assay), navigator, forceInteger);
}

/************************
 ******* Dataset ********
 ************************/

/**
 * Dataset stored as a SummarizedExperiment (or one of its subclasses) in the [**takane** format](https://github.com/ArtifactDB/takane).
 * This is intended as a virtual base class so that applications can define their own subclasses with the appropriate getter and listing methods. 
 * Subclasses should define `abbreviate()` and `serialize()` methods, as well as the static `format()` and `unserialize()` methods - 
 * see the [Dataset contract](https://github.com/kanaverse/bakana/blob/master/docs/related/custom_datasets.md) for more details.
 */
export class AbstractTakaneDataset {
    #path;
    #navigator;
    #raw_components;
    #options;

    /**
     * @param {string} path - Some kind of the path to the SummarizedExperiment.
     * The exact interpretation of this argument is left to subclasses.
     * @param {function} getter - A (possibly `async`) function that accepts a string containing the relative path to the file of interest, and returns a Uint8Array of that file's contents.
     * Each path is created by adding unix-style file separators to `path`.
     * @param {function} lister - A (possibly `async`) function that accepts a string containing the relative path to the directory of interest, and returns an array of the contents of that directory (non-recursive).
     * Each path is created by adding unix-style file separators to `path`.
     */
    constructor(path, getter, lister) {
        this.#path = path;
        this.#navigator = new TakaneNavigator(getter, lister);
        this.#raw_components = null;
        this.#options = AbstractArtifactdbDataset.defaults();
        return;
    }

    /**
     * @return {object} Default options, see {@linkcode AbstractArtifactdbDataset#setOptions setOptions} for more details.
     */
    static defaults() {
        return {
            rnaCountAssay: 0, 
            adtCountAssay: 0, 
            crisprCountAssay: 0,
            rnaExperiment: "", 
            adtExperiment: "Antibody Capture", 
            crisprExperiment: "CRISPR Guide Capture",
            primaryRnaFeatureIdColumn: null, 
            primaryAdtFeatureIdColumn: null,
            primaryCrisprFeatureIdColumn: null 
        };
    }

    /**
     * @return {object} Object containing all options used for loading.
     */
    options() {
        return { ...(this.#options) };
    }

    /**
     * @param {object} options - Optional parameters that affect {@linkcode AbstractArtifactdbDataset#load load} (but not {@linkcode AbstractArtifactdbDataset#summary summary}).
     * @param {string|number} [options.rnaCountAssay] - Name or index of the assay containing the RNA count matrix.
     * @param {string|number} [options.adtCountAssay] - Name or index of the assay containing the ADT count matrix.
     * @param {string|number} [options.crisprCountAssay] - Name or index of the assay containing the CRISPR count matrix.
     * @param {?(string|number)} [options.rnaExperiment] - Name or index of the alternative experiment containing gene expression data.
     *
     * If `i` is `null` or invalid (e.g., out of range index, unavailable name), it is ignored and no RNA data is assumed to be present.
     * If `i` is an empty string, the main experiment is assumed to contain the gene expression data.
     * @param {?(string|number)} [options.adtExperiment] - Name or index of the alternative experiment containing ADT data.
     *
     * If `i` is `null` or invalid (e.g., out of range index, unavailable name), it is ignored and no ADTs are assumed to be present.
     * If `i` is an empty string, the main experiment is assumed to contain the ADT data.
     * @param {?(string|number)} [options.crisprExperiment] - Name or index of the alternative experiment containing CRISPR guide data.
     *
     * If `i` is `null` or invalid (e.g., out of range index, unavailable name), it is ignored and no CRISPR guides are assumed to be present.
     * If `i` is an empty string, the main experiment is assumed to contain the guide data.
     * @param {?(string|number)} [options.primaryRnaFeatureIdColumn] - Name or index of the column of the `features` {@linkplain external:DataFrame DataFrame} that contains the primary feature identifier for gene expression.
     *
     * If `i` is `null` or invalid (e.g., out of range index, unavailable name), it is ignored and the primary identifier is defined as the existing row names.
     * However, if no row names are present in the SummarizedExperiment, no primary identifier is defined.
     * @param {?(string|number)} [options.primaryAdtFeatureIdColumn] - Name or index of the column of the `features` {@linkplain external:DataFrame DataFrame} that contains the primary feature identifier for the ADTs.
     *
     * If `i` is `null` or invalid (e.g., out of range index, unavailable name), it is ignored and the primary identifier is defined as the existing row names.
     * However, if no row names are present in the SummarizedExperiment, no primary identifier is defined.
     * @param {?(string|number)} [options.primaryCrisprFeatureIdColumn] - Name or index of the column of the `features` {@linkplain external:DataFrame DataFrame} that contains the primary feature identifier for the CRISPR guides.
     *
     * If `i` is `null` or invalid (e.g., out of range index, unavailable name), it is ignored and the existing row names (if they exist) are used as the primary identifier.
     * However, if no row names are present in the SummarizedExperiment, no primary identifier is defined.
     */
    setOptions(options) {
        for (const [k, v] of Object.entries(options)) {
            this.#options[k] = v;
        }
    }

    /**
     * Destroy caches if present, releasing the associated memory.
     * This may be called at any time but only has an effect if `cache = true` in {@linkcode AbstractArtifactdbDataset#load load} or {@linkcode AbstractArtifactdbDataset#summary summary}.
     */
    clear() {
        this.#navigator.clear();
        this.#raw_components = null;
    }

    async #load_components() {
        if (this.#raw_components == null) {
            this.#raw_components = await load_sce_components(this.#path, this.#navigator);
        }
        return this.#raw_components;
    }

    /**
     * @param {object} [options={}] - Optional parameters.
     * @param {boolean} [options.cache=false] - Whether to cache the intermediate results for re-use in subsequent calls to any methods with a `cache` option.
     * If `true`, users should consider calling {@linkcode AbstractArtifactdbDataset#clear clear} to release the memory once this dataset instance is no longer needed.
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
        const comp = await this.#load_components();

        const features = {};
        const assays = {};
        features[comp.main.name] = comp.main.features;
        assays[comp.main.name] = comp.main.assays;
        for (const [a, v] of Object.entries(comp.alternative)) {
            features[a] = v.features;
            assays[a] = v.assays;
        }

        const output = {
            modality_features: features,
            cells: comp.cells,
            modality_assay_names: assays,
        };

        if (!cache) {
            this.clear();
        }
        return output;
    }

    #primary_mapping() {
        return {
            RNA: this.#options.primaryRnaFeatureIdColumn, 
            ADT: this.#options.primaryAdtFeatureIdColumn,
            CRISPR: this.#options.primaryCrisprFeatureIdColumn
        };
    }

    /**
     * @param {object} [options={}] - Optional parameters.
     * @param {boolean} [options.cache=false] - Whether to cache the intermediate results for re-use in subsequent calls to any methods with a `cache` option.
     * If `true`, users should consider calling {@linkcode AbstractArtifactdbDataset#clear clear} to release the memory once this dataset instance is no longer needed.
     *
     * @return {object} An object where each key is a modality name and each value is an array (usually of strings) containing the primary feature identifiers for each row in that modality.
     * The contents are the same as the `primary_ids` returned by {@linkcode AbstractArtifactdbDataset#load load} but the order of values may be different.
     *
     * @async
     */
    async previewPrimaryIds({ cache = false } = {}) {
        const comp = await this.#load_components();

        let fmapping = {
            RNA: this.#options.rnaExperiment, 
            ADT: this.#options.adtExperiment, 
            CRISPR: this.#options.crisprExperiment 
        };

        const features = {};
        features[comp.main.name] = comp.main.features;
        for (const [a, v] of Object.entries(comp.alternative)) {
            features[a] = v.features;
        }
        const preview = futils.extractRemappedPrimaryIds(features, fmapping, this.#primary_mapping());

        if (!cache) {
            this.clear();
        }
        return preview;
    }

    /**
     * @param {object} [options={}] - Optional parameters.
     * @param {boolean} [options.cache=false] - Whether to cache the intermediate results for re-use in subsequent calls to any methods with a `cache` option.
     * If `true`, users should consider calling {@linkcode AbstractArtifactdbDataset#clear clear} to release the memory once this dataset instance is no longer needed.
     *
     * @return {object} Object containing the per-feature and per-cell annotations.
     * This has the following properties:
     *
     * - `features`: an object where each key is a modality name and each value is a {@linkplain external:DataFrame DataFrame} of per-feature annotations for that modality.
     * - `cells`: a {@linkplain external:DataFrame DataFrame} containing per-cell annotations.
     * - `matrix`: a {@linkplain external:MultiMatrix MultiMatrix} containing one {@linkplain external:ScranMatrix ScranMatrix} per modality.
     * - `primary_ids`: an object where each key is a modality name and each value is an array (usually of strings) containing the primary feature identifiers for each row in that modality.
     *
     * Modality names are guaranteed to be one of `"RNA"`, `"ADT"` or `"CRISPR"`.
     * We assume that the instance already contains an appropriate mapping from the observed feature types to each expected modality,
     * either from the {@linkcode AbstractArtifactdbDataset#defaults defaults} or with {@linkcode AbstractArtifactdbDataset#setOptions setOptions}.
     *
     * @async
     */
    async load({ cache = false } = {}) {
        const comp = await this.#load_components();

        let output = { 
            matrix: new scran.MultiMatrix,
            features: {},
            cells: comp.main.cells,
        };

        let mapping = { 
            RNA: { exp: this.#options.rnaExperiment, assay: this.#options.rnaCountAssay },
            ADT: { exp: this.#options.adtExperiment, assay: this.#options.adtCountAssay },
            CRISPR: { exp: this.#options.crisprExperiment, assay: this.#options.crisprCountAssay }
        };

        const altmap = {};
        for (const alt of comp.alternative) {
            altmap[alt.name] = alt;
        }

        try {
            for (const [k, v] of Object.entries(mapping)) {
                if (v.exp === null) {
                    continue;
                }

                let info;
                if (typeof v.exp == "string") {
                    if (v.exp === comp.main.name) {
                        info = comp.main;
                    } else if (v.exp in altmap) {
                        info = altmap[name];
                    } else {
                        continue;
                    }
                } else {
                    if (v.exp < comp.alternative.length) {
                        info = comp.alternative[v.exp];
                    } else {
                        continue;
                    }
                }

                let loaded = await extract_assay(info.path, v.assay, this.#navigator, true);
                output.matrix.add(k, loaded);
                output.features[k] = info.features;
            }

            output.primary_ids = futils.extractPrimaryIds(output.features, this.#primary_mapping());

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

/***********************
 ******* Result ********
 ***********************/

/**
 * Pre-computed analysis results stored as a SummarizedExperiment object (or one of its subclasses) in the * [**takane** format](https://github.com/ArtifactDB/takane).
 * This is intended as a virtual base class; applications should define their own subclasses with appropriate getter and listing methods. 
 */
export class AbstractTakeneResult {
    #path;
    #navigator;
    #raw_components;
    #raw_other;
    #raw_
    #options;

    /**
     * @param {string} path - Some kind of the path to the SummarizedExperiment.
     * The exact interpretation of this argument is left to subclasses.
     * @param {function} getter - A (possibly `async`) function that accepts a string containing the relative path to the file of interest, and returns a Uint8Array of that file's contents.
     * Each path is created by adding unix-style file separators to `path`.
     * @param {function} lister - A (possibly `async`) function that accepts a string containing the relative path to the directory of interest, and returns an array of the contents of that directory (non-recursive).
     * Each path is created by adding unix-style file separators to `path`.
     */
    constructor(path, getter, lister) {
        this.#path = path;
        this.#navigator = new TakaneNavigator(getter, lister);
        this.#raw_components = null;
        this.#raw_other = null;
        this.#options = AbstractTakeneResult.defaults();
    }

    /**
     * @return {object} Default options, see {@linkcode AbstractTakeneResult#setOptions setOptions} for more details.
     */
    static defaults() {
        return { 
            primaryAssay: 0,
            isPrimaryNormalized: true,
            reducedDimensionNames: null
        };
    }

    /**
     * @return {object} Object containing all options used for loading.
     */
    options() {
        return { ...(this.#options) };
    }

    /**
     * @param {object} options - Optional parameters that affect {@linkcode AbstractTakeneResult#load load} (but not {@linkcode AbstractTakeneResult#summary summary}.
     * @param {object|string|number} [options.primaryAssay] - Assay containing the relevant data for each modality.
     *
     * - If a string, this is used as the name of the assay across all modalities.
     * - If a number, this is used as the index of the assay across all modalities.
     * - If any object, the key should be the name of a modality and the value may be either a string or number specifying the assay to use for that modality.
     *   Modalities absent from this object will not be loaded.
     * @param {object|boolean} [options.isPrimaryNormalized] - Whether or not the assay for a particular modality has already been normalized.
     *
     * - If a boolean, this is used to indicate normalization status of assays across all modalities.
     *   If `false`, that modality's assay is assumed to contain count data and is subjected to library size normalization. 
     * - If any object, the key should be the name of a modality and the value should be a boolean indicating whether that modality's assay has been normalized.
     *   Modalities absent from this object are assumed to have been normalized.
     * @param {?Array} [options.reducedDimensionNames] - Array of names of the reduced dimensions to load.
     * If `null`, all reduced dimensions found in the file are loaded.
     */
    setOptions(options) {
        // Cloning to avoid pass-by-reference links.
        for (const [k, v] of Object.entries(options)) {
            this.#options[k] = bioc.CLONE(v);
        }
    }

    /**
     * Destroy caches if present, releasing the associated memory.
     * This may be called at any time but only has an effect if `cache = true` in {@linkcode AbstractTakeneResult#load load} or {@linkcode AbstractTakeneResult#summary summary}.
     */
    clear() {
        this.#raw_components = null;
        this.#raw_other = null;
        this.#navigator.clear();
    }

    async #load_components() {
        if (this.#raw_components === null) {
            this.#raw_components = { "core": await load_sce_components(this.#path, this.#navigator) };

            // Loading the other metadata.
            let listing = await this.#navigator.list_files(this.#path);
            if (listing.indexOf("other_data") != -1) {
                let other_path = this.#path + "/other_data";
                this.#raw_components.other = await extract_list(other_path, this.#navigator);
            } else {
                this.#raw_components.other = {};
            }

            // Only supporting dense_arrays in the reduced dimensions for now.
            const full_meta = await this.#navigator.get_object_file(this.#path);
            let rd_names = [];
            if ("single_cell_experiment" in full_meta) {
                const all_rd_names = await this.#navigator.get_json_file(this.#path + "/reduced_dimensions/names.json");
                for (const [i, r] of all_rd_names.entries()) {
                    const red_path = this.#path + "/reduced_dimensions/" + String(i);
                    let red_meta = await this.#navigator.get_object_file(red_path);
                    if (red_meta["type"] == "dense_array") {
                        rd_names.push({ "name": r, "path": red_path });
                    }
                }
            }
            this.#raw_components.reduced_dimensions = rd_names;
        }

        return this.#raw_components;
    }

    /**
     * @param {object} [options={}] - Optional parameters.
     * @param {boolean} [options.cache=false] - Whether to cache the results for re-use in subsequent calls to this method or {@linkcode AbstractTakeneResult#load load}.
     * If `true`, users should consider calling {@linkcode AbstractTakeneResult#clear clear} to release the memory once this dataset instance is no longer needed.
     * 
     * @return {object} Object containing the per-feature and per-cell annotations.
     * This has the following properties:
     *
     * - `modality_features`: an object where each key is a modality name and each value is a {@linkplain external:DataFrame DataFrame} of per-feature annotations for that modality.
     * - `cells`: a {@linkplain external:DataFrame DataFrame} of per-cell annotations.
     * - `modality_assay_names`: an object where each key is a modality name and each value is an Array containing the names of available assays for that modality.
     *    Unnamed assays are represented as `null` names.
     * - `reduced_dimension_names`: an Array of strings containing names of dimensionality reduction results.
     * - `other_metadata`: an object containing other metadata.
     *
     * @async 
     */
    async summary({ cache = false } = {}) {
        const comp = await this.#load_components();

        const features = {};
        const assays = {};
        features[comp.core.main.name] = comp.core.main.features;
        assays[comp.core.main.name] = comp.core.main.assays;
        for (const [a, v] of Object.entries(comp.core.alternative)) {
            features[a] = v.features;
            assays[a] = v.assays;
        }

        let output = {
            modality_features: features,
            cells: comp.core.main.cells,
            modality_assay_names: assays,
            reduced_dimension_names: comp.reduced_dimensions.map(x => x.name),
            other_metadata: comp.other,
        };

        if (!cache) {
            this.clear();
        }
        return output;
    }

    /**
     * @param {object} [options={}] - Optional parameters.
     * @param {boolean} [options.cache=false] - Whether to cache the results for re-use in subsequent calls to this method or {@linkcode AbstractTakeneResult#summary summary}.
     * If `true`, users should consider calling {@linkcode AbstractTakeneResult#clear clear} to release the memory once this dataset instance is no longer needed.
     *
     * @return {object} Object containing the per-feature and per-cell annotations.
     * This has the following properties:
     *
     * - `features`: an object where each key is a modality name and each value is a {@linkplain external:DataFrame DataFrame} of per-feature annotations for that modality.
     * - `cells`: a {@linkplain external:DataFrame DataFrame} containing per-cell annotations.
     * - `matrix`: a {@linkplain external:MultiMatrix MultiMatrix} containing one {@linkplain external:ScranMatrix ScranMatrix} per modality.
     * - `reduced_dimensions`: an object containing the dimensionality reduction results.
     *   Each value is an array of arrays, where each inner array contains the coordinates for one dimension.
     * - `other_metadata`: an object containing other metadata.
     *
     * @async
     */
    async load({ cache = false } = {}) {
        const comp = await this.#load_components();

        let output = { 
            matrix: new scran.MultiMatrix,
            features: {},
            cells: this.#raw_cells,
            reduced_dimensions: {},
            other_metadata: this.#raw_other
        };

        // Fetch the reduced dimensions first.
        {
            let reddims = this.#options.reducedDimensionNames;
            if (reddims == null) {
                reddims = comp.reduced_dimensions.entries(x => x.name);
            }

            if (reddims.length > 0) {
                let redmap = {};
                for (const [i, red] of comp.reduced_dimensions.entries()) {
                    redmap[red.name] = red.path;
                }

                for (const k of reddims) {
                    if (!(k in redmap)) {
                        continue;
                    }
                    const red_path = redmap[k];

                    let red_meta = await this.#navigator.get_object_file(red_path);
                    if (red_meta["type"] != "dense_array") {
                        console.warn("reduced dimensions of type '" + red_meta["type"] + "' are not yet supported");
                        continue;
                    }

                    let red_contents = await this.#navigator.get(red_path + "/array.h5");
                    let realized = scran.realizeFile(red_contents);
                    let acquired = [];
                    try {
                        let fhandle = new scran.H5File(realized.path);
                        let dhandle = fhandle.open("dense_array", { load: true });
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
            const altmap = {};
            for (const alt of comp.alternative) {
                altmap[alt.name] = alt;
            }

            try {
                for (const [k, v] of Object.entries(this.#raw_features)) {
                    let curassay = this.#options.primaryAssay;
                    if (typeof curassay == "object") {
                        if (k in curassay) {
                            curassay = curassay[k];
                        } else {
                            continue;
                        }
                    }

                    let curnormalized = this.#options.isPrimaryNormalized;
                    if (typeof curnormalized == "object") {
                        if (k in curnormalized) {
                            curnormalized = curnormalized[k];
                        } else {
                            curnormalized = true;
                        }
                    }

                    let info;
                    if (k === comp.core.main.name) {
                        info = comp.core.main;
                    } else if (k in altmap) {
                        info = altmap[k];
                    }

                    let loaded = await extract_assay(info.path, curassay, this.#navigator, !curnormalized);
                    output.matrix.add(k, loaded);

                    if (!curnormalized) {
                        let normed = scran.logNormCounts(loaded, { allowZeros: true });
                        output.matrix.add(k, normed);
                    }

                    output.features[k] = info.features;
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
