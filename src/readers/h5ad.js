import * as scran from "scran.js";
import * as bioc from "bioconductor";
import * as eutils from "./utils/extract.js";
import * as futils from "./utils/features.js";
import * as afile from "./abstract/file.js";

/**************************
 ******* Internals ********
 **************************/

function fetch_assay_details(handle, path) {
    let available = [];
    if ("X" in handle.children) {
        available.push("X");
    } 

    if ("layers" in handle.children) {
        let lhandle = handle.open("layers");
        if (!(lhandle instanceof scran.H5Group)) {
            throw new Error("expected a 'layers' group in a H5AD file");
        }

        for (const k of Object.keys(lhandle.children)) {
            available.push("layers/" + k);
        }
    }

    if (available.length == 0) {
        throw new Error("failed to find any assay in the H5AD file");
    }

    let deets = scran.extractHDF5MatrixDetails(path, available[0]);
    return {
        names: available,
        rows: deets.rows,
        columns: deets.columns
    };
}

function load_data_frame(handle) {
    let columns = {};

    for (const [key, val] of Object.entries(handle.children)) {
        if (val == "DataSet") {
            let dhandle = handle.open(key, { load: true });
            columns[key] = dhandle.values;
        } else if (val == "Group") {
            // Factor encoding for H5AD versions >= 0.8.0.
            let subhandle = handle.open(key);
            if ("categories" in subhandle.children && "codes" in subhandle.children) {
                let current_levels = eutils.extractHDF5Strings(subhandle, "categories");
                let codes = subhandle.open("codes", { load: true }).values;
                columns[key] = bioc.SLICE(current_levels, codes);
            }
        }
    }

    // Factor encoding for H5AD versions < 0.8.0.
    if ("__categories" in handle.children && handle.children["__categories"] == "Group") {
        let chandle = handle.open("__categories");

        for (const [key, val] of Object.entries(chandle.children)) {
            if (key in columns) {
                let current_levels = eutils.extractHDF5Strings(chandle, key);
                columns[key] = bioc.SLICE(current_levels, columns[key]);
            }
        }
    }

    if (Object.keys(columns).length == 0) {
        return null;
    } else {
        let rn = null;
        if ("_index" in columns) {
            rn = columns._index;
            delete columns._index;
        }
        return new bioc.DataFrame(columns, { rowNames: rn });
    }
}

function fetch_features(handle) {
    if ("var" in handle.children && handle.children["var"] == "Group") {
        let vhandle = handle.open("var");
        return load_data_frame(vhandle);
    }
    return null;
}

function fetch_cells(handle) {
    if ("obs" in handle.children && handle.children["obs"] == "Group") {
        let ohandle = handle.open("obs");
        return load_data_frame(ohandle);
    }
    return null;
}

/************************
 ******* Dataset ********
 ************************/

/**
 * Dataset in the H5AD format.
 */
export class H5adDataset {
    #h5_file;
    #h5_path;
    #h5_flush;
    #h5_handle;

    #raw_features;
    #raw_cells;
    #assay_details;

    #options;

    #dump_summary(fun) {
        let files = [{ type: "h5", file: fun(this.#h5_file) }]; 
        let options = this.options();
        return { files, options };
    }

    /**
     * @param {SimpleFile|string|Uint8Array|File} h5File - Contents of a H5AD file.
     * On browsers, this may be a File object.
     * On Node.js, this may also be a string containing a file path.
     */
    constructor(h5File) {
        if (h5File instanceof afile.SimpleFile) {
            this.#h5_file = h5File;
        } else {
            this.#h5_file = new afile.SimpleFile(h5File);
        }

        this.#options = H5adDataset.defaults();
        this.clear();
    }

    /**
     * @return {object} Default options, see {@linkcode H5adDataset#setOptions setOptions} for more details.
     */
    static defaults() {
        return {
            countMatrixName: null, 
            featureTypeColumnName: null, 
            featureTypeRnaName: "Gene Expression", 
            featureTypeAdtName: "Antibody Capture", 
            featureTypeCrisprName: "CRISPR Guide Capture", 
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
     * @param {object} options - Optional parameters that affect {@linkcode H5adDataset#load load} (but not {@linkcode H5adDataset#summary summary}).
     * @param {?string} [options.countMatrixName] - Name of the layer containing the count matrix.
     * If `null`, the "X" dataset is used if it is present in the file, or the first available layer if no "X" dataset is present.
     * @param {?string} [options.featureTypeColumnName] - Name of the per-feature annotation column containing the feature types.
     * If `null`, no column is assumed to contain the feature types, and all features are assumed to be genes (i.e., only the RNA modality is present).
     * @param {?string} [options.featureTypeRnaName] - Name of the feature type for gene expression.
     * Alternatively `null`, to indicate that no RNA features are to be loaded.
     * @param {?string} [options.featureTypeAdtName] - Name of the feature type for ADTs.
     * Alternatively `null`, to indicate that no ADT features are to be loaded.
     * @param {?string} [options.featureTypeCrisprName] - Name of the feature type for CRISPR guides.
     * Alternatively `null`, to indicate that no guides are to be loaded.
     * @param {?(string|number)} [options.primaryRnaFeatureIdColumn] - Name or index of the column of the `features` {@linkplain external:DataFrame DataFrame} that contains the primary feature identifier for gene expression.
     *
     * If `i` is `null` or invalid (e.g., out of range index, unavailable name), it is ignored and the row names (from the `_index` group) are used as the primary identifiers.
     * If no row names are present in this situation, no primary identifier is defined.
     * @param {?(string|number)} [options.primaryAdtFeatureIdColumn] - Name or index of the column of the `features` {@linkplain external:DataFrame DataFrame} that contains the primary feature identifier for the ADTs.
     *
     * If `i` is `null` or invalid (e.g., out of range index, unavailable name), it is ignored and the row names (from the `_index` group) are used as the primary identifiers.
     * If no row names are present in this situation, no primary identifier is defined.
     * @param {?(string|number)} [options.primaryCrisprFeatureIdColumn] - Name or index of the column of the `features` {@linkplain external:DataFrame DataFrame} that contains the primary feature identifier for the CRISPR guides.
     *
     * If `i` is `null` or invalid (e.g., out of range index, unavailable name), it is ignored and the row names (from the `_index` group) are used as the primary identifiers.
     * If no row names are present in this situation, no primary identifier is defined.
     */
    setOptions(options) {
        for (const [k, v] of Object.entries(options)) {
            this.#options[k] = v;
        }
    }

    #instantiate() {
        if (this.#h5_path != null) {
            return;
        }

        let info = scran.realizeFile(this.#h5_file.content());
        this.#h5_path = info.path;
        this.#h5_flush = info.flush;
        this.#h5_handle = new scran.H5File(this.#h5_path);
    }

    /**
     * Destroy caches if present, releasing the associated memory. 
     * This may be called at any time but only has an effect if `cache = true` in {@linkcode H5adDataset#load load} or {@linkcodeH5adDataset#summary summary}.
     */
    clear() {
        if (typeof this.#h5_flush == "function") {
            this.#h5_flush();
        }
        this.#h5_flush = null;
        this.#h5_path = null;
        this.#h5_handle = null;

        this.#raw_features = null;
        this.#raw_cells = null;
        this.#assay_details = null;
    }

    /**
     * @return {string} Format of this dataset class.
     * @static
     */
    static format() {
        return "H5AD";
    }

    /**
     * @return {object} Object containing the abbreviated details of this dataset.
     */
    abbreviate() {
        return this.#dump_summary(f => { return { name: f.name(), size: f.size() }; });
    }

    #fetch_assay_details() {
        if (this.#assay_details !== null) {
            return;
        }
        this.#instantiate();
        this.#assay_details = fetch_assay_details(this.#h5_handle, this.#h5_path);
    }

    #features() {
        if (this.#raw_features !== null) {
            return;
        }
        this.#instantiate();

        let feats = fetch_features(this.#h5_handle);
        if (feats == null) {
            this.#fetch_assay_details();
            feats = new bioc.DataFrame({}, { numberOfRows: this.#assay_details.rows });
        }

        this.#raw_features = feats;
        return;
    }

    #cells() {
        if (this.#raw_cells !== null) {
            return;
        }
        this.#instantiate();

        let cells = fetch_cells(this.#h5_handle);
        if (cells === null) {
            this.#fetch_assay_details();
            cells = new bioc.DataFrame({}, { numberOfRows: this.#assay_details.columns })
        }

        this.#raw_cells = cells;
        return;
    }

    /**
     * @param {object} [options={}] - Optional parameters.
     * @param {boolean} [options.cache=false] - Whether to cache the intermediate results for re-use in subsequent calls to any methods with a `cache` option.
     * If `true`, users should consider calling {@linkcode H5adDataset#clear clear} to release the memory once this dataset instance is no longer needed.
     *
     * @return {object} Object containing the per-feature and per-cell annotations.
     * This has the following properties:
     *
     * - `all_features`: a {@linkplain external:DataFrame DataFrame} of per-feature annotations.
     * - `cells`: a {@linkplain external:DataFrame DataFrame} of per-cell annotations.
     * - `all_assay_names`: an Array of strings containing names of potential count matrices.
     */
    summary({ cache = false } = {}) {
        this.#features();
        this.#cells();
        this.#fetch_assay_details();

        let output = {
            all_features: this.#raw_features,
            cells: this.#raw_cells,
            all_assay_names: this.#assay_details.names
        };

        if (!cache) {
            this.clear();
        }
        return output;
    }

    #feature_type_mapping() {
        return {
            RNA: this.#options.featureTypeRnaName, 
            ADT: this.#options.featureTypeAdtName,
            CRISPR: this.#options.featureTypeCrisprName
        };
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
     * If `true`, users should consider calling {@linkcode H5adDataset#clear clear} to release the memory once this dataset instance is no longer needed.
     *
     * @return {object} An object where each key is a modality name and each value is an array (usually of strings) containing the primary feature identifiers for each row in that modality.
     * The contents are the same as the `primary_ids` returned by {@linkcode H5adDataset#load load} but the order of values may be different.
     */
    previewPrimaryIds({ cache = false } = {}) {
        this.#features();
        let preview = futils.extractSplitPrimaryIds(this.#raw_features, this.#options.featureTypeColumnName, this.#feature_type_mapping(), "RNA", this.#primary_mapping());
        if (!cache) {
            this.clear();
        }
        return preview;
    }

    /**
     * @param {object} [options={}] - Optional parameters.
     * @param {boolean} [options.cache=false] - Whether to cache the intermediate results for re-use in subsequent calls to any methods with a `cache` option.
     * If `true`, users should consider calling {@linkcode H5adDataset#clear clear} to release the memory once this dataset instance is no longer needed.
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
     * either from the {@linkcode H5adDataset#defaults defaults} or with {@linkcode H5adDataset#setOptions setOptions}.
     */
    load({ cache = false } = {}) {
        this.#features();
        this.#cells();
        this.#fetch_assay_details();

        let chosen_assay = this.#options.countMatrixName;
        if (chosen_assay == null) {
            chosen_assay = this.#assay_details.names[0];
        }
        let loaded = scran.initializeSparseMatrixFromHDF5(this.#h5_path, chosen_assay);

        let output = futils.splitScranMatrixAndFeatures(loaded, this.#raw_features, this.#options.featureTypeColumnName, this.#feature_type_mapping(), "RNA");
        output.cells = this.#raw_cells;

        output.primary_ids = futils.extractPrimaryIds(output.features, this.#primary_mapping());

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
     * @param {Array} files - Array of objects like that produced by {@linkcode H5adDataset#serialize serialize}.
     * @param {object} options - Object containing additional options to be passed to the constructor.
     * @return {H5adDataset} A new instance of this class.
     * @static
     */
    static async unserialize(files, options) {
        if (files.length != 1 || files[0].type != "h5") {
            throw new Error("expected exactly one file of type 'h5' for H5AD unserialization");
        }
        return new H5adDataset(files[0].file, options);
    }
}

/************************
 ******* Results ********
 ************************/

/**
 * Pre-computed analysis results in the H5AD format.
 */
export class H5adResult {
    #h5_file;
    #h5_path;
    #h5_flush;
    #h5_handle;

    #raw_features;
    #raw_cells;
    #assay_details;
    #reddim_details;

    #options;

    /**
     * @param {SimpleFile|string|Uint8Array|File} h5File - Contents of a H5AD file.
     * On browsers, this may be a File object.
     * On Node.js, this may also be a string containing a file path.
     */
    constructor(h5File) {
        if (h5File instanceof afile.SimpleFile) {
            this.#h5_file = h5File;
        } else {
            this.#h5_file = new afile.SimpleFile(h5File);
        }

        this.#options = H5adResult.defaults();
        this.clear();
    }

    /**
     * @return {object} Default options, see {@linkcode H5adResult#setOptions setOptions} for more details.
     */
    static defaults() {
        return { 
            primaryMatrixName: null, 
            isPrimaryNormalized: true,
            featureTypeColumnName: null, 
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
     * @param {object} options - Optional parameters that affect {@linkcode H5adResult#load load} (but not {@linkcode H5adResult#summary summary}).
     * @param {?string} [options.primaryMatrixName] - Name of the layer containing the primary matrix.
     * If `null`, the "X" dataset is used if it is present in the file, or the first available layer if no "X" dataset is present.
     * @param {boolean} [options.isPrimaryNormalized] - Whether the primary matrix is already normalized.
     * If `false`, it is assumed to contain count data and is subjected to library size normalization within each modality.
     * @param {?string} [options.featureTypeColumnName] - Name of the per-feature annotation column containing the feature types.
     * If `null`, no column is assumed to contain the feature types, and all features are assumed to be genes (i.e., only the RNA modality is present).
     * @param {?Array} [options.reducedDimensionNames=null] - Array of names of the reduced dimensions to load.
     * If `null`, all reduced dimensions found in the file are loaded.
     */
    setOptions(options) {
        for (const [k, v] of Object.entries(options)) {
            if (k == "reducedDimensionNames") {
                this.#options[k] = bioc.CLONE(v); // avoid pass-by-reference links.
            } else {
                this.#options[k] = v;
            }
        }
    }

    /**
     * Destroy caches if present, releasing the associated memory. 
     * This may be called at any time but only has an effect if `cache = true` in {@linkcode H5adResult#load load} or {@linkcode H5adResult#summary summary}.
     */
    clear() {
        if (typeof this.#h5_flush == "function") {
            this.#h5_flush();
        }
        this.#h5_flush = null;
        this.#h5_path = null;
        this.#h5_handle = null;

        this.#raw_features = null;
        this.#raw_cells = null;
        this.#assay_details = null;
        this.#reddim_details = null;
    }

    #instantiate() {
        if (this.#h5_path != null) {
            return;
        }

        let info = scran.realizeFile(this.#h5_file.content());
        this.#h5_path = info.path;
        this.#h5_flush = info.flush;
        this.#h5_handle = new scran.H5File(this.#h5_path);
    }

    #fetch_assay_details() {
        if (this.#assay_details !== null) {
            return;
        }
        this.#instantiate();
        this.#assay_details = fetch_assay_details(this.#h5_handle, this.#h5_path);
        return;
    }

    #features() {
        if (this.#raw_features !== null) {
            return;
        }
        this.#instantiate();

        let feats = fetch_features(this.#h5_handle);
        if (feats == null) {
            this.#fetch_assay_details();
            feats = new bioc.DataFrame({}, { numberOfRows: this.#assay_details.rows });
        }

        this.#raw_features = feats;
        return;
    }

    #cells() {
        if (this.#raw_cells !== null) {
            return;
        }
        this.#instantiate();

        let cells = fetch_cells(this.#h5_handle);
        if (cells === null) {
            this.#fetch_assay_details();
            cells = new bioc.DataFrame({}, { numberOfRows: this.#assay_details.columns })
        }

        this.#raw_cells = cells;
        return;
    }

    #fetch_reddim_details() {
        if (this.#reddim_details !== null) {
            return;
        }
        this.#instantiate();
        
        let available = [];
        if ("obsm" in this.#h5_handle.children && this.#h5_handle.children["obsm"] == "Group") {
            let ohandle = this.#h5_handle.open("obsm");
            for (const [k, v] of Object.entries(ohandle.children)) {
                if (v == "DataSet") {
                    available.push(k);
                }
            }
        }

        this.#reddim_details = { names: available };
        return;
    }

    /**
     * @param {object} [options={}] - Optional parameters.
     * @param {boolean} [options.cache=false] - Whether to cache the results for re-use in subsequent calls to this method or {@linkcode H5adResult#load load}.
     * If `true`, users should consider calling {@linkcode H5adResult#clear clear} to release the memory once this dataset instance is no longer needed.
     *
     * @return {object} Object containing the per-feature and per-cell annotations.
     * This has the following properties:
     *
     * - `all_features`: a {@linkplain external:DataFrame DataFrame} of per-feature annotations.
     * - `cells`: a {@linkplain external:DataFrame DataFrame} of per-cell annotations.
     * - `all_assay_names`: an Array of strings containing names of potential primary matrices.
     * - `reduced_dimension_names`: an Array of strings containing names of dimensionality reduction results.
     */
    summary({ cache = false } = {}) {
        this.#features();
        this.#cells();
        this.#fetch_assay_details();
        this.#fetch_reddim_details();

        let output = {
            all_features: this.#raw_features,
            cells: this.#raw_cells,
            all_assay_names: this.#assay_details.names,
            reduced_dimension_names: this.#reddim_details.names
        };

        if (!cache) {
            this.clear();
        }
        return output;
    }

    /**
     * @param {object} [options={}] - Optional parameters.
     * @param {boolean} [options.cache=false] - Whether to cache the results for re-use in subsequent calls to this method or {@linkcode H5adResult#summary summary}.
     * If `true`, users should consider calling {@linkcode H5adResult#clear clear} to release the memory once this dataset instance is no longer needed.
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
        this.#features();
        this.#cells();
        this.#fetch_assay_details();
        this.#fetch_reddim_details();

        let chosen_assay = this.#options.primaryMatrixName;
        if (chosen_assay == null) {
            chosen_assay = this.#assay_details.names[0];
        }
        let loaded = scran.initializeSparseMatrixFromHDF5(this.#h5_path, chosen_assay, { forceInteger: !this.#options.isPrimaryNormalized });
        let output = futils.splitScranMatrixAndFeatures(loaded, this.#raw_features, this.#options.featureTypeColumnName, null, "");
        output.cells = this.#raw_cells;
        delete output.row_ids;

        if (!this.#options.isPrimaryNormalized) {
            for (const mod of output.matrix.available()) {
                let mat = output.matrix.get(mod);
                output.matrix.add(mod, scran.logNormCounts(mat, { allowZeros: true }));
            }
        }

        // Loading the dimensionality reduction results.
        let chosen_reddims = this.#options.reducedDimensionNames;
        if (chosen_reddims == null) {
            chosen_reddims = this.#reddim_details.names;
        }

        let reddims = {};
        if (chosen_reddims.length) {
            let ohandle = this.#h5_handle.open("obsm");
            for (const k of chosen_reddims) {
                let loaded = ohandle.open(k, { load: true });
                let shape = loaded.shape;
                let ndims = shape[1];
                let ncells = shape[0];

                let contents = [];
                for (var d = 0; d < ndims; d++) {
                    contents.push(new Float64Array(ncells));
                }

                for (var c = 0; c < ncells; c++) {
                    for (var d = 0; d < ndims; d++) {
                        contents[d][c] = loaded.values[c * ndims + d];
                    }
                }

                reddims[k] = contents;
            }
        }
        output.reduced_dimensions = reddims;
        
        if (!cache) {
            this.clear();
        }
        return output;
    }

}

