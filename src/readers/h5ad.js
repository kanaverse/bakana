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

function fetch_features(handle) {
    if ("var" in handle.children && handle.children["var"] == "Group") {
        let vhandle = handle.open("var");
        let index = eutils.extractHDF5Strings(vhandle, "_index");

        if (index !== null) {
            // Build it piece-by-piece for a well-defined order.
            let genes = new bioc.DataFrame({ "_index": index });

            for (const [key, val] of Object.entries(vhandle.children)) {
                if (val === "DataSet" && (key.match(/name/i) || key.match(/symb/i))) {
                    let dhandle2 = vhandle.open(key);
                    if (dhandle2.type == "String") {
                        genes.$setColumn(key, dhandle2.load());
                    }
                }
            }

            return genes;
        }
    }

    return null;
}

function fetch_cells(handle) {
    let annotations = {};

    if ("obs" in handle.children && handle.children["obs"] == "Group") {
        let ohandle = handle.open("obs");

        // Maybe it has names, maybe not, who knows; let's just add what's there.
        let index = eutils.extractHDF5Strings(ohandle, "_index");
        if (index !== null) {
            annotations["_index"] = index;
        }

        for (const [key, val] of Object.entries(ohandle.children)) {
            if (val != "DataSet") {
                continue;
            }
            let dhandle = ohandle.open(key, { load: true });
            annotations[key] = dhandle.values;
        }

        if ("__categories" in ohandle.children && ohandle.children["__categories"] == "Group") {
            let chandle = ohandle.open("__categories");

            for (const [key, val] of Object.entries(chandle.children)) {
                if (key in annotations) {
                    let current_levels = eutils.extractHDF5Strings(chandle, key);
                    annotations[key] = bioc.SLICE(current_levels, annotations[key]);
                }
            }
        }
    }

    if (Object.keys(annotations).length > 0) {
        return new bioc.DataFrame(annotations);
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

    #countMatrixName;
    #featureTypeColumnName;

    #featureTypeRnaName;
    #featureTypeAdtName;
    #featureTypeCrisprName;

    #primaryRnaFeatureIdColumn;
    #primaryAdtFeatureIdColumn;
    #primaryCrisprFeatureIdColumn;

    #dump_summary(fun) {
        let files = [{ type: "h5", file: fun(this.#h5_file) }]; 
        let options = {
            countMatrixName: this.#countMatrixName,
            featureTypeColumnName: this.#featureTypeColumnName,
            featureTypeRnaName: this.#featureTypeRnaName,
            featureTypeAdtName: this.#featureTypeAdtName,
            featureTypeCrisprName: this.#featureTypeCrisprName,
            primaryRnaFeatureIdColumn: this.#primaryRnaFeatureIdColumn,
            primaryAdtFeatureIdColumn: this.#primaryAdtFeatureIdColumn,
            primaryCrisprFeatureIdColumn: this.#primaryCrisprFeatureIdColumn
        };
        return { files, options };
    }

    /**
     * @param {SimpleFile|string|Uint8Array|File} h5File - Contents of a H5AD file.
     * On browsers, this may be a File object.
     * On Node.js, this may also be a string containing a file path.
     * @param {object} [options={}] - Optional parameters.
     * @param {?string} [options.countMatrixName=null] - See {@linkcode H5adDataset#setCountMatrixName setCountMatrixName}.
     * @param {?string} [options.featureTypeColumnName=null] - See {@linkcode H5adDataset#setFeatureTypeColumnName setFeatureTypeColumnName}.
     * @param {?string} [options.featureTypeRnaName="Gene Expression"] - See {@linkcode H5adDataset#setFeatureTypeRnaName setFeatureTypeRnaName}.
     * @param {?string} [options.featureTypeAdtName="Antibody Capture"] - See {@linkcode H5adDataset#setFeatureTypeAdtName setFeatureTypeAdtName}.
     * @param {?string} [options.featureTypeCrisprName="CRISPR Guide Capture"] - See {@linkcode H5adDataset#setFeatureTypeCrisprName setFeatureTypeCrisprName}.
     * @param {?(string|number)} [options.primaryRnaFeatureIdColumn=0] - See {@linkcode H5adDataset#setPrimaryRnaFeatureIdColumn setPrimaryRnaFeatureIdColumn}.
     * @param {?(string|number)} [options.primaryAdtFeatureIdColumn=0] - See {@linkcode H5adDataset#setPrimaryAdtFeatureIdColumn setPrimaryAdtFeatureIdColumn}.
     * @param {?(string|number)} [options.primaryCrisprFeatureIdColumn=0] - See {@linkcode H5adDataset#setPrimaryCrisprFeatureIdColumn setPrimaryCrisprFeatureIdColumn}.
     */
    constructor(h5File, { 
        countMatrixName = null, 
        featureTypeColumnName = null, 
        featureTypeRnaName = "Gene Expression", 
        featureTypeAdtName = "Antibody Capture", 
        featureTypeCrisprName = "CRISPR Guide Capture", 
        primaryRnaFeatureIdColumn = 0, 
        primaryAdtFeatureIdColumn = 0,
        primaryCrisprFeatureIdColumn = 0 
    } = {}) {
        if (h5File instanceof afile.SimpleFile) {
            this.#h5_file = h5File;
        } else {
            this.#h5_file = new afile.SimpleFile(h5File);
        }

        this.#countMatrixName = countMatrixName;
        this.#featureTypeColumnName = featureTypeColumnName;
        
        this.#featureTypeRnaName = featureTypeRnaName;
        this.#featureTypeAdtName = featureTypeAdtName;
        this.#featureTypeCrisprName = featureTypeCrisprName;

        this.#primaryRnaFeatureIdColumn = primaryRnaFeatureIdColumn;
        this.#primaryAdtFeatureIdColumn = primaryAdtFeatureIdColumn;
        this.#primaryCrisprFeatureIdColumn = primaryCrisprFeatureIdColumn;

        this.clear();
    }

    /**
     * @param {?string} name - Name of the layer containing the count matrix.
     * If `null`, the "X" dataset is used if it is present in the file, or the first available layer if no "X" dataset is present.
     */
    setCountMatrixName(name) {
        this.#countMatrixName = name;
        return;
    }

    /**
     * @param {?string} name - Name of the per-feature annotation column containing the feature types.
     * If `null`, no column is assumed to contain the feature types, and all features are assumed to be genes (i.e., only the RNA modality is present).
     */
    setFeatureTypeColumnName(name) {
        this.#featureTypeColumnName = name;
        return;
    }

    /**
     * @param {?string} name - Name of the feature type for gene expression.
     * Alternatively `null`, to indicate that no RNA features are to be loaded.
     */
    setFeatureTypeRnaName(name) {
        this.#featureTypeRnaName = name;
        return;
    }

    /**
     * @param {string} name - Name of the feature type for ADTs.
     * Alternatively `null`, to indicate that no ADT features are to be loaded.
     */
    setFeatureTypeAdtName(name) {
        this.#featureTypeAdtName = name;
        return;
    }

    /**
     * @param {?string} name - Name of the feature type for CRISPR guides.
     * Alternatively `null`, to indicate that no guides are to be loaded.
     */
    setFeatureTypeCrisprName(name) {
        this.#featureTypeCrisprName = name;
        return;
    }

    /**
     * @param {?(string|number)} i - Name or index of the column of the `features` {@linkplain external:DataFrame DataFrame} that contains the primary feature identifier for gene expression.
     * If `i` is invalid (e.g., out of range index, unavailable name), it is ignored and the primary identifier is treated as undefined.
     */
    setPrimaryRnaFeatureIdColumn(i) {
        this.#primaryRnaFeatureIdColumn = i;
        return;
    }

    /**
     * @param {?(string|number)} i - Name or index of the column of the `features` {@linkplain external:DataFrame DataFrame} that contains the primary feature identifier for the ADTs.
     * If `i` is invalid (e.g., out of range index, unavailable name), it is ignored and the primary identifier is treated as undefined.
     */
    setPrimaryAdtFeatureIdColumn(i) {
        this.#primaryAdtFeatureIdColumn = i;
        return;
    }

    /**
     * @param {?(string|number)} i - Name or index of the column of the `features` {@linkplain external:DataFrame DataFrame} that contains the primary feature identifier for the CRISPR guides.
     * If `i` is invalid (e.g., out of range index, unavailable name), it is ignored and the primary identifier is treated as undefined.
     */
    setPrimaryCrisprFeatureIdColumn(i) {
        this.#primaryCrisprFeatureIdColumn = i;
        return;
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
     * @param {boolean} [options.cache=false] - Whether to cache the results for re-use in subsequent calls to this method or {@linkcode H5adDataset#load load}.
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

    /**
     * @param {object} [options={}] - Optional parameters.
     * @param {boolean} [options.cache=false] - Whether to cache the results for re-use in subsequent calls to this method or {@linkcode H5adDataset#summary summary}.
     * If `true`, users should consider calling {@linkcode H5adDataset#clear clear} to release the memory once this dataset instance is no longer needed.
     *
     * @return {object} Object containing the per-feature and per-cell annotations.
     * This has the following properties:
     *
     * - `features`: an object where each key is a modality name and each value is a {@linkplain external:DataFrame DataFrame} of per-feature annotations for that modality.
     * - `cells`: a {@linkplain external:DataFrame DataFrame} containing per-cell annotations.
     * - `matrix`: a {@linkplain external:MultiMatrix MultiMatrix} containing one {@linkplain external:ScranMatrix ScranMatrix} per modality.
     *
     * Modality names are guaranteed to be one of `"RNA"` or `"ADT"`.
     * It is assumed that an appropriate mapping from the feature types inside the `featureFile` was previously declared,
     * either in the constructor or in {@linkcode setFeatureTypeRnaName} and {@linkcode setFeatureTypeAdtName}.
     */
    load({ cache = false } = {}) {
        this.#features();
        this.#cells();
        this.#fetch_assay_details();

        let chosen_assay = this.#countMatrixName;
        if (chosen_assay == null) {
            chosen_assay = this.#assay_details.names[0];
        }
        let loaded = scran.initializeSparseMatrixFromHDF5(this.#h5_path, chosen_assay);

        let mappings = { 
            RNA: this.#featureTypeRnaName, 
            ADT: this.#featureTypeAdtName,
            CRISPR: this.#featureTypeCrisprName
        };
        let output = futils.splitScranMatrixAndFeatures(loaded, this.#raw_features, this.#featureTypeColumnName, mappings, "RNA");
        output.cells = this.#raw_cells;

        let primaries = { 
            RNA: this.#primaryRnaFeatureIdColumn, 
            ADT: this.#primaryAdtFeatureIdColumn,
            CRISPR: this.#primaryCrisprFeatureIdColumn
        };
        futils.decorateWithPrimaryIds(output.features, primaries);

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

    #primaryMatrixName;
    #isPrimaryNormalized;
    #featureTypeColumnName;
    #reducedDimensionNames;

    /**
     * @param {SimpleFile|string|Uint8Array|File} h5File - Contents of a H5AD file.
     * On browsers, this may be a File object.
     * On Node.js, this may also be a string containing a file path.
     * @param {object} [options={}] - Optional parameters.
     * @param {?string} [options.primaryMatrixName=null] - See {@linkcode H5adResult#setPrimaryMatrixName setPrimaryMatrixName}.
     * @param {boolean} [options.isPrimaryNormalized=true] - See {@linkcode H5adResult#setIsPrimaryNormalized setIsPrimaryNormalized}.
     * @param {?string} [options.featureTypeColumnName=null] - See {@linkcode H5adResult#setFeatureTypeColumnName setFeatureTypeColumnName}.
     * @param {?Array} [options.reducedDimensionNames=null] - See {@linkcode H5adResult#setReducedDimensionNames setReducedDimensionNames}.
     */
    constructor(h5File, { 
        primaryMatrixName = null, 
        isPrimaryNormalized=true,
        featureTypeColumnName = null, 
        reducedDimensionNames = null
    } = {}) {
        if (h5File instanceof afile.SimpleFile) {
            this.#h5_file = h5File;
        } else {
            this.#h5_file = new afile.SimpleFile(h5File);
        }

        this.#primaryMatrixName = primaryMatrixName;
        this.#isPrimaryNormalized = isPrimaryNormalized;
        this.#featureTypeColumnName = featureTypeColumnName;
        this.#reducedDimensionNames = reducedDimensionNames;

        this.clear();
    }

    /**
     * Destroy caches if present, releasing the associated memory. 
     * This may be called at any time but only has an effect if `cache = true` in {@linkcode H5adResult#load load} or {@linkcodeH5adResult#summary summary}.
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

    /**
     * @param {?string} name - Name of the layer containing the primary matrix.
     * If `null`, the "X" dataset is used if it is present in the file, or the first available layer if no "X" dataset is present.
     */
    setPrimaryMatrixName(name) {
        this.#primaryMatrixName = name;
        return;
    }

    /**
     * @param {boolean} normalized - Whether the primary matrix is already normalized.
     * If `false`, it is assumed to contain count data and is subjected to library size normalization within each modality.
     */
    setIsPrimaryNormalized(normalized) {
        this.#isPrimaryNormalized = normalized;
    }

    /**
     * @param {?string} name - Name of the per-feature annotation column containing the feature types.
     * If `null`, no column is assumed to contain the feature types, and all features are assumed to be genes (i.e., only the RNA modality is present).
     */
    setFeatureTypeColumnName(name) {
        this.#featureTypeColumnName = name;
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
     * @param {boolean} [options.cache=false] - Whether to cache the results for re-use in subsequent calls to this method or {@linkcode H5adDataset#load load}.
     * If `true`, users should consider calling {@linkcode H5adDataset#clear clear} to release the memory once this dataset instance is no longer needed.
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

        let chosen_assay = this.#primaryMatrixName;
        if (chosen_assay == null) {
            chosen_assay = this.#assay_details.names[0];
        }
        let loaded = scran.initializeSparseMatrixFromHDF5(this.#h5_path, chosen_assay, { forceInteger: !this.#isPrimaryNormalized });
        let output = futils.splitScranMatrixAndFeatures(loaded, this.#raw_features, this.#featureTypeColumnName, null, "");
        output.cells = this.#raw_cells;
        delete output.row_ids;

        if (!this.#isPrimaryNormalized) {
            for (const mod of output.matrix.available()) {
                let mat = output.matrix.get(mod);
                output.matrix.add(mod, scran.logNormCounts(mat, { allowZeros: true }));
            }
        }

        // Loading the dimensionality reduction results.
        let chosen_reddims = this.#reducedDimensionNames;
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

