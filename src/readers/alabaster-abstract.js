import * as scran from "scran.js";
import * as bioc from "bioconductor";
import * as eutils from "./utils/extract.js";
import * as futils from "./utils/features.js";
import * as afile from "./abstract/file.js";
import * as jsp from "jaspagate";

/**
 * Any class that satisfies the AlabasterProjectNavigator contract, so called as it is intended to "navigate" an alabaster-formatted object directory.
 * This should provide the following methods:
 * 
 * - `get(path, asBuffer)`, a (possibly async) method that accepts a string `path` containing a relative path to a file inside an object directory. 
 *   This should return a string containing a path to the file on the local filesystem, or a Uint8Array containing the contents of the file if no local filesystem exists.
 *   If `asBuffer = true`, a Uint8Array must be returned.
 * - `exists(path)`, a (possibly async) method that accepts a string `path` containing a relative path to a file inside an object directory. 
 *   This should return a boolean indicating whether `path` exists in the object directory.
 * - `clean(localPath)`, a (possibly async) method that accepts a string containing a local path returned by `get()`.
 *   It should remove any temporary file that was created by `get()` at `localPath`.
 *   This will not be called if `get()` returns a Uint8Array.
 *
 * @typedef AlabasterProjectNavigator
 */

/****************************
 *** Jaspagate interfaces ***
 ****************************/

class AlabasterH5Group extends jsp.H5Group {
    #handle;
    #flush;

    constructor(handle, flush) {
        super();
        this.#handle = handle;
        this.#flush = flush;
    }

    attributes() {
        return this.#handle.attributes;
    }

    readAttribute(attr) {
        let ares = this.#handle.readAttribute(attr);
        return { values: ares.values, shape: ares.shape };
    }

    children() {
        return Array.from(Object.keys(this.#handle.children));
    }

    open(name) {
        let out = this.#handle.open(name);
        if (out instanceof scran.H5Group) {
            return new AlabasterH5Group(out);
        } else {
            return new AlabasterH5DataSet(out);
        }
    }

    close() {}

    _flush() {
        this.#flush();
    }
}

class AlabasterH5DataSet extends jsp.H5DataSet {
    #handle;

    constructor(handle) {
        super();
        this.#handle = handle;
    }

    attributes() {
        return this.#handle.attributes;
    }

    readAttribute(attr) {
        let ares = this.#handle.readAttribute(attr);
        return { values: ares.values, shape: ares.shape };
    }

    type() {
        let type = this.#handle.type;
        if (type instanceof scran.H5StringType) {
            return "String";
        } else if (type instanceof scran.H5CompoundType) {
            return type.members;
        } else {
            return type;
        }
    }

    shape() {
        return this.#handle.shape;
    }

    values() {
        return this.#handle.values;
    }

    close() {}
}

class AlabasterGlobalsInterface extends jsp.GlobalsInterface {
    #navigator;

    constructor(navigator) {
        super();
        this.#navigator = navigator;
    }

    get(path, options = {}) {
        const { asBuffer = false } = options;
        return this.#navigator.get(path, asBuffer);
    }

    exists(path) {
        return this.#navigator.exists(path);
    }

    clean(localPath) {
        this.#navigator.clean(localPath); 
    }

    async h5open(path) {
        let realized = scran.realizeFile(await this.get(path));
        try {
            return new AlabasterH5Group(new scran.H5File(realized.path), realized.flush);
        } catch (e) {
            realized.flush();
            throw e;
        }
    }

    h5close(handle) {
        handle._flush();
    }
}

/*********************
 *** Assay readers ***
 *********************/

class MockMatrix {
    #nrow;
    #ncol;
    #path;

    constructor(nrow, ncol, path) {
        this.#nrow = nrow;
        this.#ncol = ncol;
        this.#path = path;
    }

    _bioconductor_NUMBER_OF_ROWS() {
        return this.#nrow;
    }

    _bioconductor_NUMBER_OF_COLUMNS() {
        return this.#ncol;
    }

    async realize(globals, forceInteger, forceSparse) {
        let metadata = await jsp.readObjectFile(this.#path, globals);
        if (metadata.type == "delayed_array") {
            let contents = await globals.get(jsp.joinPath(this.#path, "array.h5"));
            try {
                let realized = scran.realizeFile(contents);
                try {
                    let handle = new scran.H5File(realized.path);
                    let output = await extract_delayed(handle.open("delayed_array"), this.#path, globals, forceInteger, forceSparse);
                    if (output == null) {
                        throw new Error("currently only supporting bakana-generated log-counts for delayed arrays");
                    }
                    return output;
                } finally {
                    realized.flush();
                }
            } finally {
                await globals.clean(contents);
            }
        } else {
            return extract_matrix(this.#path, metadata, globals, forceInteger, forceSparse);
        }
    }
}

async function extract_matrix(path, metadata, globals, { forceInteger = true, forceSparse = true } = {}) {
    if (metadata.type == "compressed_sparse_matrix") {
        let contents = await globals.get(jsp.joinPath(path, "matrix.h5"));
        try {
            let realized = scran.realizeFile(contents);
            try {
                let fhandle = new scran.H5File(realized.path);
                const name = "compressed_sparse_matrix";

                let dhandle = fhandle.open(name);
                const shape = dhandle.open("shape").values; 
                const layout = dhandle.readAttribute("layout").values[0];

                let out = scran.initializeSparseMatrixFromHdf5Group(realized.path, name, shape[0], shape[1], (layout == "CSR"), { forceInteger });
                return out;
            } finally {
                realized.flush();
            }
        } finally {
            await globals.clean(contents);
        }

    } else if (metadata.type == "dense_array") {
        let contents = await globals.get(jsp.joinPath(path, "array.h5"));
        try {
            let realized = scran.realizeFile(contents);
            try {
                let fhandle = new scran.H5File(realized.path);
                const name = "dense_array";
                let dhandle = fhandle.open(name);
                let transposed = false;
                if (dhandle.attributes.indexOf("transposed") >= 0) {
                    let trans_info = dhandle.readAttribute("transposed");
                    transposed = (trans_info.values[0] != 0);
                }

                return scran.initializeMatrixFromHdf5Dataset(realized.path, name + "/data", { transposed, forceInteger, forceSparse });
            } finally {
                realized.flush();
            }
        } finally {
            await globals.clean(contents);
        }

    } else {
        throw new Error("unknown matrix type '" + metadata.type + "'");
    }
}

async function extract_delayed(handle, path, globals, forceInteger, forceSparse) {
    const dtype = handle.readAttribute("delayed_type").values[0];
    if (dtype === "operation") {
        let optype = handle.readAttribute("delayed_operation").values[0];

        if (optype === "unary arithmetic") {
            const seed = await extract_delayed(handle.open("seed"), path, globals, forceInteger, forceSparse)

            const dhandle = handle.open("value");
            let arg = dhandle.values;
            if (dhandle.shape.length == 0) {
                arg = arg[0];
            }

            return scran.delayedArithmetic(
                seed,
                handle.open("method").values[0],
                arg,
                {
                    right: handle.open("side").values[0] !== "right",
                    inPlace: true
                }
            );

        } else if (optype === "unary math") {
            const seed = await extract_delayed(handle.open("seed"), path, globals, forceInteger, forceSparse)
            const meth = handle.open("method").values[0];
            let base = null;
            if (meth == "log") {
                if ("base" in handle.children) {
                    base = handle.open("base").values[0];
                }
            }
            return scran.delayedMath(
                seed,
                meth,
                {
                    logBase: base,
                    inPlace: true
                }
            );

        } else if (optype == "transpose") {
            const seed = await extract_delayed(handle.open("seed"), path, globals, forceInteger, forceSparse)
            const perm = handle.open("permutation").values;
            if (perm[0] == 1 && perm[1] == 0) {
                return scran.transpose(seed, { inPlace: true });
            } else if (perm[0] == 0 && perm[1] == 1) {
                return seed;
            } else {
                throw new Error("invalid permutation for transposition operation at '" + path + "'");
            }

        } else if (optype == "subset") {
            let mat = await extract_delayed(handle.open("seed"), path, globals, forceInteger, forceSparse)
            const ihandle = handle.open("index");
            if ("0" in ihandle.children) {
                mat = scran.subsetRows(mat, ihandle.open("0").values, { inPlace: true });
            }
            if ("1" in ihandle.children) {
                mat = scran.subsetColumns(mat, ihandle.open("1").values, { inPlace: true });
            }
            return mat;

        } else if (optype == "combine") {
            const shandle = handle.open("seeds");
            let seeds = [];
            try {
                const nchildren = Object.keys(shandle.childen).length;
                for (var c = 0; c < nchildren; c++) {
                    seeds.push(await extract_delayed(shandle.open(String(c)), path, globals, forceInteger, forceSparse));
                }
                if (handle.open("along").values[0] == 0) {
                    return scran.rbind(seeds);
                } else {
                    return scran.cbind(seeds);
                }
            } finally {
                for (const s of seeds) {
                    scran.free(s);
                }
            }

        } else {
            throw new Error("unsupported delayed operation '" + optype + "'");
        }

    } else if (dtype === "array") {
        let atype = handle.readAttribute("delayed_array").values[0];

        if (atype === "custom takane seed array") {
            let index = handle.open("index").values[0];
            let seed_path = jsp.joinPath(path, "seeds", String(index));
            let seed_metadata = await jsp.readObjectFile(seed_path, globals);
            let mat;
            let output;
            try {
                return await extract_matrix(seed_path, seed_metadata, globals, forceInteger, forceSparse); 
            } finally {
                scran.free(mat);
            }

        } else if (atype == "dense array") {
            let is_native = handle.open("native").values[0] != 0;
            return scran.initializeMatrixFromHdf5Dataset(handle.file, handle.name + "/data", { transposed: !is_native, forceInteger, forceSparse });

        } else if (atype == "sparse matrix") {
            const shape = handle.open("shape").values; 
            const is_csr = handle.open("by_column").values[0] == 0;
            return scran.initializeSparseMatrixFromHdf5Group(handle.file, handle.name, shape[0], shape[1], is_csr, { forceInteger });

        } else {
            throw new Error("unsupported delayed array '" + atype + "'");
        }

    } else {
        throw new Error("unsupported delayed type '" + dtype + "'");
    }

    return output;
}

function readMockAssay(nrow, ncol, path, metadata, globals, options) {
    return new MockMatrix(nrow, ncol, path);
}

function readMockReducedDimension(ncol, path, metadata, globals, options) {
    return new MockMatrix(ncol, 2, path);
}

/*************************
 *** Utility functions ***
 *************************/

function apply_over_experiments(se, fun) {
    let main_experiment_name = "";
    let is_sce = se instanceof bioc.SingleCellExperiment;
    if (is_sce && se.mainExperimentName() !== null) {
        main_experiment_name = se.mainExperimentName();
    }

    let output = {};
    output[main_experiment_name] = fun(se);
    if (is_sce) {
        for (const alt of se.alternativeExperimentNames()) {
            output[alt] = fun(se.alternativeExperiment(alt));
        }
    }
    return output;
}

function extract_all_features(se) {
    return apply_over_experiments(se, x => x.rowData());
}

function extract_all_assay_names(se) {
    return apply_over_experiments(se, x => x.assayNames());
}

function simplify_List_columns(df) { // avoid the hassle of dealing with List compatibility problems in the rest of bakana.
    for (const k of df.columnNames()) {
        let col = df.column(k);
        if (col instanceof bioc.List) {
            df.setColumn(k, col.toArray(), { inPlace: true });
        }
    }
    return null;
}

/************************
 ******* Dataset ********
 ************************/

/**
 * Dataset stored as a SummarizedExperiment in the **alabaster** format.
 * This is intended as a virtual base class; applications should define subclasses that are tied to a specific {@linkplain AlabasterProjectNavigator} class.
 * Subclasses should define `abbreviate()` and `serialize()` methods, as well as the static `format()` and `unserialize()` methods - 
 * see the [Dataset contract](https://github.com/LTLA/bakana/blob/master/docs/related/custom_readers.md) for more details.
 */
export class AbstractAlabasterDataset {
    #navigator;
    #raw_se;
    #options;

    /**
     * @param {AlabasterProjectNavigator} navigator - A navigator object that describes how to obtain files from the alabaster-formatted object directory.
     */
    constructor(navigator) {
        this.#navigator = navigator;
        this.#options = AbstractAlabasterDataset.defaults();
        this.#raw_se = null;
    }

    /**
     * @return {object} Default options, see {@linkcode AbstractAlabasterDataset#setOptions setOptions} for more details.
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
     * @param {object} options - Optional parameters that affect {@linkcode AbstractAlabasterDataset#load load} (but not {@linkcode AbstractAlabasterDataset#summary summary}).
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
     * This may be called at any time but only has an effect if `cache = true` in {@linkcode AbstractAlabasterDataset#load load} or {@linkcode AbstractAlabasterDataset#summary summary}.
     */
    clear() {
        this.#raw_se = null;
    }

    #create_globals() { 
        return new AlabasterGlobalsInterface(this.#navigator);
    }

    async #populate() {
        if (this.#raw_se === null) {
            this.#raw_se = await jsp.readObject(
                ".",
                null, 
                this.#create_globals(),
                {
                    DataFrame_readNested: false,
                    DataFrame_readMetadata: false,
                    SummarizedExperiment_readAssay: readMockAssay,
                    SummarizedExperiment_readMetadata: false,
                    SingleCellExperiment_readReducedDimension: false
                }
            );
            simplify_List_columns(this.#raw_se.columnData());
            apply_over_experiments(this.#raw_se, y => simplify_List_columns(y.rowData()));
        }
    }

    /**
     * @param {object} [options={}] - Optional parameters.
     * @param {boolean} [options.cache=false] - Whether to cache the intermediate results for re-use in subsequent calls to any methods with a `cache` option.
     * If `true`, users should consider calling {@linkcode AbstractAlabasterDataset#clear clear} to release the memory once this dataset instance is no longer needed.
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
        await this.#populate();

        let output = {
            modality_features: extract_all_features(this.#raw_se),
            cells: this.#raw_se.columnData(),
            modality_assay_names: extract_all_assay_names(this.#raw_se)
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
     * If `true`, users should consider calling {@linkcode AbstractAlabasterDataset#clear clear} to release the memory once this dataset instance is no longer needed.
     *
     * @return {object} An object where each key is a modality name and each value is an array (usually of strings) containing the primary feature identifiers for each row in that modality.
     * The contents are the same as the `primary_ids` returned by {@linkcode AbstractAlabasterDataset#load load} but the order of values may be different.
     *
     * @async
     */
    async previewPrimaryIds({ cache = false } = {}) {
        await this.#populate();

        let fmapping = {
            RNA: this.#options.rnaExperiment, 
            ADT: this.#options.adtExperiment, 
            CRISPR: this.#options.crisprExperiment 
        };

        let raw_features = extract_all_features(this.#raw_se);
        let preview = futils.extractRemappedPrimaryIds(raw_features, fmapping, this.#primary_mapping());

        if (!cache) {
            this.clear();
        }
        return preview;
    }

    /**
     * @param {object} [options={}] - Optional parameters.
     * @param {boolean} [options.cache=false] - Whether to cache the intermediate results for re-use in subsequent calls to any methods with a `cache` option.
     * If `true`, users should consider calling {@linkcode AbstractAlabasterDataset#clear clear} to release the memory once this dataset instance is no longer needed.
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
     * either from the {@linkcode AbstractAlabasterDataset#defaults defaults} or with {@linkcode AbstractAlabasterDataset#setOptions setOptions}.
     *
     * @async
     */
    async load({ cache = false } = {}) {
        await this.#populate();

        let output = { 
            matrix: new scran.MultiMatrix,
            features: {},
            cells: this.#raw_se.columnData()
        };

        let mapping = { 
            RNA: { exp: this.#options.rnaExperiment, assay: this.#options.rnaCountAssay },
            ADT: { exp: this.#options.adtExperiment, assay: this.#options.adtCountAssay },
            CRISPR: { exp: this.#options.crisprExperiment, assay: this.#options.crisprCountAssay }
        };

        let experiments_by_name = apply_over_experiments(this.#raw_se, x => x);
        let num_alts = 0;
        if (this.#raw_se instanceof bioc.SingleCellExperiment) {
            num_alts = this.#raw_se.alternativeExperimentNames().length;
        }

        try {
            for (const [k, v] of Object.entries(mapping)) {
                if (v.exp === null) {
                    continue;
                }

                let chosen_se;
                if (typeof v.exp == "string") {
                    if (!(v.exp in experiments_by_name)) {
                        continue;
                    }
                    chosen_se = experiments_by_name[v.exp];
                } else {
                    if (v.exp >= num_alts) {
                        continue;
                    }
                    chosen_se = this.#raw_se.alternativeExperiment(v.exp);
                }

                let loaded = await chosen_se.assay(v.assay).realize(this.#create_globals(), /* forceInteger = */ true, /* forceSparse = */ true);
                output.matrix.add(k, loaded);
                output.features[k] = chosen_se.rowData();
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
 * Pre-computed analysis results stored as a SummarizedExperiment object (or one of its subclasses) in the **ArtifactDB** format.
 * This is intended as a virtual base class; applications should define subclasses that are tied to a specific {@linkplain AlabasterProjectNavigator} class.
 */
export class AbstractAlabasterResult {
    #navigator;
    #raw_se;
    #options;

    /**
     * @param {AlabasterNavigator} navigator - A navigator object that describes how to obtain files from an alabaster-formatted object directory.
     */
    constructor(navigator) {
        this.#navigator = navigator;
        this.#options = AbstractAlabasterResult.defaults();
        this.#raw_se = null;
    }

    /**
     * @return {object} Default options, see {@linkcode AbstractAlabasterResults#setOptions setOptions} for more details.
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
     * @param {object} options - Optional parameters that affect {@linkcode AbstractAlabasterResult#load load} (but not {@linkcode AbstractAlabasterResult#summary summary}.
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
     * This may be called at any time but only has an effect if `cache = true` in {@linkcode AbstractAlabasterResult#load load} or {@linkcode AbstractAlabasterResult#summary summary}.
     */
    clear() {
        this.#raw_se = null;
    }

    #create_globals() { 
        return new AlabasterGlobalsInterface(this.#navigator);
    }

    async #populate() {
        if (this.#raw_se === null) {
            this.#raw_se = await jsp.readObject(
                ".",
                null, 
                this.#create_globals(),
                {
                    DataFrame_readNested: false,
                    DataFrame_readMetadata: false,
                    SummarizedExperiment_readAssay: readMockAssay,
                    SummarizedExperiment_readMetadata: false,
                    SingleCellExperiment_readReducedDimension: readMockReducedDimension
                }
            );
            simplify_List_columns(this.#raw_se.columnData());
            apply_over_experiments(this.#raw_se, y => simplify_List_columns(y.rowData()));
        }
    }

    /**
     * @param {object} [options={}] - Optional parameters.
     * @param {boolean} [options.cache=false] - Whether to cache the results for re-use in subsequent calls to this method or {@linkcode AbstractAlabasterResult#load load}.
     * If `true`, users should consider calling {@linkcode AbstractAlabasterResult#clear clear} to release the memory once this dataset instance is no longer needed.
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
        await this.#populate();

        let output = {
            modality_features: extract_all_features(this.#raw_se),
            cells: this.#raw_se.columnData(),
            modality_assay_names: extract_all_assay_names(this.#raw_se),
            reduced_dimension_names: []
        };

        if (this.#raw_se instanceof bioc.SingleCellExperiment) {
            output.reduced_dimension_names = this.#raw_se.reducedDimensionNames();
        }

        if (!cache) {
            this.clear();
        }
        return output;
    }

    /**
     * @param {object} [options={}] - Optional parameters.
     * @param {boolean} [options.cache=false] - Whether to cache the results for re-use in subsequent calls to this method or {@linkcode AbstractAlabasterResult#summary summary}.
     * If `true`, users should consider calling {@linkcode AbstractAlabasterResult#clear clear} to release the memory once this dataset instance is no longer needed.
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
        await this.#populate();

        let output = { 
            matrix: new scran.MultiMatrix,
            features: {},
            cells: this.#raw_se.columnData(),
            reduced_dimensions: {}
        };

        if (this.#raw_se instanceof bioc.SingleCellExperiment) {
            let chosen_rd = this.#options.reducedDimensionNames;
            if (chosen_rd === null) {
                chosen_rd = this.#raw_se.reducedDimensionNames();
            }
            for (const k of chosen_rd) {
                let current = await this.#raw_se.reducedDimension(k).realize(this.#create_globals(), /* forceInteger = */ false, /* forceSparse = */ false);
                try {
                    let collected = [];
                    let ncol = current.numberOfColumns();
                    for (var c = 0; c < ncol; c++) {
                        collected.push(current.column(c));
                    }
                    output.reduced_dimensions[k] = collected;
                } finally {
                    scran.free(current);
                }
            }
        }

        // Now fetching the assay matrix.
        const experiments_by_name = apply_over_experiments(this.#raw_se, x => x);
        try {
            for (const [name, chosen_se] of Object.entries(experiments_by_name)) {
                let curassay = this.#options.primaryAssay;
                if (typeof curassay == "object") {
                    if (name in curassay) {
                        curassay = curassay[name];
                    } else {
                        continue;
                    }
                }

                let curnormalized = this.#options.isPrimaryNormalized;
                if (typeof curnormalized == "object") {
                    if (name in curnormalized) {
                        curnormalized = curnormalized[name];
                    } else {
                        curnormalized = true;
                    }
                }

                let loaded = await chosen_se.assay(curassay).realize(this.#create_globals(), /* forceInteger= */ !curnormalized, /* forceSparse = */ true);
                output.matrix.add(name, loaded);

                if (!curnormalized) {
                    let normed = scran.normalizeCounts(loaded, { allowZeros: true });
                    output.matrix.add(name, normed);
                }

                output.features[name] = chosen_se.rowData();
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
