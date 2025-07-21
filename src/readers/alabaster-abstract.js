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
 * - `file(p)`, a (possibly async) method that accepts a string `p` containing a relative path inside an object directory and returns the contents of the file at that path.
 *   The return value should typically be a Uint8Array; on Node.js, methods may alternatively return a string containing a path to the file on the local file system.
 *
 * Optionally, the AlabasterProjectNavigator class may implement a `clear()` method to remove any cached content.
 * This will be called by {@linkcode AbstractAlabasterDataset#clear AbstractAlabasterDataset.clear} and  {@linkcode AbstractAlabasterResult#clear AbstractAlabasterResult.clear}.
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
    constructor(handle) {
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

class AlabasterFsInterface extends jsp.GlobalFsInterface {
    #navigator;

    constructor(navigator) {
        this.#navigator = navigator;
    }

    get(path, options = {}) {
        const { asBuffer = false } = options;
        return this.#navigator.file(path, asBuffer);
    }

    exists(path) {
        return this.#navigator.exists(path);
    }

    clean(path) {
        this.#navigator.clean(path); 
    }
}

class AlabasterH5Interface extends jsp.GlobalH5Interface {
    open(x, options) {
        let realized = scran.realizeFile(x);
        let output = new AlabasterH5Group(new scran.H5File(realized.path), realized.flush);
        return output;
    }

    close(handle) {
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

    numberOfRows() {
        return this.#nrow;
    }

    numberOfColumns() {
        return this.#ncol;
    }

    realize(globals, forceInteger) {
        let metadata = jsp.readObjectFile(this.#path + "/OBJECT", globals);
        if (metadata.type == "delayed_array") {
            let contents = await globals.fs.get(path + "/array.h5");
            try {
                let handle = await globals.h5.open(contents);
                try {
                    let output = extract_logcounts(handle, path, globals);
                    if (output == null) {
                        throw new Error("currently only supporting bakana-generated log-counts for delayed arrays");
                    }
                    return output;
                } finally {
                    await realized.close();
                }
            } finally {
                await globals.fs.clean(contents);
            }
        } else {
            return extract_matrix(path, metadata, globals, forceInteger);
        }
    }
}

async extract_matrix(path, metadata, globals, forceInteger) {
    if (metadata.type == "compressed_sparse_matrix") {
        let contents = globals.fs.get(path + "/matrix.h5");
        try {
            let realized = scran.realizeFile(contents);
            try {
                return scran.initializeSparseMatrixFromHdf5(realized.path, name, { forceInteger });
            } finally {
                realized.flush();
            }
        } finally {
            globals.fs.clean(contents);
        }

    } else if (metadata.type == "dense_array") {
        let contents = globals.fs.get(path + "/array.h5");
        try {
            let realized = scran.realizeFile(contents);
            try {
                return scran.initializeDenseMatrixFromHdf5(realized.path, name, { forceInteger });
            } finally {
                realized.flush();
            }
        } finally {
            globals.fs.clean(contents);
        }

    } else {
        throw new Error("unknown matrix type '" + metadata.type + "'");
    }
    
    return null;
}

// This specifically loads the log-counts created by the dumper.
// At some point we may provide more general support for delayed operations, but this would require appropriate support in scran.js itself.
async function extract_logcounts(handle, path, navigator) {
    if (handle.readAttribute("delayed_type").values[0] !== "operation") {
        return null;
    }
    if (handle.readAttribute("delayed_operation").values[0] !== "unary arithmetic") {
        return null;
    }
    if (Math.abs(handle.open("value").values[0] - Math.log(2)) > 0.00000001) {
        return null;
    }
    if (handle.open("method").values[0] !== "/") {
        return null;
    }
    if (handle.open("side").values[0] !== "right") {
        return null;
    }

    let ghandle2 = handle.open("seed");
    if (ghandle2.readAttribute("delayed_type").values[0] !== "operation") {
        return null;
    }
    if (ghandle2.readAttribute("delayed_operation").values[0] !== "unary math") {
        return null;
    }
    if (ghandle2.open("method").values[0] !== "log1p") {
        return null;
    }

    let ghandle3 = ghandle2.open("seed");
    if (ghandle3.readAttribute("delayed_type").values[0] !== "operation") {
        return null;
    }
    if (ghandle3.readAttribute("delayed_operation").values[0] !== "unary arithmetic") {
        return null;
    }
    if (ghandle3.open("method").values[0] !== "/") {
        return null;
    }
    if (ghandle3.open("side").values[0] !== "right") {
        return null;
    }
    if (ghandle3.open("along").values[0] !== 1) {
        return null;
    }
    let sf = ghandle3.open("value").values;

    let ahandle = ghandle3.open("seed");
    if (ahandle.readAttribute("delayed_type").values[0] !== "array") {
        return null;
    }
    if (ahandle.readAttribute("delayed_array").values[0] !== "custom takane seed array") {
        return null;
    }
    let index = ahandle.open("index").values[0];

    let mat;
    let output;
    try {
        mat = await extract_matrix(path + "/seeds/" + String(index), navigator, false); // don't force it to be integer, but we don't mind if it is.
        output = scran.logNormCounts(mat, { sizeFactors: sf, center: false });
    } finally {
        scran.free(mat);
    }

    return output;
}

function readMockAssay(nrow, ncol, path, metadata, globals, options) {
    return new MockMatrix(nrow, ncol, path);
}

/*************************
 *** Utility functions ***
 *************************/

function extract_from_experiments(se, fun) {
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
    return extract_from_experiments(se, x => x.rowData());
}

function extract_all_assay_names(se) {
    return extract_from_experiments(se, x => x.assayNames());
}

/************************
 ******* Dataset ********
 ************************/

/**
 * Dataset stored as a SummarizedExperiment in the **alabaster** format.
 * This is intended as a virtual base class; applications should define subclasses that are tied to a specific {@linkplain ArtifactdbProjectNavigator} class.
 * Subclasses should define `abbreviate()` and `serialize()` methods, as well as the static `format()` and `unserialize()` methods - 
 * see the [Dataset contract](https://github.com/LTLA/bakana/blob/master/docs/related/custom_readers.md) for more details.
 */
export class AbstractArtifactdbDataset {
    #path;
    #navigator;
    #se;
    #options;

    /**
     * @param {string} path - Path to the SummarizedExperiment in the ArtifactDB project directory.
     * @param {ArtifactdbProjectNavigator} navigator - A navigator object that describes how to obtain the various assets from the project directory containing `path`.
     */
    constructor(path, navigator) {
        this.#path = path;
        this.#navigator = navigator;
        this.#options = AbstractArtifactdbDataset.defaults();
        this.#raw_se = null;
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
        this.#se = null;
    }

    async #populate() {
        if (this.#se === null) {
            this.#se = await jsp.readObject(
                "",
                null, 
                {
                    fs: new AlabasterFsInterface(this.#navigator),
                    h5: new AlabasterH5Interface
                },
                {
                    DataFrame_readNested: false,
                    SummarizedExperiment_readAssay: readMockAssay,
                    SingleCellExperiment_readReducedDimension: false
                }
            );
        }
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
        await this.#populate();

        let output = {
            modality_features: extract_all_features(this.#se),
            cells: this.#se.columnData(),
            modality_assay_names: extract_all_assay_names(this.#se)
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
        await this.#populate();

        let fmapping = {
            RNA: this.#options.rnaExperiment, 
            ADT: this.#options.adtExperiment, 
            CRISPR: this.#options.crisprExperiment 
        };

        let raw_features = extract_all_features(this.#se);
        let preview = futils.extractRemappedPrimaryIds(raw_features, fmapping, this.#primary_mapping());

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
        await this.#features();
        await this.#cells();

        let output = { 
            matrix: new scran.MultiMatrix,
            features: {},
            cells: this.#raw_cells
        };

        let mapping = { 
            RNA: { exp: this.#options.rnaExperiment, assay: this.#options.rnaCountAssay },
            ADT: { exp: this.#options.adtExperiment, assay: this.#options.adtCountAssay },
            CRISPR: { exp: this.#options.crisprExperiment, assay: this.#options.crisprCountAssay }
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
                        meta = await this.#navigator.metadata(altmap[v.exp]);
                    }
                } else {
                    if (v.exp >= alts.length) {
                        continue;
                    }
                    name = alts[v.exp].name;
                    meta = await this.#navigator.metadata(alts[v.exp].resource.path);
                }

                let loaded = await extract_assay(meta, v.assay, this.#navigator, true);
                output.matrix.add(k, loaded);
                output.features[k] = this.#raw_features[name]; 
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
 * This is intended as a virtual base class; applications should define subclasses that are tied to a specific {@linkplain ArtifactdbProjectNavigator} class.
 */
export class AbstractArtifactdbResult {
    #path;
    #navigator;

    #raw_features;
    #raw_cells;
    #raw_other;

    #options;

    /**
     * @param {string} path - Path to the SummarizedExperiment in the ArtifactDB project directory.
     * @param {ArtifactdbProjectNavigator} navigator - A navigator object that describes how to obtain the various assets from the project directory containing `path`.
     */
    constructor(path, navigator) {
        this.#path = path;
        this.#navigator = new MetadataCacheWrapper(navigator);
        this.#options = AbstractArtifactdbResult.defaults();

        // Don't call clear() here, see comments above in the Dataset constructor.
        this.#reset_local_caches();
    }

    /**
     * @return {object} Default options, see {@linkcode AbstractArtifactdbResults#setOptions setOptions} for more details.
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
     * @param {object} options - Optional parameters that affect {@linkcode AbstractArtifactdbResult#load load} (but not {@linkcode AbstractArtifactdbResult#summary summary}.
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

    #reset_local_caches() {
        this.#raw_features = null;
        this.#raw_cells = null;
        this.#raw_other = null;
    }

    /**
     * Destroy caches if present, releasing the associated memory.
     * This may be called at any time but only has an effect if `cache = true` in {@linkcode AbstractArtifactdbResult#load load} or {@linkcode AbstractArtifactdbResult#summary summary}.
     */
    clear() {
        this.#reset_local_caches();
        this.#navigator.clear();
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

    async #other() {
        if (this.#raw_other !== null) {
            return;
        }

        let full_meta = await this.#navigator.metadata(this.#path);
        if ("other_data" in full_meta.summarized_experiment) {
            let other_path = full_meta.summarized_experiment.other_data.resource.path;
            this.#raw_other = await extract_other_data(other_path, this.#navigator);
        } else {
            this.#raw_other = {};
        }
        return;
    }

    async #get_all_reddim_names(rd_meta, store) {
        for (const red of rd_meta) {
            let redmeta = await this.#navigator.metadata(red.resource.path);
            if (redmeta["$schema"].startsWith("hdf5_dense_array/") && redmeta.array.dimensions.length == 2) {
                store.push(red.name);
            }
        }
    }

    /**
     * @param {object} [options={}] - Optional parameters.
     * @param {boolean} [options.cache=false] - Whether to cache the results for re-use in subsequent calls to this method or {@linkcode AbstractArtifactdbResult#load load}.
     * If `true`, users should consider calling {@linkcode AbstractArtifactdbResult#clear clear} to release the memory once this dataset instance is no longer needed.
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
        await this.#features();
        await this.#cells();
        await this.#other();

        let output = {
            modality_features: this.#raw_features,
            cells: this.#raw_cells,
            modality_assay_names: await extract_all_assay_names(this.#path, this.#navigator),
            reduced_dimension_names: [],
            other_metadata: this.#raw_other
        };

        let full_meta = await this.#navigator.metadata(this.#path);
        if ("single_cell_experiment" in full_meta) {
            let reddim_meta = full_meta.single_cell_experiment.reduced_dimensions;
            await this.#get_all_reddim_names(reddim_meta, output.reduced_dimension_names);
        }

        if (!cache) {
            this.clear();
        }
        return output;
    }

    /**
     * @param {object} [options={}] - Optional parameters.
     * @param {boolean} [options.cache=false] - Whether to cache the results for re-use in subsequent calls to this method or {@linkcode AbstractArtifactdbResult#summary summary}.
     * If `true`, users should consider calling {@linkcode AbstractArtifactdbResult#clear clear} to release the memory once this dataset instance is no longer needed.
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
        await this.#features();
        await this.#cells();
        await this.#other();

        let full_meta = await this.#navigator.metadata(this.#path);

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
            let reddim_meta = full_meta.single_cell_experiment.reduced_dimensions;

            if (reddims == null) {
                reddims = [];
                await this.#get_all_reddim_names(reddim_meta, reddims);
            }

            if (reddims.length > 0) {
                let redmap = {};
                for (const red of reddim_meta) {
                    redmap[red.name] = red.resource.path;
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

                    let meta;
                    if (k === "") {
                        meta = full_meta;
                    } else {
                        meta = await this.#navigator.metadata(altmap[k]);
                    }

                    let loaded = await extract_assay(meta, curassay, this.#navigator, !curnormalized);
                    output.matrix.add(k, loaded);

                    if (!curnormalized) {
                        let normed = scran.logNormCounts(loaded, { allowZeros: true });
                        output.matrix.add(k, normed);
                    }

                    output.features[k] = this.#raw_features[k];
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
