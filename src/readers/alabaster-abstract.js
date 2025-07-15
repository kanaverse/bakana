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

// This specifically loads the log-counts created by the dumper.
async function extract_logcounts(handle, navigator) {
    if (handle.readAttribute("delayed_type").values[0] !== "operation") {
        return null;
    }
    if (handle.readAttribute("delayed_operation").values[0] !== "unary arithmetic") {
        return null;
    }
    if (Math.abs(handle.open("value", { load: true }).values[0] - Math.log(2)) > 0.00000001) {
        return null;
    }
    if (handle.open("method", { load: true }).values[0] !== "/") {
        return null;
    }
    if (handle.open("side", { load: true }).values[0] !== "right") {
        return null;
    }

    let ghandle2 = handle.open("seed");
    if (ghandle2.readAttribute("delayed_type").values[0] !== "operation") {
        return null;
    }
    if (ghandle2.readAttribute("delayed_operation").values[0] !== "unary math") {
        return null;
    }
    if (ghandle2.open("method", { load: true }).values[0] !== "log1p") {
        return null;
    }

    let ghandle3 = ghandle2.open("seed");
    if (ghandle3.readAttribute("delayed_type").values[0] !== "operation") {
        return null;
    }
    if (ghandle3.readAttribute("delayed_operation").values[0] !== "unary arithmetic") {
        return null;
    }
    if (ghandle3.open("method", { load: true }).values[0] !== "/") {
        return null;
    }
    if (ghandle3.open("side", { load: true }).values[0] !== "right") {
        return null;
    }
    if (ghandle3.open("along", { load: true }).values[0] !== 1) {
        return null;
    }
    let sf = ghandle3.open("value", { load: true }).values;

    let ahandle = ghandle3.open("seed");
    if (ahandle.readAttribute("delayed_type").values[0] !== "array") {
        return null;
    }
    if (ahandle.readAttribute("delayed_array").values[0] !== "custom alabaster local array") {
        return null;
    }
    let path = ahandle.open("path", { load: true }).values[0];

    let mat;
    let output;
    try {
        mat = await extract_assay_raw(path, navigator, false); // don't force it to be integer, but we don't mind if it is.
        output = scran.logNormCounts(mat, { sizeFactors: sf, center: false });
    } finally {
        scran.free(mat);
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
            output = scran.initializeSparseMatrixFromHdf5(stuff.path, name, { forceInteger });
        } finally {
            stuff.flush();
        }

    } else if (assmeta["$schema"].startsWith("hdf5_delayed_array/")) {
        let stuff = scran.realizeFile(contents);
        try {
            let fhandle = new scran.H5File(stuff.path);
            let ghandle = fhandle.open(assmeta.hdf5_delayed_array.group);

            // TODO: replace with calls to chihaya.js.
            output = await extract_logcounts(ghandle, navigator);
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

/*************************************
 ******* Jaspalite interfaces ********
 *************************************/

class ScranH5Group extends jsp.H5Group {
    constructor(handle) {
        super();
        this.#handle = handle;
    }

    attributes() {
        return this.#handle.attributes;
    }

    readAttribute(name) {
        return this.#handle.readAttribute(name);
    }

    writeAttribute(attr, type, shape, data, options = {}) {
        this.#handle.writeAttribute(attr, type, shape, data, options);
    }

    open(name) {
        let child = this.#handle.open(name);
        if (child instanceof scran.Group) {
            return new ScranH5Group(child);
        } else if (child instanceof scran.Dataset) {
            return new ScranH5DataSet(child);
        } else {
            throw new Error("unknown object type");
        }
    }

    createGroup(name) {
        return new ScranH5Group(this.#handle.createGroup(name));
    }

    createDataSet(name, type, shape, options = {}) {
        let handle = this.#handle.createDataSet(name, type, shape, options);
        if (!("returnHandle" in options) || options.returnHandle) {
            return new ScranH5DataSet(handle);
        }
    }
}

class ScranH5DataSet extends jsp.H5DataSet {
    constructor(handle) {
        super();
        this.#handle = handle;
    }

    attributes() {
        return this.#handle.attributes;
    }

    readAttribute(name) {
        return this.#handle.readAttribute(name);
    }

    writeAttribute(attr, type, shape, data, options = {}) {
        this.#handle.writeAttribute(attr, type, shape, data, options);
    }

    type() {
        return this.#handle.type;
    }

    shape() {
        return this.#handle.shape;
    }

    values() {
        return this.#handle.values;
    }

    write(data) {
        this.#handle.write(data);
    }

    close() {
        this.#handle;
    }
}

// As we'll be using the metadata often, we cache it at this level. This
// removes the burden of caching on the implementation of the navigator. 
class MetadataCacheWrapper {
    #navigator;
    #metadata_cache;

    constructor(nav) {
        this.#navigator = nav;
        this.#metadata_cache = {};
    }

    clear() {
        this.#metadata_cache = {};
        if ("clear" in this.#navigator) {
            this.#navigator.clear();
        }
    }

    async metadata(path) {
        if (path in this.#metadata_cache) {
            return this.#metadata_cache[path];
        } else {
            let content = await this.#navigator.metadata(path);
            this.#metadata_cache[path] = content;
            return content;
        }
    }

    file(path) {
        return this.#navigator.file(path);
    }
};

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

    #raw_features;
    #raw_cells;

    #options;

    /**
     * @param {string} path - Path to the SummarizedExperiment in the ArtifactDB project directory.
     * @param {ArtifactdbProjectNavigator} navigator - A navigator object that describes how to obtain the various assets from the project directory containing `path`.
     */
    constructor(path, navigator) {
        this.#path = path;
        this.#navigator = new MetadataCacheWrapper(navigator);
        this.#options = AbstractArtifactdbDataset.defaults();

        // Don't call this.clear() here. We don't want to clear the navigator's
        // cache at this point, as the navigator might contain some cached
        // values when passed to the constructor. We should respect any caches
        // until we're specifically told to discard it with clear() or cache =
        // false in load() or summary().
        this.#reset_local_caches();
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

    #reset_local_caches() {
        this.#raw_features = null;
        this.#raw_cells = null;
    }

    /**
     * Destroy caches if present, releasing the associated memory.
     * This may be called at any time but only has an effect if `cache = true` in {@linkcode AbstractArtifactdbDataset#load load} or {@linkcode AbstractArtifactdbDataset#summary summary}.
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
        await this.#features();

        let fmapping = {
            RNA: this.#options.rnaExperiment, 
            ADT: this.#options.adtExperiment, 
            CRISPR: this.#options.crisprExperiment 
        };

        let preview = futils.extractRemappedPrimaryIds(this.#raw_features, fmapping, this.#primary_mapping());

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

export const AlabasterSummarizedExperimentDatasetBase = AbstractArtifactdbDataset;

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

export const AlabasterSummarizedExperimentResultBase = AbstractArtifactdbResult;
