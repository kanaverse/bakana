import * as scran from "scran.js";
import * as bioc from "bioconductor";
import * as afile from "./abstract/file.js";
import * as eutils from "./utils/extract.js";
import * as futils from "./utils/features.js";

/**
 * Dataset in the 10X HDF5 feature-barcode matrix format, see [here](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/output/matrices) for details.
 */
export class TenxHdf5Dataset { 
    #h5_file;
    #h5_path;
    #h5_flush;

    #raw_features;
    #raw_cells;
    #raw_shape;

    #options;

    #dump_summary(fun) {
        let files = [{ type: "h5", file: fun(this.#h5_file) }];
        let options = this.options();
        return { files, options };
    }

    /**
     * @param {SimpleFile|string|Uint8Array|File} h5File - Contents of a HDF5 file in the 10X feature-barcode format.
     * On browsers, this may be a File object.
     * On Node.js, this may also be a string containing a file path.
     */
    constructor(h5File) {
        if (h5File instanceof afile.SimpleFile) {
            this.#h5_file = h5File;
        } else {
            this.#h5_file = new afile.SimpleFile(h5File);
        }

        this.#options = TenxHdf5Dataset.defaults();
        this.clear();
    }

    /**
     * @return {object} Default options, see {@linkcode TenxHdf5Dataset#setOptions setOptions} for more details.
     */
    static defaults() {
        return {
            featureTypeRnaName: "Gene Expression", 
            featureTypeAdtName: "Antibody Capture", 
            featureTypeCrisprName: "CRISPR Guide Capture", 
            primaryRnaFeatureIdColumn: 0, 
            primaryAdtFeatureIdColumn: 0,
            primaryCrisprFeatureIdColumn: 0
        };
    }

    /**
     * @return {object} Object containing all options used for loading.
     */
    options() {
        return { ...(this.#options) };
    }

    /**
     * @param {object} options - Optional parameters that affect {@linkcode TenxHdf5Dataset#load load} (but not {@linkcode TenxHdf5Dataset#summary summary}).
     * @param {?string} [options.featureTypeRnaName] - Name of the feature type for gene expression.
     * If `null` or the string is not present among the feature types, no RNA features are to be loaded.
     *
     * If no feature type information is available in the dataset, all features are considered to be genes by default.
     * This behavior can also be explicitly requested by setting this argument to the only non-`null` value among all `featureType*Name` parameters.
     * @param {?string} [options.featureTypeAdtName] - Name of the feature type for ADTs.
     * If `null` or the string is not present among the feature types, no ADT features are to be loaded.
     *
     * If no feature type information is available in the dataset and this argument is set to the only non-`null` value among all `featureType*Name` parameters, all features are considered to be ADTs.
     * @param {?string} [options.featureTypeCrisprName] - Name of the feature type for CRISPR guides.
     * If `null` or the string is not present among the feature types, no guides are to be loaded.
     * 
     * If no feature type information is available in the dataset and this argument is set to the only non-`null` value among all `featureType*Name` parameters, all features are considered to be guides.
     * @param {string|number} [options.primaryRnaFeatureIdColumn] - Name or index of the column of the `features` {@linkplain external:DataFrame DataFrame} that contains the primary feature identifier for gene expression.
     * If `i` is invalid (e.g., out of range index, unavailable name), it is ignored and the primary identifier is treated as undefined.
     * @param {string|number} [options.primaryAdtFeatureIdColumn] - Name or index of the column of the `features` {@linkplain external:DataFrame DataFrame} that contains the primary feature identifier for the ADTs.
     * If `i` is invalid (e.g., out of range index, unavailable name), it is ignored and the primary identifier is treated as undefined.
     * @param {string|number} [options.primaryCrisprFeatureIdColumn] - Name or index of the column of the `features` {@linkplain external:DataFrame DataFrame} that contains the primary feature identifier for the CRISPR guides.
     * If `i` is invalid (e.g., out of range index, unavailable name), it is ignored and the primary identifier is treated as undefined.
     */
    setOptions(options) {
        for (const [k, v] of Object.entries(options)) {
            this.#options[k] = v;
        }
    }

    #instantiate() {
        if (this.#h5_path !== null) {
            return;
        }

        let info = scran.realizeFile(this.#h5_file.content());
        this.#h5_path = info.path;
        this.#h5_flush = info.flush;
    }

    /**
     * Destroy caches if present, releasing the associated memory.
     * This may be called at any time but only has an effect if `cache = true` in {@linkcode TenxHdf5Dataset#load load} or {@linkcodeTenxHdf5Dataset#summary summary}. 
     */
    clear() {
        if (typeof this.#h5_flush == "function") {
            this.#h5_flush();
        }
        this.#h5_flush = null;
        this.#h5_path = null;

        this.#raw_features = null;
        this.#raw_cells = null;
    }

    /**
     * @return {string} Format of this dataset class.
     * @static
     */
    static format() {
        return "10X";
    }

    /**
     * @return {object} Object containing the abbreviated details of this dataset,
     * in a form that can be cheaply stringified.
     */
    abbreviate() {
        return this.#dump_summary(f => { return { name: f.name(), size: f.size() }; });
    }

    #features() {
        if (this.#raw_features !== null) {
            return;
        }

        this.#instantiate();
        let handle = new scran.H5File(this.#h5_path);
        if (!("matrix" in handle.children) || handle.children["matrix"] != "Group") {
            throw new Error("expected a 'matrix' group at the top level of the file");
        }
        let mhandle = handle.open("matrix");

        if (!("features" in mhandle.children) || mhandle.children["features"] != "Group") {
            throw new Error("expected a 'matrix/features' group containing the feature annotation");
        }
        let fhandle = mhandle.open("features");

        let ids = eutils.extractHdf5Strings(fhandle, "id");
        if (ids == null) {
            throw new Error("expected a 'matrix/features/id' string dataset containing the feature IDs");
        }
        let feats = new bioc.DataFrame({ id: ids }); // build it piece-by-piece for a well-defined ordering.

        let names = eutils.extractHdf5Strings(fhandle, "name");
        if (names !== null) {
            feats.$setColumn("name", names);
        }

        let ftype = eutils.extractHdf5Strings(fhandle, "feature_type");
        if (ftype !== null) {
            feats.$setColumn("type", ftype);
        }

        this.#raw_features = feats;
        return;
    }

    #cells() {
        if (this.#raw_cells !== null) {
            return;
        }

        this.#instantiate();

        let fhandle = new scran.H5File(this.#h5_path);
        let dhandle = fhandle.open("matrix");
        let shandle = dhandle.open("shape");
        let shape = shandle.values;
        if (shape.length != 2 || !shape.every(x => (typeof x === "number" && x >= 0 && Number.isInteger(x)))) {
            throw new Error("expected 'shape' to contain 2 non-negative integers");
        }
        this.#raw_shape = shandle.values;

        this.#raw_cells = new bioc.DataFrame({}, { numberOfRows: this.#raw_shape[1] });
    }

    /**
     * @param {object} [options={}] - Optional parameters.
     * @param {boolean} [options.cache=false] - Whether to cache the intermediate results for re-use in subsequent calls to any methods with a `cache` option.
     * If `true`, users should consider calling {@linkcode TenxHdf5Dataset#clear clear} to release the memory once this dataset instance is no longer needed.
     *
     * @return {object} Object containing the per-feature and per-cell annotations.
     * This has the following properties:
     *
     * - `modality_features`: an object where each key is a modality name and each value is a {@linkplain external:DataFrame DataFrame} of per-feature annotations for that modality.
     *   Unlike {@linkcode TenxMatrixMarketDataset#load load}, modality names are arbitrary.
     * - `cells`: a {@linkplain external:DataFrame DataFrame} of per-cell annotations.
     */
    summary({ cache = false } = {}) {
        this.#features();
        this.#cells();

        let output = {
            "modality_features": futils.reportFeatures(this.#raw_features, "type"),
            "cells": this.#raw_cells
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
     * If `true`, users should consider calling {@linkcode TenxHdf5Dataset#clear clear} to release the memory once this dataset instance is no longer needed.
     *
     * @return {object} An object where each key is a modality name and each value is an array (usually of strings) containing the primary feature identifiers for each row in that modality.
     * The contents are the same as the `primary_ids` returned by {@linkcode TenxHdf5Dataset#load load} but the order of values may be different.
     */
    previewPrimaryIds({ cache = false } = {}) {
        this.#features();
        let preview = futils.extractSplitPrimaryIds(this.#raw_features, "type", this.#feature_type_mapping(), "RNA", this.#primary_mapping());
        if (!cache) {
            this.clear();
        }
        return preview;
    }

    /**
     * @param {object} [options={}] - Optional parameters.
     * @param {boolean} [options.cache=false] - Whether to cache the intermediate results for re-use in subsequent calls to any methods with a `cache` option.
     * If `true`, users should consider calling {@linkcode TenxHdf5Dataset#clear clear} to release the memory once this dataset instance is no longer needed.
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
     * either from the {@linkcode TenxHdf5Dataset#defaults defaults} or with {@linkcode TenxHdf5Dataset#setOptions setOptions}.
     *
     * If the feature annotation lacks information about the feature types, it is assumed that all features are genes, i.e., only the RNA modality is present.
     */
    load({ cache = false } = {}) {
        this.#features();
        this.#cells();

        let loaded = scran.initializeSparseMatrixFromHdf5Group(
            this.#h5_path,
            "matrix",
            this.#raw_shape[0],
            this.#raw_shape[1],
            /* byRow = */ false,
            { forceInteger: true, layered: true }
        ); // collection gets handled inside splitScranMatrixAndFeatures.

        let output = futils.splitScranMatrixAndFeatures(loaded, this.#raw_features, "type", this.#feature_type_mapping(), "RNA");
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
     * @param {Array} files - Array of objects like that produced by {@linkcode TenxHdf5Dataset#serialize serialize}.
     * @param {object} options - Object containing additional options to be passed to the constructor.
     * @return {TenxHdf5Dataset} A new instance of this class.
     * @static
     */
    static async unserialize(files, options) {
        if (files.length != 1 || files[0].type != "h5") {
            throw new Error("expected exactly one file of type 'h5' for 10X HDF5 unserialization");
        }
        let output = new TenxHdf5Dataset(files[0].file);
        output.setOptions(output);
        return output;
    }
}
