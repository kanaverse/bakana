import * as scran from "scran.js";
import * as bioc from "bioconductor";
import * as afile from "./abstract/file.js";
import * as eutils from "./utils/extract.js";
import * as futils from "./utils/features.js";

/**
 * Dataset in the 10X Matrix Market format, see [here](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/advanced/matrices) for details.
 */
export class TenxMatrixMarketDataset {
    #matrix_file;
    #feature_file;
    #barcode_file;

    #dimensions;
    #raw_features;
    #raw_cells;

    #options;

    #dump_summary(fun) {
        let files = [{ type: "mtx", file: fun(this.#matrix_file) }];

        if (this.#feature_file !== null) {
            files.push({ type: "genes", file: fun(this.#feature_file) });
        }

        if (this.#barcode_file !== null) {
            files.push({ type: "annotations", file: fun(this.#barcode_file) });
        }

        let options = this.options(); 
        return { files, options };
    }

    /**
     * @param {SimpleFile|string|Uint8Array|File} matrixFile - A Matrix Market file.
     * On browsers, this may be a File object.
     * On Node.js, this may also be a string containing a file path.
     * @param {?(SimpleFile|string|Uint8Array|File)} featureFile - Contents of a feature annotation file.
     * If `null`, it is assumed that no file was available.
     * @param {?(SimpleFile|string|Uint8Array|File)} barcodeFile - Contents of a barcode annotation file.
     * If `null`, it is assumed that no file was available.
     */
    constructor(matrixFile, featureFile, barcodeFile, { 
        featureTypeRnaName = "Gene Expression", 
        featureTypeAdtName = "Antibody Capture", 
        featureTypeCrisprName = "CRISPR Guide Capture", 
        primaryRnaFeatureIdColumn = 0, 
        primaryAdtFeatureIdColumn = 0,
        primaryCrisprFeatureIdColumn = 0
    } = {}) {
        if (matrixFile instanceof afile.SimpleFile) {
            this.#matrix_file = matrixFile;
        } else {
            this.#matrix_file = new afile.SimpleFile(matrixFile);
        }

        if (featureFile instanceof afile.SimpleFile || featureFile == null) {
            this.#feature_file = featureFile;
        } else {
            this.#feature_file = new afile.SimpleFile(featureFile);
        }

        if (barcodeFile instanceof afile.SimpleFile || barcodeFile == null) {
            this.#barcode_file = barcodeFile;
        } else {
            this.#barcode_file = new afile.SimpleFile(barcodeFile);
        }

        this.#options = TenxMatrixMarketDataset.defaults();
        this.clear();
    }

    /**
     * @return {object} Default options, see {@linkcode TenxMatrixMarketDataset#setOptions setOptions} for more details.
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
     * @param {object} options - Optional parameters that affect {@linkcode TenxMatrixMarketDataset#load load} (but not {@linkcode TenxMatrixMarketDataset#summary summary}.
     * @param {?string} [options.featureTypeRnaName] - Name of the feature type for gene expression.
     * Alternatively `null`, to indicate that no RNA features are to be loaded.
     * @param {?string} [options.featureTypeAdtName] - Name of the feature type for ADTs.
     * Alternatively `null`, to indicate that no ADT features are to be loaded.
     * @param {?string} [options.featureTypeCrisprName] - Name of the feature type for CRISPR guides.
     * Alternatively `null`, to indicate that no guides are to be loaded.
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

    /**
     * Destroy caches if present, releasing the associated memory.
     * This may be called at any time but only has an effect if `cache = true` in {@linkcode TenxMatrixMarketDataset#load load} or {@linkcodeTenxMatrixMarketDataset#summary summary}. 
     */
    clear() {
        this.#dimensions = null;
        this.#raw_features = null;
        this.#raw_cells = null;
    }

    /**
     * @return {string} Format of this dataset class.
     * @static
     */
    static format() {
        return "MatrixMarket";
    }

    /**
     * @return {object} Object containing the abbreviated details of this dataset,
     * in a form that can be cheaply stringified.
     */
    abbreviate(args) {
        return this.#dump_summary(f => { return { name: f.name(), size: f.size() }; });
    }

    #fetch_dimensions() {
        if (this.#dimensions !== null) {
            return;
        }
        var is_gz = this.#matrix_file.name().endsWith(".gz");
        let headers = scran.extractMatrixMarketDimensions(this.#matrix_file.content(), { "compressed": is_gz });
        this.#dimensions = [headers.rows, headers.columns];
    }

    async #features() {
        if (this.#raw_features !== null) {
            return;
        }

        this.#fetch_dimensions();
        let NR = this.#dimensions[0];
        if (this.#feature_file == null) {
            this.#raw_features = new bioc.DataFrame({}, { numberOfRows: NR });
            return;
        }

        let fname = this.#feature_file.name();
        var is_gz = fname.endsWith(".gz");
        let parsed = await eutils.readTable2(this.#feature_file.content(), { compression: (is_gz ? "gz" : "none") });

        if (parsed.length == NR + 1) {
            // If it seems to have a header, we just use that directly.
            let output = {};
            let headers = parsed.shift();
            headers.forEach((x, i) => {
                output[x] = parsed.map(y => y[i]);
            });
            this.#raw_features = output;
            return;
        }

        // Otherwise, we assume it's standard 10X CellRanger output, without a header.
        if (parsed.length !== NR) {
            throw new Error("number of matrix rows is not equal to the number of rows in '" + fname + "'");
        } 

        var ids = [], symb = [];
        parsed.forEach(x => {
            ids.push(x[0]);
            symb.push(x[1]);
        });

        let output = new bioc.DataFrame({}, { numberOfRows: NR }); // build it piece-by-piece for a well-defined order.
        output.$setColumn("id", ids);
        output.$setColumn("name", symb);

        if (parsed[0].length > 2) {
            let types = [];
            parsed.forEach(x => { types.push(x[2]); });
            output.$setColumn("type", types);
        }

        this.#raw_features = output;
        return;
    }

    async #cells() {
        if (this.#raw_cells !== null) {
            return;
        }

        this.#fetch_dimensions();
        if (this.#barcode_file == null) {
            this.#raw_cells = new bioc.DataFrame({}, { numberOfRows: this.#dimensions[1] });
            return;
        }

        let bname = this.#barcode_file.name();
        var is_gz = bname.endsWith(".gz");
        let parsed = await eutils.readTable2(this.#barcode_file.content(), { compression: (is_gz ? "gz" : "none") });

        // Check if a header is present or not. Standard 10X output doesn't have a 
        // header but we'd like to support some kind of customization.
        let diff = this.#dimensions[1] - parsed.length;
        let headers;
        if (diff == 0) {
            headers = parsed[0]; // whatever, just using the first row. Hope it's unique enough!
        } else if (diff == -1) {
            headers = parsed.shift();
        } else {
            throw new Error("number of matrix columns is not equal to the number of rows in '" + bname + "'");
        }

        let annotations = {}
        headers.forEach((x, i) => {
            annotations[x] = parsed.map(y => y[i]);
        });

        for (const [k, v] of Object.entries(annotations)) {
            let conv = eutils.promoteToNumber(v);
            if (conv !== null) {
                annotations[k] = conv;
            }
        }

        this.#raw_cells = new bioc.DataFrame(annotations);
        return;
    }

    /**
     * @param {object} [options={}] - Optional parameters.
     * @param {boolean} [options.cache=false] - Whether to cache the results for re-use in subsequent calls to this method or {@linkcode TenxMatrixMarketDataset#load load}.
     * If `true`, users should consider calling {@linkcode TenxMatrixMarketDataset#clear clear} to release the memory once this dataset instance is no longer needed.
     *
     * @return {object} Object containing the per-feature and per-cell annotations.
     * This has the following properties:
     *
     * - `modality_features`: an object where each key is a modality name and each value is a {@linkplain external:DataFrame DataFrame} of per-feature annotations for that modality.
     *   Unlike {@linkcode TenxMatrixMarketDataset#load load}, modality names are arbitrary.
     * - `cells`: a {@linkplain external:DataFrame DataFrame} of per-cell annotations.
     *
     * @async
     */
    async summary({ cache = false } = {}) {
        await this.#features();
        await this.#cells();

        let output = {
            "modality_features": futils.reportFeatures(this.#raw_features, "type"),
            "cells": this.#raw_cells
        };

        if (!cache) {
            this.clear();
        }
        return output;
    }

    /**
     * @param {object} [options={}] - Optional parameters.
     * @param {boolean} [options.cache=false] - Whether to cache the results for re-use in subsequent calls to this method or {@linkcode TenxMatrixMarketDataset#summary summary}.
     * If `true`, users should consider calling {@linkcode TenxMatrixMarketDataset#clear clear} to release the memory once this dataset instance is no longer needed.
     *
     * @return {object} Object containing the per-feature and per-cell annotations.
     * This has the following properties:
     *
     * - `features`: an object where each key is a modality name and each value is a {@linkplain external:DataFrame DataFrame} of per-feature annotations for that modality.
     * - `cells`: a {@linkplain external:DataFrame DataFrame} containing per-cell annotations.
     * - `matrix`: a {@linkplain external:MultiMatrix MultiMatrix} containing one {@linkplain external:ScranMatrix ScranMatrix} per modality.
     * - `primary_ids`: an object where each key is a modality name and each value is an array of strings containing the feature identifiers for each row in that modality.
     *
     * Modality names are guaranteed to be one of `"RNA"` or `"ADT"`.
     * It is assumed that an appropriate mapping from the feature types inside the `featureFile` was previously declared,
     * either in the constructor or in {@linkcode setFeatureTypeRnaName} and {@linkcode setFeatureTypeAdtName}.
     *
     * If the feature annotation file is absent or if it lacks a third column for the feature types,
     * `load()` will assume that all features are genes (i.e., only the RNA modality is present).
     *
     * @async
     */
    async load({ cache = false } = {}) {
        await this.#features();
        await this.#cells();

        var is_gz = this.#matrix_file.name().endsWith(".gz");
        let loaded = scran.initializeSparseMatrixFromMatrixMarket(this.#matrix_file.content(), { "compressed": is_gz });

        let mappings = { 
            RNA: this.#options.featureTypeRnaName, 
            ADT: this.#options.featureTypeAdtName,
            CRISPR: this.#options.featureTypeCrisprName
        };
        let output = futils.splitScranMatrixAndFeatures(loaded, this.#raw_features, "type", mappings, "RNA"); 
        output.cells = this.#raw_cells;

        let primaries = { 
            RNA: this.#options.primaryRnaFeatureIdColumn, 
            ADT: this.#options.primaryAdtFeatureIdColumn,
            CRISPR: this.#options.primaryCrisprFeatureIdColumn
        };
        output.primary_ids = futils.extractPrimaryIds(output.features, primaries);

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
    async serialize() {
        return this.#dump_summary(f => f);
    }

    /**
     * @param {Array} files - Array of objects like that produced by {@linkcode TenxMatrixMarketDataset#serialize serialize}.
     * @param {object} options - Object containing additional options to be passed to the constructor.
     * @return {TenxMatrixMarketDataset} A new instance of this class.
     * @static
     */
    static async unserialize(files, options) {
        let args = {};
        for (const x of files) {
            if (x.type in args) {
                throw new Error("duplicate file of type '" + x.type + "' detected during MatrixMarket unserialization");
            }
            args[x.type] = x.file;
        }

        if (!("mtx" in args)) {
            throw new Error("expected file of type 'mtx' for during MatrixMarket unserialization");
        }

        let feat = null;
        if ("genes" in args) {
            feat = args.genes;
        }

        let barcode = null;
        if ("annotations" in args) {
            barcode = args.annotations;
        }

        return new TenxMatrixMarketDataset(args.mtx, feat, barcode, options);
    }
}
