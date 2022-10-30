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

    /**
     * @param {SimpleFile|string|Uint8Array|File} matrixFile - A Matrix Market file.
     * On browsers, this may be a File object.
     * On Node.js, this may also be a string containing a file path.
     * @param {?(SimpleFile|string|Uint8Array|File)} featureFile - Contents of a feature annotation file.
     * If `null`, it is assumed that no file was available.
     * @param {?(SimpleFile|string|Uint8Array|File)} barcodeFile - Contents of a barcode annotation file.
     * If `null`, it is assumed that no file was available.
     */
    constructor(matrixFile, featureFile, barcodeFile) {
        super();

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

        this.#dimensions = null;
        this.#raw_features = null;
        this.#raw_cells = null;
    }

    static format() {
        return "MatrixMarket";
    }

    abbreviate(args) {
        var formatted = {
            "matrix": {
                name: this.#matrix_file.name(),
                size: this.#matrix_file.size()
            }
        };

        if (this.#feature_file !== null) {
            formatted.features = {
                name: this.#feature_file.name(),
                size: this.#feature_file.size()
            };
        }

        if (this.#barcode_file !== null) {
            formatted.barcodes = {
                name: this.#barcode_file.name(),
                size: this.#barcode_file.size()
            };
        }

        return formatted;
    }

    #fetch_dimensions() {
        if (this.#dimensions !== null) {
            return;
        }
        var is_gz = this.#matrix_file.name().endsWith(".gz");
        let headers = scran.extractMatrixMarketDimensions(this.#matrix_file.content(), { "compressed": is_gz });
        this.#dimensions = [headers.rows, headers.columns];
    }

    #features() {
        if (this.#raw_features !== null) {
            return;
        }

        this.#fetch_dimensions();
        let NR = this.#dimensions[0];
        if (this.#feature_file == null) {
            this.#raw_features = new bioc.DataFrame({}, { numberOfRows: NR });
            return;
        }

        const content = new Uint8Array(this.#feature_file.buffer());
        let fname = this.#feature_file.name();
        var is_gz = fname.endsWith(".gz");
        let parsed = eutils.readTable(content, { compression: (is_gz ? "gz" : "none") });

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
        let output = { "id": ids, "symbol": symb };

        if (parsed[0].length >= 3) {
            let types = [];
            parsed.forEach(x => { types.push(x[2]); });
            output.type = types;
        }

        this.#raw_features = new bioc.DataFrame(output);
        return;
    }

    #cells() {
        if (this.#raw_cells !== null) {
            return;
        }

        this.#fetch_dimensions();
        if (this.#barcode_file == null) {
            this.#raw_cells = new bioc.DataFrame({}, { numberOfRows: this.#dimensions[1] });
            return;
        }

        const content = new Uint8Array(this.#barcode_file.buffer());
        let bname = this.#barcode_file.name();
        var is_gz = bname.endsWith(".gz");
        let parsed = eutils.readTable(content, { compression: (is_gz ? "gz" : "none") });

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

    annotations() {
        this.#features();
        this.#cells();

        let anno = {};
        for (const k of this.#raw_cells.columnNames()) {
            anno[k] = eutils.summarizeArray(this.#raw_cells.column(k));
        }

        return {
            "features": futils.reportFeatures(this.#raw_features, "type"),
            "cells": {
                "number": this.#raw_cells.numberOfRows(),
                "summary": anno
            }
        };
    }

    load() {
        this.#features();
        this.#cells();

        var is_gz = this.#matrix_file.name().endsWith(".gz");
        let loaded = scran.initializeSparseMatrixFromMatrixMarket(this.#matrix_file.content(), { "compressed": is_gz });
        let output = futils.splitScranMatrixAndFeatures(loaded, this.#raw_features, "type");
        output.cells = bioc.CLONE(this.#raw_cells);
        return output;
    }

    // These 'type's are largely retained for back-compatibility.
    async serialize() {
        let files = [ { type: "mtx", file: this.#matrix_file } ];

        if (this.#feature_file !== null) {
            files.push({ type: "genes", file: this.#feature_file });
        }

        if (this.#barcode_file !== null) {
            files.push({ type: "annotations", file: this.#barcode_file });
        }

        return files;
    }

    static async unserialize(files) {
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

        return new TenxMatrixMarketDataset(args.mtx, feat, barcode);
    }
}
