import * as scran from "scran.js";
import * as bioc from "bioconductor";
import * as eutils from "./utils/extract.js";
import * as futils from "./utils/features.js";
import * as afile from "./abstract/file.js";

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
    #raw_cell_factors;

    #chosen_assay;
    #assay_details;

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

        this.clear();
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
     * This may be called at any time but only has an effect if `cache = true` in {@linkcode H5adDataset#load load} or {@linkcodeH5adDataset#annotations annotations}.
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
        this.#raw_cell_factors = null;

        this.#chosen_assay = null;
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
        return { 
            "h5": {
                "name": this.#h5_file.name(),
                "size": this.#h5_file.size()
            }
        };
    }

    #choose_assay() {
        if (this.#chosen_assay !== null) {
            return;
        }

        if ("X" in this.#h5_handle.children) {
            this.#chosen_assay = "X";
            return;
        } 

        if (!("layers" in this.#h5_handle.children)) {
            throw new Error("expected 'X' or 'layers' in a H5AD file");
        }

        let lhandle = this.#h5_handle.open("layers");
        if (!(lhandle instanceof scran.H5Group)) {
            throw new Error("expected a 'layers' group in a H5AD file");
        }

        for (const k of Object.keys(lhandle.children)) {
            if (k.match(/^count/i)) {
                this.#chosen_assay = "layers/" + k;
                return;
            }
        }

        throw new Error("failed to find any count-like layer in a H5AD file");
    }

    #fetch_assay_details() {
        if (this.#assay_details !== null) {
            return;
        }

        this.#choose_assay();
        this.#assay_details = scran.extractHDF5MatrixDetails(this.#h5_path, this.#chosen_assay);
        return;
    }

    #features() {
        if (this.#raw_features !== null) {
            return;
        }

        this.#instantiate();
        let handle = this.#h5_handle;
        if ("var" in handle.children && handle.children["var"] == "Group") {
            let vhandle = handle.open("var");
            let index = eutils.extractHDF5Strings(vhandle, "_index");
            if (index !== null) {
                let genes = { "_index": index };

                for (const [key, val] of Object.entries(vhandle.children)) {
                    if (val === "DataSet" && (key.match(/name/i) || key.match(/symb/i))) {
                        let dhandle2 = vhandle.open(key);
                        if (dhandle2.type == "String") {
                            genes[key] = dhandle2.load();
                        }
                    }
                }

                this.#raw_features = new bioc.DataFrame(genes);
                return;
            }
        }

        this.#fetch_assay_details();
        this.#raw_features = new bioc.DataFrame({}, { numberOfRows: this.#assay_details.rows });
        return;
    }

    #cells() {
        if (this.#raw_cells !== null) {
            return;
        }

        this.#instantiate();
        let handle = this.#h5_handle;
        let annotations = {};
        let factors = {};

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
                        factors[key] = eutils.extractHDF5Strings(chandle, key);
                    }
                }
            }
        }

        if (Object.keys(annotations).length == 0) {
            this.#fetch_assay_details();
            this.#raw_cells = new bioc.DataFrame({}, { numberOfRows: this.#assay_details.columns })
        } else {
            this.#raw_cells = new bioc.DataFrame(annotations);
        }
        this.#raw_cell_factors = factors;
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
     * - `features`: an object where each key is a modality name and each value is a {@linkplain external:DataFrame DataFrame} of per-feature annotations for that modality.
     * - `cells`: an object containing:
     *   - `number`: the number of cells in this dataset.
     *   - `summary`: an object where each key is the name of a per-cell annotation field and its value is a summary of that annotation field,
     *     following the same structure as returned by {@linkcode summarizeArray}.
     */
    annotations({ cache = false } = {}) {
        this.#features();
        this.#cells();

        let summaries = {};
        for (const k of this.#raw_cells.columnNames()) {
            let x;
            if (k in this.#raw_cell_factors) {
                x = this.#raw_cell_factors[k];
            } else {
                x = this.#raw_cells.column(k);
            }
            summaries[k] = eutils.summarizeArray(x);
        }

        let output = {
            features: { "": futils.cloneCached(this.#raw_features, cache) },
            cells: {
                number: this.#raw_cells.numberOfColumns(),
                summary: summaries
            }
        };

        if (!cache) {
            this.clear();
        }
        return output;
    }

    /**
     * @param {object} [options={}] - Optional parameters.
     * @param {boolean} [options.cache=false] - Whether to cache the results for re-use in subsequent calls to this method or {@linkcode H5adDataset#annotations annotations}.
     * If `true`, users should consider calling {@linkcode H5adDataset#clear clear} to release the memory once this dataset instance is no longer needed.
     *
     * @return {object} Object containing the per-feature and per-cell annotations.
     * This has the following properties:
     *
     * - `features`: an object where each key is a modality name and each value is a {@linkplain external:DataFrame DataFrame} of per-feature annotations for that modality.
     * - `cells`: a {@linkplain external:DataFrame DataFrame} containing per-cell annotations.
     * - `matrix`: a {@linkplain external:MultiMatrix MultiMatrix} containing one {@linkplain external:ScranMatrix ScranMatrix} per modality.
     * - `row_ids`: an object where each key is a modality name and each value is an integer array containing the feature identifiers for each row in that modality.
     */
    load({ cache = false } = {}) {
        this.#features();
        this.#cells();
        this.#choose_assay();

        let remap = (df, k, v) => {
            let cats = this.#raw_cell_factors[k];
            let temp = new Array(v.length);
            v.forEach((x, i) => {
                temp[i] = cats[x]
            }); 
            df.$setColumn(k, temp);
        };

        let cells;
        let src = this.#raw_cells;
        if (cache) {
            cells = new bioc.DataFrame({}, { numberOfRows: src.numberOfRows(), metadata: bioc.CLONE(src.metadata()) });
            for (const k of src.columnNames()) {
                let v = src.column(k);
                if (k in this.#raw_cell_factors) {
                    remap(cells, k, v);
                } else {
                    cells.$setColumn(k, bioc.CLONE(v)); // Making explicit copies to avoid problems with pass-by-reference.
                }
            }
        } else {
            cells = src;
            for (const k of cells.columnNames()) {
                if (k in this.#raw_cell_factors) {
                    remap(cells, k, cells.column(k));
                }
            }
        }

        let loaded = scran.initializeSparseMatrixFromHDF5(this.#h5_path, this.#chosen_assay);
        let output = futils.splitScranMatrixAndFeatures(loaded, this.#raw_features, null, cache);
        output.cells = cells;

        if (!cache) {
            this.clear();
        }
        return output;
    }

    /**
     * @return {Array} Array of objects representing the files used in this dataset.
     * Each object corresponds to a single file and contains:
     * - `type`: a string denoting the type.
     * - `file`: a {@linkplain SimpleFile} object representing the file contents.
     */
    serialize() {
        return [ { type: "h5", file: this.#h5_file } ];
    }

    /**
     * @param {Array} files - Array of objects like that produced by {@linkcode H5adDataset#serialize serialize}.
     * @return {H5adDataset} A new instance of this class.
     * @static
     */
    static async unserialize(files) {
        if (files.length != 1 || files[0].type != "h5") {
            throw new Error("expected exactly one file of type 'h5' for H5AD unserialization");
        }
        return new H5adDataset(files[0].file);
    }
}
