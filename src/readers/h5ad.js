import * as scran from "scran.js";
import { Dataset } from "./base.js";
import * as eutils from "./utils/extract.js";
import * as futils from "./utils/features.js";
import * as afile from "../abstract/file.js";

/**
 * Dataset in the H5AD format, see [here](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/advanced/matrices) for details.
 * @augments Dataset
 */
export class H5adDataset extends Dataset {
    #h5_file;
    #h5_path;
    #h5_flush;

    #check_features;
    #raw_features;
    #check_cells;
    #raw_cells;
    #raw_cell_factors;

    #handle;
    #counts_details;

    /**
     * @param {SimpleFile|string|Uint8Array|File} h5File - Contents of a HDF5 file.
     * On browsers, this may be a File object.
     * On Node.js, this may also be a string containing a file path.
     */
    constructor(h5_file) {
        super();

        if (h5_file instanceof afile.SimpleFile) {
            this.#h5_file = h5_file;
        } else {
            this.#h5_file = new afile.SimpleFile(h5);
        }

        this.#h5_path = info.path;
        this.#h5_flush = info.flush;

        this.#check_features = false;
        this.#check_cells = false;

        try {
            this.#handle = new scran.H5File(this.#h5_path);
        } catch (e) {
            this.#h5_flush();
            throw e;
        }

        this.#counts_details = null;
    }

    static format() {
        return "H5AD";
    }

    abbreviate() {
        return { 
            "h5": {
                "name": this.#h5_file.name(),
                "size": this.#h5_file.size()
            }
        };
    }

    #fetch_counts_details() {
        if (this.#counts_details !== null) {
            return;
        }

        if ("X" in this.#handle.children) {
            let shape = this.#handle.open("X").shape;
            this.#counts_details = {
                name: "X",
                cells: shape[0],
                features: shape[1]
            };
            return;
        } 

        if (!("layers" in this.#handle.children)) {
            throw new Error("expected 'X' or 'layers' in a H5AD file");
        }

        let lhandle = this.#handle.open("layers");
        if (!(lhandle instanceof scran.H5Group)) {
            throw new Error("expected a 'layers' group in a H5AD file");
        }

        for (const k of Object.keys(lhandle.children)) {
            if (k.match(/^count/i)) {
                let shape = lhandle.open(k);
                this.#counts_details = {
                    name: "layers/" + k,
                    cells: shape[0],
                    features: shape[1]
                };
                return;
            }
        }

        throw new Error("failed to find any count-like layer in a H5AD file");
    }

    #features() {
        if (this.#check_features) {
            return;
        }
        this.#check_features = true;

        let handle = this.#handle;
        if ("var" in handle.children && handle.children["var"] == "Group") {
            let vhandle = handle.open("var");
            let index = eutils.extractHDF5Strings(vhandle, "_index");
            if (index !== null) {
                genes = { "_index": index };

                for (const [key, val] of Object.entries(vhandle.children)) {
                    if (val === "DataSet" && (key.match(/name/i) || key.match(/symb/i))) {
                        let dhandle2 = vhandle.open(key);
                        if (dhandle2.type == "String") {
                            genes[key] = dhandle2.load();
                        }
                    }
                }
                
                this.#raw_features = genes;
                return;
            }
        }

        this.#fetch_counts_details();
        let NR = this.#counts_details.features;
        let ids = [];
        for (var i = 0; i < NR; i++) {
            ids.push("Feature " + String(i));
        }
        this.#raw_features = { id: id };
        return;
    }

    #cells() {
        if (this.#check_cells) {
            return;
        }
        this.#check_cells = true;

        let handle = this.#handle;
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

        this.#raw_cells = annotations;
        this.#raw_cell_factors = factors;
        return;
    }

    annotations() {
        this.#features();
        this.#cells();

        let summaries = {};
        for (const [k, v] of Object.entries(this.#raw_cells)) {
            if (k in this.#raw_cell_factors) {
                summaries[k] = eutils.summarizeArray(this.#raw_cell_factors[k]);
            } else {
                summaries[k] = eutils.summarizeArray(v);
            }
        }

        return {
            features: { "": scran.cloneArrayCollection(this.#raw_features) },
            cells: summaries
        };
    }

    load() {
        this.#features();
        this.#cells();
        this.#fetch_counts_details();

        let cells = {};
        for (const [k, v] of Object.entries(this.#raw_cells)) {
            if (!(k in this.#raw_cell_factors)) {
                cells[k] = v.slice();
                continue;
            }

            let cats = this.#raw_cell_factors[k];
            let temp = new Array(v.length);
            v.forEach((x, i) => {
                temp[i] = cats[x]
            }); 
            cells[k] = temp;
        }

        let loaded = scran.initializeSparseMatrixFromHDF5(tmppath, this.#counts_details.name);
        let output = futils.splitScranMatrixAndFeatures(loaded, this.#raw_features, null);
        output.cells = scran.cloneArrayCollection(this.#raw_cells);
        return output;
    }

    free() {
        this.#h5_flush();
    }

    serialize() {
        return [ { type: "h5", file: this.#h5_file } ];
    }

    static async unserialize(files) {
        if (files.length != 1 || files[0].type != "h5") {
            throw new Error("expected exactly one file of type 'h5' for H5AD unserialization");
        }
        return new H5adDataset(files[0].file);
    }
}
