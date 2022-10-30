import * as scran from "scran.js";
import * as bioc from "bioconductor";
import * as eutils from "./utils/extract.js";
import * as futils from "./utils/features.js";
import * as afile from "./abstract/file.js";

/**
 * Dataset in the H5AD format, see [here](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/advanced/matrices) for details.
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
     * @param {SimpleFile|string|Uint8Array|File} h5File - Contents of a HDF5 file.
     * On browsers, this may be a File object.
     * On Node.js, this may also be a string containing a file path.
     */
    constructor(h5File) {
        super();

        if (h5File instanceof afile.SimpleFile) {
            this.#h5_file = h5File;
        } else {
            this.#h5_file = new afile.SimpleFile(h5File);
        }

        let info = scran.realizeFile(this.#h5_file.content());
        this.#h5_path = info.path;
        this.#h5_flush = info.flush;

        try {
            this.#h5_handle = new scran.H5File(this.#h5_path);
        } catch (e) {
            this.#h5_flush();
            throw e;
        }

        this.#raw_features = null;
        this.#raw_cells = null;
        this.#chosen_assay = null;
        this.#assay_details = null;
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

    annotations() {
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
        
        return {
            features: { "": bioc.CLONE(this.#raw_features) },
            cells: {
                number: this.#raw_cells.numberOfColumns(),
                summary: summaries
            }
        };
    }

    load() {
        this.#features();
        this.#cells();
        this.#choose_assay();

        let cells = new bioc.DataFrame({}, { numberOfRows: this.#raw_cells.numberOfRows(), metadata: bioc.CLONE(this.#raw_cells.metadata()) });
        for (const k of this.#raw_cells.columnNames()) {
            let v = this.#raw_cells.column(k);
            if (!(k in this.#raw_cell_factors)) {
                cells.$setColumn(k, bioc.CLONE(v));
            } else {
                let cats = this.#raw_cell_factors[k];
                let temp = new Array(v.length);
                v.forEach((x, i) => {
                    temp[i] = cats[x]
                }); 
                cells.$setColumn(k, temp);
            }
        }

        let loaded = scran.initializeSparseMatrixFromHDF5(this.#h5_path, this.#chosen_assay);
        let output = futils.splitScranMatrixAndFeatures(loaded, this.#raw_features, null);
        output.cells = cells;
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
