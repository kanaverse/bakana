import * as scran from "scran.js";
import * as utils from "./utils/general.js";
import * as iutils from "./utils/inputs.js";

/**
 * This step handles the loading of the input count matrices into memory.
 * This wraps various matrix initialization functions in [**scran.js**](https://github.com/jkanche/scran.js),
 * depending on the format of the supplied matrices.
 *
 * Methods not documented here are not part of the stable API and should not be used by applications.
 * @hideconstructor
 */
export class InputsState {
    #parameters;
    #cache;
    #abbreviated;

    constructor(parameters = null, cache = null) {
        this.#parameters = (parameters === null ? {} : parameters);
        this.#cache = (cache === null ? {} : cache);
        this.#abbreviated = {};
        this.changed = false;
        return;
    }

    free() {
        utils.freeCache(this.#cache.matrix);
        utils.freeCache(this.#cache.block_ids);

        let alt = this.#cache.alternatives;
        if (alt) {
            for (const v of Object.values(alt)) {
                utils.freeCache(v.matrix);
            }
        }

        return;
    }

    /***************************
     ******** Getters **********
     ***************************/

    fetchCountMatrix() {
        return this.#cache.matrix;
    }

    fetchGenes() {
        return this.#cache.genes;
    }

    fetchGeneTypes() {
        return this.#cache.gene_types;
    }

    hasAdt() {
        return this.#cache.adt_key !== null;
    }

    fetchAdtCountMatrix() {
        return this.#cache.alternatives[this.#cache.adt_key].matrix;
    }

    fetchAdtGenes() {
        return this.#cache.alternatives[this.#cache.adt_key].genes;
    }
  
    /**
     * Fetch an annotation for all cells in the dataset.
     * This considers all cells in the dataset before QC filtering - 
     * see {@linkcode QualityControlState#fetchFilteredAnnotations QualityControlState.fetchFilteredAnnotations} for an alternative.
     *
     * @param {string} col - Name of the annotation field of interest.
     *
     * @return An object containing the requested annotations.
     * This will contain:
     * - `type`: string specifying the the type of annotations.
     *   This can be either `"factor"` or `"array"`.
     *
     * For `type: "factor"`, the object will additionally contain:
     * - `levels`: an array of strings containing the unique factor levels.
     * - `index`: an Int32Array of length equal to the number of cells, referencing entries in `levels`.
     * 
     * For `type: "array"`, the object will additionally contain:
     * - `values`: a TypedArray of values of length equal to the number of cells.
     */
    fetchAnnotations(col) {
        let annots = this.#cache.annotations;
        let size = this.#cache.matrix.numberOfColumns();

        if (annots === null || !(col in annots)) {
            throw new Error(`${col} does not exist in the column annotations`);
        }

        let current = annots[col];

        // i.e., is it already a factor? In which case, we make a copy of its contents.
        // This ensures we don't have any accidental writes or transfers. 
        if (ArrayBuffer.isView(current)) {
            return {
                "type": "array",
                "values": current.slice()
            };
        } else if (utils.isObject(current)) {
            if (!("type" in current) || current.type != "factor") {
                throw new Error("annotation column should have 'type: \"factor\"' if it is an object");
            }
            return { 
                type: current.type,
                index: current.index.slice(),
                levels: current.levels.slice()
            };
        }

        let uniq_vals = [];
        let uniq_map = {};
        let indices = new Int32Array(size);
        annots[col].map((x, i) => {
            if (!(x in uniq_map)) {
                uniq_map[x] = uniq_vals.length;
                uniq_vals.push(x);
            }
            indices[i] = uniq_map[x];
        });

        return {
            "type": "factor",
            "index": indices,
            "levels": uniq_vals
        };
    }

    fetchBlock() {
        return this.#cache.block_ids;
    }

    fetchBlockLevels() {
        return this.#cache.block_levels;
    }

    /***************************
     ******** Compute **********
     ***************************/

    /**
     * This method should not be called directly by users, but is instead invoked by {@linkcode runAnalysis}.
     * `matrices` is taken from the argument of the same name in {@linkcode runAnalysis},
     * while `sample_factor` is taken from the property of the same name in the `inputs` property of the `parameters`.
     *
     * @param {object} matrices - An object containing data for one or more count matrices.
     * Each property corresponds to a single matrix and should contain a `format` string property along with any number of File objects (browser) or file paths (Node.js).
     * See the description of the argument of the same name in {@linkcode runAnalysis}.
     * @param {?string} sample_factor - Name of the column of the cell annotations specifying the sample of origin for each cell.
     * This is only used if a single count matrix is supplied.
     *
     * If `null`, all cells are assumed to originate from the same sample.
     *
     * @return The object is updated with the new results.
     * A promise is returned that resolves to `null` once input loading is complete - this should be resolved before any downstream steps are run.
     */
    async compute(matrices, sample_factor) {
        // Don't bother proceeding with any of the below
        // if we're operating from a reloaded state.
        if (matrices === null) {
            this.changed = false;
            return;
        }

        let entries = Object.entries(matrices);
        let tmp_abbreviated = {};
        for (const [key, val] of entries) {
            let namespace = iutils.chooseReader(val.format);
            tmp_abbreviated[key] = namespace.abbreviate(val);
        }

        if (!utils.changedParameters(tmp_abbreviated, this.#abbreviated) && this.#parameters.sample_factor === sample_factor) {
            this.changed = false;
            return;
        }

        let new_matrices = {};
        for (const [key, val] of entries) {
            let namespace = iutils.chooseReader(val.format);
            new_matrices[key] = new namespace.Reader(val);
        }

        this.free();
        this.#cache = await process_and_cache(new_matrices, sample_factor);

        this.#abbreviated = tmp_abbreviated;
        this.#parameters.matrices = new_matrices;
        this.#parameters.sample_factor = sample_factor;

        this.changed = true;

        return null;
    }

    /***************************
     ******** Results **********
     ***************************/

    /**
     * Obtain a summary of the state, typically for display on a UI like **kana**.
     *
     * @return An object containing:
     *
     * - `dimensions`: an object containing `num_genes` and `num_cells`, the number of genes and cells respectively.
     *   For multiple matrices, the number of cells is the total number across all matrices.
     * - `genes`: an object containing the per-gene annotation.
     *   Each property is an array of length equal to the number of genes, usually containing strings with gene identifiers or symbols.
     *   Property names are arbitrary.
     * - (optional) `annotations`: an array of strings containing the names of available cell annotation fields.
     */
    summary() {
        var output = {
            "dimensions": {
                "num_genes": this.#cache.matrix.numberOfRows(),
                "num_cells": this.#cache.matrix.numberOfColumns()
            },
            "genes": { ...(this.#cache.genes) }
        };
        if (this.#cache.annotations !== null) {
            output.annotations = Object.keys(this.#cache.annotations);
        }
        return output;
    }

    /*************************
     ******** Saving *********
     *************************/

    async serialize(handle, embeddedSaver) {
        let ghandle = handle.createGroup("inputs");

        let multifile = false;
        {
            let phandle = ghandle.createGroup("parameters");

            let formats = [];
            let names = [];
            let numbers = [];
            let fihandle = phandle.createGroup("files");
            let sofar = 0;

            for (const [key, val] of Object.entries(this.#parameters.matrices)) {
                formats.push(val.format());
                names.push(key);

                let files = await val.serialize(embeddedSaver);
                numbers.push(files.length);

                for (const obj of files) {
                    let curhandle = fihandle.createGroup(String(sofar));
                    curhandle.writeDataSet("type", "String", [], obj.type);
                    curhandle.writeDataSet("name", "String", [], obj.name);

                    if (typeof obj.id == "string") {
                        curhandle.writeDataSet("id", "String", [], obj.id);
                    } else if (typeof obj.offset == "number" && typeof obj.size == "number") {
                        curhandle.writeDataSet("offset", "Uint32", [], obj.offset);
                        curhandle.writeDataSet("size", "Uint32", [], obj.size);
                    } else {
                        throw new Error("object should contain either an 'id' string or 'offset' and 'size' numbers"); 
                    }

                    sofar++;
                }
            }

            if (formats.length > 1) {
                multifile = true;
                phandle.writeDataSet("format", "String", null, formats);
                phandle.writeDataSet("sample_groups", "Int32", null, numbers);
                phandle.writeDataSet("sample_names", "String", null, names);
            } else {
                phandle.writeDataSet("format", "String", [], formats[0]);
                if (this.#parameters.sample_factor !== null) {
                    phandle.writeDataSet("sample_factor", "String", [], this.#parameters.sample_factor);
                }
            }
        }

        {
            let rhandle = ghandle.createGroup("results");
            rhandle.writeDataSet("dimensions", "Int32", null, [this.#cache.matrix.numberOfRows(), this.#cache.matrix.numberOfColumns()]);

            // For diagnostic purposes, we store the number of samples;
            // this may not be captured by the parameters if we're dealing
            // with a sample_factor from a single file.
            if (this.#cache.block_levels !== null) {
                rhandle.writeDataSet("num_samples", "Int32", [], this.#cache.block_levels.length); 
            }

            rhandle.writeDataSet("indices", "Int32", null, this.#cache.matrix.identities());
        }

        return;
    }
}

/**************************
 ******* Internals ********
 **************************/

async function process_datasets(matrices, sample_factor) {
    // Loading all of the individual matrices. 
    let datasets = {};
    try {
        for (const [key, val] of Object.entries(matrices)) {
            // Too much hassle to convert this into a Promise.all(), because we
            // need to make sure it gets freed properly on failure.
            datasets[key] = await val.load();
        }
    } catch (e) {
        // If any one fails, we free the rest.
        for (const [key, val] of Object.entries(datasets)){
            utils.freeCache(val.matrix);
        }
        throw e;
    }

    // Ensure we have a reproducible order; otherwise the batch
    // order becomes dependent on the JS engine's ordering.
    let dkeys = Object.keys(datasets);
    dkeys.sort();

    if (dkeys.length == 1) {
        let current = datasets[dkeys[0]];
        let blocks = null;
        let block_levels = null;

        if (sample_factor !== null) {
            // Single matrix with a batch factor.
            try {
                let anno_batch = current.annotations[sample_factor];

                let anno_levels = null;
                if (utils.isObject(anno_batch)) {
                    anno_levels = anno_batch.levels;
                    anno_batch = anno_batch.index;
                }

                let ncols = current.matrix.numberOfColumns();
                if (anno_batch.length != ncols) {
                    throw new Error("length of sample factor '" + sample_factor + "' should be equal to the number of cells"); 
                }

                blocks = scran.createInt32WasmArray(ncols);
                block_levels = [];
                let block_arr = blocks.array();

                let uvals = {};
                anno_batch.forEach((x, i) => {
                    if (!(x in uvals)) {
                        uvals[x] = block_levels.length;
                        block_levels.push(x);
                    }
                    block_arr[i] = uvals[x];
                });

                if (anno_levels !== null) {
                    block_levels = Array.from(block_levels).map(i => anno_levels[i]);
                }

            } catch (e) {
                utils.freeCache(blocks);
                utils.freeCache(current.matrix);
                throw e;
            }
        }

        current.block_ids = blocks;
        current.block_levels = block_levels;
        return current;

    } else {
        // Multiple matrices, each of which represents a batch.
        let output = {}
        let blocks;

        try {
            // Identify the gene columns to use
            let genes = {};
            for (const k of dkeys) {
                genes[k] = datasets[k].genes;
                if (genes[k] === null) {
                    throw new Error("no gene annotations found in matrix '" + k + "'");
                }
            }

            let result = iutils.commonFeatureTypes(genes);
            if (result.best_type === null) {
                throw new Error("no common feature types available across all matrices");
            }
            let best_fields = result.best_fields;

            let gnames = [];
            let mats = [];
            let total = 0;
            for (const k of dkeys) {
                let current = datasets[k];
                gnames.push(current.genes[best_fields[k]]);
                mats.push(current.matrix);
                total += current.matrix.numberOfColumns();
            }

            blocks = scran.createInt32WasmArray(total);
            let barr = blocks.array();
            let nice_barr = new Array(total);
            let sofar = 0;
            for (var i = 0; i < dkeys.length; i++) {
                let old = sofar;
                let current = datasets[dkeys[i]];
                sofar += current.matrix.numberOfColumns();
                barr.fill(i, old, sofar);
                nice_barr.fill(dkeys[i], old, sofar);
            }
            output.block_ids = blocks;
            output.block_levels = dkeys;

            let merged = scran.cbindWithNames(mats, gnames);
            output.matrix = merged.matrix;

            // Extracting gene information from the first object. We won't make
            // any attempt at merging and deduplication across objects.
            let included = new Set(merged.indices);
            output.genes = {};
            let first_genes = datasets[dkeys[0]].genes;
            for (const [key, val] of Object.entries(first_genes)) {
                output.genes[key] = val.filter((x, i) => included.has(i));
            }

            // Get all annotations keys across datasets; we then concatenate
            // columns with the same name, or we just fill them with missings.
            let ckeys = new Set();
            for (const d of dkeys) {
                let current = datasets[d];
                if (current.annotations !== null) {
                    for (const a of Object.keys(current.annotations)) {
                        ckeys.add(a);
                    }
                }
            }
            let anno_keys = [...ckeys];

            let combined_annotations = {};
            for (const i of anno_keys) {
                let current_combined = [];
                for (const d of dkeys) {
                    let current = datasets[d];
                    let x;
                    if (current.annotations && current.annotations[i]) {
                        x = current.annotations[i];
                        if (!(x instanceof Array)) {
                            x = Array.from(x);
                        }
                    } else {
                        x = new Array(current.matrix.numberOfColumns());
                    }
                    current_combined = current_combined.concat(x);
                }
                combined_annotations[i] = current_combined;
            }

            output.annotations = combined_annotations;
            output.annotations["__batch__"] = nice_barr;

        } catch (e) {
            utils.freeCache(blocks);
            throw e;

        } finally {
            // Once the merged dataset is created, the individual dataset
            // matrices are no longer useful, so we need to delete them anyway.
            for (const v of Object.values(datasets)) {
                utils.freeCache(v.matrix);
            }
        }

        return output;
    }
}

async function process_and_cache(new_matrices, sample_factor) {
    let cache = await process_datasets(new_matrices, sample_factor);

    var gene_info_type = {};
    var gene_info = cache.genes;
    for (const [key, val] of Object.entries(gene_info)) {
        gene_info_type[key] = scran.guessFeatures(val);
    }
    cache.gene_types = gene_info_type;

    cache.adt_key = null;
    if ("alternatives" in cache) {
        for (const k of Object.keys(cache.alternatives)) {
            if (k.match(/antibody/i) || k.match(/adt/i)) {
                cache.adt_key = k;
                break;
            }
        }
    }

    return cache;
}

/**************************
 ******** Loading *********
 **************************/

export async function unserialize(handle, embeddedLoader) {
    let ghandle = handle.open("inputs");
    let phandle = ghandle.open("parameters");

    // Extracting the files.
    let fihandle = phandle.open("files");
    let kids = fihandle.children;
    let all_files = new Array(kids.length);

    for (const x of Object.keys(kids)) {
        let current = fihandle.open(x);

        let curfile = {};
        for (const field of ["type", "name"]) {
            let dhandle = current.open(field, { load: true });
            curfile[field] = dhandle.values[0];
        }

        if ("id" in current.children) {
            curfile.id = current.open("id", { load: true }).values[0];
        } else {
            for (const field of ["offset", "size"]) {
                curfile[field] = current.open(field, { load: true }).values[0];
            }
        }

        let idx = Number(x);
        all_files[idx] = curfile;
    }

    // Extracting the format and organizing the files.
    let parameters = { matrices: {}, sample_factor: null };
    let fohandle = phandle.open("format", { load: true });
    let solofile = (fohandle.shape.length == 0);

    if (solofile) {
        let format = fohandle.values[0];
        let namespace = iutils.chooseReader(format);
        parameters.matrices["default"] = await namespace.unserialize(all_files, embeddedLoader);
        if ("sample_factor" in phandle.children) {
            parameters.sample_factor = phandle.open("sample_factor", { load: true }).values[0];
        } else {
            parameters.sample_factor = null;
        }

    } else {
        let formats = fohandle.values;
        let sample_names = phandle.open("sample_names", { load: true }).values;
        let sample_groups = phandle.open("sample_groups", { load: true }).values;

        let sofar = 0;
        for (var i = 0; i < formats.length; i++) {
            let start = sofar;
            sofar += sample_groups[i];
            let curfiles = all_files.slice(start, sofar);
            let namespace = iutils.chooseReader(formats[i]);
            parameters.matrices[sample_names[i]] = await namespace.unserialize(curfiles, embeddedLoader);
        }
    }

    // Loading matrix data.
    let cache = await process_and_cache(parameters.matrices, parameters.sample_factor);

    // We need to do something if the permutation is not the same.
    let rhandle = ghandle.open("results");
    let permuter;

    let perm = null;
    if (solofile) {
        if ("permutation" in rhandle.children) {
            // v1.0-v1.1
            let dhandle = rhandle.open("permutation", { load: true });
            let ids = new Int32Array(dhandle.values.length);
            dhandle.values.forEach((x, i) => { ids[x] = i; });
            perm = scran.updateRowIdentities(cache.matrix, ids);
        } else if ("identities" in rhandle.children) {
            // v1.2+
            let dhandle = rhandle.open("identities", { load: true });
            perm = scran.updateRowIdentities(cache.matrix, dhandle.values);
        } else {
            // Otherwise, we're dealing with v0 states. We'll just
            // assume it was the same, I guess. Should be fine as we didn't change
            // the permutation code in v0.
        }
    } else {
        let old_ids;
        if ("indices" in rhandle.children) {
            // v1.1
            old_ids = rhandle.open("indices", { load: true }).values;
        } else {
            // v1.2+
            old_ids = rhandle.open("identities", { load: true }).values;
        }
        perm = scran.updateRowIdentities(cache.matrix, old_ids);
    }

    if (perm === null) {
        permuter = (x) => {}
    } else {
        permuter = (x) => {
            let copy = x.slice();
            x.forEach((y, i) => {
                x[i] = copy[perm[i]];
            });
        };
    }

    return { 
        state: new InputsState(parameters, cache),
        parameters: { sample_factor: parameters.sample_factor }, // only returning the sample factor - we don't pass the files back. 
        permuter: permuter
    }
}
