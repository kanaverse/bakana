import * as scran from "scran.js";
import * as utils from "./utils/general.js";
import * as iutils from "./utils/inputs.js";

/**
 * This step handles the loading of the input count matrices into memory.
 * This wraps various matrix initialization functions in [**scran.js**](https://github.com/jkanche/scran.js),
 * depending on the format of the supplied matrices.
 *
 * The parameters in {@linkcode runAnalysis} should be an object containing:
 *
 * - `sample_factor`: a string specifying the column of the cell annotations specifying the sample of origin for each cell.
 *   This is only used if a single count matrix is supplied.
 *   If `null`, all cells are assumed to originate from the same sample.
 *
 * Calling the **`results()`** method for the relevant state instance will return an object containing:
 *
 * - `dimensions`: an object containing `num_genes` and `num_cells`, the number of genes and cells respectively.
 *   For multiple matrices, the number of cells is the total number across all matrices.
 * - `genes`: an object containing the per-gene annotation.
 *   Each property is an array of length equal to the number of genes, usually containing strings with gene identifiers or symbols.
 *   Property names are arbitrary.
 * - (optional) `annotations`: an array of strings containing the names of available cell annotation fields.
 *
 * If `annotations` is present, users can also call `fetchAnnotations(col)` on the state instance.
 *
 * - `col` should be a string containing the name of a cell annotation field.
 * - The return value is either:
 *   - An array of length equal to the number of cells in the original dataset.
 *   - An object representing a factor.
 *     This should contain `factor`, an array of strings containing the unique factor levels;
 *     and `index`, an array of length equal to the number of cells, referencing entries in `factor`. 
 * - Note that the return value refers to the cells prior to any QC filtering.
 *   Actual use will require further subsetting based on the discard vector from {@linkcode quality_control}.
 * 
 * Methods not documented here are not part of the stable API and should not be used by applications.
 *
 * @namespace inputs
 */

export class State {
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

    fetchAnnotations(col) {
        let annots = this.#cache.annotations;
        let size = this.#cache.matrix.numberOfColumns();

        if (!(col in annots)) {
            throw new Error(`${col} does not exist in the column annotations`);
        }

        // i.e., is it already a factor?
        if (utils.isObject(annots[col]) && "type" in annots[col]) {
            return annots[col];
        }

        let uvals = {};
        let uTypedAray = new Uint8Array(size);
        annots[col].map((x, i) => {
            if (!(x in uvals)) {
                uvals[x] = Object.keys(uvals).length;
            }

            uTypedAray[i] = uvals[x];
        });

        return {
            "index": Object.keys(uvals),
            "factor": uTypedAray
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

    compute(files, sample_factor) {
        // Don't bother proceeding with any of the below
        // if we're operating from a reloaded state.
        if (files === null) {
            this.changed = false;
            return;
        }

        let entries = Object.entries(files);
        let tmp_abbreviated = {};
        for (const [key, val] of entries) {
            let namespace = iutils.chooseReader(val.format);
            tmp_abbreviated[key] = namespace.abbreviate(val);
        }

        if (!utils.changedParameters(tmp_abbreviated, this.#abbreviated) && this.#parameters.sample_factor === sample_factor) {
            this.changed = false;
            return;
        }

        let new_files = {};
        for (const [key, val] of entries) {
            let namespace = iutils.chooseReader(val.format);
            new_files[key] = new namespace.Reader(val);
        }

        utils.freeCache(this.#cache.matrix);
        utils.freeCache(this.#cache.block_ids);
        this.#cache = process_and_cache(new_files, sample_factor);

        this.#abbreviated = tmp_abbreviated;
        this.#parameters.files = new_files;
        this.#parameters.sample_factor = sample_factor;

        this.changed = true;

        return;
    }

    /***************************
     ******** Results **********
     ***************************/

    results() {
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

            for (const [key, val] of Object.entries(this.#parameters.files)) {
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

            if (multifile) {
                rhandle.writeDataSet("indices", "Int32", null, this.#cache.indices);
            } else {
                rhandle.writeDataSet("permutation", "Int32", null, this.#cache.matrix.permutation({ copy: "view" }));
            }
        }

        return;
    }
}

/**************************
 ******* Internals ********
 **************************/

function dummy_genes(numberOfRows) {
    let genes = []
    for (let i = 0; i < numberOfRows; i++) {
        genes.push(`Gene ${i + 1}`);
    }
    return { "id": genes };
}

function process_datasets(files, sample_factor) {
    // Loading all of the individual matrices.
    let datasets = {};
    try {
        for (const [key, val] of Object.entries(files)) {
            let current = val.load();

            if (!("genes" in current)) {
                current.genes = dummy_genes(current.matrix.numberOfRows());
            } 
            if (current.matrix.isPermuted()) {
                scran.permuteFeatures(current.matrix, current.genes);
            }

            datasets[key] = current;
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

        if (sample_factor && sample_factor !== null) {
            // Single matrix with a batch factor.
            try {
                let anno_batch = current.annotations[sample_factor];
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

            // Storing the identities of the genes in terms of the
            // original row indices of the first matrix.
            let firstperm = mats[0].permutation({ restore: false });
            output.indices = merged.indices.map(i => firstperm[i]);

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

function process_and_cache(new_files, sample_factor) {
    let cache = process_datasets(new_files, sample_factor);

    var gene_info_type = {};
    var gene_info = cache.genes;
    for (const [key, val] of Object.entries(gene_info)) {
        gene_info_type[key] = scran.guessFeatures(val);
    }
    cache.gene_types = gene_info_type;

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
    let parameters = { files: {}, sample_factor: null };
    let fohandle = phandle.open("format", { load: true });
    let solofile = (fohandle.shape.length == 0);

    if (solofile) {
        let format = fohandle.values[0];
        let namespace = iutils.chooseReader(format);
        parameters.files["default"] = await namespace.unserialize(all_files, embeddedLoader);
        if ("sample_factor" in phandle.children) {
            parameters.sample_factor = phandle.open("sample_factor", { load: true }).values[0];
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
            parameters.files[sample_names[i]] = await namespace.unserialize(curfiles, embeddedLoader);
        }
    }

    // Loading matrix data.
    let cache = process_and_cache(parameters.files, parameters.sample_factor);

    // We need to do something if the permutation is not the same.
    let rhandle = ghandle.open("results");
    let permuter;
    if (solofile) {
        let perm = null;
        if ("permutation" in rhandle.children) {
            let dhandle = rhandle.open("permutation", { load: true });
            perm = scran.updatePermutation(cache.matrix, dhandle.values);
        } else {
            // Otherwise, we're dealing with v0 states. We'll just
            // assume it was the same, I guess. Should be fine as we didn't change
            // the permutation code in v0.
        }

        if (perm !== null) {
            // Adding a permuter function for all per-gene vectors.
            permuter = (x) => {
                let temp = x.slice();
                x.forEach((y, i) => {
                    temp[i] = x[perm[i]];
                });
                x.set(temp);
                return;
            };
        } else {
            permuter = (x) => { };
        }
    } else {
        let old_indices = rhandle.open("indices", { load: true }).values;
        let new_indices = cache.indices;
        if (old_indices.length != new_indices.length) {
            throw new Error("old and new indices must have the same length for results to be interpretable");
        }

        let same = true;
        for (var i = 0; i < old_indices.length; i++) {
            if (old_indices[i] != new_indices[i]) {
                same = false
                break
            }
        }

        if (same) {
            permuter = (x) => {};
        } else {
            let remap = {};
            old_indices.forEach((x, i) => {
                remap[x] = i;
            });
            let perm = new_indices.map(o => {
                if (!(o in remap) || remap[o] === null) {
                    throw new Error("old and new indices should contain the same row identities");
                }
                let pos = remap[o];
                remap[o] = null;
                return pos;
            });
            permuter = (x) => {
                return x.map((y, i) => x[perm[i]]);
            };
        }
    }

    return { 
        state: new State(parameters, cache),
        parameters: { sample_factor: parameters.sample_factor }, // only returning the sample factor - we don't pass the files back. 
        permuter: permuter
    }
}