import * as scran from "scran.js";
import * as utils from "./utils/general.js";
import * as iutils from "../readers/index.js";
export const step_name = "inputs";

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

    constructor(parameters = null, cache = null, abbreviated = null) {
        this.#parameters = (parameters === null ? {} : parameters);
        this.#cache = (cache === null ? {} : cache);
        this.#abbreviated = (abbreviated === null ? {} : abbreviated);
        this.changed = false;
        return;
    }

    free() {
        utils.freeCache(this.#cache.block_ids);
        utils.freeCache(this.#cache.matrix);
    }

    /***************************
     ******** Getters **********
     ***************************/

    listAvailableTypes() {
        return this.#cache.matrix.available();
    }

    hasAvailable(type) {
        return this.#cache.matrix.has(type);
    }

    fetchCountMatrix({ type = "RNA" } = {}) {
        return this.#cache.matrix.get(type);
    }

    fetchGenes({ type = "RNA" } = {}) {
        return this.#cache.genes[type];
    }

    fetchGeneTypes() {
        return this.#cache.gene_types;
    }

    /**
     * Fetch an annotation for all cells in the dataset.
     * This considers all cells in the dataset before QC filtering - 
     * see {@linkcode QualityControlState#fetchFilteredAnnotations QualityControlState.fetchFilteredAnnotations} for an alternative.
     *
     * @param {string} col - Name of the annotation field of interest.
     *
     * @return {Array|TypedArray} Array of length equal to the total number of cells, containing the requested annotations.
     */
    fetchAnnotations(col) {
        let annots = this.#cache.annotations;
        let size = this.#cache.matrix.numberOfColumns();

        if (annots === null || !(col in annots)) {
            throw new Error(`${col} does not exist in the column annotations`);
        }

        // Make a copy, avoid accidental writes or transfers. 
        return annots[col].slice();
    }

    fetchBlock() {
        return this.#cache.block_ids;
    }

    fetchBlockLevels() {
        return this.#cache.block_levels;
    }

    fetchParameters({ subset = true } = {}) {
        // Cloning the parameters to avoid pass-by-reference behavior affecting the
        // InputsState object. We don't pass the files back here.
        let output = { ...this.#parameters };

        let replacement = null;
        if (subset && output.subset !== null) {
            replacement = { ...output.subset };
            for (const k of ["indices", "values"]) {
                if (k in replacement) {
                    replacement[k] = replacement[k].slice();
                }
            }
        }
        output.subset = replacement;

        return output;
    }

    /***************************
     ******** Compute **********
     ***************************/

    static defaults() {
        return {
            sample_factor: null,
            subset: null
        };
    }

    static areSubsetsEqual(sub1, sub2) {
        let isund1 = (typeof sub1 == "undefined");
        let isund2 = (typeof sub2 == "undefined");
        if (isund1 !== isund2) {
            return false;
        } else if (isund1) {
            return true;
        }

        let isnull1 = (sub1 === null);
        let isnull2 = (sub2 === null);
        if (isnull1 !== isnull2) {
            return false;
        } else if (isnull1) {
            return true;
        }

        let has_indices1 = ("indices" in sub1);
        let has_indices2 = ("indices" in sub2);
        if (has_indices1 != has_indices2) {
            return false;
        } else if (!has_indices1) {
            return !utils.changedParameters(sub1, sub2);
        }

        // Checking everything other than 'indices'; don't try
        // to stringify potentially long arrays.
        let rest1 = { ...sub1, indices: null };
        let rest2 = { ...sub2, indices: null };
        if (utils.changedParameters(rest1, rest2)) {
            return false;
        }

        let n = sub1.indices.length;
        if (n != sub2.indices.length) {
            return false;
        }

        for (var i = 0; i < n; i++) {
            if (sub1.indices[i] != sub2.indices[i]) {
                return false;
            }
        }

        return true;
    }

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
     * @param {?subset} subset - Object describing if any pre-analysis subsetting should be applied.
     * This may contain either:
     *
     * - `indices`, an array of integers that specifies indices of the columns to retain in the loaded dataset.
     *   For multiple `matrices`, the loaded dataset is created by sorting the individual matrices by the keys of `matrices` before combining them by column.
     * - `field`, a string specifying a field of the column annotation; and `values`, an array of allowed values for that annotation.
     *   Cells are only retained if they are associated with any of the allowable values for that annotation field.
     *
     * If both are present in `subset`, `indices` takes precedence.
     * If `subset` is `null`, no subsetting is performed and all cells are used in the downstream analysis.
     *
     * @return The object is updated with the new results.
     * A promise is returned that resolves to `null` once input loading is complete - this should be resolved before any downstream steps are run.
     */
    async compute(matrices, sample_factor, subset) {
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

        if (!utils.changedParameters(tmp_abbreviated, this.#abbreviated) && this.#parameters.sample_factor === sample_factor && InputsState.areSubsetsEqual(this.#parameters.subset, subset)) {
            this.changed = false;
            return;
        }

        let new_matrices = {};
        for (const [key, val] of entries) {
            let namespace = iutils.chooseReader(val.format);
            new_matrices[key] = new namespace.Reader(val);
        }

        this.free();
        this.#cache = await process_and_cache(new_matrices, sample_factor, subset);

        this.#abbreviated = tmp_abbreviated;
        this.#parameters.matrices = new_matrices;
        this.#parameters.sample_factor = sample_factor;
        this.#parameters.subset = subset;

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
        let ngenes = {};
        for (const a of this.#cache.matrix.available()) {
            ngenes[a] = this.#cache.matrix.get(a).numberOfRows();
        }
        
        var output = {
            "num_cells": this.#cache.matrix.numberOfColumns(),
            "num_genes": ngenes,
            "genes": { ...(this.#cache.genes) }
        };
        if (this.#cache.annotations !== null) {
            output.annotations = Object.keys(this.#cache.annotations);
        }
        return output;
    }

    /**************************
     ******** Subset **********
     **************************/

    subsetCells(subset) {
        utils.checkIndices(subset, this.#cache.matrix.numberOfColumns());

        // Pass by reference of the cache's internal contents is probably okay
        // here, given that everything is readonly anyway.
        let new_cache = { ...this.#cache };
        subset_datasets(new_cache, subset, { freeExisting: false });

        // Faking the original indices for caching.
        if (this.#parameters.subset !== null) {
            let old_keep = harvest_subset_indices(this.#parameters.subset, null, this.#cache.annotations);
            subset = subset.map(i => old_keep[i]);
        }

        let new_params = this.fetchParameters({ subset: false });
        new_params.subset = { indices: subset };

        return new InputsState(new_params, new_cache, { ...this.#abbreviated });
    }

    /*************************
     ******** Saving *********
     *************************/

    async serialize(handle, embeddedSaver) {
        let ghandle = handle.createGroup("inputs");

        let multifile = false;
        {
            let phandle = ghandle.createGroup("parameters");

            // Make sure we're in sorted order, for consistency with
            // how the merge is done.
            let names = Object.keys(this.#parameters.matrices);
            names.sort();

            let formats = [];
            let numbers = [];
            let fihandle = phandle.createGroup("files");
            let sofar = 0;

            for (const key of names) {
                let val = this.#parameters.matrices[key];
                formats.push(val.format());

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

            if (this.#parameters.subset !== null) {
                let shandle = phandle.createGroup("subset");
                let schandle = shandle.createGroup("cells");

                if ("indices" in this.#parameters.subset) {
                    schandle.writeDataSet("indices", "Int32", null, this.#parameters.subset.indices);
                } else if ("field" in this.#parameters.subset) {
                    schandle.writeDataSet("field", "String", [], this.#parameters.subset.field);
                    schandle.writeDataSet("values", "String", null, this.#parameters.subset.values);
                } else {
                    throw new Error("unrecognized specification for 'subset'");
                }
            }
        }

        {
            let rhandle = ghandle.createGroup("results");
            rhandle.writeDataSet("num_cells", "Int32", [], this.#cache.matrix.numberOfColumns());

            let fhandle = rhandle.createGroup("num_features");
            for (const a of this.#cache.matrix.available()) {
                fhandle.writeDataSet(a, "Int32", [], this.#cache.matrix.get(a).numberOfRows());
            }

            // For diagnostic purposes, we store the number of samples;
            // this may not be captured by the parameters if we're dealing
            // with a sample_factor from a single file.
            if (this.#cache.block_levels !== null) {
                rhandle.writeDataSet("num_samples", "Int32", [], this.#cache.block_levels.length); 
            }

            // Looping through all available matrices.
            let ihandle = rhandle.createGroup("identities");
            for (const a of this.#cache.matrix.available()) {
                ihandle.writeDataSet(a, "Int32", null, this.#cache.matrix.get(a).identities());
            }
        }

        return;
    }
}

/**************************
 ******* Internals ********
 **************************/

// Exported for testing only.
export function commonFeatureTypes(genes) {
    let scores = {
        "symbol-mouse": [],
        "symbol-human": [],
        "ensembl-mouse": [],
        "ensembl-human": []
    };
    let fields = JSON.parse(JSON.stringify(scores));

    let names = Object.keys(genes);
    for (const name of names) {
        let curgenes = genes[name];

        let best_scores = {};
        let best_fields = {};
        for (const [k, v] of Object.entries(curgenes)) {
            let fscore = scran.guessFeatures(v);
            let curname = fscore.type + "-" + fscore.species;
            if (!(curname in best_scores) || fscore.confidence > best_scores[curname]) {
                best_scores[curname] = fscore.confidence;
                best_fields[curname] = k;
            }
        }

        for (const [k, v] of Object.entries(best_fields)) {
            fields[k].push(v);
            scores[k].push(best_scores[k]);
        }
    }

    let best_score = -1000;
    let best_type = null;

    for (const [k, v] of Object.entries(scores)) {
        if (v.length == names.length) { // skipping if not represented in all entries.
            let nscore = v.reduce((a, b) => a * b);
            if (nscore > best_score) {
                best_score = nscore;
                best_type = k;
            }
        }
    }

    let best_fields = {};
    let best_features = null;

    if (best_type !== null) {
        let best_type_cols = fields[best_type];
        let best_features_sub = best_type.split("-");
        best_features = {
            type: best_features_sub[0],
            species: best_features_sub[1]
        };
        for (var i = 0; i < names.length; i++) {
            best_fields[names[i]] = best_type_cols[i];
        }
    }

    return {
        "best_type": best_features,
        "best_fields": best_fields
    };
}

function bind_single_modality(dkeys, datasets, type) {
    let output = {};

    try {
        // Identify the gene columns to use.
        let genes = {};
        for (const k of dkeys) {
            genes[k] = datasets[k].genes[type];
            if (genes[k] === null) {
                throw new Error("no gene annotations found in matrix '" + k + "'");
            }
        }

        let result = commonFeatureTypes(genes);
        if (result.best_type === null) {
            throw new Error("no common feature types available across all matrices");
        }
        let best_fields = result.best_fields;

        let gnames = [];
        let mats = [];
        for (const k of dkeys) {
            gnames.push(genes[k][best_fields[k]]);
            mats.push(datasets[k].matrix.get(type));
        }

        let merged = scran.cbindWithNames(mats, gnames);
        output.matrix = merged.matrix;

        // Extracting gene information from the first object. We won't make
        // any attempt at merging and deduplication across objects.
        let first_genes = genes[dkeys[0]];
        output.genes = scran.subsetArrayCollection(first_genes, merged.indices);

    } catch (e) {
        utils.freeCache(output.matrix);
        throw e;
    }

    return output;
}

function bind_datasets(dkeys, datasets) {
    // Checking which feature types are available across all datasets.
    let available = null;
    for (const k of dkeys) {
        if (available === null) {
            available = datasets[k].matrix.available();
        } else {
            let present = new Set(datasets[k].matrix.available());
            available = available.filter(x => present.has(x));
        }
    }

    let blocks;
    let output = { 
        matrix: new scran.MultiMatrix, 
        genes: {} 
    };

    try {
        for (const a of available) {
            let current = bind_single_modality(dkeys, datasets, a);
            output.matrix.add(a, current.matrix);
            output.genes[a] = current.genes;
        }

        // Get all annotations keys across datasets; we then concatenate
        // columns with the same name, or we just fill them with missings.
        let lengths = [];
        let annos = [];
        for (const d of dkeys) {
            let current = datasets[d];
            if (current.annotations !== null) {
                annos.push(current.annotations);
            } else {
                annos.push({});
            }
            lengths.push(current.matrix.numberOfColumns());
        }
        output.annotations = scran.combineArrayCollections(annos, { lengths: lengths });

        // Generating a block vector.
        let ncells = new Array(dkeys.length);
        dkeys.forEach((x, i) => { ncells[i] = datasets[x].matrix.numberOfColumns(); });
        blocks = scran.createBlock(ncells);
        output.block_ids = blocks;
        output.block_levels = dkeys;

        let nice_barr = new Array(blocks.length);
        blocks.forEach((x, i) => { nice_barr[i] = dkeys[x]; })
        output.annotations["__batch__"] = nice_barr;

    } catch (e) {
        utils.freeCache(blocks);
        utils.freeCache(output.matrix);
        throw e;
    } 

    return output;
}

async function load_datasets(matrices, sample_factor) {
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
                if (anno_batch.length != current.matrix.numberOfColumns()) {
                    throw new Error("length of sample factor '" + sample_factor + "' should be equal to the number of cells"); 
                }
                let converted = scran.convertBlock(anno_batch);
                blocks = converted.ids;
                block_levels = converted.levels;
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
        let output;
        try {
            output = bind_datasets(dkeys, datasets);
        } finally {
            // No need to hold references to the individual matrices
            // once the full matrix is loaded.
            for (const [k, v] of Object.entries(datasets)) {
                utils.freeCache(v.matrix);
            }
        }
        return output;
    }
}

function harvest_subset_indices(subset, total_ncols, annotations) {
    let keep = null;

    if (subset !== null) {
        // Applying the subset.
        if ("indices" in subset) {
            keep = subset.indices;
            if (total_ncols !== null) {
                utils.checkIndices(keep, total_ncols);
            }

        } else if ("field" in subset) {
            if (!(subset.field in annotations)) {
                throw new Error("failed to find 'subset.field' in the column annotations");
            }
            let anno = annotations[subset.field];

            keep = [];
            let allowed = new Set(subset.values);
            anno.forEach((x, i) => {
                if (allowed.has(x)) {
                    keep.push(i);
                }
            });

        } else {
            throw new Error("unrecognized specification for 'subset'");
        }
    }

    return keep;
}

function subset_datasets(cache, keep, { freeExisting = true } = {}) {
    if (keep === null) {
        return;
    }

    let replacement = new scran.MultiMatrix;
    try {
        for (const key of cache.matrix.available()) {
            let current = cache.matrix.get(key);
            replacement.add(key, scran.subsetColumns(current, keep));
        }
        if (freeExisting) {
            cache.matrix.free();
        }
        cache.matrix = replacement;
    } catch (e) {
        replacement.free();
        throw e;
    }

    if (cache.annotations !== null) {
        cache.annotations = scran.subsetArrayCollection(cache.annotations, keep);
    }

    if (cache.block_ids !== null) {
        let subblock = scran.subsetBlock(cache.block_ids, keep);
        try {
            let dropped = scran.dropUnusedBlock(subblock);
            cache.block_levels = dropped.map(x => cache.block_levels[x]);
            if (freeExisting) {
                cache.block_ids.free();
            }
            cache.block_ids = subblock;
        } catch (e) {
            subblock.free();
            throw e;
        }
    }
}

async function process_and_cache(new_matrices, sample_factor, subset) {
    let cache = {};

    try {
        cache = await load_datasets(new_matrices, sample_factor);
        let keep = harvest_subset_indices(subset, cache.matrix.numberOfColumns(), cache.annotations);
        subset_datasets(cache, keep);

        var gene_info_type = {};
        var gene_info = cache.genes["RNA"];
        for (const [key, val] of Object.entries(gene_info)) {
            gene_info_type[key] = scran.guessFeatures(val);
        }
        cache.gene_types = gene_info_type;

    } catch (e) {
        utils.freeCache(cache.matrix);
        utils.freeCache(cache.block_ids);
        throw e;
    }

    return cache;
}

/**************************
 ******** Loading *********
 **************************/

function createPermuter(perm) {
    return x => {
        let copy = x.slice();
        x.forEach((y, i) => {
            x[i] = copy[perm[i]];
        });
    };
}

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

    // Figuring out the subset.
    let subset = null;
    if ("subset" in phandle.children) {
        let shandle = phandle.open("subset");
        subset = {};

        if ("cells" in shandle.children) {
            let schandle = shandle.open("cells");
            if ("indices" in schandle.children) {
                subset.indices = schandle.open("indices", { load: true }).values;
            } else if ("field" in schandle.children) {
                subset.field = schandle.open("field", { load: true }).values[0];
                subset.values = schandle.open("values", { load: true }).values;
            } else {
                throw new Error("unrecognized specification for 'subset'");
            }
        }
    }

    parameters.subset = subset;

    // Loading matrix data.
    let cache = await process_and_cache(parameters.matrices, parameters.sample_factor, parameters.subset);

    // We need to do something if the permutation is not the same.
    let rhandle = ghandle.open("results");

    let perm = {};
    if (solofile) {
        if ("permutation" in rhandle.children) {
            // v1.0-v1.1
            let dhandle = rhandle.open("permutation", { load: true });
            let ids = new Int32Array(dhandle.values.length);
            dhandle.values.forEach((x, i) => { ids[x] = i; });
            perm.RNA = scran.updateRowIdentities(cache.matrix.get("RNA"), ids);
        } else if ("identities" in rhandle.children) {
            if (rhandle.children["identities"] == "DataSet") {
                // v1.2
                let dhandle = rhandle.open("identities", { load: true });
                perm.RNA = scran.updateRowIdentities(cache.matrix.get("RNA"), dhandle.values);
            } else {
                // v2.0
                let ihandle = rhandle.open("identities");
                for (const a of Object.keys(ihandle.children)) {
                    if (cache.matrix.has(a)) {
                        let dhandle = ihandle.open(a, { load: true });
                        perm[a] = scran.updateRowIdentities(cache.matrix.get(a), dhandle.values);
                    }
                }
            }
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

            let ref = cache.matrix.get("RNA").identities().sort();
            let old_ids2 = old_ids.slice().sort();
            for (var i = 0; i < old_ids2.length; i++) {
                if (ref[i] != old_ids2[i]) {
                    console.log([i, ref[i], old_ids2[i]]);
                    break;
                }
            }
            perm.RNA = scran.updateRowIdentities(cache.matrix.get("RNA"), old_ids);
        } else {
            if (rhandle.children["identities"] == "DataSet") {
                // v1.2+
                old_ids = rhandle.open("identities", { load: true }).values;
                perm.RNA = scran.updateRowIdentities(cache.matrix.get("RNA"), old_ids);
            } else {
                // v2.0
                let ihandle = rhandle.open("identities");
                for (const a of Object.keys(ihandle.children)) {
                    if (cache.matrix.has(a)) {
                        let dhandle = ihandle.open(a, { load: true });
                        perm[a] = scran.updateRowIdentities(cache.matrix.get(a), dhandle.values);
                    }
                }
            }
        }
    }

    let permuters = {};
    for (const a of cache.matrix.available()) {
        if (a in perm && perm[a] !== null) {
            permuters[a] = createPermuter(perm[a]); 
        } else {
            permuters[a] = x => {};
        }
    }

    /*
     * We could try to construct 'abbreviated', but there isn't really
     * any point because callers are expected to set 'matrices = null'
     * in their calls to 'compute()' on an unserialized analysis, so 
     * any setting of '#abbreviated' wouldn't even get used.
     */

    return { 
        state: new InputsState(parameters, cache),
        permuters: permuters
    };
}
