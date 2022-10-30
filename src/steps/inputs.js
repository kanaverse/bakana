import * as scran from "scran.js";
import * as bioc from "bioconductor";
import * as utils from "./utils/general.js";
import * as iutils from "../readers/index.js";
export const step_name = "inputs";

const RAW_SUBSET_OVERRIDE = "raw_subset_indices";

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
        utils.freeCache(this.#cache.matrix);
        utils.freeCache(this.#cache.raw_matrix);
        utils.freeCache(this.#cache.block_ids);
        utils.freeCache(this.#cache.raw_block_ids);
        utils.freeCache(this.#cache.multi_block_ids);
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
        // Just returning the entire object for the time being,
        // until downstream code is adjusted for the use of DFs.
        return this.#cache.genes[type]._columns;
    }

    fetchRowIds({ type = "RNA" } = {}) {
        return this.#cache.row_ids[type];
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
        if (annots === null || !annots.hasColumn(col)) {
            throw new Error(`${col} does not exist in the column annotations`);
        }

        // Make a copy, avoid accidental writes or transfers. 
        return bioc.CLONE(annots.column(col));
    }

    fetchBlock() {
        return this.#cache.block_ids;
    }

    fetchBlockLevels() {
        return this.#cache.block_levels;
    }

    fetchParameters() {
        // Cloning the parameters to avoid pass-by-reference behavior affecting the
        // InputsState object. We don't pass the files back here.
        let output = { ...this.#parameters };
        output.subset = this.constructor.#cloneSubset(output.subset);
        return output;
    }

    fetchDirectSubset() {
        if (RAW_SUBSET_OVERRIDE in this.#cache) {
            return this.#cache[RAW_SUBSET_OVERRIDE].slice();
        } else {
            return null;
        }
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

    /**
     * This method should not be called directly by users, but is instead invoked by {@linkcode runAnalysis}.
     * `datasets` is taken from the `matrices` argument in {@linkcode runAnalysis},
     * while `sample_factor` is taken from the property of the same name in the `inputs` property of the `parameters`.
     *
     * @param {object} datasets - An object containing data for one or more datasets.
     * Each property corresponds to a single dataset and its value should be a {@linkplain Dataset} object.
     * See the description of the argument of the same name in {@linkcode runAnalysis}.
     * @param {?string} sample_factor - Name of the column of the cell annotations specifying the sample of origin for each cell.
     * This is only used if a single count matrix is supplied.
     *
     * If `null`, all cells are assumed to originate from the same sample.
     * @param {?subset} subset - Object describing if any pre-analysis subsetting should be applied.
     * This should contain `field`, a string specifying a field of the column annotation.
     *
     * - For categorical variables, the object should also contain `values`, an array of allowed values for that annotation.
     *   Cells are only retained if they are associated with any of the allowable values for that annotation field.
     * - For continuous variables, the object should also contain `ranges`, an array of arrays.
     *   Each inner array should contain two numbers defining the start and end of a range.
     *   Ranges should be sorted and non-overlapping (boundaries excepted).
     *
     * If `subset` is `null`, no subsetting is performed and all cells are used in the downstream analysis.
     *
     * @return The object is updated with the new results.
     * A promise is returned that resolves to `null` once input loading is complete - this should be resolved before any downstream steps are run.
     */
    async compute(datasets, sample_factor, subset) {
        this.changed = false;

        // Don't bother proceeding with any of the below
        // if we're operating from a reloaded state.
        if (datasets !== null) {
            let tmp_abbreviated = {};
            for (const [key, val] of Object.entries(datasets)) {
                tmp_abbreviated[key] = { format: val.constructor.format(), details: val.abbreviate() };
            }

            if (utils.changedParameters(tmp_abbreviated, this.#abbreviated)) {
                await load_and_cache(datasets, this.#cache);
                this.#abbreviated = tmp_abbreviated;
                this.#cache.datasets = datasets;
                this.changed = true;
            }
        }

        if (this.changed || this.#parameters.sample_factor !== sample_factor) {
            block_and_cache(sample_factor, this.#cache);
            this.#parameters.sample_factor = sample_factor;
            this.changed = true;
        }

        if (this.changed || (!(RAW_SUBSET_OVERRIDE in this.#cache) && utils.changedParameters(subset, this.#parameters.subset))) {
            subset_and_cache(subset, this.#cache);
            this.#parameters.subset = this.constructor.#cloneSubset(subset);
            this.changed = true;
        }

        return null;
    }

    /******************************
     ******** Subsetting **********
     ******************************/

    static #cloneSubset(subset) {
        // We use a dedicated cloning function to handle Infs,
        // as these get converted to nulls by the JSON stringify.
        if (subset == null) {
            return subset;
        }

        let clone = { ...subset };
        if ("values" in clone) {
            clone.values = clone.values.slice();
        }

        if ("ranges" in clone) {
            clone.ranges = clone.ranges.map(x => x.slice());
        }

        return clone;
    }

    /**
     * Undo the effect of subsetting on an array of indices.
     *
     * @param {Array|TypedArray} indices - Array of column indices to the subsetted matrix.
     * Note that this will be modified in-place.
     *
     * @return Entries of `indices` are replaced with indices to the pre-subsetted matrix.
     */
    undoSubset(indices) {
        // Setting the subset to null, if the parameter-level subset hasn't
        // been set yet. This is because we might get indirectly called via
        // setDirectSubset() before compute() has been run.
        let subset = null;
        if ("subset" in this.#parameters) {
            subset = this.#parameters.subset;
        }

        let keep = harvest_subset_indices(subset, this.#cache);
        if (keep !== null) {
            indices.forEach((x, i) => { indices[i] = keep[x] });
        }
    }

    #configureIndices(indices, copy, onOriginal) {
        // scran.js's subset functions will pick up out-of-range indices.
        utils.checkIndices(indices, null);

        // We make a copy here to take ownership of the underlying memory,
        // otherwise any edits in the caller would mutate the new InputsState's
        // indices by reference.
        if (copy) {
            indices = indices.slice();
        }

        if (!onOriginal) {
            this.undoSubset(indices);
        }

        return indices;
    }

    /**
     * Unlike most of the other methods, `setDirectSubset` can be called on an InputsState before {@linkcode InputsState#compute compute}.
     * This means that a user can create the state object from {@linkcode createAnalysis},
     * specify a subset of cells via `setDirectSubset` on the (currently empty) InputsState object in `inputs`,
     * and then call {@linkcode runAnalysis} to execute an analysis on the desired subset of cells.
     * 
     * @param {TypedArray|Array} indices - Array containing the indices for the desired subset of cells.
     * This should be sorted and non-duplicate.
     * Any existing subset in this object will be overridden by `indices`.
     * @param {object} [options] - Optional parameters.
     * @param {boolean} [options.copy=true] - Whether to make a copy of `indices` before storing it inside the returned state object.
     * If `false`, it is assumed that the caller makes no further use of the passed `indices`.
     * @param {boolean} [options.onOriginal=false] - Whether `indices` contains indices on the original dataset or on the dataset in `state`.
     * This distinction is only relevant if the current InputsState object already contains a specified subset.
     * If `false`, the `indices` are assumed to refer to the already-subsetted dataset that exists in `state`;
     * if `true`, the `indices` are assumed to refer to the original dataset from which the subset in `state` was created.
     *
     * @return The dataset in this InputsState object is subsetted to the desired `indices`.
     */
    setDirectSubset(indices, { copy = true, onOriginal = false } = {}) {
        if (indices !== null) {
            this.#cache[RAW_SUBSET_OVERRIDE] = this.#configureIndices(indices, copy, onOriginal);
        } else {
            delete this.#cache[RAW_SUBSET_OVERRIDE];            
        }

        // If it's already got a matrix entry, we re-run it.
        if ("matrix" in this.#cache) {
            subset_and_cache(this.#parameters.subset, this.#cache);
            this.changed = true;
        }
    }

    createDirectSubset(indices, { copy = true, onOriginal = false } = {}) {
        let new_cache = {};
        new_cache[RAW_SUBSET_OVERRIDE] = this.#configureIndices(indices, copy, onOriginal);

        // Making explicit clones to take ownership.
        new_cache.raw_matrix = this.#cache.raw_matrix.clone();
        for (const x of [ "multi_block_ids", "raw_block_ids" ]) {
            if (x in this.#cache) {
                if (this.#cache[x] === null) {
                    new_cache[x] = null;
                } else {
                    new_cache[x] = this.#cache[x].clone();
                }
            }
        }

        // These can probably be copied directly, given that they are always
        // replaced wholesale in the various *_and_cache functions, rather than
        // being modified in-place.
        for (const x of [ "raw_annotations", "genes", "gene_types", "multi_block_levels", "raw_block_levels" ]) {
            if (x in this.#cache) {
                new_cache[x] = this.#cache[x];
            }
        }

        subset_and_cache(null, new_cache);

        let new_params = this.fetchParameters();
        new_params.subset = null;

        return new InputsState(new_params, new_cache, this.#abbreviated);
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
     *   For multiple datasets, the number of cells is defined as the sum across all datasets.
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

        let gene_info = {};
        for (const [k, v] of Object.entries(this.#cache.genes)) {
            let info = {};
            for (const c of v.columnNames()) {
                info[c] = bioc.CLONE(v.column(c));
            }
            gene_info[k] = info;
        }
        
        var output = {
            "num_cells": this.#cache.matrix.numberOfColumns(),
            "num_genes": ngenes,
            "genes": gene_info 
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

            // Make sure we're in sorted order, for consistency with
            // how the merge is done.
            let names = Object.keys(this.#cache.datasets);
            names.sort();

            let formats = [];
            let numbers = [];
            let fihandle = phandle.createGroup("files");
            let sofar = 0;

            for (const key of names) {
                let val = this.#cache.datasets[key];
                formats.push(val.constructor.format());

                let files = await val.serialize();
                numbers.push(files.length);

                for (const obj of files) {
                    let curhandle = fihandle.createGroup(String(sofar));
                    curhandle.writeDataSet("type", "String", [], obj.type);
                    curhandle.writeDataSet("name", "String", [], obj.file.name());

                    if (embeddedSaver === null) {
                        if (file2link == null) {
                            throw new Error("no valid linking function from 'setCreateLink'");
                        }
                        let id = await file2link(obj.file.content());
                        curhandle.writeDataSet("id", "String", [], id);
                    } else {
                        let saved = await embeddedSaver(obj.file.content(), obj.file.size());
                        if (!(typeof saved.offset == "number") || !(typeof saved.size == "number")) {
                            throw new Error("output of 'embeddedSaver' should contain 'offset' and 'size' numbers");
                        }
                        curhandle.writeDataSet("offset", "Uint32", [], saved.offset);
                        curhandle.writeDataSet("size", "Uint32", [], saved.size);
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

            if (this.#parameters.subset !== null || RAW_SUBSET_OVERRIDE in this.#cache) {
                let shandle = phandle.createGroup("subset");
                let schandle = shandle.createGroup("cells");

                if (RAW_SUBSET_OVERRIDE in this.#cache) {
                    schandle.writeDataSet("indices", "Int32", null, this.#cache[RAW_SUBSET_OVERRIDE]);
                } else if ("field" in this.#parameters.subset) {
                    schandle.writeDataSet("field", "String", [], this.#parameters.subset.field);

                    if ("values" in this.#parameters.subset) {
                        schandle.writeDataSet("values", "String", null, this.#parameters.subset.values);
                    } else {
                        let raw_ranges = this.#parameters.subset.ranges;
                        let ranges = [].concat(...raw_ranges);
                        check_subset_ranges(ranges);
                        schandle.writeDataSet("ranges", "Float64", [ranges.length/2, 2], ranges);
                    }
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

            // Looping through all available modalities.
            let ihandle = rhandle.createGroup("identities");
            for (const a of this.#cache.matrix.available()) {
                ihandle.writeDataSet(a, "Int32", null, this.#cache.row_ids[a]);
            }
        }

        return;
    }
}

/************************************
 ******* Internals - guessing *******
 ************************************/

// TODO: remove all of this and require explicit mention of which features are which.

// Exported for testing only.
export function commonFeatureTypes(genes) {
    let scores = {
        "symbol-mouse": [],
        "symbol-human": [],
        "ensembl-mouse": [],
        "ensembl-human": []
    };
    let fields = bioc.CLONE(scores);

    for (var i = 0; i < genes.length; i++) {
        let curgenes = genes[i];
        let best_scores = {};
        let best_fields = {};

        for (const k of curgenes.columnNames()) {
            let v = curgenes.column(k);
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
        if (v.length == genes.length) { // skipping if not represented in all entries.
            let nscore = v.reduce((a, b) => a * b);
            if (nscore > best_score) {
                best_score = nscore;
                best_type = k;
            }
        }
    }

    let best_fields = [];
    let best_features = null;

    if (best_type !== null) {
        let best_type_cols = fields[best_type];
        let best_features_sub = best_type.split("-");
        best_features = {
            type: best_features_sub[0],
            species: best_features_sub[1]
        };
        for (var i = 0; i < genes.length; i++) {
            best_fields.push(best_type_cols[i]);
        }
    }

    return {
        "best_type": best_features,
        "best_fields": best_fields
    };
}

export function guessDefaultModalities(loaded) {
    // Checking which modalities are available across all datasets.
    let available = {};
    for (var i = 0; i < loaded.length; i++) {
        let current = Object.keys(loaded[i].features);

        for (const modality of current) {
            let chosen = null;
            if (modality == "") {
                chosen = "RNA";
            } else if (modality.match(/gene expression/i)) {
                chosen = "RNA";
            } else if (modality.match(/^rna/i)) {
                chosen = "RNA";
            } else if (modality.match(/^gex/i)) {
                chosen = "RNA";
            } else if (modality.match(/^adt/i)) {
                chosen = "ADT";
            } else if (modality.match(/antibody capture/i)) {
                chosen = "ADT";
            }

            if (chosen !== null) {
                if (!(chosen in available)) {
                    available[chosen] = new Array(loaded.length);
                    available[chosen].fill(null);
                }
                if (available[chosen][i] == null) {
                    available[chosen][i] = modality;
                }
            }
        }
    }

    let choices = {};
    for (const [k, v] of Object.entries(available)) {
        let failed = false;
        for (const v2 of v) {
            if (v2 == null) {
                failed = true;
                break;
            }
        }
        if (!failed) {
            choices[k] = v;
        }
    }

    return choices;
}

/************************************
 ******* Internals - loading ********
 ************************************/

function bind_single_modality(modality, mapping, loaded) {
    let output = {};

    try {
        let genes = [];
        for (var i = 0; i < mapping.length; i++) {
            let mod_name = mapping[i];
            genes.push(loaded[i].features[mod_name]);
        }

        let result = commonFeatureTypes(genes);
        if (result.best_type === null) {
            throw new Error("no common feature types available across all datasets for '" + modality + "'");
        }

        let gnames = [];
        let mats = [];
        for (var i = 0; i < mapping.length; i++) {
            gnames.push(genes[i].column(result.best_fields[i]));
            mats.push(loaded[i].matrix.get(mapping[i]));
        }

        let merged = scran.cbindWithNames(mats, gnames);
        output.matrix = merged.matrix;

        // Extracting gene information from the first object. We won't make
        // any attempt at merging and deduplication across objects.
        output.features = bioc.SLICE(genes[0], merged.indices);
        output.row_ids = bioc.SLICE(loaded[0].row_ids[mapping[0]], merged.indices);

    } catch (e) {
        utils.freeCache(output.matrix);
        throw e;
    }

    return output;
}

function bind_datasets(names, loaded) {
    let common_modes = guessDefaultModalities(loaded);
    if (Object.keys(common_modes).length == 0) {
        throw new Error("failed to guess common modalities across all datasets");
    }

    let blocks;
    let output = { 
        matrix: new scran.MultiMatrix, 
        features: {},
        row_ids: {}
    };

    try {
        for (const [k, v] of Object.entries(common_modes)) {
            let current = bind_single_modality(k, v, loaded);
            output.matrix.add(k, current.matrix);
            output.features[k] = current.features;
            output.row_ids[k] = current.row_ids;
        }

        let annos = loaded.map(x => x.cells);
        output.cells = bioc.flexibleCombineRows(annos);

        // Generating a block vector.
        let ncells = new Array(loaded.length);
        loaded.forEach((x, i) => { ncells[i] = x.matrix.numberOfColumns(); });
        blocks = scran.createBlock(ncells);
        output.block_ids = blocks;
        output.block_levels = names;

        let nice_barr = new Array(blocks.length);
        blocks.forEach((x, i) => { nice_barr[i] = names[x]; })
        output.cells.$setColumn("__batch__", nice_barr);

    } catch (e) {
        utils.freeCache(blocks);
        utils.freeCache(output.matrix);
        throw e;
    } 

    return output;
}

function rename_dataset(single) {
    let common_modes = guessDefaultModalities([single]);
    if (Object.keys(common_modes).length == 0) {
        throw new Error("failed to guess common modalities across all datasets");
    }

    let output = { 
        matrix: new scran.MultiMatrix, 
        features: {},
        row_ids: {}
    };

    try {
        for (const [k, v] of Object.entries(common_modes)) {
            let current = v[0];
            output.matrix.add(k, single.matrix.get(current));
            output.features[k] = single.features[current];
            output.row_ids[k] = single.row_ids[current];
        }
    } catch (e) {
        scran.free(output.matrix);
        throw e;
    }

    output.cells = single.cells;
    output.block_ids = null;
    output.block_levels = null;

    return output;
}

async function load_datasets(datasets) {
    // Ensure we have a reproducible order; otherwise the batch
    // order becomes dependent on the JS engine's ordering.
    let names = Object.keys(datasets);
    names.sort();

    let loaded = [];
    try {
        for (const key of names) {
            // Too much hassle to convert this into a Promise.all(), because we
            // need to make sure it gets freed properly on failure.
            loaded.push(await datasets[key].load());
        }
    } catch (e) {
        // If any one fails, we free the rest.
        for (const x of loaded) {
            scran.free(x.matrix);
        }
        throw e;
    }

    let output;
    if (names.length == 1) {
        try {
            output = rename_dataset(loaded[0]);
        } catch (e) {
            scran.free(loaded[0].matrix);
            throw e;
        }
    } else {
        try {
            output = bind_datasets(names, loaded);
        } finally {
            // No need to hold references to the individual matrices once the
            // binding is complete, so we release them.
             for (const x of loaded) {
                scran.free(x.matrix);
            }
        }
    }

    return output;
}

/******************************************
 ******* Internals - miscellaneous ********
 ******************************************/

function harvest_subset_indices(subset, cache) {
    if (RAW_SUBSET_OVERRIDE in cache) {
        return cache[RAW_SUBSET_OVERRIDE];
    }

    if (subset == null) {
        return null;
    }

    if (!cache.raw_annotations.hasColumn(subset.field)) {
        throw new Error("failed to find '" + subset.field + "' in the column annotations");
    }

    let anno = cache.raw_annotations.column(subset.field);
    let keep = [];

    if ("values" in subset) {
        let allowed = new Set(subset.values);
        anno.forEach((x, i) => {
            if (allowed.has(x)) {
                keep.push(i);
            }
        });
    } else {
        // Check each entry to see whether it belongs to the range.
        // This is cheaper than sorting anything, assuming there 
        // aren't that many ranges.
        anno.forEach((x, i) => {
            for (const r of subset.ranges) {
                if (x >= r[0] && x <= r[1]) {
                    keep.push(i);
                    return;
                }
            }
        });
    }

    return keep;
}

function check_subset_ranges(ranges) { 
    if (ranges.length % 2 !== 0) {
        throw new Error("'ranges' should have two columns in 'subset'");
    }
    for (var i = 1; i < ranges.length; i++) {
        if (ranges[i] < ranges[i-1]) {
            throw new Error("'ranges' should be sorted in increasing order");
        }
    }
}

/************************************
 ******* Internals - caching ********
 ************************************/

async function load_and_cache(new_datasets, cache) {
    utils.freeCache(cache.raw_matrix);
    utils.freeCache(cache.matrix); // freeing this as well, to release all references and potentially release memory.
    utils.freeCache(cache.multi_block_ids);

    let res = await load_datasets(new_datasets);
    cache.raw_matrix = res.matrix;
    cache.row_ids = res.row_ids;
    cache.raw_annotations = res.cells;
    cache.multi_block_ids = res.block_ids;
    cache.multi_block_levels = res.block_levels;

    cache.genes = res.features;
    var gene_info_type = {};
    var gene_info = cache.genes["RNA"];
    for (const k of gene_info.columnNames()) {
        let v = gene_info.column(k);
        gene_info_type[k] = scran.guessFeatures(v);
    }
    cache.gene_types = gene_info_type;
}

function block_and_cache(sample_factor, cache) {
    utils.freeCache(cache.raw_block_ids);

    let blocks = null;
    let block_levels = null;

    if (sample_factor !== null) {
        // Single matrix with a batch factor.
        try {
            let anno_batch = cache.raw_annotations.column(sample_factor);
            if (anno_batch.length != cache.raw_matrix.numberOfColumns()) {
                throw new Error("length of sample factor '" + sample_factor + "' should be equal to the number of cells"); 
            }
            let converted = scran.convertBlock(anno_batch);
            blocks = converted.ids;
            block_levels = converted.levels;
        } catch (e) {
            utils.freeCache(blocks);
            throw e;
        }
    } else {
        if (cache.multi_block_ids !== null) { 
            // Creating a view so that freeing of this object is a no-op.
            // We're downstream of load_and_cache so any freeing of
            // multi_block_ids would require block_and_cache to rerun
            // anyway, so we don't have to worry about invalidation.
            blocks = cache.multi_block_ids.view();
        } else {
            blocks = null;
        }
        block_levels = cache.multi_block_levels;
    }

    cache.raw_block_ids = blocks;
    cache.raw_block_levels = block_levels;
}

function subset_and_cache(subset, cache) {
    utils.freeCache(cache.matrix);
    utils.freeCache(cache.block_ids);

    let keep = harvest_subset_indices(subset, cache);

    let new_annotations;
    let new_matrix;
    let new_block_ids;
    let new_block_levels;

    try {
        if (keep === null) {
            new_annotations = cache.raw_annotations;

            // Need to make a clone so that it can be freed independently of the original.
            // This is cheap as only the shared pointer is cloned, not the underlying data.
            new_matrix = cache.raw_matrix.clone();

            if (cache.raw_block_ids !== null) {
                // A view also works, given that we're downstream of the generating
                // process for raw_block_ids and thus our lifetime is always tied to it.
                new_block_ids = cache.raw_block_ids.view();
                new_block_levels = cache.raw_block_levels;
            } else {
                new_block_ids = null;
                new_block_levels = null;
            }

        } else {
            new_annotations = bioc.SLICE(cache.raw_annotations, keep);

            if (cache.raw_block_ids !== null) {
                new_block_ids = scran.subsetBlock(cache.raw_block_ids, keep);
                let dropped = scran.dropUnusedBlock(new_block_ids);
                new_block_levels = dropped.map(x => cache.raw_block_levels[x]);
            } else {
                new_block_ids = null;
                new_block_levels = null;
            }

            new_matrix = new scran.MultiMatrix;
            for (const key of cache.raw_matrix.available()) {
                let current = cache.raw_matrix.get(key);
                new_matrix.add(key, scran.subsetColumns(current, keep));
            }
        }

    } catch (e) {
        utils.freeCache(new_matrix);
        utils.freeCache(new_block_ids);
        throw e;
    }

    cache.annotations = new_annotations;
    cache.block_levels = new_block_levels;
    cache.block_ids = new_block_ids;
    cache.matrix = new_matrix;
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

async function unserialize_Dataset(format, all_files, loader) {
    if (!(format in iutils.availableReaders)) {
        throw new Error("unknown format '" + format + "' during unserialization");
    }
    let cls = iutils.availableReaders[format];

    let handles = [];
    for (const f of all_files) {
        let b;
        if (loader == null) {
            if (link2file == null) {
                throw new Error("no valid linking function from 'setResolveLink'");
            }
            b = await link2file(f.id);
        } else {
            b = await loader(f.offset, f.size);
        }
        let handle = new iutils.SimpleFile(b, { name: f.name }) 
        handles.push({ type: f.type, file: handle });
    }

    return await cls.unserialize(handles);
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
    let readers = {};
    let parameters = { sample_factor: null };
    let fohandle = phandle.open("format", { load: true });
    let solofile = (fohandle.shape.length == 0);

    if (solofile) {
        let format = fohandle.values[0];
        readers["default"] = await unserialize_Dataset(format, all_files, embeddedLoader);
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
            let format = formats[i];
            readers[sample_names[i]] = await unserialize_Dataset(format, curfiles, embeddedLoader);
        }
    }

    // Figuring out the subset.
    let subset = null;
    let raw_indices = null;
    if ("subset" in phandle.children) {
        let shandle = phandle.open("subset");

        if ("cells" in shandle.children) {
            let schandle = shandle.open("cells");
            if ("indices" in schandle.children) {
                raw_indices = schandle.open("indices", { load: true }).values;
            } else if ("field" in schandle.children) {
                subset = { field: schandle.open("field", { load: true }).values[0] };

                if ("values" in schandle.children) {
                    subset.values = schandle.open("values", { load: true }).values;
                } else {
                    let ranges = schandle.open("ranges", { load: true }).values;
                    check_subset_ranges(ranges);
                    let reranges = [];
                    for (var i = 0; i < ranges.length/2; i++) {
                        reranges.push([ ranges[2*i], ranges[2*i + 1] ]);
                    }
                    subset.ranges = reranges;
                }
            } else {
                throw new Error("unrecognized specification for 'subset'");
            }
        }
    }

    parameters.subset = subset;

    // Loading matrix data.
    let cache = { readers: readers };
    await load_and_cache(readers, cache);
    block_and_cache(parameters.sample_factor, cache);

    if (raw_indices !== null) {
        cache[RAW_SUBSET_OVERRIDE] = raw_indices;
    }
    subset_and_cache(parameters.subset, cache);

    // We need to do something if the permutation is not the same.
    let rhandle = ghandle.open("results");

    let perm = {};
    if (solofile) {
        if ("permutation" in rhandle.children) {
            // v1.0-v1.1
            let dhandle = rhandle.open("permutation", { load: true });
            let ids = new Int32Array(dhandle.values.length);
            dhandle.values.forEach((x, i) => { ids[x] = i; });
            perm.RNA = scran.updateRowIdentities(cache.row_ids["RNA"], ids);
        } else if ("identities" in rhandle.children) {
            if (rhandle.children["identities"] == "DataSet") {
                // v1.2
                let dhandle = rhandle.open("identities", { load: true });
                perm.RNA = scran.updateRowIdentities(cache.row_ids["RNA"], dhandle.values);
            } else {
                // v2.0
                let ihandle = rhandle.open("identities");
                for (const a of Object.keys(ihandle.children)) {
                    if (cache.matrix.has(a)) {
                        let dhandle = ihandle.open(a, { load: true });
                        perm[a] = scran.updateRowIdentities(cache.row_ids[a], dhandle.values);
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

            let ref = cache.row_ids["RNA"].slice().sort();
            let old_ids2 = old_ids.slice().sort();
            for (var i = 0; i < old_ids2.length; i++) {
                if (ref[i] != old_ids2[i]) {
                    console.log([i, ref[i], old_ids2[i]]);
                    break;
                }
            }
            perm.RNA = scran.updateRowIdentities(cache.row_ids["RNA"], old_ids);
        } else {
            if (rhandle.children["identities"] == "DataSet") {
                // v1.2+
                old_ids = rhandle.open("identities", { load: true }).values;
                perm.RNA = scran.updateRowIdentities(cache.row_ids["RNA"], old_ids);
            } else {
                // v2.0
                let ihandle = rhandle.open("identities");
                for (const a of Object.keys(ihandle.children)) {
                    if (cache.matrix.has(a)) {
                        let dhandle = ihandle.open(a, { load: true });
                        perm[a] = scran.updateRowIdentities(cache.row_ids[a], dhandle.values);
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
     * any point because callers are expected to set 'datasets = null'
     * in their calls to 'compute()' on an unserialized analysis, so 
     * any setting of '#abbreviated' wouldn't even get used.
     */

    return { 
        state: new InputsState(parameters, cache),
        permuters: permuters
    };
}

/**************************
 ******** Linking *********
 **************************/

var file2link = null;
var link2file = null;

/**
 * Specify a function to create links for data files.
 * By default, this only affects files for the MatrixMarket, H5AD and 10X formats, and is only used when linking is requested.
 *
 * @param {function} fun - Function that returns a linking idenfier to a data file.
 * The function should accept the following arguments:
 *
 * - A string specifying the type of the file, e.g., `"mtx"`, `"h5"`.
 * - A string containing the name of the file.
 * - An ArrayBuffer containing the file content (for browsers) or a string containing the file path (for Node.js),
 *
 * The function is expected to return a string containing some unique identifier to the file.
 * This is most typically used to register the file with some user-specified database system for later retrieval.
 *
 * @return `fun` is set as the global link creator for this step. 
 * The _previous_ value of the creator is returned.
 */
export function setCreateLink(fun) {
    let previous = file2link;
    file2link = fun;
    return previous;
}

/**
 * Specify a function to resolve links for data files.
 * By default, this only affects files for the MatrixMarket, H5AD and 10X formats, and is only used when links are detected.
 *
 * @param {function} fun - Function that accepts a string containing a linking idenfier and returns an ArrayBuffer containing the file content (for browsers) or a string containing the file path (for Node.js),
 * This is most typically used to retrieve a file from some user-specified database system.
 *
 * @return `fun` is set as the global resolver for this step. 
 * The _previous_ value of the resolver is returned.
 */
export function setResolveLink(fun) {
    let previous = link2file;
    link2file = fun;
    return previous;
}
