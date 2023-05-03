import * as scran from "scran.js";
import * as bioc from "bioconductor";
import * as utils from "./utils/general.js";
import * as iutils from "../readers/index.js";
export const step_name = "inputs";

const RAW_SUBSET_OVERRIDE = "raw_subset_indices";

/**
 * This step handles the loading of all datasets into memory.
 * This wraps various matrix initialization functions in [**scran.js**](https://github.com/kanaverse/scran.js),
 * depending on the format of the supplied datasets.
 *
 * Methods not documented here are not part of the stable API and should not be used by applications.
 * @hideconstructor
 */
export class InputsState {
    #parameters;
    #cache;
    #abbreviated;
    #preserve_dataset_cache;

    constructor(parameters = null, cache = null, abbreviated = null) {
        this.#parameters = (parameters === null ? {} : parameters);
        this.#cache = (cache === null ? {} : cache);
        this.#abbreviated = (abbreviated === null ? {} : abbreviated);
        this.#preserve_dataset_cache = false;
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

    /**
     * @return {external:MultiMatrix} A {@linkplain external:MultiMatrix MultiMatrix} object containing counts for one or more modalities.
     * Each modality is represented by a separate count matrix, where each row of the matrix represents a feature of that modality.
     * All matrices have the same number and ordering of cells in their columns.
     */
    fetchCountMatrix() {
        return this.#cache.matrix;
    }

    /**
     * @return {object} Object where each key is the name of a modality and each value is a {@linkplain external:DataFrame DataFrame}.
     * Each row of the DataFrame corresponds to a feature in that modality 
     * (i.e., a row in the corresponding matrix from {@linkcode InputsState#fetchCountMatrix fetchCountMatrix})
     * and each column represents a per-feature annotation field.
     */
    fetchFeatureAnnotations() {
        return this.#cache.genes;
    }

    /**
     * @return {object} Object where each key is the name of a modality and each value is an Int32Array.
     * Each entry of an Int32Array specifies the identity of the corresponding row of its count matrix from {@linkcode InputsState#fetchCountMatrix fetchCountMatrix}.
     */
    fetchRowIds() {
        return this.#cache.row_ids;
    }

    /**
     * @return {external:DataFrame} {@linkplain external:DataFrame DataFrame} containing per-cell annotations.
     * Each row of the DataFrame corresponds to a cell in {@linkcode InputsState#fetchCountMatrix fetchCountMatrix},
     * and each column represents a per-cell annotation field.
     *
     * Note that this considers all cells in the dataset before QC filtering - 
     * see {@linkcode QualityControlState#applyFilter QualityControlState.applyFilter} to obtain a filtered version of each column.
     */
    fetchCellAnnotations() {
        return this.#cache.annotations;
    }

    /**
     * @return {?Int32Array} Array of length equal to the number of cells in the dataset,
     * identifying the block to which each cell is assigned.
     * Alternatively `null`, if no blocking is performed.
     */
    fetchBlock() {
        return this.#cache.block_ids;
    }

    /**
     * @return {?Array} Array of names of the blocks, or `null` if no blocking is performed.
     */
    fetchBlockLevels() {
        return this.#cache.block_levels;
    }

    /**
     * @return {object} Object containing the parameters.
     */
    fetchParameters() {
        // Cloning the parameters to avoid pass-by-reference behavior affecting the
        // InputsState object. We don't pass the files back here.
        let output = { ...this.#parameters };
        output.subset = this.constructor.#cloneSubset(output.subset);
        return output;
    }

    fetchDatasets() {
        return this.#cache.datasets;
    }

    /**
     * @param {object} [options={}] - Optional parameters.
     * @param {boolean} [options.copy=true] - Whether to return a copy of the subsets to avoid pass-by-reference behaviors.
     *
     * @return {?Int32Array} Array containing the indices to use for direct subsetting -
     * see {@linkcode InputsState#setDirectSubset setDirectSubset} for more information.
     * Alternatively `null`, if direct subsetting is not performed.
     */
    fetchDirectSubset({ copy = true } = {}) {
        if (RAW_SUBSET_OVERRIDE in this.#cache) {
            let candidate = this.#cache[RAW_SUBSET_OVERRIDE];
            return (copy ? candidate.slice() : candidate);
        } else {
            return null;
        }
    }

    guessRnaFeatureTypes() {
        if (!("RNA" in this.#cache.genes)) {
            return null;
        }

        if (!("inferred_rna_types" in this.#cache)) {
            this.#cache.inferred_rna_types = utils.guessFeatureTypes(this.#cache.genes["RNA"]);
        }

        return this.#cache.inferred_rna_types;
    }

    /***************************
     ******** Compute **********
     ***************************/

    static defaults() {
        return {
            block_factor: null,
            subset: null
        };
    }

    /**
     * Allow each {@linkplain Dataset} reader (i.e., the `datasets` in {@linkcode InputsState#compute compute}) to cache any intermediate results during loading.
     * By default, this is disabled as caching increases memory usage of the analysis without any major runtime improvements to `compute` when the `datasets` do not change.
     *
     * Setting `cache = true` is only useful if the instances in `datasets` are to be re-used outside of **bakana**, or if they are to be re-used in `compute()` in different combinations. 
     * In such cases, there may be a runtime improvement that warrants the increase in memory usage.
     * If caching is used, the user is responsible for releasing cached resources via each instance's `clear()` method once they are no longer needed.
     *
     * @param {boolean} cache - Whether to allow {@linkplain Dataset} instances to cache their results.
     */
    enableDatasetCache(cache) {
        this.#preserve_dataset_cache = cache;
        return;
    }

    /**
     * This method should not be called directly by users, but is instead invoked by {@linkcode runAnalysis}.
     *
     * @param {object} datasets - An object containing data for one or more datasets.
     * Each property corresponds to a single dataset and its value should satisfy the {@linkplain Dataset} contract.
     * See the description of the argument of the same name in {@linkcode runAnalysis}.
     * @param {object} parameters - Parameter object, equivalent to the `inputs` property of the `parameters` in {@linkcode runAnalysis}.
     * @param {?string} parameters.block_factor - Name of the column of the cell annotations specifying the sample of origin for each cell.
     * This is only used if a single count matrix is supplied.
     *
     * If `null`, all cells are assumed to originate from the same sample.
     *
     * If any entries of the blocking factor are invalid (i.e., `null`), they are removed from downstream analyses.
     * @param {?subset} parameters.subset - Object describing if any pre-analysis subsetting should be applied.
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
    async compute(datasets, parameters) {
        let { block_factor, subset } = parameters;
        this.changed = false;

        // Don't bother proceeding with any of the below
        // if we're operating from a reloaded state.
        if (datasets !== null) {
            let tmp_abbreviated = {};
            for (const [key, val] of Object.entries(datasets)) {
                tmp_abbreviated[key] = { format: val.constructor.format(), details: val.abbreviate() };
            }

            if (utils.changedParameters(tmp_abbreviated, this.#abbreviated)) {
                await load_and_cache(datasets, this.#cache, this.#preserve_dataset_cache);
                this.#abbreviated = tmp_abbreviated;
                this.#cache.datasets = { ...datasets }; // making a deep-ish copy to avoid pass-by-reference links.
                delete this.#cache.inferred_rna_types;
                this.changed = true;
            }
        }

        if (this.changed || this.#parameters.block_factor !== block_factor) {
            block_and_cache(block_factor, this.#cache);
            this.#parameters.block_factor = block_factor;
            this.changed = true;
        }

        // final condition handles loss of 'matrix' when setDirectSubset() is called.
        if (this.changed || (!(RAW_SUBSET_OVERRIDE in this.#cache) && utils.changedParameters(subset, this.#parameters.subset)) || !("matrix" in this.#cache)) { 
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
     * This works with all subset specifications, e.g., via `parameters.subset` in {@linkcode InputsState#compute compute}, with {@linkcode InputsState#setDirectSubset setDirectSubset},
     * or even from the implicit subsetting when the factor specified in `parameters.block` contains `null` entries.
     *
     * @param {Array|TypedArray} indices - Array of column indices to the subsetted matrix.
     * Note that this will be modified in-place.
     *
     * @return Entries of `indices` are replaced with indices to the pre-subsetted matrix.
     */
    undoSubset(indices) {
        if ("matrix" in this.#cache) {
            let max_index = this.fetchCountMatrix().numberOfColumns();
            for (const x of indices) {
                if (x < 0 || x >= max_index) {
                    throw new Error("entries of 'indices' should be less than the number of cells in the dataset");
                }
            }
        }

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

        // Flag that it needs to be rerun.
        scran.free(this.#cache.matrix);
        delete this.#cache.matrix;
    }

    createDirectSubset(indices, { copy = true, onOriginal = false } = {}) {
        let new_cache = {};
        new_cache[RAW_SUBSET_OVERRIDE] = this.#configureIndices(indices, copy, onOriginal);

        // Need to manually copy everything in 'this.#cache' that is set in
        // load_and_cache or block_and_cache.

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
        for (const x of [ "row_ids", "raw_annotations", "genes", "multi_block_levels", "raw_block_levels" ]) {
            if (x in this.#cache) {
                new_cache[x] = this.#cache[x];
            }
        }

        subset_and_cache(null, new_cache);

        let new_params = this.fetchParameters();
        new_params.subset = null;

        return new InputsState(new_params, new_cache, this.#abbreviated);
    }
}

/************************************
 ******* Internals - loading ********
 ************************************/

const known_modalities = [ "RNA", "ADT", "CRISPR" ];

function bind_single_modality(modality, loaded) {
    let output = {};

    try {
        let gnames = [];
        let mats = [];
        for (var i = 0; i < loaded.length; i++) {
            mats.push(loaded[i].matrix.get(modality));

            let primary_id = loaded[i].primary_ids[modality];
            if (primary_id == null) {
                throw new Error("modality '" + modality + "' lacks a primary identifier for dataset " + String(i));
            }
            gnames.push(primary_id);
        }

        let merged = scran.cbindWithNames(mats, gnames);
        output.matrix = merged.matrix;

        // Extracting gene information from the first object. We won't make
        // any attempt at merging and deduplication across objects.
        output.features = bioc.SLICE(loaded[0].features[modality], merged.indices);
        output.row_ids = bioc.SLICE(loaded[0].row_ids[modality], merged.indices);

    } catch (e) {
        utils.freeCache(output.matrix);
        throw e;
    }

    return output;
}

function bind_datasets(names, loaded) {
    let common_modes = [];
    for (const mod of known_modalities) {
        let okay = true;
        for (const l of loaded) {
            if (!l.matrix.has(mod)) {
                okay = false;
                break;
            }
        }
        if (okay) {
            common_modes.push(mod);
        }
    }

    if (common_modes.length == 0) {
        throw new Error("failed to find common modalities across all datasets");
    }

    let blocks;
    let output = { 
        matrix: new scran.MultiMatrix, 
        features: {},
        row_ids: {}
    };

    try {
        for (const k of common_modes) {
            let current = bind_single_modality(k, loaded);
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
    let modalities = single.matrix.available();
    if (modalities.length == 0) {
        throw new Error("");
    }

    let output = { 
        matrix: new scran.MultiMatrix, 
        features: {},
        row_ids: {}
    };

    try {
        for (const k of known_modalities) {
            if (!single.matrix.has(k)) {
                continue;
            }

            output.matrix.add(k, single.matrix.get(k));
            output.features[k] = single.features[k];
            output.row_ids[k] = single.row_ids[k];
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

async function load_datasets(datasets, preserve_dataset_cache) {
    // Ensure we have a reproducible order; otherwise the batch
    // order becomes dependent on the JS engine's ordering.
    let names = Object.keys(datasets);
    names.sort();

    let loaded = [];
    try {
        for (const key of names) {
            // Too much hassle to convert this into a Promise.all(), because we
            // need to make sure it gets freed properly on failure.
            loaded.push(await datasets[key].load({ cache: preserve_dataset_cache }));
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

const invalid_block_id = -1;

function harvest_subset_indices(subset, cache) {
    let keep;

    if (RAW_SUBSET_OVERRIDE in cache) {
        keep = cache[RAW_SUBSET_OVERRIDE];
    } else if (subset == null) {
        keep = null;
    } else {
        if (!cache.raw_annotations.hasColumn(subset.field)) {
            throw new Error("failed to find '" + subset.field + "' in the column annotations");
        }

        let anno = cache.raw_annotations.column(subset.field);
        keep = [];

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
    }

    // Filter out invalid block IDs. Note that this might get called
    // before compute() is run (via undoSubset), so we need to protect
    // against the case where the raw_block_ids has not been set yet.
    if ("raw_block_ids" in cache && cache.raw_block_ids !== null) {
        let bids = cache.raw_block_ids.array();

        let keep2 = [];
        if (keep !== null) {
            for (const i of keep) {
                if (bids[i] !== invalid_block_id) {
                    keep2.push(i);
                }
            }
        } else {
            for (var i = 0; i < bids.length; i++) {
                if (bids[i] !== invalid_block_id) {
                    keep2.push(i);
                }
            }
        }
        keep = keep2;
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

async function load_and_cache(new_datasets, cache, preserve_dataset_cache) {
    utils.freeCache(cache.raw_matrix);
    utils.freeCache(cache.matrix); // freeing this as well, to release all references and potentially release memory.
    utils.freeCache(cache.multi_block_ids);

    let res = await load_datasets(new_datasets, preserve_dataset_cache);
    cache.raw_matrix = res.matrix;
    cache.row_ids = res.row_ids;
    cache.raw_annotations = res.cells;
    cache.multi_block_ids = res.block_ids;
    cache.multi_block_levels = res.block_levels;
    cache.genes = res.features;
}

function block_and_cache(block_factor, cache) {
    utils.freeCache(cache.raw_block_ids);

    let blocks = null;
    let block_levels = null;

    if (block_factor !== null) {
        // Single matrix with a batch factor.
        try {
            let anno_batch = cache.raw_annotations.column(block_factor);
            if (anno_batch.length != cache.raw_matrix.numberOfColumns()) {
                throw new Error("length of blocking factor '" + block_factor + "' should be equal to the number of cells"); 
            }
            let converted = scran.factorize(anno_batch, { action: "none", placeholder: invalid_block_id });
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

export function updateRowIdentities(current, old) {
    if (current.length == old.length) {
        if (current.every((x, i) => x == old[i])) {
            return y => y;
        }
    }

    let mapping = {};
    old.forEach((x, i) => {
        mapping[x] = i;
    });

    let perm = new Int32Array(current.length);
    current.forEach((x, i) => {
        if (x in mapping) {
            perm[i] = mapping[x];
        } else {
            perm[i] = -1;
        }
    });

    return y => {
        let copy = new y.constructor(perm.length);
        perm.forEach((i, j) => {
            if (i == -1) {
                if (copy instanceof Array) {
                    copy[j] = null;
                } else if (copy instanceof Float64Array) {
                    copy[j] = Number.NaN;
                } else {
                    copy[j] = -1; // dunno what else to do here.
                }
            } else {
                copy[j] = y[i];
            }
        });
        return copy;
    };
}

function extract_serialized_files(handle) {
    let kids = handle.children;
    let all_files = new Array(kids.length);

    for (const x of Object.keys(kids)) {
        let current = handle.open(x);

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

    return all_files;
}

async function unserialize_Dataset(format, all_files, all_options, loader) {
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

    return await cls.unserialize(handles, all_options);
}

export async function unserialize(handle, embeddedLoader) {
    let ghandle = handle.open("inputs");
    let phandle = ghandle.open("parameters");

    let readers = {};
    let parameters = { block_factor: null };
    let solofile = false; // legacy argument.

    if ("datasets" in phandle.children) {
        let dhandle = phandle.open("datasets");
        for (const k of Object.keys(dhandle.children)) {
            let curdhandle = dhandle.open(k);
            let format = curdhandle.open("format", { load: true }).values[0];
            let name = curdhandle.open("name", { load: true }).values[0];

            let fihandle = curdhandle.open("files");
            let curfiles = extract_serialized_files(fihandle);

            let options = {};
            if ("options" in curdhandle.children) {
                options = JSON.parse(curdhandle.open("options", { load: true }).values[0]);
            }

            readers[name] = await unserialize_Dataset(format, curfiles, options, embeddedLoader);
        }

        if ("block_factor" in phandle.children) {
            parameters.block_factor = phandle.open("block_factor", { load: true }).values[0];
        }

    } else {
        // Extracting the files.
        let fihandle = phandle.open("files");
        let all_files = extract_serialized_files(fihandle);

        // Extracting the format and organizing the files.
        let fohandle = phandle.open("format", { load: true });
        solofile = (fohandle.shape.length == 0);

        if (solofile) {
            let format = fohandle.values[0];
            readers["default"] = await unserialize_Dataset(format, all_files, {}, embeddedLoader);
            if ("sample_factor" in phandle.children) {
                parameters.block_factor = phandle.open("sample_factor", { load: true }).values[0];
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
                readers[sample_names[i]] = await unserialize_Dataset(format, curfiles, {}, embeddedLoader);
            }
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
    block_and_cache(parameters.block_factor, cache);

    if (raw_indices !== null) {
        cache[RAW_SUBSET_OVERRIDE] = raw_indices;
    }
    subset_and_cache(parameters.subset, cache);

    // We need to do something if the permutation is not the same.
    let rhandle = ghandle.open("results");
    let perm = {};

    if ("feature_identities" in rhandle.children) { 
        // v3.0
        let ihandle = rhandle.open("feature_identities");
        for (const a of Object.keys(ihandle.children)) {
            if (cache.matrix.has(a)) {
                let dhandle = ihandle.open(a, { load: true });
                perm[a] = updateRowIdentities(cache.row_ids[a], dhandle.values);
            }
        }
    } else {
        if (solofile) {
            if ("permutation" in rhandle.children) {
                // v1.0-v1.1
                let dhandle = rhandle.open("permutation", { load: true });
                let ids = new Int32Array(dhandle.values.length);
                dhandle.values.forEach((x, i) => { ids[x] = i; });
                perm.RNA = updateRowIdentities(cache.row_ids["RNA"], ids);
            } else if ("identities" in rhandle.children) {
                if (rhandle.children["identities"] == "DataSet") {
                    // v1.2
                    let dhandle = rhandle.open("identities", { load: true });
                    perm.RNA = updateRowIdentities(cache.row_ids["RNA"], dhandle.values);
                } else {
                    // v2.0
                    let ihandle = rhandle.open("identities");
                    for (const a of Object.keys(ihandle.children)) {
                        if (cache.matrix.has(a)) {
                            let dhandle = ihandle.open(a, { load: true });
                            perm[a] = updateRowIdentities(cache.row_ids[a], dhandle.values);
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
                perm.RNA = updateRowIdentities(cache.row_ids["RNA"], old_ids);
            } else {
                if (rhandle.children["identities"] == "DataSet") {
                    // v1.2+
                    old_ids = rhandle.open("identities", { load: true }).values;
                    perm.RNA = updateRowIdentities(cache.row_ids["RNA"], old_ids);
                } else {
                    // v2.0
                    let ihandle = rhandle.open("identities");
                    for (const a of Object.keys(ihandle.children)) {
                        if (cache.matrix.has(a)) {
                            let dhandle = ihandle.open(a, { load: true });
                            perm[a] = updateRowIdentities(cache.row_ids[a], dhandle.values);
                        }
                    }
                }
            }
        }
    }

    // Any missing modalities, for whatever reason.
    for (const a of cache.matrix.available()) {
        if (!(a in perm)) {
            perm[a] = y => y;
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
        permuters: perm
    };
}

/**************************
 ******** Linking *********
 **************************/

var file2link = null;
var link2file = null;

/**
 * Specify a function to create links for data files.
 *
 * @param {function} fun - Function that accepts:
 *
 * - `format`: the string containing the format of the dataset that owns the file.
 * - `file`: a {@linkplain SimpleFile} representing the file contents.
 *
 * It should return a string containing some unique identifier to the file.
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
 *
 * @param {function} fun - Function that accepts a string containing a linking idenfier and returns any value that can be used in the {@linkplain SimpleFile} constructor
 * i.e., a Uint8Array, File (on browser) or string containing a file path (on Node.js).
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
