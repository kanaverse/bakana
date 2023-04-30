import * as scran from "scran.js";
import * as bioc from "bioconductor";
import * as utils from "./utils/general.js";
import * as rutils from "../readers/index.js";
import * as inputs_module from "./inputs.js";
import * as norm_module from "./rna_normalization.js";

const baseUrl = "https://github.com/kanaverse/singlepp-references/releases/download/2023-04-28";

export const step_name = "cell_labelling";

/************************************
 ****** Internals for loading *******
 ************************************/

var download_fun  = utils.defaultDownload;

function set_download(fun) {
    let previous = download_fun;
    download_fun = fun;
    return previous;
}

async function acquire_file(name, suffix) {
    let full = name + "_" + suffix;
    let b = await download_fun(baseUrl + "/" + full);
    return new rutils.SimpleFile(b, { name: full })
}

const all_loaded = {};

function flush_prepared(cache) {
    if ("prepared" in cache) {
        for (const v of Object.values(cache.prepared)) {
            v.built.raw.free();
        }
        delete cache.prepared;
    }
}

async function process_genes(file) {
    let gene_lines = await rutils.readLines2(file.content(), { compression: "gz" }); // gene names
    let acquired = [];

    for (const x of gene_lines) {
        let val = null;
        if (x !== "") {
            val = x.split("\t");
            if (val.length == 1) {
                val = val[0];
            }
        }
        acquired.push(val);
    }

    return acquired;
}

async function load_reference(name, gene_id_type) {
    let gene_suffix = "genes_" + gene_id_type.toLowerCase() + ".csv.gz";

    if (name in all_loaded) {
        let output = all_loaded[name];
        let known_genes = output.genes;
        if (!(gene_id_type in known_genes)) {
            known_genes[gene_id_type] = await process_genes(await acquire_file(name, gene_suffix));
        }
        return output;
    }

    const suffixes = [ 
        "labels_fine.csv.gz",
        "label_names_fine.csv.gz",
        "markers_fine.gmt.gz",
        "matrix.csv.gz",
        gene_suffix
    ];

    let contents = await Promise.all(suffixes.map(x => acquire_file(name, x)));

    let loaded;
    let stored;
    try {
        loaded = scran.loadLabelledReferenceFromBuffers(
            contents[3].buffer(), // rank matrix
            contents[2].buffer(), // markers
            contents[0].buffer()  // label per sample
        );

        let labels = await rutils.readLines2(contents[1].content(), { compression: "gz" }); // full label names
        stored = {
            "raw": loaded, 
            "labels": labels,
            "genes": {}
        };

        stored.genes[gene_id_type] = await process_genes(contents[4]);
        all_loaded[name] = stored;

    } catch (e) {
        utils.freeCache(loaded);
        throw e;
    }

    return stored;
}

function flush_loaded() {
    for (const [k, v] of Object.entries(all_loaded)) {
        v.raw.free();
        delete all_loaded[k];
    }
}

/*************************************
 ****** Internals for building *******
 *************************************/

const available_references = {
    "9606": [ "BlueprintEncode", "DatabaseImmuneCellExpression", "HumanPrimaryCellAtlas", "MonacoImmune", "NovershternHematopoietic" ],
    "10090": [ "ImmGen", "MouseRNAseq" ]
};

function internal_build_reference(name, gene_ids, gene_id_type) {
    let built;
    let output;
    try {
        let current = all_loaded[name];
        let loaded = current.raw;

        if (!(gene_id_type in current.genes)) {
            throw new Error("unknown gene type '" + gene_id_type + "'");
        }
        let chosen_ids = current.genes[gene_id_type];

        built = scran.buildLabelledReference(gene_ids, loaded, chosen_ids); 
        output = {
            "loaded": current,
            "built": {
                "features": chosen_ids,
                "raw": built
            }
        };

    } catch (e) {
        utils.freeCache(built);
        throw e;
    }

    return output;
}

async function build_reference(cache, references, automatic, species, gene_id_column, gene_id_type, old_parameters, annofun, guessfun) {
    if (
        automatic !== old_parameters.automatic ||
        utils.changedParameters(references, old_parameters.references) ||
        (
            !automatic &&
            (
                species !== old_parameters.species ||
                gene_id_column !== old_parameters.gene_id_column ||
                gene_id_type !== old_parameters.gene_id_type
            )
        )
    ) {
        let species2 = species;
        let gene_id_column2 = gene_id_column;
        let gene_id_type2 = gene_id_type;

        if (automatic) {
            let auto = CellLabellingState.configureFeatureParameters(guessfun());
            species2 = auto.species;
            gene_id_column2 = auto.gene_id_column;
            gene_id_type2 = auto.gene_id_type;
        }

        let allowable = new Set;
        for (const s of species2) {
            if (s in available_references) {
                available_references[s].forEach(x => { allowable.add(x); });
            }
        }

        // Building each individual reference.
        let feats = annofun();
        let gene_ids = (gene_id_column2 == null ? feats.rowNames() : feats.column(gene_id_column2));
        cache.gene_ids = gene_ids;

        let valid = {};
        if (gene_ids !== null) {
            for (const ref of references) {
                if (allowable.has(ref)) {
                    await load_reference(ref, gene_id_type2);
                    valid[ref] = internal_build_reference(ref, gene_ids, gene_id_type2);
                }
            }
        }

        flush_prepared(cache);
        cache.prepared = valid;

        // Building an integrated reference, if necessary.
        let used_refs = Object.keys(valid);
        if (used_refs.length > 1) {
            let arr = Object.values(valid);
            let loaded = arr.map(x => x.loaded.raw);
            let feats = arr.map(x => x.built.features);
            let built = arr.map(x => x.built.raw);

            utils.freeCache(cache.integrated);
            cache.integrated = scran.integrateLabelledReferences(gene_ids, loaded, feats, built);
        } else {
            utils.freeCache(cache.integrated);
            delete cache.integrated;
        }
        cache.used_refs = used_refs;

       return true;
    }

    return false;
}

function transplant_parameters(references, automatic, species, gene_id_column, gene_id_type, parameters) {
    parameters.references = bioc.CLONE(references); // make a copy to avoid pass-by-reference behavior.
    parameters.automatic = automatic;
    parameters.species = bioc.CLONE(species);
    parameters.gene_id_column = gene_id_column;
    parameters.gene_id_type = gene_id_type;
}

/************************************
 ****** Internals for compute *******
 ************************************/

function transform_results(names, results, assigned) {
    let nclusters = results.numberOfCells();
    let ntargets = names.length;
    let output = new Array(nclusters);

    for (var r = 0; r < nclusters; r++) {
        let all_scores = {};
        let cscores = results.scoresForCell(r);
        for (var l = 0; l < ntargets; l++) {
            all_scores[names[l]] = cscores[l];
        }
        output[r] = { best: names[assigned[r]], all: all_scores };
    }

    return output;
}

function assign_labels(x, cache) {
    let matrix = x;
    let temp_cluster_means;
    let temp_matrix;

    // Converting marker results into means.
    if (x instanceof scran.ScoreMarkersResults) {
        let ngroups = x.numberOfGroups();
        let ngenes = (cache.gene_ids !== null ? cache.gene_ids.length : null);

        if (cache.gene_ids === null) {
            matrix = null;                
        } else {
            let ngenes = cache.gene_ids.length;

            // Creating a column-major array of mean vectors for each cluster.
            temp_cluster_means = scran.createFloat64WasmArray(ngroups * ngenes);
            for (var g = 0; g < ngroups; g++) {
                let means = x.means(g, { copy: false }); // Warning: direct view in wasm space - be careful.
                let cluster_array = temp_cluster_means.array();
                cluster_array.set(means, g * ngenes);
            }

            temp_matrix = scran.ScranMatrix.createDenseMatrix(ngenes, ngroups, temp_cluster_means, { columnMajor: true, copy: false });
            matrix = temp_matrix;
        }
    }

    // Running classifications; this is a no-op if gene_ids = null as 'valid' should be empty.
    let valid = cache.prepared;
    let results = { per_reference: {} };
    let raw = {};
    for (const [key, ref] of Object.entries(valid)) {
        let current = scran.labelCells(matrix, ref.built.raw);
        raw[key] = current;
        results.per_reference[key] = transform_results(ref.loaded.labels, current, current.predictedLabels({ copy: false }));
    }

    if ("integrated" in cache) {
        let single_results = [];
        for (const key of cache.used_refs) { // enforce correct order.
            single_results.push(raw[key]);
        }

        let current = scran.integrateCellLabels(matrix, single_results, cache.integrated);
        results.integrated = transform_results(cache.used_refs, current, current.predictedReferences({ copy: false }));
        current.free();
    }

    for (const v of Object.values(raw)) {
        v.free();
    }
    utils.freeCache(temp_matrix);
    utils.freeCache(temp_cluster_means);

    return results;
}

/********************
 ****** State *******
 ********************/

/**
 * Cell labelling involves assigning cell type labels to clusters using the [**SingleR** algorithm](https://github.com/LTLA/CppSingleR),
 * based on [pre-formatted reference expression profiles](https://github.com/clusterfork/singlepp-references).
 * This wraps [`labelCells`](https://kanaverse.github.io/scran.js/global.html#labelCells)
 * and related functions from [**scran.js**](https://github.com/kanaverse/scran.js).
 *
 * In theory, we could do this at the single-cell level, but we use clusters instead to expedite the computation and simplify interpretation.
 * If multiple references are requested, we will use each for assignment before attempting to choose the best label for each cluster across references.
 *
 * Methods not documented here are not part of the stable API and should not be used by applications.
 * @hideconstructor
 */
export class CellLabellingState {
    #inputs;
    #normalized;
    #parameters;
    #cache;

    constructor(inputs, normalized, parameters = null, cache = null) {
        if (!(inputs instanceof inputs_module.InputsState)) {
            throw new Error("'inputs' should be a State object from './inputs.js'");
        }
        this.#inputs = inputs;

        if (!(normalized instanceof norm_module.RnaNormalizationState)) {
            throw new Error("'normalized' should be a RnaNormalizationState object from './rna_normalization.js'");
        }
        this.#normalized = normalized;

        this.#parameters = (parameters === null ? {} : parameters);
        this.#cache = (cache === null ? {} : cache);
        this.changed = false;
    }

    free() {
        flush_prepared(this.#cache);
    }

    /***************************
     ******** Getters **********
     ***************************/

    valid() {
        let mat = this.#inputs.fetchCountMatrix();
        return mat.has("RNA");
    }

    /**
     * @return {object} Object containing the parameters.
     */
    fetchParameters() {
        // Avoid any pass-by-reference activity.
        let out = { ...this.#parameters };
        out.references = bioc.CLONE(out.references);
        out.species = bioc.CLONE(out.species);
        return out;
    }

    /**
     * @return {object} Object where each key is the name of a reference and each value is the number of shared features between the test and reference daatasets.
     */
    fetchNumberOfSharedFeatures() {
        let output = {};
        for (const key of this.#cache.used_refs) {
            output[key] = this.#cache.prepared[key].built.raw.sharedFeatures();
        }
        return output;
    }

    /****************************
     ******** Defaults **********
     ****************************/

    /**
     * @return {object} Default parameters that may be modified and fed into {@linkcode CellLabellingCore#compute compute}.
     */
    static defaults() {
        return {
            references: null,
            automatic: true,
            species: [],
            gene_id_column: null,
            gene_id_type: "ENSEMBL"
        };
    }

    static configureFeatureParameters(guesses) {
        let best_key = null;
        let best = { type: "symbol", species: "human", confidence: 0 };

        if ("row_names" in guesses) {
            let val = guesses.row_names;
            if (val.confidence > best.confidence) {
                best = val;
            }
        }

        for (const [key, val] of Object.entries(guesses.columns)) {
           if (val.confidence > best.confidence) {
                best = val;
                best_key = key;
            }
        }

        return {
            gene_id_column: best_key,
            species: [best.species],
            gene_id_type: best.type.toUpperCase()
        };
    }

    /***************************
     ******** Remotes **********
     ***************************/

    /**
     * Available references for each species.
     * Each key is a taxonomy ID and each value is an array of strings containing the names of references for that species.
     * @type {object}
     */
    static availableReferences = available_references;

    /**
     * Flush all cached references.
     *
     * By default, {@linkcode CellLabellingState#compute compute} will cache the loaded references in a global cache for re-use across {@linkplain CellLabellingState} instances.
     * These cached references are not tied to any single instance and will not be removed by garbage collectors or by {@linkcode freeAnalysis}.
     * Rather, this function should be called to release the relevant memory.
     */
    static flush() {
        flush_loaded();
        return;
    }

    /**
     * Specify a function to download references for the cell labelling step.
     *
     * @param {function} fun - Function that accepts a single string containing a URL and returns any value that can be used in the {@linkplain SimpleFile} constructor.
     * This is most typically a Uint8Array of that URL's contents, but it can also be a path to a locally cached file on Node.js.
     *
     * @return `fun` is set as the global downloader for this step. 
     * The _previous_ value of the downloader is returned.
     */
    static setDownload(fun) {
        return set_download(fun);
    }

    /***************************
     ******** Compute **********
     ***************************/

    /**
     * This method should not be called directly by users, but is instead invoked by {@linkcode runAnalysis}.
     *
     * @param {object} parameters - Parameter object, equivalent to the `cell_labelling` property of the `parameters` of {@linkcode runAnalysis}.
     * @param {?Array} parameters.references - Array of strings specifying the names of the reference datasets, see {@linkcode CellLabellingState.availableReferences availableReferences} for more details.
     * If `null`, all reference datasets from all species are used.
     * @param {boolean} parameters.automatic - Automatically choose feature-based parameters based on the feature annotation for the RNA modality.
     * If `true`, the column of the annotation that best matches human/mouse Ensembl/symbols is identified and used to set `species`, `gene_id_column` and `gene_id_type`.
     * @param {Array} parameters.species - Array of strings specifying zero, one or more species involved in this dataset.
     * Each entry should be a taxonomy ID (e.g. `"9606"`, `"10090"`) as specified in {@linkcode CellLabellingState.availableReferences availableReferences}.
     * This is used internally to filter `references` to the entries relevant to these species. 
     * Ignored if `automatic = true`.
     * @param {?(string|number)} parameters.gene_id_column - Name or index of the column of the RNA entry of {@linkcode InputsState#fetchFeatureAnnotations InputsState.fetchFeatureAnnotations} containing the identity of each gene. 
     * If `null`, identifiers are taken from the row names.
     * Ignored if `automatic = true`.
     * @param {string} parameters.gene_id_type - Type of feature identifier in `gene_id_column`.
     * This should be one of `"ENSEMBL"`, `"SYMBOL"` or `"ENTREZ"`
     * Ignored if `automatic = true`.
     *
     * @return The object is updated with the new results.
     * @async
     */
    async compute(parameters) {
        let references;
        let automatic;
        let species;
        let gene_id_column;
        let gene_id_type;

        if ("references" in parameters) {
            references = parameters.references;
            automatic = parameters.automatic;
            species = parameters.species;
            gene_id_column = parameters.gene_id_column;
            gene_id_type = parameters.gene_id_type;
        } else {
            references = null;
            automatic = true;
            let def = CellLabellingState.defaults();
            species = def.species;
            gene_id_column = def.gene_id_column;
            gene_id_type = def.gene_id_type;
        }

        this.changed = false;

        if (this.valid()) {
            this.changed = await build_reference(
                this.#cache, 
                references, 
                automatic, 
                species, 
                gene_id_column, 
                gene_id_type, 
                this.#parameters, 
                () => this.#inputs.fetchFeatureAnnotations()["RNA"],
                () => this.#inputs.guessRnaFeatureTypes()
            );
        }

        transplant_parameters(
            references, 
            automatic, 
            species, 
            gene_id_column, 
            gene_id_type, 
            this.#parameters
        );
    }

    /**
     * @param {ScranMatrix|ScoreMarkersResults} x - A matrix of (normalized or unnormalized) expression values, with genes in rows and cells/clusters in columns.
     * Alternatively, an object containing marker results, e.g., as computed by {@linkcode MarkerDetectionState}. 
     * 
     * @return {object} Object containing:
     *
     * - `per_reference`: an object where each key is the name of a reference dataset and its value is an array.
     *   This inner array is of length equal to the number of columns of `x` (if matrix) or groups in `x` (if marker results).
     *   Each entry is an object containing `best`, the name of the best label assigned to a column/group in this reference;
     *   and `all`, an object where each key is a label in this reference dataset and its value is the score for assigning that label to this column/group.
     * - (optional) `integrated`: an array of length equal to the number of columns/groups.
     *   Each entry is an object containing `best`, the name of the best reference for this column/group;
     *   and `all`, an object where each key is the name of a reference dataset and its value is the score for this column/group.
     *   This property is only reported if multiple references are used.
     */
    computeLabels(x) {
        return assign_labels(x, this.#cache);
    }
}

/**************************
 ******** Loading *********
 **************************/

export function unserialize(handle, inputs, markers) {
    let parameters =  {
        mouse_references: [],
        human_references: []
    };
    let cache = { results: {} };

    // Protect against old analysis states that don't have cell_labelling.
    if ("cell_labelling" in handle.children) {
        let ghandle = handle.open("cell_labelling");
        
        {
            let phandle = ghandle.open("parameters");
            parameters.mouse_references = phandle.open("mouse_references", { load: true }).values;
            parameters.human_references = phandle.open("human_references", { load: true }).values;
        }

        {
            let rhandle = ghandle.open("results");

            if ("per_reference" in rhandle.children) {
                let perhandle = rhandle.open("per_reference");
                for (const key of Object.keys(perhandle.children)) {
                    cache.results[key] = perhandle.open(key, { load: true }).values;
                }
                if ("integrated" in rhandle.children) {
                    cache.integrated_results = rhandle.open("integrated", { load: true }).values;
                }
            }
        }
    }

    return new CellLabellingState(inputs, markers, parameters, cache);
}
