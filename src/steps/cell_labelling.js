import * as scran from "scran.js";
import * as utils from "./utils/general.js";
import * as rutils from "../readers/index.js";
import * as inputs_module from "./inputs.js";
import * as markers_module from "./marker_detection.js";

var downloadFun = async (url) => {
    let resp = await fetch(url);
    if (!resp.ok) {
        throw new Error("failed to fetch content at " + url + "(" + resp.status + ")");
    }
    return new Uint8Array(await resp.arrayBuffer());
};

const baseUrl = "https://github.com/LTLA/singlepp-references/releases/download/v2.0.0";

// Loaded references are constant, independent of the dataset;
// so we can keep these as globals for re-use across States.
const all_loaded = {};

export const step_name = "cell_labelling";

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
    #markers;
    #parameters;
    #cache;

    constructor(inputs, markers, parameters = null, cache = null) {
        if (!(inputs instanceof inputs_module.InputsState)) {
            throw new Error("'inputs' should be a State object from './inputs.js'");
        }
        this.#inputs = inputs;

        if (!(markers instanceof markers_module.MarkerDetectionState)) {
            throw new Error("'markers' should be a State object from './marker_detection.js'");
        }
        this.#markers = markers;

        this.#parameters = (parameters === null ? {} : parameters);
        this.#cache = (cache === null ? {} : cache);
        this.changed = false;
    }

    #flush_prepared() {
        if ("prepared" in this.#cache) {
            for (const v of Object.values(this.#cache.prepared)) {
                v.built.raw.free();
            }
            delete this.#cache.prepared;
        }
    }

    free() {
        utils.freeCache(this.#cache.buffer);
        this.#flush_prepared();
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
        return { ...this.#parameters };
    }

    /**
     * @return {object} An object containing:
     *
     * - `per_reference`: an object where keys are the reference names and the values are arrays of strings.
     *   Each array is of length equal to the number of clusters and contains the cell type classification for each cluster.
     * - `integrated`: an array of length equal to the number of clusters.
     *   Each element is a string specifying the name of the reference with the best label for each cluster.
     *   Only available if multiple references are requested.
     *
     * This is available after running {@linkcode CellLabellingState#compute compute}.
     */
    fetchResults() {
        // No real need to clone these, they're string arrays
        // so they can't be transferred anyway.
        let perref = {};
        for (const [key, val] of Object.entries(this.#cache.results)) {
            perref[key] = val;
        }

        let output = { "per_reference": perref };
        if ("integrated_results" in this.#cache) {
            output.integrated = this.#cache.integrated_results;
        }

        return output;
    }

    /**
     * @return {object} Object where each key is the name of a reference and each value is the number of shared features between the test and reference daatasets.
     */
    fetchNumberOfSharedFeatures() {
        let output = {};
        for (const key of Object.keys(this.#cache.results)) {
            output[key] = this.#cache.prepared[key].built.raw.sharedFeatures();
        }
        return output;
    }

    /****************************
     ******** Defaults **********
     ****************************/

    static defaults() {
        return {
            references: [],
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

    /**
     * Available references for each species.
     * Each key is a taxonomy ID and each value is an array of strings containing the names of references for that species.
     * @type {object}
     */
    static availableReferences = {
        "9606": [ "BlueprintEncode", "DatabaseImmuneCellExpression", "HumanPrimaryCellAtlas", "MonacoImmune", "NovershternHematopoietic" ],
        "10090": [ "ImmGen", "MouseRNAseq" ]
    };

    /***************************
     ******** Compute **********
     ***************************/

    async #load_reference(name) {
        if (name in all_loaded) {
            return;
        }

        const suffixes = [ 
            "genes.csv.gz",
            "labels_fine.csv.gz",
            "label_names_fine.csv.gz",
            "markers_fine.gmt.gz",
            "matrix.csv.gz"
        ];

        let contents = await Promise.all(
            suffixes.map(
                async suffix => {
                    let full = name + "_" + suffix;
                    let b = await downloadFun(baseUrl + "/" + full);
                    return new rutils.SimpleFile(b, { name: full })
                }
            )
        );

        let loaded;
        try {
            loaded = scran.loadLabelledReferenceFromBuffers(
                contents[4].buffer(), // rank matrix
                contents[3].buffer(), // markers
                contents[1].buffer()  // label per sample
            );

            let gene_lines = await rutils.readLines2(contents[0].content(), { compression: "gz" }); // gene names
            let ensembl = [];
            let symbol = [];
            let entrez = [];
            let empty2null = x => (x == "" ? null : x);

            gene_lines.forEach(x => {
                let fields = x.split(",");
                ensembl.push(empty2null(fields[0]));
                symbol.push(empty2null(fields[1]));
                entrez.push(empty2null(fields[2]));
            });

            let labels = await rutils.readLines2(contents[2].content(), { compression: "gz" }); // full label names
            all_loaded[name] = { 
                "raw": loaded, 
                "genes": {
                    "ENSEMBL": ensembl,
                    "SYMBOL": symbol,
                    "ENTREZ": entrez
                },
                "labels": labels
            };

        } catch (e) {
            utils.freeCache(loaded);
            throw e;
        }
    }

    #build_reference(name, gene_ids, gene_id_type) {
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

    /**
     * This method should not be called directly by users, but is instead invoked by {@linkcode runAnalysis}.
     *
     * @param {object} parameters - Parameter object, equivalent to the `cell_labelling` property of the `parameters` of {@linkcode runAnalysis}.
     * @param {Array} parameters.references - Array of strings specifying the names of the reference datasets, see {@linkcode CellLabellingState.availableReferences availableReferences} for more details.
     * @param {boolean} parameters.automatic - Automatically choose feature-based parameters based on the feature annotation for the RNA modality.
     * If `true`, the column of the annotation that best matches human/mouse Ensembl/symbols is identified and used to set `species`, `gene_id_column` and `gene_id_type`.
     * @param {Array} parameters.species - Array of strings specifying zero, one or more species involved in this dataset.
     * Each entry can either be the common name (e.g., `"mouse"`, `"human"`) or a taxonomy ID (e.g. 9606, 10090).
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
            references = [ ...(parameters.human_references), ...(parameters.mouse_references) ];
            automatic = true;
            let def = CellLabellingState.defaults();
            species = def.species;
            gene_id_column = def.gene_id_column;
            gene_id_type = def.gene_id_type;
        }

        this.changed = false;

        if (this.valid()) {
            // Gathering the references.
            if (
                this.#inputs.changed ||
                automatic !== this.#parameters.automatic ||
                utils.changedParameters(references, this.#parameters.references) ||
                (
                    !automatic &&
                    (
                        species !== this.#parameters.species ||
                        gene_id_column !== this.#parameters.gene_id_column ||
                        gene_id_type !== this.#parameters.gene_id_type
                    )
                )
            ) {
                let species2 = species;
                let gene_id_column2 = gene_id_column;
                let gene_id_type2 = gene_id_type;

                if (automatic) {
                    let guesses = this.#inputs.guessRnaFeatureTypes();
                    let auto = CellLabellingState.configureFeatureParameters(guesses);
                    species2 = auto.species;
                    gene_id_column2 = auto.gene_id_column;
                    gene_id_type2 = auto.gene_id_type;
                }

                let allowable = new Set;
                for (const s0 of species2) {
                    let s = utils.toTaxonomy(s0);
                    if (s in CellLabellingState.availableReferences) {
                        CellLabellingState.availableReferences[s].forEach(x => { allowable.add(x); });
                    }
                }

                // Building each individual reference.
                let feats = this.#inputs.fetchFeatureAnnotations()["RNA"];
                let gene_ids = (gene_id_column2 == null ? feats.rowNames() : feats.column(gene_id_column2));
                this.#cache.gene_ids = gene_ids;

                let valid = {};
                if (gene_ids !== null) {
                    for (const ref of references) {
                        if (allowable.has(ref)) {
                            await this.#load_reference(ref);
                            valid[ref] = this.#build_reference(ref, gene_ids, gene_id_type2);
                        }
                    }
                }

                this.#flush_prepared();
                this.#cache.prepared = valid;

                // Building an integrated reference, if necessary.
                let used_refs = Object.keys(valid);
                if (used_refs.length > 1) {
                    let arr = Object.values(valid);
                    let loaded = arr.map(x => x.loaded.raw);
                    let feats = arr.map(x => x.built.features);
                    let built = arr.map(x => x.built.raw);

                    utils.freeCache(this.#cache.integrated);
                    this.#cache.integrated = scran.integrateLabelledReferences(gene_ids, loaded, feats, built);
                } else {
                    utils.freeCache(this.#cache.integrated);
                    delete this.#cache.integrated;
                }
                this.#cache.used_refs = used_refs;

                this.changed = true;
            }

            let marker_results = this.#markers.fetchResults()["RNA"];
            let ngroups = marker_results.numberOfGroups();
            let ngenes = (this.#cache.gene_ids !== null ? this.#cache.gene_ids.length : null);
            let cluster_means = this.#cache.buffer;

            if (this.#markers.changed) {
                if (ngenes !== null) {
                    // Creating a column-major array of mean vectors for each cluster.
                    cluster_means = utils.allocateCachedArray(ngroups * ngenes, "Float64Array", this.#cache);
                    for (var g = 0; g < ngroups; g++) {
                        let means = marker_results.means(g, { copy: false }); // Warning: direct view in wasm space - be careful.
                        let cluster_array = cluster_means.array();
                        cluster_array.set(means, g * ngenes);
                    }
                }
                this.changed = true;
            }

            if (this.changed) {
                // Running classifications on the cluster means. This is a
                // no-op if gene_ids = null as 'valid' should be empty.
                let valid = this.#cache.prepared;

                this.#cache.results = {};
                for (const [key, ref] of Object.entries(valid)) {
                    let output = scran.labelCells(cluster_means, ref.built.raw, { numberOfFeatures: ngenes, numberOfCells: ngroups });
                    let labels = [];
                    for (const o of output) {
                        labels.push(ref.loaded.labels[o]);
                    }
                    this.#cache.results[key] = labels;
                }

                // Performing additional integration, if necessary. 
                if ("integrated" in this.#cache) {
                    let results = [];
                    for (const key of this.#cache.used_refs) {
                        results.push(this.#cache.results[key]);
                    }

                    let out = scran.integrateCellLabels(cluster_means, results, this.#cache.integrated, { numberOfFeatures: ngenes, numberOfCells: ngroups });
                    let as_names = [];
                    out.forEach(i => {
                        as_names.push(this.#cache.used_refs[i]);
                    });
                    this.#cache.integrated_results = as_names;
                } else {
                    delete this.#cache.integrated_results;
                }
            }
        } else {
            this.#cache.results = {};
            delete this.#cache.integrated_results;
        }

        this.#parameters.references = references;
        this.#parameters.automatic = automatic;
        this.#parameters.species = species;
        this.#parameters.gene_id_column = gene_id_column;
        this.#parameters.gene_id_type = gene_id_type;

        return;
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

/**************************
 ******** Setters *********
 **************************/

/**
 * Specify a function to download references for the cell labelling step.
 *
 * @param {function} fun - Function that accepts a single string containing a URL and returns any value that can be used in the {@linkplain SimpleFile} constructor.
 * This is most typically a Uint8Array of that URL's contents, but it can also be a path to a locally cached file on Node.js.
 *
 * @return `fun` is set as the global downloader for this step. 
 * The _previous_ value of the downloader is returned.
 */
export function setCellLabellingDownload(fun) {
    let previous = downloadFun;
    downloadFun = fun;
    return previous;
}
