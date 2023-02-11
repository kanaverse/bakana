import * as scran from "scran.js";
import * as bioc from "bioconductor";
import * as utils from "./utils/general.js";
import * as mutils from "./utils/markers.js";
import * as rutils from "../readers/index.js";
import * as inputs_module from "./inputs.js";
import * as filter_module from "./cell_filtering.js";
import * as norm_module from "./rna_normalization.js";
import * as markers_module from "./marker_detection.js";

var downloadFun = async (url) => {
    let resp = await fetch(url);
    if (!resp.ok) {
        throw new Error("failed to fetch content at " + url + "(" + resp.status + ")");
    }
    return new Uint8Array(await resp.arrayBuffer());
};

const base = "https://github.com/LTLA/kana-feature-sets/releases/download";

export const step_name = "feature_set_enrichment";

export class FeatureSetEnrichmentState {
    #inputs;
    #filter;
    #normalized;
    #markers;
    #parameters;
    #cache;

    constructor(inputs, filter, normalized, markers, parameters = null, cache = null) {
        if (!(inputs instanceof inputs_module.InputsState)) {
            throw new Error("'inputs' should be a State object from './inputs.js'");
        }
        this.#inputs = inputs;

        if (!(filter instanceof filter_module.CellFilteringState)) {
            throw new Error("'filter' should be a CellFilteringState object");
        }
        this.#filter = filter;

        if (!(normalized instanceof norm_module.RnaNormalizationState)) {
            throw new Error("'normalized' should be a RnaNormalizationState object from './rna_normalization.js'");
        }
        this.#normalized = normalized;

        if (!(markers instanceof markers_module.MarkerDetectionState)) {
            throw new Error("'markers' should be a State object from './marker_detection.js'");
        }
        this.#markers = markers;

        this.#parameters = (parameters === null ? {} : parameters);
        this.#cache = (cache === null ? {} : cache);

        for (const key of [ "loaded", "mapped", "filtered", "adhoc_results" ]) {
            this.#cache[key] = {};
        }

        this.changed = false;
    }

    free() {
        this.#cache = {};
        return; // nothing extra to free here.
    }

    /***************************
     ******** Getters **********
     ***************************/

    static defaults() {
        return {
            collections: [],
            dataset_id_column: "id", 
            reference_id_column: "ENSEMBL", 
            top_markers: 100
        };
    }

    valid() {
        let mat = this.#inputs.fetchCountMatrix();
        return mat.has("RNA");
    }

    async fetchCollectionDetails() {
        let p = this.#parameters;
        let collected = await this.#prepare_collections(p.collections, p.dataset_id_column, p.reference_id_column);
        let output = {};
        for (const [k, v] of Object.entries(collected)) {
            output[k] = { names: v.names, descriptions: v.descriptions, sizes: v.sizes, universe: v.indices.length };
        }
        return output;
    }

    async fetchGroupResults(group, effect_size, summary) {
        let key = String(group) + ":" + effect_size + ":" + summary;

        if (!(key in this.#cache.adhoc_results)) {
            let res = this.#markers.fetchResults()["RNA"];
            let p = this.#parameters;
            let collections = await this.#prepare_collections(p.collections, p.dataset_id_column, p.reference_id_column);

            if (effect_size == "delta_detected") {
                effect_size = "deltaDetected";
            }
            let min_threshold = effect_size == "auc" ? 0.5 : 0;
            let use_largest = effect_size !== "min_rank"; 
            let sumidx = mutils.summaries2int[summary];

            let output = {};
            for (const [name, info] of Object.entries(collections)) {
                let stats = res[effect_size](group, { summary: sumidx, copy: false });
                let curstats = bioc.SLICE(stats, info.indices);
                let threshold = scran.computeTopThreshold(curstats, p.top_markers, { largest: use_largest });
                let in_set = [];

                if (use_largest) {
                    if (threshold < min_threshold) {
                        threshold = min_threshold;
                    }
                    curstats.forEach((x, i) => {
                        if (x >= threshold) {
                            in_set.push(i);
                        }
                    });
                } else {
                    curstats.forEach((x, i) => {
                        if (x <= threshold) {
                            in_set.push(i);
                        }
                    });
                }

                let computed = scran.testFeatureSetEnrichment(in_set, info.sets, curstats.length);
                output[name] = { 
                    counts: computed.count, 
                    pvalues: computed.pvalue,
                    num_markers: in_set.length
                };
            }

            this.#cache.adhoc_results[key] = output;
        }

        return this.#cache.adhoc_results[key];
    }

    async fetchPerCellScores(collection, set_index) {
        let p = this.#parameters;
        let collected = await this.#prepare_collections(p.collections, p.dataset_id_column, p.reference_id_column);
        let current = collected[collection];
        let set = current.sets[set_index];
        let indices = current.indices;
        // console.log(bioc.SLICE(this.#inputs.fetchFeatureAnnotations().RNA.column("id"), bioc.SLICE(indices, set).sort()));

        let mat = this.#normalized.fetchNormalizedMatrix();
        let features = utils.allocateCachedArray(mat.numberOfRows(), "Uint8Array", this.#cache, "set_buffer");
        features.fill(0);
        let farr = features.array();
        set.forEach(x => { farr[indices[x]] = 1; }); 

        return scran.scoreFeatureSet(mat, features, { block: this.#filter.fetchFilteredBlock() });
    }

    fetchParameters() {
        return { ...this.#parameters };
    }

    /****************************************
     ******** Preparing references **********
     ****************************************/

    async #load_collection(name) {
        if (!(name in this.#cache.loaded)) {
            let suffixes = [
                "features.csv.gz",
                "sets.txt.gz"
            ];

            let contents = await Promise.all(
                suffixes.map(
                    async suffix => {
                        let full = name + "_" + suffix;
                        let b = await downloadFun(base + "/" + name + "/" + full);
                        return new rutils.SimpleFile(b, { name: full })
                    }
                )
            );

            let gene_info = {};
            {
                let genes = await rutils.readTable2(contents[0].content(), { delim: "," });
                let headers = genes.shift();
                for (const x of headers) {
                    gene_info[x] = [];
                }
                for (const line of genes) {
                    for (var i = 0; i < headers.length; i++) {
                        gene_info[headers[i]].push(line[i]);
                    }
                }
            }

            let set_members = [];
            let set_name = [];
            let set_description = [];
            {
                let features = await rutils.readLines2(contents[1].content());
                for (const line of features) {
                    let values = line.split("\t");
                    set_name.push(values[0]);
                    set_description.push(values[1]);

                    let last = Number(values[2]);
                    let members = [last];
                    for (var i = 3; i < values.length; i++) {
                        let latest = Number(values[i]) + last;
                        members.push(latest);
                        last = latest;
                    }
                    set_members.push(members);
                }
            }

            this.#cache.loaded[name] = {
                features: gene_info,
                members: set_members,
                names: set_name,
                descriptions: set_description
            };
        }

        return this.#cache.loaded[name];
    }

    #remap_collection(name, data_id, ref_id) {
        if (!(name in this.#cache.mapped)) {
            let feats = this.#inputs.fetchFeatureAnnotations()["RNA"];
            if (!feats.hasColumn(data_id)) {
                throw new Error("no column '" + data_id + "' in the feature annotations");
            }
            
            let loaded = this.#cache.loaded[name];
            if (!(ref_id in loaded.features)){
                throw new Error("no column '" + ref_id + "' in the feature set annotations");
            }

            let output = scran.remapFeatureSets(feats.column(data_id), loaded.features[ref_id], loaded.members);
            let sizes = new Int32Array(output.sets.length);
            output.sets.forEach((x, i) => { sizes[i] = x.length });
            output.sizes = sizes;

            this.#cache.mapped[name] = output;
        }

        return this.#cache.mapped[name];
    }

    async #prepare_collections(collections, dataset_id_column, reference_id_column) {
        let collected = {};
        for (const x of collections) {
            let loaded = await this.#load_collection(x);
            let mapped = this.#remap_collection(x, dataset_id_column, reference_id_column);

            collected[x] = {
                names: loaded.names,
                descriptions: loaded.descriptions,
                indices: mapped.target_indices,
                sets: mapped.sets,
                sizes: mapped.sizes
            };
        }
        return collected;
    }

    /***************************
     ******** Compute **********
     ***************************/

    async compute(collections, dataset_id_column, reference_id_column, top_markers) {
        this.changed = false;
        if (this.#inputs.changed) {
            this.#cache.mapped = {};
            this.#cache.filtered = {};
            this.changed = true;
        }

        if (!this.valid()) {
            return;
        }

        let collected = await this.#prepare_collections(collections, dataset_id_column, reference_id_column);

        if (utils.changedParameters(this.#parameters.collections, collections) ||
            this.#parameters.dataset_id_column !== dataset_id_column ||
            this.#parameters.reference_id_column !== reference_id_column)
        {
            this.#parameters.collections = collections;
            this.#parameters.dataset_id_column = dataset_id_column;
            this.#parameters.reference_id_column = reference_id_column;
            this.changed = true;
        }

        if (this.changed || this.#markers.changed || top_markers !== this.#parameters.top_markers) {
            this.#cache.adhoc_results = {};
            this.#parameters.top_markers = top_markers;
            this.changed = true;
        }
    }

    /*************************
     ******** Saving *********
     *************************/

    async serialize(handle) {
        let ghandle = handle.createGroup("feature_set_enrichment");
        
        {
            let phandle = ghandle.createGroup("parameters");
            phandle.writeDataSet("feature_sets", "String", null, this.#parameters.feature_sets);
            phandle.writeDataSet("dataset_id_column", "String", null, this.#parameters.dataset_id_column);
            phandle.writeDataSet("reference_id_column", "String", null, this.#parameters.reference_id_column);
            phandle.writeDataSet("top_markers", "Int32", null, this.#parameters.top_markers);
        }

        {
            let rhandle = ghandle.createGroup("results");
        }

        return;
    }
};

/**************************
 ******** Loading *********
 **************************/

export function unserialize(handle, inputs, filter, normalized, markers) {
    let parameters = {};
    let cache = {};

    // Protect against old analysis states that don't have cell_labelling.
    if ("feature_set_enrichment" in handle.children) {
        let ghandle = handle.open("feature_set_enrichment");
        
        {
            let phandle = ghandle.open("parameters");
            parameters.feature_sets = phandle.open("feature_sets", { load: true }).values;
            for (const k of [ "dataset_id_column", "reference_id_column", "top_markers" ]) {
                parameters[k] = phandle.open(k, { load: true }).values[0];
            }
        }
    }

    return new FeatureSetEnrichmentState(inputs, filter, normalized, markers, parameters, cache);
}

/**************************
 ******** Setters *********
 **************************/

/**
 * Specify a function to download feature sets.
 *
 * @param {function} fun - Function that accepts a single string containing a URL and returns any value that can be used in the {@linkplain SimpleFile} constructor.
 * This is most typically a Uint8Array of that URL's contents, but it can also be a path to a locally cached file on Node.js.
 *
 * @return `fun` is set as the global downloader for this step.
 * The _previous_ value of the downloader is returned.
 */
export function setFeatureSetEnrichmentDownload(fun) {
    let previous = downloadFun;
    downloadFun = fun;
    return previous;
}
