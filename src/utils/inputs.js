import * as scran from "scran.js";
import * as TENxReader from "./../readers/10x.js";
import * as H5ADReader from "./../readers/h5ad.js";
import * as MtxReader from "./../readers/mtx.js";

export function commonFeatureTypes(genes) {
    // Choosing the features to intersect.
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

export function chooseReader(format) {
    if (!(format in availableReaders)) {
        throw "unknown matrix format '" + format + "'";
    }
    return availableReaders[format];
}

/**
 * List of available readers.
 * Each key specifies a matrix format; the corresponding value should be an ES6 module containing methods to interpret that format.
 */
export var availableReaders = {
    "MatrixMarket": MtxReader,
    "10X": TENxReader,
    "H5AD": H5ADReader
};
