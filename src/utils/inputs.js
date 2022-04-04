import * as scran from "scran.js";
import * as TENxReader from "./../readers/10x.js";
import * as H5ADReader from "./../readers/h5ad.js";
import * as MtxReader from "./../readers/mtx.js";

export function intersectGenes(genes) {
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

    // Now actually performing an intersection.
    let intersection = [];
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
            let name = names[i];
            let best = best_type_cols[i];
            let curgenes = genes[name][best];

            if (i == 0) {
                intersection = curgenes;
            } else {
                let dset = new Set(curgenes);
                intersection = intersection.filter(n => dset.has(n));
            }
            best_fields[name] = best;
        }
    }

    return { 
        "intersection": intersection, 
        "best_type": best_features,
        "best_fields": best_fields
    };
}

export function chooseNamespace(format) {
    let namespace;
    switch (format) {
        case "MatrixMarket":
            namespace = MtxReader;
            break;
        case "10X":
            namespace = TENxReader;
            break;
        case "H5AD":
            namespace = H5ADReader;
            break;
        default:
            throw "unknown matrix format '" + format + "'";
    }
    return namespace;
}
