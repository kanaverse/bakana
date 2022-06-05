import * as utils from "./utils/general.js";
import * as iutils from "./utils/inputs.js";
import * as f from "./abstract/file.js";

/**
 * Perform preflight validation of the annotations in the input matrices.
 *
 * This is usually done in regards to the consistency of gene annotations for multiple matrices.
 * If multiple matrices are present, an error is raised if any matrix does not contain gene annotations;
 * or there is no common species and feature type for annotations across all matrices;
 * or the intersection of genes is empty.
 *
 * @param {object} matrices - An object where each property is itself an object representing a single input matrix.
 * See the argument of the same name in {@linkcode runAnalysis} for more details.
 *
 * @return A promise that resolves to an object containing preflight check information:
 * - `annotations`: an object where each property corresponds to an input matrix in `matrices`.
 *   Each element is an array containing the names of the cell annotation fields in the corresponding matrix;
 *   or `null`, if that matrix does not contain any annotations.
 * - `features`: an object where each key is the name of a modality (e.g., `"RNA"`, `"ADT"`) and the property is another object.
 *   The latter contains:
 *   - `common`: an integer containing the number of common genes across all matrices.
 *      This may be `null` if `matrices` has only a single entry that lacks gene annotation.
 *   - `fields`: an object where each property corresponds to an input matrix in `matrices`.
 *      Each property is a string containing the name of the gene annotation field that was chosen for the intersection in the corresponding matrix.
 *      This may be `null` if `matrices` has only a single entry that lacks gene annotation.
 */
export async function validateAnnotations(matrices) {
    let mkeys = Object.keys(matrices);
    let multi = mkeys.length > 1;

    let promises = [];
    for (const key of mkeys) {
        let val = matrices[key];
        let namespace = iutils.chooseReader(val.format);
        promises.push(namespace.preflight(val));
    }
    let collected = await Promise.all(promises);

    let genes = {};
    let annotations = {};
    for (const [i, key] of Object.entries(mkeys)) {
        let stuff = collected[i];
        if (stuff.genes !== null) {
            genes[key] = stuff.genes;
        } else if (multi) {
            throw new Error("cannot find gene annotations for matrix '" + key + "'");
        }
        annotations[key] = stuff.annotations;
    }

    // Find the intersection of all modalities.
    let modalities = null;
    for (const [k, v] of Object.entries(genes)) {
        if (modalities == null) {
            modalities = new Set(Object.keys(v));
        } else {
            let alt = Object.keys(v).filter(x => modalities.has(x));
            modalities = new Set(alt);
        }
    }

    // For each modality, intersect the features.
    let feature_info = {};
    
    if (modalities !== null) {
        for (const m of modalities) {
            feature_info[m] = {};

            let genes2 = {};
            for (const [k, v] of Object.entries(genes)) {
                genes2[k] = v[m];
            }
            let results = iutils.commonFeatureTypes(genes2);
            if (results.best_type === null) {
                throw new Error("cannot find common feature types across all matrices");
            }
            feature_info[m].fields = results.best_fields;

            if (multi) {
                let intersection = null;
                for (const [k, v] of Object.entries(results.best_fields)) {
                    let curgenes = genes2[k][v];
                    if (intersection === null) {
                        intersection = curgenes;
                    } else {
                        let dset = new Set(curgenes);
                        intersection = intersection.filter(n => dset.has(n));
                    }
                }
                feature_info[m].common = intersection.length;
            } else {
                feature_info[m].common = Object.values(Object.values(genes2)[0])[0].length;
            }
        }
    } else {
        feature_info["RNA"] = { common: null, fields: null };
    }

    return { 
        annotations: annotations,
        features: feature_info
    };
}
