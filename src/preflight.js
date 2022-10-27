import * as rutils from "./readers/index.js";
import * as iutils from "./steps/inputs.js";
import * as f from "./abstract/file.js";

/**
 * Perform preflight validation of the annotations in the input datasets.
 *
 * This is usually done in regards to the consistency of gene annotations for multiple datasets.
 * If multiple datasets are present, an error is raised if any matrix does not contain gene annotations;
 * or there is no common species and feature type for annotations across all datasets;
 * or the intersection of genes is empty.
 *
 * @param {object} datasets - An object where each key is the name of a dataset and each property is a {@linkplain Dataset}.
 * See the argument of the same name in {@linkcode runAnalysis} for more details.
 *
 * @return {object} An object containing preflight check information for the `annotations` and `features`.
 *
 * `annotations` is an object where each property corresponds to an input matrix in `datasets`.
 * Each element is itself an object describing the cell annotation fields in the corresponding matrix.
 * Each key of the inner object contains the name of an annotation field, while each value is an object summarizing that annotation.
 * Each annotation summary contains:
 *   - a `type` property, indicating whether the annotation is `"categorical"` or "`continuous"`.
 *   - for categorical annotations, a `values` array containing the unique values.
 *     This may be truncated for brevity, in which case the `truncated` property is `true`.
 *   - for continuous annotations, the `min` and `max` properties containing the minimum and maximum values.
 *
 * `features` is an object where each key is the name of a modality (e.g., `"RNA"`, `"ADT"`) and the property is another object.
 * Each inner object contains:
 * - `common`: an integer containing the number of common genes across all datasets.
 *    This may be `null` if `datasets` has only a single entry that lacks gene annotation.
 * - `fields`: an object where each property corresponds to an input matrix in `datasets`.
 *    Each property is a string containing the name of the gene annotation field that was chosen for the intersection in the corresponding matrix.
 *    This may be `null` if `datasets` has only a single entry that lacks gene annotation.
 *
 * @async
 */
export async function validateAnnotations(datasets) {
    let mkeys = Object.keys(datasets);
    let multi = mkeys.length > 1;

    let promises = [];
    for (const key of mkeys) {
        promises.push(datasets[key].annotations());
    }
    let collected = await Promise.all(promises);

    let modal_mappings = iutils.guessDefaultModalities(collected);
    if (Object.keys(modal_mappings).length == 0) {
        throw new Error("failed to find any common modalities");
    }

    let annotations = {};
    for (const [i, key] of Object.entries(mkeys)) {
        annotations[key] = collected[i].cells;
    }

    // For each modality, intersect the features.
    let feature_info = {};
    
    for (const [m, chosen] of Object.entries(modal_mappings)) {
        feature_info[m] = {};

        let genes2 = [];
        for (var i = 0; i < collected.length; i++) {
            genes2.push(collected[i].features[chosen[i]]);
        }

        let results = iutils.commonFeatureTypes(genes2);
        if (results.best_type === null) {
            throw new Error("cannot find common feature types across all datasets");
        }
        feature_info[m].fields = {};
        for (var i = 0; i < collected.length; i++) {
            feature_info[m].fields[mkeys[i]] = results.best_fields[i];
        }

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

    return { 
        annotations: annotations,
        features: feature_info
    };
}
