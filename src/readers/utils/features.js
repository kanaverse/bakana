import * as scran from "scran.js";

export function createMockIds(n) {
    let ids = [];
    for (var i = 0; i < n; i++) {
        ids.push("Feature " + String(i));
    }
    return ids;
}

const default_modality = "";

function create_solo_default_object(value) {
    let output = {};
    output[default_modality] = value;
    return output;
}

export function reportFeatures(rawFeatures, typeField) {
    if (typeField in rawFeatures) {
        let by_type = scran.splitByFactor(rawFeatures[typeField]);
        let copy = { ...rawFeatures };
        delete copy[typeField];
        return scran.splitArrayCollection(copy, by_type);
    } else {
        // Cloning this instance to avoid complications if the caller modifies the return value.
        return create_solo_default_object(scran.cloneArrayCollection(rawFeatures));
    }
}

export function splitScranMatrixAndFeatures(loaded, rawFeatures, typeField) {
    let output = { matrix: new scran.MultiMatrix };

    try {
        let out_mat = loaded.matrix;
        let out_ids = loaded.row_ids;
        output.matrix.add(default_modality, out_mat);

        let current_features;
        if (out_ids !== null) {
            current_features = scran.subsetArrayCollection(rawFeatures, out_ids);
        } else {
            current_features = scran.cloneArrayCollection(rawFeatures);
            out_ids = new Int32Array(out_mat.numberOfRows());
            out_ids.forEach((x, i) => { out_ids[i] = i });
        }

        if (typeField !== null && typeField in current_features) {
            let by_type = scran.splitByFactor(current_features[typeField]);
            let type_keys = Object.keys(by_type);

            if (type_keys.length > 1) {
                let replacement = new scran.MultiMatrix({ store: scran.splitRows(out_mat, by_type) });
                scran.free(output.matrix);
                output.matrix = replacement;
            } else {
                output.matrix.rename(default_modality, type_keys[0]);
            }

            delete current_features[typeField];
            output.features = scran.splitArrayCollection(current_features, by_type);
            output.row_ids = scran.splitArray(out_ids, by_type);

        } else {
            output.row_ids = create_solo_default_object(out_ids);
            output.features = create_solo_default_object(current_features);
        }
    } catch (e) {
        scran.free(output.matrix);
        throw e;
    }

    return output;
}
