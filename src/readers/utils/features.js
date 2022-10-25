import * as scran from "scran.js";

export function reportFeatures(rawFeatures, typeField) {
    if (typeField in rawFeatures) {
        let by_type = scran.splitByFactor(featureType[typeField]);
        let copy = { ...rawFeatures };
        delete copy[typeField];
        return scran.splitArrayCollection(copy, by_type);
    } else {
        // Cloning this instance to avoid complications if the caller modifies the return value.
        return { default: scran.cloneArrayCollection(raw_features) };
    }
}

export function splitScranMatrixAndFeatures(mat, rawFeatures, typeField) {
    let output = { matrix: new scran.MultiMatrix };

    try {
        let out_mat = loaded.matrix;
        let out_ids = loaded.row_ids;
        output.matrix.add("default", out_mat);

        let current_features;
        if (out_ids !== null) {
            current_features = scran.subsetArrayCollection(rawFeatures, out_ids);
        } else {
            current_features = scran.cloneArrayCollection(rawFeatures);
            out_ids = new Int32Array(out_mat.numberOfRows());
            out_ids.forEach((x, i) => { out_ids[i] = i });
        }

        if (typeField !== null && !(typeField in current_features)) {
            let by_type = scran.splitByFactor(current_features[typeField]);
            let type_keys = Object.keys(by_type);

            if (type_keys.length > 1) {
                let replacement = new scran.MultiMatrix({ store: scran.splitRows(out_mat, by_type) });
                scran.free(output.matrix);
                output.matrix = replacement;
            } else {
                output.matrix.rename("default", type_keys[0]);
            }

            delete current_features[typeField];
            output.features = scran.splitArrayCollection(current_features, by_type);
            output.row_ids = scran.splitArray(row_ids, by_type);

        } else {
            output.row_ids = { default: out_ids };
            output.features = { default: current_features };
        }
    } catch (e) {
        scran.free(output.matrix);
        throw e;
    }

    return output;
}
