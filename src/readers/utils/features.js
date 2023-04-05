import * as scran from "scran.js";
import * as bioc from "bioconductor";

function create_solo_default_object(value, modality) {
    let output = {};
    output[modality] = value;
    return output;
}

export function reportFeatures(rawFeatures, typeField) {
    if (rawFeatures.hasColumn(typeField)) {
        let by_type = bioc.presplitFactor(rawFeatures.column(typeField));
        let copy = bioc.CLONE(rawFeatures, { deepCopy: false }); // SPLIT will make a copy anyway.
        copy.$removeColumn(typeField);
        return bioc.SPLIT(copy, by_type);
    } else {
        return create_solo_default_object(rawFeatures, "");
    }
}

function is_subset_noop(indices, full_length) {
    if (indices.length != full_length) {
        return false;
    }
    for (var i = 0; i < full_length; i++) {
        if (i !== indices[i]) {
            return false;
        }
    }
    return true;
}

function renameByModality(input, featureTypeMapping) {
    let output = {};
    for (const [k, v] of Object.entries(featureTypeMapping)) {
        if (v !== null && v in input) {
            output[k] = input[v];
        }
    }
    return output;
}

function splitByModality(features, typeField, featureTypeMapping) {
    let by_type = bioc.presplitFactor(features.column(typeField));
    if (featureTypeMapping === null) {
        return by_type;
    }
    return renameByModality(by_type, featureTypeMapping);
}

export function extractSplitPrimaryIds(features, typeField, featureTypeMapping, featureTypeDefault, primary) {
    if (typeField !== null && features.hasColumn(typeField)) {
        let by_type = splitByModality(features, typeField, featureTypeMapping);
        for (const [k, v] of Object.entries(by_type)) {
            let col = extractPrimaryIdColumn(k, features, primary);
            by_type[k] = bioc.SLICE(col, v);
        }
        return by_type;
    }

    let output = {};
    output[featureTypeDefault] = extractPrimaryIdColumn(featureTypeDefault, features, primary);
    return output;
}

export function splitScranMatrixAndFeatures(loaded, rawFeatures, typeField, featureTypeMapping, featureTypeDefault) {
    let output = { matrix: new scran.MultiMatrix };

    try {
        let out_mat = loaded.matrix;
        let out_ids = loaded.row_ids;
        output.matrix.add("", out_mat);

        let current_features;
        if (out_ids !== null) {
            current_features = bioc.SLICE(rawFeatures, out_ids);
        } else {
            current_features = bioc.CLONE(rawFeatures, { deepCopy: false }); // because we're deleting a column.
            out_ids = new Int32Array(out_mat.numberOfRows());
            out_ids.forEach((x, i) => { out_ids[i] = i });
        }

        if (typeField !== null && current_features.hasColumn(typeField)) {
            let by_type = splitByModality(current_features, typeField, featureTypeMapping);
            let type_keys = Object.keys(by_type);
            let skip_subset = is_subset_noop(type_keys[0], out_mat.numberOfRows());

            if (type_keys.length > 1 || !skip_subset) {
                let replacement = new scran.MultiMatrix({ store: scran.splitRows(out_mat, by_type) });
                scran.free(output.matrix);
                output.matrix = replacement;
            } else {
                output.matrix.rename("", type_keys[0]);
            }

            delete current_features[typeField];
            output.features = bioc.SPLIT(current_features, by_type);
            output.row_ids = bioc.SPLIT(out_ids, by_type);

        } else {
            output.matrix.rename("", featureTypeDefault);
            output.row_ids = create_solo_default_object(out_ids, featureTypeDefault);
            output.features = create_solo_default_object(current_features, featureTypeDefault);
        }
    } catch (e) {
        scran.free(output.matrix);
        throw e;
    }

    return output;
}

function extractPrimaryIdColumn(modality, modality_features, primary) {
    if (!(modality in primary)) {
        throw new Error("modality '" + modality + "' has no primary key identifier");  
    }
    let id = primary[modality];

    if ((typeof id == "string" && modality_features.hasColumn(id)) || (typeof id == "number" && id < modality_features.numberOfColumns())) {
        return modality_features.column(id);
    } 

    return modality_features.rowNames();
}

export function extractPrimaryIds(features, primary) {
    let output = {};
    for (const [k, v] of Object.entries(features)) {
        output[k] = extractPrimaryIdColumn(k, v, primary);
    }
    return output;
}

export function extractRemappedPrimaryIds(features, featureTypeMapping, primary) {
    let renamed = renameByModality(features, featureTypeMapping);
    return extractPrimaryIds(renamed, primary);
}
