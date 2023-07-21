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
        let copy = rawFeatures.removeColumn(typeField);
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

function findUnnamedDefault(featureTypeMapping, featureTypeDefault) {
    let found = null;
    let multiple = false;
    for (const [k, v] of Object.entries(featureTypeMapping)) {
        if (v !== null) {
            if (found !== null) {
                multiple = true;
            }
            found = k;
        }
    }

    if (found === null || multiple) {
        return featureTypeDefault;
    } else {
        return found;
    }
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

    // Seeing if any featureTypeMapping is set to the unnamed string.
    let new_default = findUnnamedDefault(featureTypeMapping, featureTypeDefault);
    let output = {};
    output[new_default] = extractPrimaryIdColumn(new_default, features, primary);
    return output;
}

export function splitScranMatrixAndFeatures(loaded, rawFeatures, typeField, featureTypeMapping, featureTypeDefault) {
    let output = { matrix: new scran.MultiMatrix };

    try {
        output.matrix.add("", loaded);

        let current_features = bioc.CLONE(rawFeatures, { deepCopy: false }); // because we're deleting a column.
        if (typeField !== null && current_features.hasColumn(typeField)) {
            let by_type = splitByModality(current_features, typeField, featureTypeMapping);
            let type_keys = Object.keys(by_type);
            let skip_subset = is_subset_noop(type_keys[0], out_mat.numberOfRows());

            if (type_keys.length > 1 || !skip_subset) {
                let replacement = new scran.MultiMatrix({ store: scran.splitRows(loaded, by_type) });
                scran.free(output.matrix);
                output.matrix = replacement;
            } else {
                output.matrix.rename("", type_keys[0]);
            }

            delete current_features[typeField];
            output.features = bioc.SPLIT(current_features, by_type);

        } else {
            output.matrix.rename("", featureTypeDefault);
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
