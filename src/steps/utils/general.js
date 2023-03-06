import * as scran from "scran.js";

export function freeCache(object) {
    // Just an alias for back-compatibility.
    scran.free(object);
    return;
}

function changedParametersIllegal(x, y, xskip, yskip) {
    // Failing if this is a TypedArray or ArrayBuffer;
    // we shouldn't be seeing these things here anyway.
    if (!xskip) {
        if (x instanceof ArrayBuffer || ArrayBuffer.isView(x)) {
            throw new Error("parameters cannot contain ArrayBuffers or their views");
        }
    }
    if (!yskip) {
        if (y instanceof ArrayBuffer || ArrayBuffer.isView(y)) {
            throw new Error("parameters cannot contain ArrayBuffers or their views");
        }
    }
}

export function changedParameters(x, y) {
    if (typeof x != typeof y) {
        changedParametersIllegal(x, y, false, false);
        return true;
    } else if (typeof x != "object") {
        return x != y;
    }

    //Handling nulls (which are objects).
    let xnull = x === null;
    let ynull = y === null;
    if (xnull !== ynull) {
        changedParametersIllegal(x, y, xnull, ynull);
        return true;
    } else if (xnull) {
        return false;
    }

    // Handling arrays (which are also objects).
    let xarr = x instanceof Array;
    let yarr = y instanceof Array;
    if (xarr != yarr) {
        changedParametersIllegal(x, y, xarr, yarr);
        return true;
    } else if (xarr) {
        if (x.length != y.length) {
            return true;
        }

        for (var i = 0; i < x.length; i++) {
            if (changedParameters(x[i], y[i])) {
                return true;
            }
        }

        return false;
    }

    changedParametersIllegal(x, y, false, false);
    
    // Now actually handling objects. We don't 
    // worry about the order of the keys here.
    let xkeys = Object.keys(x);
    let ykeys = Object.keys(y);
    if (xkeys.length != ykeys.length) {
        return true;
    }

    xkeys.sort();
    ykeys.sort();
    for (var i = 0; i < xkeys.length; i++) {
        if (xkeys[i] != ykeys[i]) {
            return true;
        }
    }

    for (const k of xkeys) {
        if (changedParameters(x[k], y[k])) {
            return true;
        }
    }

    return false;
}

export function allocateCachedArray(size, type, cache, name = "buffer") {
    var reallocate = true;
    if (name in cache) {
        var candidate = cache[name];

        // Views also trigger reallocation, because it is assumed that the
        // caller of this function does not own the view, but downstream
        // uses of the array will involve writing to it.
        if (candidate.size != size || candidate.constructor.className != type || candidate.owner !== null) { 
            candidate.free();
        } else {
            reallocate = false;
        }
    }
  
    if (reallocate) {
        switch (type) {
            case "Uint8Array":
                cache[name] = scran.createUint8WasmArray(size);
                break;
            case "Int32Array":
                cache[name] = scran.createInt32WasmArray(size);
                break;
            case "Float64Array":
                cache[name] = scran.createFloat64WasmArray(size);
                break;
            default:
                // We only ever use one of the three above types in our 
                // internal data stores, so no need to go all-out here.
                throw "allocating '" + type + "' not yet supported";
        }
    }

    return cache[name];
}

export function findValidUpstreamStates(states, msg) {
    let to_use = [];
    for (const [k, v] of Object.entries(states)) {
        if (v.valid()) {
            to_use.push(k);
        }
    }
    if (to_use.length == 0) {
        throw new Error("expected at least one valid upstream " + msg + " state");
    }
    return to_use;
}

export function checkIndices(indices, max) {
    if (max !== null) {
        for (const i of indices) {
            if (i < 0 || i >= max) {
                throw new Error("subset indices are out of range");
            }
        }
    }

    for (var i = 1; i < indices.length; i++) {
        if (indices[i] <= indices[i-1]) {
            throw new Error("subset indices must be sorted and unique");
        }
    }
}

export async function defaultDownload(url) {
    let resp = await fetch(url);
    if (!resp.ok) {
        throw new Error("failed to fetch content at " + url + "(" + resp.status + ")");
    }
    return new Uint8Array(await resp.arrayBuffer());
}

export function guessFeatureTypes(genes) {
    let output = { columns: {} };

    let rn = genes.rowNames();
    if (rn !== null) {
        output.row_names = scran.guessFeatures(rn, { forceTaxonomy: true });
    }

    for (const key of genes.columnNames()) {
        let curcol = genes.column(key);
        if (curcol instanceof Array) {
            output.columns[key] = scran.guessFeatures(genes.column(key), { forceTaxonomy: true });
        }
    }

    return output;
}
