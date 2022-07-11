import * as scran from "scran.js";

export function mimicGetter(value, copy) {
    // Inheritance seems to be namespaced by module,
    // so we can't use instanceof.
    if ("className" in value.constructor && value.constructor.className.endsWith("WasmArray")) { 
        if (copy == "view" || copy == "hdf5") {
            return value.view();
        } else if (copy) {
            return value.slice();
        } else {
            return value.array();
        }
    } else {
        if (copy === true) {
            return value.slice();
        } else if (copy == "view") {
            // If the caller actually wanted a WasmArray, they would
            // have generated a WasmArray during the unserialization.
            throw new Error("'copy: \"view\"' not supported for mimics");
        } else {
            // Includes copy = "hdf5", where a TypedArray or WasmArray can be used.
            return value;
        }
    }
}

export function freeCache(object) {
    // Just an alias for simplicity.
    scran.safeFree(object);
    return;
}

export function changedParameters(x, y) {
    return JSON.stringify(x) != JSON.stringify(y);
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

export function isObject(object) {
    return typeof object === 'object' && Array.isArray(object) === false && ArrayBuffer.isView(object) === false;
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
