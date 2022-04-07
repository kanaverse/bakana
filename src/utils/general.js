import * as scran from "scran.js";
import * as wa from "wasmarrays.js";

export function mimicGetter(value, copy) {
    if (value instanceof wa.WasmArray) {
        if (copy == "view") {
            return value.view();
        } else if (copy) {
            return value.slice();
        } else {
            return value.array();
        }
    } else {
        if (copy === true) {
            return value.slice();
        } else {
            // Includes copy = "view"; we just provide a no-copy and assume
            // that, if the caller actually wanted a WasmArray, they would
            // have generated a WasmArray during the unserialization.
            return value;
        }
    }
}

export function freeCache(object) {
    if (object !== undefined && object !== null) {
        object.free();
    }
    return;
}

export function changedParameters(x, y) {
    return JSON.stringify(x) != JSON.stringify(y);
}

export function allocateCachedArray(size, type, cache, name = "buffer") {
    var reallocate = true;
    if (name in cache) {
        var candidate = cache[name];
        if (candidate.size != size || candidate.constructor.className != type) {
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
    return typeof object === 'object' && Array.isArray(object) === false;
}