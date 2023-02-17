import * as pako from "pako";
import * as wa from "wasmarrays.js";

function dump_internal(x) {
    let output;

    if (x instanceof Array) {
        output = { "type": "list", "values": [] };

        if (x.length) {
            let all_strings = true;
            let all_bools = true;
            let all_numbers = true;
            for (const e of x) {
                if (e !== null) {
                    if (typeof e !== "string") {
                        all_strings = false;
                    }
                    if (typeof e !== "boolean") {
                        all_bools = false;
                    }
                    if (typeof e !== "number") {
                        all_numbers = false;
                    }
                }
            }

            if (all_strings) {
                output.type = "string";
                output.values = x;
            } else if (all_bools) {
                output.type = "boolean";
                output.values = x;
            } else if (all_numbers) {
                output.type = "number";
                output.values = x;
            } else {
                for (const e of x) {
                    output.values.push(dump_internal(e));
                }
            }
        }

    } else if (x.constructor === Object) {
        output = { "type": "list", "values": [], "names": [] };
        for (const [k, v] of Object.entries(x)) {
            output.names.push(k);
            output.values.push(dump_internal(v));
        }

    } else if (x instanceof Int32Array) {
        output = { "type": "integer", "values": Array.from(x) }

    } else if (x instanceof wa.Int32WasmArray) {
        output = { "type": "integer", "values": Array.from(x.array()) }

    } else if (x instanceof Float64Array) {
        output = { "type": "number", "values": Array.from(x) }

    } else if (x instanceof wa.Float64WasmArray) {
        output = { "type": "number", "values": Array.from(x.array()) }

    } else if (typeof x == "number") {
        output = { "type": "number", "values": [x] };

    } else if (typeof x == "string") {
        output = { "type": "string", "values": [x] };

    } else if (typeof x == "boolean") {
        output = { "type": "boolean", "values": [x] };

    } else {
        throw new Error("don't know how to save entry of type '" + typeof x + "'");
    }

    return output;
}

export function dumpList(x, path) {
    let values = dump_internal(x);
    let encoded = JSON.stringify(values, null, 2) + "\n"; // add trailing newline to terminate file.
    let contents = pako.gzip(encoded);
    return {
        metadata: {
            "$schema": "json_simple_list/v1.json",
            "path": path + "/simple.json.gz",
            "simple_list": {
                "children": []
            },
            "json_simple_list": {
                "compression": "gzip"
            }
        },
        contents: contents
    };
}
