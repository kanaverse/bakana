import * as sutils from "./utils/serialize.js";
import * as scran from "scran.js";
export { FORMAT_VERSION } from "./utils/serialize.js";

export function parseKanaFileInternal(input, statePath) {
    return sutils.parseKanaFileFromBuffer(input, statePath);
}
