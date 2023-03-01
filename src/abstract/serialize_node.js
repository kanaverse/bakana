import * as sutils from "./utils/serialize.js";
import * as fs from "fs";
import * as Path from "path";
import * as v0 from "./utils/legacy_v0.js";
import * as pako from "pako";
import * as os from "os";
export { FORMAT_VERSION } from "./utils/serialize.js";

export async function parseKanaFileInternal(input, statePath, { stageDir = null } = {}) {
    if (input instanceof Uint8Array) {
        return sutils.parseKanaFileFromBuffer(input, statePath);
    }

    if (stageDir === null) {
        stageDir = fs.mkdtempSync(Path.join(os.tmpdir(), "kana-"));
    }

    let fd = fs.openSync(input);
    let prebuffer = new Uint8Array(24);
    fs.readSync(fd, prebuffer, 0, prebuffer.length, null);
    fs.closeSync(fd);

    let parsed = sutils.parsePreamble(prebuffer);
    let state_len = parsed.state;
    let delta = parsed.offset + state_len;

    if (parsed.version < 1000000) {
        let fd = fs.openSync(input);
        let statebuffer = new Uint8Array(state_len);
        fs.readSync(fd, statebuffer, 0, state_len, parsed.offset);
        fs.closeSync(fd);

        var contents = pako.ungzip(statebuffer, { "to": "string" });
        let state = JSON.parse(contents);
        v0.convertFromVersion0(state, statePath);

    } else {
        // Piping it into the output file.
        let istream = fs.createReadStream(input, { start: parsed.offset, end: delta - 1 });
        let ostream = fs.createWriteStream(statePath);
        let piped = istream.pipe(ostream);

        await new Promise((resolve, reject) => {
            piped.on("finish", () => resolve(null));
            piped.on("error", e => reject(e));
        });
    }

    let output = { 
        version: parsed.version,
        embedded: parsed.embedded,
        load: null
    };

    if (parsed.embedded) {
        // Safest to just reopen the damn file and write it to the location.
        // However, if we already rewrote it, then we skip the process.
        output.load = async (offset, size) => {
            let opath = Path.join(stageDir, String(offset));
            if (!fs.existsSync(opath)) {
                let istream = fs.createReadStream(input, { start: delta + offset, end: delta + offset + size - 1});
                let ostream = fs.createWriteStream(opath);
                let piped = istream.pipe(ostream);
                
                await new Promise((resolve, reject) => {
                    piped.on("finish", () => resolve(null));
                    piped.on("error", e => reject(e));
                });
            }
            return opath;
        };
    }

    return output;
}
