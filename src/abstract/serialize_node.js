import * as sutils from "../utils/serialize.js";
import * as fs from "fs";
import * as Path from "path";

/**
 * This contains a function to create and load a kana file with Node.
 */
export async function createKanaFileInternal(statePath, inputFiles, { outputPath = null } = {}) {
    if (outputPath === null) {
        let dir = fs.mkdtempSync("kana-");
        outputPath = Path.join(dir, "analysis.kana");
    }

    let stream = fs.createWriteStream(outputPath, { flags: 'w' });

    let embedded = (inputFiles !== null);
    let preamble = sutils.createPreamble(embedded, fs.statSync(statePath).size);
    stream.write(Buffer.from(preamble));

    let stateStream = fs.createReadStream(statePath);
    let piped = stateStream.pipe(stream, { end: false });

    await new Promise((resolve, reject) => {
        piped.on("unpipe", () => resolve(true));
        piped.on("error", e => reject(e));
    });

    if (embedded) {
        for (const ipath of inputFiles) {
            let istream = fs.createReadStream(ipath);
            let piped = istream.pipe(stream, { end: false });

            await new Promise((resolve, reject) => {
                piped.on("unpipe", () => resolve(true));
                piped.on("error", e => reject(e));
            });
        }
    }

    stream.end();

    return new Promise((resolve, reject) => {
        stream.on("finish", () => resolve(outputPath));
        stream.on("error", (e) => reject(e));
    });
}

export function parseKanaFileInternal(input, statePath, { stageDir = null } = {}) {
    if (stageDir === null) {
        stageDir = fs.mkdtempSync("kana-");
    }

    let fd = fs.openSync(input);
    let prebuffer = new Uint8Array(24);
    fs.readSync(fd, prebuffer, 0, prebuffer.length, null);

    // TODO: if I had more Node skills, we could stream it into the file rather 
    // than loading it into memory before dumping it out again. Oh well.
    let parsed = sutils.parsePreamble(prebuffer.buffer);
    let state_len = parsed.state;
    let statebuffer = new Uint8Array(state_len);
    fs.readSync(fd, statebuffer, 0, state_len, null);
    fs.writeFileSync(statePath, statebuffer);

    fs.closeSync(fd);

    let delta = parsed.offset + state_len;

    if (parsed.embedded) {
        // Safest to just reopen the damn file and write it to the location.
        // However, if we already rewrote it, then we skip the process.
        return (offset, size) => {
            let opath = Path.join(stageDir, String(offset));
            if (!fs.existsSync(opath)) {
                let contents = new Uint8Array(size);
                let handle = fs.openSync(input, 'r');
                fs.readSync(handle, contents, 0, size, delta + offset);
                fs.writeFileSync(opath, contents);
            }
            return opath;
        };
    } else {
        return null;
    }
}
