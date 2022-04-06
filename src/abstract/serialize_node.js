import * as sutils from "../utils/serialize.js";
import * as fs from "fs";
import * as Path from "path";

/**
 * This contains a function to create and load a kana file with Node.
 */
export async function createKanaFileInternal(path, statePath, inputFiles) {
    let stream = fs.createWriteStream(path, { flags: 'a' });

    let embedded = (inputFiles === null):
    let preamble = sutils.createPreamble(embedded, fs.statSync(statePath).size);
    stream.write(Buffer.from(preamble));

    let stateStream = fs.createReadStream(statePath);
    stateStream.pipe(stream);

    if (embedded) {
        for (const ipath of inputFiles) {
            let istream = fs.createReadStream(ipath);
            istream.pipe(stream);
        }
    }

    return new Promise((resolve, reject) => {
        stream.on("finish", () => resolve(true));
        stream.on("error", (e) => reject(e));
        stream.end();
    });
}

export async function parseKanaFileInternal(input, statePath, { inputDir: null } = {}) {
    if (inputDir === null) {
        inputDir = fs.mkdtempSync("kana-");
    }

    let fd = fs.openSync(path);
    let prebuffer = new Uint8Array(24);
    fs.readSync(fd, prebuffer, 0, prebuffer.length, null);

    // TODO: if I had more Node skills, we could stream it into the file rather 
    // than loading it into memory before dumping it out again. Oh well.
    let parsed = sutils.parsePreamble(prebuffer.buffer);
    let state_len = parsed.state;
    let statebuffer = new Uint8Array(state_len);
    fs.readSync(fd, statebuffer, 0, state_len, null);
    fs.writeFileSync(statePath, statebuffer);

    fs.closeSync();

    if (parsed.embedded) {
        // Safest to just reopen the damn file and write it to the location.
        return (offset, size) => {
            let opath = Path.join(inputDir, offset);
            if (fs.existsFileSync(opath)) {
                let contents = new Uint8Array(size);
                let handle = fs.openSync(path);
                fs.readSync(handle, contents, 0, size, offset);
                fs.writeFileSync(opath, contents);
            }
            return opath;
        };
    } else {
        return null;
    }
}
