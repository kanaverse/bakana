import * as fs from "fs";
import * as path from "path";
import * as crypto from "crypto";
import { md5 } from 'hash-wasm';

export async function attachMd5sums(files) {
    for (const x of files) {
        if (!("contents" in x)) {
            continue;
        }

        if (x.contents instanceof Uint8Array) {
            // Not technically necessary as we could continue using crypto
            // here. However, as the browser requires hash-wasm, we might as
            // well use it here to double-check that it works properly.
            x.metadata.md5sum = await md5(x.contents);
        } else {
            x.metadata.md5sum = await (new Promise((resolve, reject) => {
                let hash = crypto.createHash('md5');
                const rs = fs.createReadStream(x.contents)
                rs.on('error', reject)
                rs.on('data', chunk => hash.update(chunk))
                rs.on('end', () => resolve(hash.digest('hex')))
            }));
        }
    }
}

export async function realizeDirectory(files, directory, path_) {
    let base_path = path.join(directory, path_)
    if (fs.existsSync(base_path)) {
        fs.rmSync(base_path, { recursive: true, force: true });
    }

    for (const x of files) {
        let full_path = path.join(directory, x.metadata.path);
        let dir_path = path.dirname(full_path);
        fs.mkdirSync(dir_path, { recursive: true });

        let suffix;
        if ("contents" in x) {
            suffix = ".json";

            if (x.contents instanceof Uint8Array) {
                fs.writeFileSync(full_path, x.contents);
            } else {
                fs.renameSync(x.contents, full_path);
            }
            x.contents = full_path;

        } else if (x.metadata["$schema"].startsWith("redirection/")) {
            suffix = ".json";
        } else {
            suffix = "";
        }

        // Add a trailing newline to avoid no-newline warnings. 
        fs.writeFileSync(full_path + suffix, JSON.stringify(x.metadata, null, 2) + "\n");
    }

    return;
}
