import * as fs from "fs";
import * as path from "path";
import * as crypto from "crypto";
import { md5 } from 'hash-wasm';

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

export function loadFilePath(p) {
    return fs.readFileSync(p);
}
