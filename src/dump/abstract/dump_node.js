import * as fs from "fs";
import * as path from "path";

export function fsexists() {
    return true;
}

export function read(dir, path, asBuffer) {
    let loc = path.join(dir, path);
    if (asBuffer) {
        return new Uint8Array(fs.readFileSync(loc));
    } else {
        return loc;
    }
}

export function write(dir, path, x) {
    fs.writeFileSync(path.join(dir, path), x);
}

export function mkdir(dir, path) {
    fs.mkdirSync(path.join(dir, path), { recursive: true });
}

export function copy(dir, from, to) {
    fs.copyFileSync(path.join(dir, from), path.join(dir, to));
}
