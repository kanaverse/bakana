import * as fs from "fs";
import * as pp from "path";

export function fsexists() {
    return true;
}

export function read(dir, path, asBuffer) {
    let loc = pp.join(dir, path);
    if (asBuffer) {
        return new Uint8Array(fs.readFileSync(loc));
    } else {
        return loc;
    }
}

export function write(dir, path, x) {
    fs.writeFileSync(pp.join(dir, path), x);
}

export function mkdir(dir, path) {
    fs.mkdirSync(pp.join(dir, path), { recursive: true });
}

export function copy(dir, from, to) {
    fs.copyFileSync(pp.join(dir, from), pp.join(dir, to));
}
