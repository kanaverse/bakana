import * as fs from "fs";
import * as path from "path";

export function initialize(host) {
    return host;
}

export function read(host, path, asBuffer) {
    let loc = path.join(host, path);
    if (asBuffer) {
        return new Uint8Array(fs.readFileSync(loc));
    } else {
        return loc;
    }
}

export function write(host, path, x) {
    fs.writeFileSync(path.join(host, path), x);
}

export function mkdir(host, path) {
    fs.mkdirSync(path.join(host, path), { recursive: true });
}

export function copy(host, from, to) {
    fs.copyFileSync(path.join(host, from), path.join(host, to));
}

export function join(host, path) {
    return path.join(host, path);
}
