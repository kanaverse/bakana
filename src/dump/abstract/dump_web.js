export function initialize(host) {
    return {};
}

export function read(host, path, asBuffer) {
    return host[path];
}

export function write(host, path, x) {
    host[path] = x;
}

export function mkdir(host, path) {
    host[path] = null;
}

export function copy(host, from, to) {
    host[to] = host[from];
}
