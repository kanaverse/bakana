import * as pako from "pako";

export function extractHDF5Strings(handle, name) {
    if (!(name in handle.children)) {
        return null;
    }

    if (handle.children[name] !== "DataSet") {
        return null;
    }

    let content = handle.open(name);
    if (content.type !== "String") {
        return null;
    }

    return content.load();
}

export function readTextLines(buffer, compression = "gz") {
    let txt = buffer;
    if (compression == "gz") {
        txt = pako.ungzip(buffer);
    }

    const dec = new TextDecoder();
    let decoded = dec.decode(txt);

    let lines = decoded.split("\n");
    if (lines.length > 0 && lines[lines.length - 1] == "") { // ignoring the trailing newline.
        lines.pop();
    }

    return lines;    
}

export function readDSVFromBuffer(content, fname, delim = "\t") {
    var ext = fname.name.split('.').pop();
    let lines = readTextLines(content, ext);
    lines.forEach((x, i) => {
        lines[i] = x.split(delim);
    });
    return lines;
}

export function generateRandomName(prefix = "", suffix = "") {
    return prefix + String(Number(new Date())) + suffix
}
