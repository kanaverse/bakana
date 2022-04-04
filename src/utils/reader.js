import * as pako from "pako";
import * as afile from "../abstract/file.js";

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

var cache = {
    file2link: null, 
    link2file: null
};

export function setCreateLink(fun) {
    cache.file2link = fun;
    return;
}

export function setResolveLink(fun) {
    cache.link2file = fun;
    return;
}

export async function standardSerialize(details, type, embeddedSaver) {
    let output = { 
        "type": type,
        "name": details.name 
    };
    let serialized = details.content.serialized();

    if (embeddedSaver !== null) {
        let eout = await embeddedSaver(serialized, details.content.size());
        output.offset = eout.offset;
        output.size = eout.size;
    } else {
        let fun = cache.file2link;
        if (fun === null) {
            throw new Error("link-creating function has not been set by 'setCreateLink'");
        }
        output.id = await fun(serialized);
    }

    return output;
}

export async function standardUnserialize(details, embeddedLoader) {
    let output = { name: details.name };

    if ("id" in details) {
        let fun = cache.link2file;
        if (fun === null) {
            throw new Error("link-resolving function has not been set by 'setResolveLink'");
        }
        output.content = new afile.LoadedFile(await fun(details.id));
    } else {
        output.content = new afile.LoadedFile(await embeddedLoader(details.offset, details.size));
    }

    return output;
}

export function formatFile(file, sizeOnly) {
    let output = { "name": afile.rawName(file) };
    if (sizeOnly) {
        output.size = afile.rawSize(file);
    } else {
        output.content = new afile.LoadedFile(file);
    }
    return output;
}
