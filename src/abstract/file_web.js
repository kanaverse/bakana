/**
 * 'obj' is expected to be a File object.
 */

export function size(obj) {
    return obj.size;
}

export function buffer(obj) {
    let reader = new FileReaderSync();
    return reader.readAsArrayBuffer(obj);
}
