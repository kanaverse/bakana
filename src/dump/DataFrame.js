import * as wa from "wasmarrays.js";

// Monkey-patching these methods so that we can use these WasmArrays as columns in a bioc.DataFrame.
wa.Uint8WasmArray.prototype._bioconductor_LENGTH = function() { return this.length; };
wa.Int32WasmArray.prototype._bioconductor_LENGTH = function() { return this.length; };
wa.Float64WasmArray.prototype._bioconductor_LENGTH = function() { return this.length; };
