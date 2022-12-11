import * as scran from "scran.js";
import * as utils from "./general.js";

export class PcaStateBase {}

export class PcaMimic { 
    constructor(pcs, var_exp) {
        this.var_exp = var_exp;
        try {
            this.pcs = scran.createFloat64WasmArray(pcs.length);
            this.pcs.set(pcs);
        } catch (e) {
            utils.freeCache(this.pcs);
            throw e;
        }
    }

    principalComponents({ copy }) {
        return utils.mimicGetter(this.pcs, copy);
    }

    numberOfCells() {
        return this.pcs.length / this.numberOfPCs();
    }

    numberOfPCs() {
        return this.var_exp.length;
    }

    varianceExplained({ copy = true } = {}) {
        return utils.mimicGetter(this.var_exp, copy);
    }

    totalVariance () {
        return 1;
    }

    free() {
        this.pcs.free();
    }
}
