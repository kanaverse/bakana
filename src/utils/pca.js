import * as scran from "scran.js";
import * as utils from "./general.js";

export function formatPCs(pcs) {
    return {
        "pcs": pcs.principalComponents({ copy: "view" }),
        "num_pcs": pcs.numberOfPCs(),
        "num_obs": pcs.numberOfCells()
    };
}

export function formatSummary(pcs) {
    var var_exp = pcs.varianceExplained();
    var total_var = pcs.totalVariance();
    var_exp.forEach((x, i) => {
        var_exp[i] = x/total_var;
    });
    return { 
        "var_exp": var_exp 
    };
}

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
