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
