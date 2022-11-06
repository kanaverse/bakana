import * as inputs from "./../../src/steps/inputs.js";

test("modality guessing works as expected", () => {
    {
        let output = inputs.guessDefaultModalities([ { features: { "gene expression": null } } ]);
        expect(output.RNA).toEqual(["gene expression"]);
        expect("ADT" in output).toBe(false);
    }

    {
        let output = inputs.guessDefaultModalities([ { features: { "RNA": null, "Antibody Capture": null } } ]);
        expect(output.RNA).toEqual(["RNA"]);
        expect(output.ADT).toEqual(["Antibody Capture"]);
    }

    // works with multiple samples.
    {
        let output = inputs.guessDefaultModalities([ 
            { features: { "RNA": null, "Antibody Capture": null } },
            { features: { "gex": null, "hto": null } } 
        ]);
        expect(output.RNA).toEqual(["RNA", "gex"]);
        expect(output.ADT).toEqual(["Antibody Capture", "hto"]);
    }

    // Drops incomplete modalities.
    {
        let output = inputs.guessDefaultModalities([ 
            { features: { "RNA": null, "Antibody Capture": null } },
            { features: { "gex": null } }
        ]);
        expect(output.RNA).toEqual(["RNA", "gex"]);
        expect("ADT" in output).toBe(false);
    }

    // Pattern priority is respected.
    {
        let output = inputs.guessDefaultModalities([ 
            { features: { "RNA": null, "": null, "ADT": null, "HTO": null } },
            { features: { "": null, "RNA": null, "HTO": null, "ADT": null } }
        ]);
        expect(output.RNA).toEqual(["RNA", "RNA"]);
        expect(output.ADT).toEqual(["ADT", "ADT"]);
    }
})
