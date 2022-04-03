base=https://github.com/clusterfork/singlepp-references/releases/download

curl -L ${base}/hs-latest/BlueprintEncode_genes.csv.gz > BlueprintEncode_genes.csv.gz
curl -L ${base}/hs-latest/BlueprintEncode_matrix.csv.gz > BlueprintEncode_matrix.csv.gz
curl -L ${base}/hs-latest/BlueprintEncode_label_names_fine.csv.gz > BlueprintEncode_label_names_fine.csv.gz
curl -L ${base}/hs-latest/BlueprintEncode_labels_fine.csv.gz > BlueprintEncode_labels_fine.csv.gz
curl -L ${base}/hs-latest/BlueprintEncode_markers_fine.gmt.gz > BlueprintEncode_markers_fine.gmt.gz

curl -L ${mbase}/mm-latest/ImmGen_genes.csv.gz > ImmGen_genes.csv.gz
curl -L ${mbase}/mm-latest/ImmGen_matrix.csv.gz > ImmGen_matrix.csv.gz
curl -L ${mbase}/mm-latest/ImmGen_label_names_fine.csv.gz > ImmGen_label_names_fine.csv.gz
curl -L ${mbase}/mm-latest/ImmGen_labels_fine.csv.gz > ImmGen_labels_fine.csv.gz
curl -L ${mbase}/mm-latest/ImmGen_markers_fine.gmt.gz > ImmGen_markers_fine.gmt.gz
