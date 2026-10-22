#

We work an a diamond blastp run on _Thermococcus kodakarensis_ against the proteins of MetaCyc.

We first extract the best hits:

```bash
sort --key=1 --field-separator=$'\t' "blastp.tsv" > "./sorted_blastp.tsv"
awk -f "diamond_blastp_best_hit.awk" "./sorted_blastp.tsv" > "./best_hit.blastp.tsv"
```

To get a set of "confident" enzyme annotation:

```bash
confidence_min_identity=60
less_confidence_min_identity=40
awk -v confidence_min_identity="$confidence_min_identity" -v less_confidence_min_identity="$less_confidence_min_identity" ' /^[^#]/ {
    identity_percentage=$3
    if (identity_percentage > confidence_min_identity) {
        print $0 > "confident.blastp.tsv"
    } else if (identity_percentage > less_confidence_min_identity) {
        print $0 > "less_confident.blastp.tsv"
    }
} ' ./best_hit.blastp.tsv
```

Then we reuse a code sample from the `pan2met-wf` workflow, to get MetaCyc reaction identifiers:

```bash
monomer_to_reactions="/home/sortion/data/KEGG_MetaCyc_crossref/2025-11-06/metacyc_monomers_to_reactions.tsv"
awk -F $'\t' -v OFS=$'\t' \
        -f "hash_join.awk" \
        -v key1=1 -v value1=2 \
        -v key2=2 -v value2=1 \
            "$monomer_to_reactions" \
         "confident.blastp.tsv" > "confident_enzyme_reaction.tsv"
```


We create a .gt file from the PPanGGOLiN HDF5 file:
```bash
ppanggolin write_pangenome -p pangenome.h5 --gt --output ./output/ -f
```
