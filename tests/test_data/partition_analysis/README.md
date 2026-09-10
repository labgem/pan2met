Test cli:

``` console
python3 -m pan2met.analysis.partition --pangenome /home/sortion/Documents/projects/analysis/NILSmetabolism/results/ppanggolin/pangenome.h5 --reactions ./tests/test-data/partition_analysis/reactions.list --reaction-enzyme ./tests/test-data/partition_analysis/reactions_enzyme_genes.tsv --output /tmp/test_analysis_reaction_partition.tsv
diff /tmp/test_analysis_reaction_partition.tsv ./tests/test-data/partition_analysis/expected_result.tsv

```
