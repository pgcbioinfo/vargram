## Sequence files (`sequences/`)

Except for the reference, the 80 BA.1 (`sc2_BA1_n80.fasta`) and 80 BA.2 (`sc2_BA2_n80.fasta`) samples were sequenced at the Philippine Genome Center with funding support from the Epidemiology Bureau of the Department of Health. These sequences were initially shared publicly through the [EpiCoV database of the Global Initiative for 221 Sharing All Influenza Data (GISAID)](https://gisaid.org/).

The reference sequence (`sc2_wuhan_2019.fasta`) was downloaded from [Nextclade](https://github.com/nextstrain/nextclade_data/tree/master/data/nextstrain/sars-cov-2/wuhan-hu-1/orfs) along with the genome annotation (`sc2.gff`).

## Analysis files (`analysis/`)

[Analysis output files](https://docs.nextstrain.org/projects/nextclade/en/stable/user/output-files/04-results-tsv.html) of the 80 BA.1 (`BA1_analysis_cli.tsv`) and 80 BA.2 (`BA2_analysis_cli.tsv`) sequences were obtained by running Nextclade CLI v3.9.1 and the above-mentioned reference and genome annotation files as input. The analysis files were then concatenated to form `omicron_analysis_cli.tsv`. The file `XBB_analysis_web.tsv` was downloaded from Nextclade Web using the example XBB sequences as input.

## Key files (`keys/`)

The key files `BA1_key.csv` and `BA2_key.csv` were generated from sets of BA.1-assigned and BA.2-assigned sequences downloaded from GISAID. The mutations here largely match the designated mutations of these lineages on [outbreak.info](https://outbreak.info).

## Summary files (`summary/`)

The summary files were generated from the same 160 BA.1 and BA.2 sequences (see [documentation](https://pgcbioinfo.github.io/vargram/mutation_profile/#getting-summary-data)). The file `omicron_summary_data_t1.csv` is created from the two Omicron FASTA files but with the threshold set to one.

## Metadata (`metadata/`)

The sole metadata file `voc_metadata.csv` consists of metadata information from a subset of SARS-CoV-2 sequences from the Philippines.