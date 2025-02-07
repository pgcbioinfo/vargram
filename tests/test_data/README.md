## Sequence files

Except for the reference, the 80 BA.1 (`sc2_BA1_n80.fasta`) and 80 BA.2 (`sc2_BA2_n80.fasta`) samples were sequenced at the Philippine Genome Center with funding support from the Epidemiology Bureau of the Department of Health. These sequences were initially shared publicly through the [EpiCoV database of the Global Initiative for 221 Sharing All Influenza Data (GISAID)](https://gisaid.org/).

The reference sequence (`sc2_wuhan_2019.fasta`) was downloaded from [Nextclade](https://github.com/nextstrain/nextclade_data/tree/master/data/nextstrain/sars-cov-2/wuhan-hu-1/orfs) along with the genome annotation (`sc2.gff`).

## Analysis files

[Analysis output files](https://docs.nextstrain.org/projects/nextclade/en/stable/user/output-files/04-results-tsv.html) of the 80 BA.1 (`BA1_analysis_cli.tsv`) and 80 BA.2 (`BA2_analysis_cli.tsv`) sequences were obtained by running Nextclade CLI v3.9.1 and the above-mentioned reference and genome annotation files as input. The analysis files were then concatenated to form `omicron_analysis_cli.tsv`. The file `XBB_analysis_web.tsv` was downloaded from Nextclade Web using the example XBB sequences as input.
