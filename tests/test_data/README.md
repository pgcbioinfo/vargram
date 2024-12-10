## Sequence files

Except for the reference, the 50 AY.1 (`sc2_AY1_n50.fasta`) and 50 P.3 (`sc2_P3_n50.fasta`) samples were sequenced at the Philippine Genome Center with funding support from the Epidemiology Bureau of the Department of Health. These sequences were initially shared publicly through the [EpiCoV database of the Global Initiative for 221 Sharing All Influenza Data (GISAID)](https://gisaid.org/).

The reference sequence (`sc2_wuhan_2019.fasta`) was downloaded from [Nextclade](https://github.com/nextstrain/nextclade_data/tree/master/data/nextstrain/sars-cov-2/wuhan-hu-1/orfs) along with the genome annotation (`sc2.gff`).

## Analysis files

[Analysis output files](https://docs.nextstrain.org/projects/nextclade/en/stable/user/output-files/04-results-tsv.html) of the 50 AY.1 (`AY1_analysis_cli.tsv`) and 50 P.3 (`P3_analysis_cli.tsv`) were produced by running Nextclade CLI v3.9.1 and the abovementioned FASTA files as input. These files were concatenated to form `batches_analysis_cli.tsv`.

The remaining files were produced during testing (`test_profile_input.py`).

