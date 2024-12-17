![Docs Build Status](https://github.com/cjpalpallatoc/test-repo/actions/workflows/docs.yml/badge.svg?)
![Tests Status](https://github.com/cjpalpallatoc/test-repo/actions/workflows/tests.yml/badge.svg?)

<div style="text-align: center;">
    <img src="https://github.com/pgcbioinfo/vargram/blob/main/docs/assets/images/vargram_header.png?raw=True" alt="VARGRAM Header" />
</div>
<h1 style="text-align:center;">VARGRAM (Visual ARrays for GRaphical Analysis of Mutations)</h1>


🧬 VARGRAM is a Python package that makes it easy to generate insightful figures for genomic surveillance, born out of our experience during the COVID-19 pandemic. With the latest update, VARGRAM can be used to generate mutation profiles straight from sequence files by hooking into existing tools such as Nextclade. The figures can be easily customized within a Python script or Jupyter notebook using a declarative syntax.

🔥 We are actively developing VARGRAM into a full visualization library for common use cases in molecular epidemiology. More modules will be added in the coming months. If you have a feature request or find a bug, please [submit an issue](https://github.com/pgcbioinfo/vargram/issues). 

## Documentation

Full installation instructions and tutorials are available on [the VARGRAM documentation website](https://pgcbioinfo.github.io/vargram/).

## Installation

Install with [pip](https://pip.pypa.io/en/stable/):
```bash
pip install vargram
``` 
Python version ≥3.11 is required.

VARGRAM relies on [Nextclade](https://clades.nextstrain.org/) to perform mutation calling when sequence files are provided. Make sure to [download the Nextclade CLI](https://docs.nextstrain.org/projects/nextclade/en/stable/user/nextclade-cli/installation/index.html) and add it to the path. (Alternatively, you may provide Nextclade's analysis CSV output directly and VARGRAM can still produce a mutation profile without Nextclade installed.)

## Quickstart Guide

To produce a mutation profile, VARGRAM requires a single FASTA file (or a directory of FASTA files) of samples, a FASTA file for the reference, and a genome annotation file following the [GFF3](https://docs.nextstrain.org/projects/nextclade/en/stable/user/input-files/03-genome-annotation.html) format.

A mutation profile can be generated in just four lines of code:
```python
from vargram import vargram # Importing the package

vg = vargram(seq='path/to/covid_samples/', # Provide sample sequences
            ref='path/to/covid_reference.fa', # Provide reference sequence
            gene='path/to/covid_annotation.gff') # Provide genome annotation
vg.profile() # Tell VARGRAM you want to create a mutation profile
vg.show() # And show the resulting figure
```

Alternatively, you can simply provide a CSV file. For example, you can upload your sequences to the [Nextclade web app](https://clades.nextstrain.org/) and download the analysis CSV output. VARGRAM recognizes this output and can process it:
```python
from vargram import vargram

vg = vargram(data='path/to/nextclade_analysis.csv') # Provide Nextclade analysis file
vg.profile()
vg.show()
```
Calling the mutation profile this way does not require Nextclade CLI to be installed.

## Sample Output

A sample mutation profile is shown below:
<div style="text-align: center;">
    <img src="https://github.com/pgcbioinfo/vargram/blob/main/docs/assets/images/readme_profile.png?raw=True" alt="mutation profile" />
</div>

Note that by default, VARGRAM favors placing genes with the most number of mutations first. Thus, the S gene above is shown at the very top. You may wish to force VARGRAM to show the genes in order based on the start position. You can do so by setting `vg.profile(order=True)`. The mutation profile above will then look like the following:

<div style="text-align: center;">
    <img src="https://github.com/pgcbioinfo/vargram/blob/main/docs/assets/images/ordered_profile.png?raw=True" alt="mutation profile with genes ordered" />
</div>