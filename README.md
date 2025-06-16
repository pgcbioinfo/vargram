![Docs Build Status](https://github.com/pgcbioinfo/vargram/actions/workflows/docs.yml/badge.svg?)
![Tests Status](https://github.com/pgcbioinfo/vargram/actions/workflows/tests.yml/badge.svg?)

<div style="text-align: center;">
    <img src="https://github.com/pgcbioinfo/vargram/blob/main/docs/assets/images/vargram_header.png?raw=True" alt="VARGRAM Header" />
</div>
<h1 style="text-align:center;">VARGRAM (Visual ARrays for GRaphical Analysis of Mutations)</h1>


🧬 VARGRAM is a Python package that makes it easy to generate insightful figures for genomic surveillance, born out of our experience during the COVID-19 pandemic. Currently, VARGRAM supports generating mutation profiles straight from sequence files by hooking into existing tools such as Nextclade. The figures can be easily customized within a Python script or Jupyter notebook using a declarative syntax.

🔥 We are actively developing VARGRAM into a full visualization library for common use cases in molecular epidemiology. More modules will be added in the coming months. If you have a feature request or find a bug, please [submit an issue](https://github.com/pgcbioinfo/vargram/issues). 

## Documentation

Full installation instructions and tutorials are available on [the VARGRAM documentation website](https://pgcbioinfo.github.io/vargram/).

## Installation

Install with [pip](https://pip.pypa.io/en/stable/):
```bash
pip install vargram
``` 
Python version ≥3.11 is required.

VARGRAM relies on [Nextclade](https://clades.nextstrain.org/) to perform mutation calling when sequence files are provided. Make sure to [download the Nextclade CLI](https://docs.nextstrain.org/projects/nextclade/en/stable/user/nextclade-cli/installation/index.html) and [add it to the path](https://pgcbioinfo.github.io/vargram/install_nextclade/#adding-an-executable-to-the-path). You may also just provide Nextclade's analysis CSV output directly and VARGRAM can still produce a mutation profile without Nextclade installed.

## Quickstart Guide

To produce a mutation profile, VARGRAM requires a single FASTA file (or a directory of FASTA files) of samples, a FASTA file for the reference, and a genome annotation file following the [GFF3](https://docs.nextstrain.org/projects/nextclade/en/stable/user/input-files/03-genome-annotation.html) format.

A mutation profile can be generated in just four lines of code:
```python
from vargram import vargram # Importing the package

vg = vargram(seq='path/to/<samples-directory>', # Provide sample sequences
            ref='path/to/<reference.fa>', # Provide reference sequence
            gene='path/to/<annotation.gff>') # Provide genome annotation
vg.profile() # Tell VARGRAM you want to create a mutation profile
vg.show() # And show the resulting figure
```

Alternatively, you can simply provide a CSV file. For example, you can upload your sequences to the [Nextclade web app](https://clades.nextstrain.org/) and download the analysis CSV output. VARGRAM recognizes this output and can process it:
```python
from vargram import vargram

vg = vargram(data='path/to/<nextclade_analysis.csv>') # Provide Nextclade analysis file
vg.profile()
vg.show()
```
Calling the mutation profile this way does not require Nextclade CLI to be installed.

## Sample Output

Install VARGRAM and try out the following snippet, which will download test data for you. Nextclade CLI does not need to be installed for the following example:
```python
# Import main VARGRAM module and module to download external data
from vargram import vargram 
from vargram.data import example

# Download test data into test_data directory
example.get('test_data') 

# Generate the mutation profile
vg = vargram(data='test_data/analysis/omicron_analysis_cli.tsv') # Provide data
vg.profile() # Tell VARGRAM you want to create a mutation profile
vg.show() # Show the figure
vg.save("default_profile.png", dpi=300) # Save the figure
```
This will produce the following figure:
<div style="text-align: center;">
    <img src="https://github.com/pgcbioinfo/vargram/blob/main/docs/assets/images/default_profile.png?raw=True" alt="mutation profile" />
</div>

Note that by default, VARGRAM favors placing genes with the most number of mutations first. The figure can be customized to show genes by their start position, to force a horizontal layout and other options:
```python
vg = vargram(data='test_data/analysis/omicron_analysis_cli.tsv', # Provide data
            gene='test_data/sc2.gff') # Provide annotation file
vg.profile(threshold=5, # Set minimum count for a mutation to be included
        ytype='counts') # Set y-axis to show raw count
vg.aes(stack_title='Region', # Change batch legend title
    stack_label=['Foreign', 'Local'], # Change batch names
    stack_color=['#009193', '#E33E84'], # Change batch bar colors
    order=True, # Order the genes based on the annotation file
    flat=True) # Force a horizontal layout
vg.key('test_data/keys/BA1_key.csv', label='BA.1') # Show key mutations of BA.1
vg.key('test_data/keys/BA2_key.csv', label='BA.2') # Show key mutations of BA.2
vg.show() # Show the figure
```
This results to the following figure:
<div style="text-align: center;">
    <img src="https://github.com/pgcbioinfo/vargram/blob/main/docs/assets/images/ordered_flat_profile.png?raw=True" alt="mutation profile with genes ordered" />
</div>