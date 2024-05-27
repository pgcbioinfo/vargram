"""Module to perform auxiliary tasks for VARGRAM class methods."""

import pandas as pd
import os

def read_comma_or_tab(file):
    root, extension = os.path.splitext(file)

    if extension == '.csv':
        return pd.read_csv(file)
    elif extension == '.tsv':
        return pd.read_csv(file, delimiter='\t')
    else:
        raise ValueError("Incorrect file format. Expecting CSV or TSV file.")