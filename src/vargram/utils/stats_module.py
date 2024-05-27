"""Module to generate a dataframe of unique mutations and their counts."""

from numpy import split
import pandas as pd
import re
import warnings

def parse_mutation(aa_mutation, part):
    """Gets the gene or position from the mutation. Can also remove the gene prefix of the mutation.

    Args:
        aa_mutation (string): The amino acid mutation as given by the nomenclature of Nextclade, e.g. ORF1b:G662S.
        part (string): The part (gene name or position) of the mutation to be retrieved. Part can also indicate removal.

    Returns:
        string: The gene name or position on which the mutation occurs. Or the mutation without the gene prefix.
    """
    
    # Removing prefix
    if part == 'gene_removal':
        pattern = r'^[^:]+:'
        stripped = re.sub(pattern, '', aa_mutation)
        return stripped

    # Patterns for the gene name and position
    if part == 'gene':
        pattern = r'^([^:]+):'
    
    if part == 'position':
        pattern = r':.*?([0-9]+)'

    # Obtaining the match
    match = re.search(pattern, aa_mutation)
    if not match:
        print(part, aa_mutation)
        raise ValueError("Failed to parse mutation for gene or position.")
    retrieved  = match.group(1)

    return retrieved

def get_mutation_type(mutation):
    """Creates a new column for mutation type.
    
    Args:
        mutation (string): The gene-stripped mutation name.

    Returns:
        string: The type of mutation.
    """

    if '-' in mutation:
        return 'del'
    elif ':' in mutation:
        return 'in'
    else:
        return 'sub'

def process_nextread(nextread_output):
    """Gets the unique mutations and their individual counts.

    Args:
        nextread_subset (pandas.DataFrame): A subset (for multiple batches) of the DataFrame produced by nextread().
    
    Returns:
        pandas.DataFrame: A DataFrame of unique mutations and their counts.
    """
    # Nextclade columns
    aa_sub = 'aaSubstitutions'
    aa_del = 'aaDeletions'
    aa_ins = 'aaInsertions'


    # Compiling all amino acid mutations (substitution, deletion, insertion) into a single column
    single_column = nextread_output.copy()

    single_column[aa_sub] = single_column[aa_sub].fillna('')
    single_column[aa_del] = single_column[aa_del].fillna('')
    single_column[aa_ins] = single_column[aa_ins].fillna('')

    single_column['mutation'] = single_column.apply(lambda row: ','.join(filter(None, [row[aa_sub], row[aa_del], row[aa_ins]])), axis=1)
    single_column['mutation'] = single_column['mutation'].apply(lambda x: x.split(',')).copy()

    # Generating a row per mutation
    # Each row therefore may not be unique
    exploded = single_column.explode('mutation')


    # Creating new column for gene and position
    exploded['gene'] = exploded['mutation'].apply(lambda x: parse_mutation(x, 'gene'))
    exploded['position'] = exploded['mutation'].apply(lambda x: int(parse_mutation(x, 'position')))

    # Removing gene prefix of mutations
    exploded['mutation'] = exploded['mutation'].apply(lambda x: parse_mutation(x, 'gene_removal'))

    # Adding mutation type
    exploded['type'] = exploded['mutation'].apply(lambda x: get_mutation_type(x))

    # Sorting by position
    exploded.sort_values(by='position', inplace=True)
    exploded.reset_index(drop=True, inplace=True)

    # Rearranging columns
    processed_nextread = exploded[['batch', 'gene', 'position', 'mutation', 'type']]

    return processed_nextread