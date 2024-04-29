import pandas as pd

def pivot_batch_counts(batch_counts, batches, yaxis):
    """Generates a summary count of mutations.
    
    Args:
        batch_counts (pandas.DataFrame): A DataFrame of mutation counts with labeled batches.
    
    Returns:
        pivoted_batch_counts (pandas.DataFrame): Pivoted batch counts with unique mutations and separate columns for the batches.
    """

    pivoted_batch_counts = batch_counts.pivot_table(index=['gene', 'position', 'mutation'], columns='batch', values='count', fill_value=0).reset_index()
    pivoted_batch_counts = pivoted_batch_counts.rename_axis(None, axis=1)

    if yaxis == 'Weights':
        for batch in batches:
            batch_sum = pivoted_batch_counts[batch].sum()
            pivoted_batch_counts[batch] = pivoted_batch_counts[batch] * 100 / batch_sum
    
    pivoted_batch_counts['sum'] = pivoted_batch_counts.drop(columns=['gene', 'position', 'mutation']).sum(axis=1)

    return pivoted_batch_counts

def add_key_mutations(pivoted_batch_counts, key_lineage, key_name):
    """Add key mutations from reference lineages
    
    Args:
        pivoted_batch_counts (pandas.DataFrame):  Pivoted batch counts with unique mutations and separate columns for the batches.
        key_lineage (pandas.DataFrame): The DataFrame containing key mutations per gene.
        key_name (str): The name of the key lineage.

    Returns:
        keyful_pivoted_counts (pandas.DataFrame): The DataFrame containing summary mutation counts per batch including key mutations.
    """
    key_lineage[key_name] = 1
    keyful_pivoted_counts = pd.merge(pivoted_batch_counts, key_lineage, on=['gene', 'position', 'mutation'], how='outer')
    keyful_pivoted_counts.fillna(0, inplace=True) # batch and sum columns are filled with 0
    keyful_pivoted_counts.sort_values(by='position', inplace=True)
    keyful_pivoted_counts.reset_index(drop=True, inplace=True)

    return keyful_pivoted_counts




    









    

