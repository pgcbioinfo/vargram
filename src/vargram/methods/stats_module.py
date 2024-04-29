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

def count(nextread_subset):
    """Gets the unique mutations and their individual counts.

    Args:
        nextread_subset (pandas.DataFrame): A subset (for multiple batches) of the DataFrame produced by nextread().
    
    Returns:
        pandas.DataFrame: A DataFrame of unique mutations and their counts.
    """
    c_nextread_subset = nextread_subset.copy()

    # Compiling all amino acid mutations (substitution, deletion, insertion) into a single column
    c_nextread_subset['aaSubstitutions'].fillna('', inplace=True)
    c_nextread_subset['aaDeletions'].fillna('', inplace=True)
    c_nextread_subset['aaInsertions'].fillna('', inplace=True)

    c_nextread_subset['mutations'] = c_nextread_subset.apply(lambda row: ','.join(filter(None, [row['aaSubstitutions'], row['aaDeletions'], row['aaInsertions']])), axis=1)
    aa_combined = c_nextread_subset[['mutations']]
    aa_combined = aa_combined['mutations'].apply(lambda x: x.split(','))

    # Generating a row per mutation
    # Each row therefore may not be unique
    aa_exploded = aa_combined.explode('mutations')

    # Getting counts per unique mutation
    aa_counted = aa_exploded.value_counts()

    # Compiling unique mutations and their counts
    aa = pd.DataFrame({'mutation': aa_counted.index, 'count': aa_counted.values})

    # Creating new column for gene and position
    aa['gene'] = aa['mutation'].apply(lambda x: parse_mutation(x, 'gene'))
    aa['position'] = aa['mutation'].apply(lambda x: int(parse_mutation(x, 'position')))

    # Removing gene prefix of mutations
    aa['mutation'] = aa['mutation'].apply(lambda x: parse_mutation(x, 'gene_removal'))

    # Adding mutation type
    aa['type'] = aa['mutation'].apply(lambda x: get_mutation_type(x))

    # Sorting by position
    aa.sort_values(by='position', inplace=True)
    aa.reset_index(drop=True, inplace=True)

    # Rearranging columns
    aa = aa[['gene', 'position', 'mutation', 'type', 'count']]

    return aa

def count_per_batch(nextread_out, threshold = 10):
    """Gets the individual mutation counts per batch.

    Args:
        nextread_out (pandas.DataFrame): The DataFrame output of nextread().
        threshold (integer): The minimum no. of mutation counts for a mutation to be kept.

    Returns:
        pandas.DataFrame: A DataFrame of mutation counts with labeled batches.
    """

    # Getting unique batch names
    batches = nextread_out['batch'].unique()

    # Getting mutation counts per batch
    outputs = []
    for batch in batches:
        # Getting counts per batch
        batch_df = nextread_out[nextread_out['batch'] == batch]
        batch_mutations = count(batch_df)
        batch_mutations.insert(0, 'batch', batch)

        #Appending output
        outputs.append(batch_mutations)
    
    aa_batched = pd.concat(outputs, ignore_index=True)  
    
    # Applying threshold
    aa_batched = aa_batched[aa_batched['count'] > threshold]
    aa_batched = aa_batched.reset_index(drop=True)

    if aa_batched.empty:
        raise ValueError("Mutation count DataFrame is empty. Possibly no mutation met count threshold.", UserWarning)

    return aa_batched


def count_per_gene(geneless_counts):
    """Gets the mutation counts per gene

    Args:
        geneless_counts (pandas.DataFrame): The DataFrame containing mutation counts per batch.

    Returns:
        geneful_counts (pandas.DataFrame): The DataFrame containing genes and their no. of unique mutations.
    """
    
    genes = geneless_counts['gene'].unique().tolist()

    unique_counts = []
    for gene in genes:
        genes_df = geneless_counts[geneless_counts['gene'] == gene]
        unique_counts.append(len(genes_df['mutation'].unique()))

    geneful_counts = pd.DataFrame({'gene': genes, 'count': unique_counts})
    geneful_counts.sort_values(by='count', ascending=False, inplace=True)
    geneful_counts.reset_index(drop=True, inplace=True)

    return geneful_counts

def add_key_mutations(keyless_counts, key_lineage, key_name):
    """Add key mutations from reference lineages
    
    Args:
        keyless_counts (pandas.DataFrame): The DataFrame containing mutation counts per batch.
        key_lineage (pandas.DataFrame): The DataFrame containing key mutations per gene.

    Returns:
        keyful_counts (pandas.DataFrame): The DataFrame containing mutation counts per batch including key mutations.
    """
    # Add mutation type column and count to key lineage dataframe
    key_lineage['type'] = key_lineage['mutation'].apply(lambda x: get_mutation_type(x))
    key_lineage['count'] = 0

    # Add key column to keyless dataframe
    keyless_counts['key'] = pd.NA

    # 
    key_mutations = key_lineage[['gene', 'position', 'mutation', 'type']].values.tolist()
    missing_mutations = []

    # Get unique mutations for the first batch
    keyless_batches = keyless_counts['batch'].unique()
    first_batch = keyless_batches[0]
    keyful_batch = keyless_counts[keyless_counts['batch'] == first_batch]
    batch_mutations = keyful_batch[['gene', 'mutation']].values.tolist()
        
    # Finding missing key mutations for the current batch
    missing = [mutation for mutation in key_mutations if mutation not in batch_mutations]
        
    # Extending missing_mutations list with missing combinations along with other columns
    missing_mutations.extend([(first_batch, mutation[0], mutation[1], mutation[2], mutation[3], 0, key_name) for mutation in missing])

    # Creating a DataFrame from missing_mutations
    missing_counts = pd.DataFrame(missing_mutations, columns=['batch', 'gene', 'position', 'mutation', 'type', 'count', 'key'])

    # Concatenate keyless_counts and missing_counts
    keyful_counts = pd.concat([keyless_counts, missing_counts], ignore_index=True)
    keyful_counts.sort_values(by=['batch', 'position'], inplace=True)

    return keyful_counts


















    



    




