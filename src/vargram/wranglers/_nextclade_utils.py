"""Module to perform auxiliary tasks for nextread function."""

import re
import subprocess
import os

def check_file_extension(file_path, valid_extensions):
    """Checks for the validity of the file extension.
    
    Args: 
        file_path (string): File path to be checked.
        valid_extensions (list): List of valid extensions to be checked against.

    Returns:
        file_extension (string): Extension of the file path provided.
        validity (bool): Truth value of validity of the file extension.
    """

    file_extension = file_path.lower().split('.')[-1]

    if file_extension in valid_extensions:
        return file_extension, True
    else:
        return file_extension, False

def check_reference(ref):
    """Checks whether provided reference name is in Nextclade.
    
    Args:
        ref (string): The dataset name.

    Returns:
        None.
    """

    dataset_names_only = subprocess.check_output(['nextclade', 'dataset', 'list', '--only-names']).decode(encoding='utf-8').split()  
    dataset_list_full = subprocess.check_output(['nextclade', 'dataset', 'list']).decode(encoding='utf-8')
    shortcut_parentheses = r'\(shortcuts:(.*?)\)' # getting all shortcuts enclosed in parentheses
    shortcut_quotes = r'"(.*?)"' # getting the individual shortcuts enclosed in double quotes
    dataset_shortcuts = re.findall(shortcut_parentheses, dataset_list_full)
    dataset_shortcuts = ' '.join(dataset_shortcuts)
    dataset_shortcuts = re.findall(shortcut_quotes, dataset_shortcuts)

    if ref not in dataset_names_only and ref not in dataset_shortcuts:
        try:
            raise ValueError("Nextclade reference name not recognized. See valid names and shortcuts below.")
        except ValueError as e:
            dataset_list_command = ['nextclade', 'dataset', 'list']
            subprocess.run(dataset_list_command, check=True)
    return None
    
def input_checker(kwargs):
    """Takes the keyword arguments of nextread() and checks for errors.

    Args:
        kwarg (dict): Set of keyword arguments passed to nextread().

    Returns:
        None.
    """

    expected_keys = {"seq", "ref", "gene"}#, "meta"}

    # Raise error for unknown keys
    unexpected_keys = set(kwargs.keys()) - expected_keys
    if unexpected_keys:
        raise TypeError(f"Unexpected keyword(s): {', '.join(unexpected_keys)}")
    
    # Raise error when "seq" is provided but "ref" isn't 
    if "seq" in kwargs and "ref" not in kwargs:
        raise ValueError("Path of reference sequence is not provided.")
        
    # Raise error when "ref" is provided but "seq" isn't
    if "ref" in kwargs and "seq" not in kwargs:
        raise ValueError("Path of sequences is not provided.")  

    # Raise error when "gene" is not provided and "ref" is a file
    if "gene" not in kwargs and os.path.isfile(kwargs["ref"]):
        raise ValueError("Genome annotation is not provided.")   
        
    for key in kwargs.keys():
        valid_extensions = {"seq": ["fa", "fasta"], 
                                "ref": ["fa", "fasta"],
                                "gene": ["gff", "gff3"],
                                "meta": ["csv", "tsv"]}
        file_type = {"seq": "sequence",
                         "ref": "reference",
                         "gene": "gene annotation",
                         "meta": "metadata"}
            
        if key == "seq" and not os.path.isfile(kwargs["seq"]):
            continue

        if key == "ref" and not os.path.isfile(kwargs["ref"]):
            check_reference(kwargs["ref"]) # Check if provided reference name is recognized instead
            continue

        file_extension, validity = check_file_extension(kwargs[key], valid_extensions[key])
        if not validity:
            raise ValueError(f"Unsupported {file_type[key]} file format: {file_extension}")
        
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