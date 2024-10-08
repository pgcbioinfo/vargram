"""Module to perform auxiliary tasks for nextclade function."""

import re
import subprocess
import os


def check_file_extension(file_path, valid_extensions):
    """Checks for the validity of the file extension.
    
    Parameters
    ----------
    file_path : str
        File path to be checked.
    valid_extensions : list
        List of valid extensions to be checked against.

    Returns
    -------
    file_extension : str 
        Extension of the file path provided.
    validity : bool
        Truth value of validity of the file extension.

    """
    file_extension = file_path.lower().split('.')[-1]
    if file_extension in valid_extensions:
        return file_extension, True
    else:
        return file_extension, False


def check_reference(ref):
    """Checks whether provided reference name is in Nextclade.
    
    Parameters
    ----------
    ref : str 
        The dataset name.

    Returns
    -------
    None

    Raises
    ------
    ValueError
        If Nextclade reference dataset name is not recognized.

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
            name_not_recognized = ("Nextclade reference name not recognized. " 
                                    "See valid names and shortcuts below.")
            raise ValueError(name_not_recognized)
        except ValueError as e:
            dataset_list_command = ['nextclade', 'dataset', 'list']
            subprocess.run(dataset_list_command, check=True)
    return None


def input_checker(kwargs):
    """Takes the keyword arguments of nextclade() and checks for errors.

    Parameters
    ----------
    kwarg : dict 
        Set of keyword arguments passed to nextclade().

    Returns
    -------
    None

    Raises
    ------
    ValueError
        If 'seq' and 'ref' are not provided together.
        If 'gene' is not provided when 'ref' is not a Nextclade reference dataset name.
        If file extension is not expected.
    TypeError
        If an unexpected keyword is provided.

    """
    # Raise error for unknown keys
    expected_keys = {"seq", "ref", "gene"}#, "meta"}
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
    
    # Checking extensions
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

    Parameters
    ----------
    aa_mutation : str
        The amino acid mutation as given by the nomenclature of Nextclade,
        e.g. ORF1b:G662S.
    part : str
        The part (gene name or position) of the mutation to be retrieved. 
        Part can also indicate removal.

    Returns
    -------
    str 
        The gene name or position on which the mutation occurs. 
        Or the mutation without the gene prefix.
    
    Raises
    ------
    ValueError
        If position or gene is not parsed.

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
        pattern = r'(\d+)'

    # Obtaining the match
    match = re.findall(pattern, aa_mutation)
    if not match:
        raise ValueError(f"Failed to parse mutation for gene or position: {aa_mutation}.")
    retrieved  = match[-1]
    return retrieved


def get_mutation_type(mutation):
    """Creates a new column for mutation type.
    
    Parameters
    ----------
    mutation : str
        The gene-stripped mutation name.

    Returns
    -------
    str
        The type of mutation.

    """
    if '-' in mutation:
        return 'del'
    elif ':' in mutation:
        return 'in'
    else:
        return 'sub'


def process_nextclade(nextclade_output):
    """Gets the unique mutations and their individual counts.

    Parameters
    ----------
    nextclade_output : pandas.DataFrame
        A subset (for multiple batches) of the DataFrame produced by nextclade().
    
    Returns
    -------
    pandas.DataFrame
        A DataFrame of unique mutations and their counts.

    """
    # Nextclade columns
    aa_sub = 'aaSubstitutions'
    aa_del = 'aaDeletions'
    aa_ins = 'aaInsertions'

    # Compiling all amino acid mutations (substitution, deletion, insertion) into a single column
    single_column = nextclade_output.copy()
    single_column[aa_sub] = single_column[aa_sub].fillna('')
    single_column[aa_del] = single_column[aa_del].fillna('')
    single_column[aa_ins] = single_column[aa_ins].fillna('')
    single_column['mutation'] = single_column.apply(lambda row: ','.join(filter(None, [row[aa_sub], row[aa_del], row[aa_ins]])), axis=1)
    single_column['mutation'] = single_column['mutation'].apply(lambda x: x.split(',')).copy()

    # Generating a row per mutation
    # Each row therefore may not be unique
    exploded = single_column.explode('mutation')

    # Creating gene column and removing gene prefix of mutations
    exploded['gene'] = exploded['mutation'].apply(lambda x: parse_mutation(x, 'gene'))
    exploded['mutation'] = exploded['mutation'].apply(lambda x: parse_mutation(x, 'gene_removal'))

    # Rearranging
    exploded.reset_index(drop=True, inplace=True)
    processed_nextclade = exploded[['batch', 'gene', 'mutation']]
    return processed_nextclade