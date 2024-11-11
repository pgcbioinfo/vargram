"""Main module for reading user input data and processing with Nextclade"""

# Types of sequences input
# 1. A single FASTA containing multiple sequences (1 batch)
# 2. A directory containing N FASTA files, each containing multiple sequences (N batches)

from ._nextclade_utils import input_checker
from ._nextclade_cli import create_command, capture_output
import os
import pandas as pd
import tempfile
import shutil


def nextclade(**kwargs):
    """Takes input sequence FASTA files and transforms it 
    into a DataFrame to be plotted by VARGRAM.

    Parameters
    ----------
    seq : str
        FASTA file path of the sequences.
    ref : str
        FASTA file path of the reference sequence.
    gene : str
        GFF3 file path of the genome annotation.
    
    Returns
    -------
    None

    Raises
    ------
    ValueError
        If Nextclade analysis dataframe is empty.

    """
    input_checker(kwargs)
    try:
        # Creating secure temporary directory to store Nextclade analysis output file
        secure_analysis_dir = tempfile.mkdtemp(prefix="secure_analysis_dir")

        # Creating secure temporary directory to store Nextclade 
        # reference and genome annotation if Nextclade reference name is provided
        secure_ref_dir = tempfile.mkdtemp(prefix="secure_ref_dir")

        if os.path.isdir(kwargs["seq"]): # Case 1: A directory of FASTA files is provided
            batches = os.listdir(kwargs["seq"])
            try: # Removing nuissance folder attribute file created in macOS folders 
                batches.remove('.DS_Store')  
            except ValueError:
                pass

            # Getting Nextclade analysis output per FASTA batch and concatenating 
            # the dataframes
            seq_dir = kwargs["seq"]
            outputs = []
            kwargs_mod = kwargs
            for batch in batches:
                # Getting Nexctlade output
                kwargs_mod["seq"] = os.path.join(seq_dir,batch)
                nextclade_command, gene_path = create_command(input = kwargs_mod, secure_analysis_dir = secure_analysis_dir, secure_ref_dir = secure_ref_dir) 
                out = capture_output(nextclade_command)

                # Appending output
                batch_name = os.path.splitext(batch)[0]
                out.insert(0, 'batch', batch_name)
                outputs.append(out)
            nextclade_output = pd.concat(outputs, ignore_index=True)   
        else: # Case 2: One FASTA file provided
            nextclade_command, gene_path = create_command(input = kwargs, secure_analysis_dir = secure_analysis_dir, secure_ref_dir = secure_ref_dir) 
            nextclade_output = capture_output(nextclade_command)

            # Add batch name
            batch_name = os.path.basename(kwargs['seq'])
            batch_name = os.path.splitext(batch_name)[0]
            nextclade_output.insert(0, 'batch', batch_name)

        # Sorting by batch name and seq name:
        nextclade_output.sort_values(by=['batch', 'seqName'], inplace=True)
        nextclade_output.reset_index(drop=True, inplace=True)

        # Getting annotation
        gff_columns = ["seqname", "source", "feature", "start", "end", "score",
                       "strand", "frame", "attribute"]
        annotation = pd.read_csv(gene_path, sep="\t", comment="#", 
                                 header=None, names=gff_columns)

    # Remove created directories
    finally:
        if os.path.exists(secure_analysis_dir):
            shutil.rmtree(secure_analysis_dir)
        if os.path.exists(secure_ref_dir):
            shutil.rmtree(secure_ref_dir)

    # Removing sequences with warnings or errors
    nextclade_output.drop(nextclade_output.dropna(subset=['warnings']).index)
    nextclade_output.drop(nextclade_output.dropna(subset=['errors']).index)
    if nextclade_output.empty:
        raise ValueError("Nextclade analysis DataFrame is empty.")
    return nextclade_output, annotation