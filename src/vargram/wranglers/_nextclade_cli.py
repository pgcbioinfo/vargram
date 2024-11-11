"""Module that creates Nextclade command and captures Nextclade output."""

import subprocess
import os
import pandas as pd


def create_command(**kwargs):
    """Creates Nextclade command based on user input.
    
    Parameters
    ----------
    input : dict 
        Dictionary containing the user input to nextclade().
    secure_analysis_dir : str
        Secure temporary directory to store Nexclade analysis TSV output.
    secure_ref_dir : str
        Secure temporary directory to store reference sequence from Nextclade.

    Returns
    -------
    list
        A list whose elements constitute the Nextclade command.

    """

    # Getting user input
    input = kwargs["input"]
    seq = input["seq"]
    ref = input["ref"]   

    # Getting nextread-created directories
    secure_analysis_dir = kwargs["secure_analysis_dir"]
    secure_ref_dir = kwargs["secure_ref_dir"]

    # Creating command based on user input
    if os.path.isfile(ref): # Reference FASTA is provided
        nextclade_command = [
            "nextclade",
            "run",
            "-r",
            ref,
            "-t",
            os.path.join(secure_analysis_dir, 'analysis.tsv')
        ]
        if "gene" in input: # Gene annotation is provided
            gene = input["gene"]
            nextclade_command = nextclade_command + ["-m", gene]
        nextclade_command = nextclade_command + [seq]
        gene_path = gene
    else: # Name of Nextclade reference is provided
        secure_ref_dir = kwargs["secure_analysis_dir"]
        # Get dataset first and then run the command
        gene_path = os.path.join(secure_ref_dir, 'genome_annotation.gff3')
        nextclade_command = [
            [
                "nextclade",
                "dataset",
                "get",
                "-n",
                ref, 
                "-o",
                secure_ref_dir
            ],
            ["nextclade",
                "run",
                "-r",
                os.path.join(secure_ref_dir, 'reference.fasta'),
                "-m",
                gene_path,
                "-t", 
                os.path.join(secure_analysis_dir, 'analysis.tsv'),
                seq
            ]
        ]
    return nextclade_command, gene_path


def capture_output(command):
    """Runs Nextclade CLI and captures the output.

    Parameters
    ----------
    command : str
        The Nextclade CLI command.

    Returns
    -------
    pandas.DataFrame
        A DataFrame of the Nextclade analysis TSV output.

    """
    try:
        if len(command) != 2: # Reference FASTA is provided
            subprocess.run(command, check=True)
            analysis_arg_ind = command.index("-t") + 1
            analysis_file_path = command[analysis_arg_ind]
        else: # Name of Nextclade reference is provided
            analysis_arg_ind = command[1].index("-t") + 1
            analysis_file_path = command[1][analysis_arg_ind]
            for line in command:
                subprocess.run(line, check=True)

        # Reading Nextclade analysis TSV output
        analysis_dataframe = pd.read_csv(analysis_file_path, delimiter='\t')
        # Removing index column
        analysis_dataframe.drop('index', axis=1, inplace=True)
        return analysis_dataframe
    except subprocess.CalledProcessError as e:
        print(f"Error running Nextclade: {e}")
    except FileNotFoundError:
        print("Nextclade executable not found. Make sure it is included in the system $PATH.")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")