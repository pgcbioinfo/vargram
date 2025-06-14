"""Module to fetch Nextclade data."""

import subprocess
import os

def run_nextclade(command,  write=False):
    """Run Nextclade CLI.
    
    Parameters
    ----------
    command : list 
        Directory where test data will be downloaded.
    write : bool, default:False
        Determines whether to print output of command.
    
    Returns
    -------
    None
    
    """
 
    try:
        result = subprocess.run(command, check=True, capture_output=True, text=True)
        if write:
            print(result.stdout)
    except subprocess.CalledProcessError as e:
        print(f"Error running Nextclade: {e}")
    except FileNotFoundError:
        print("Nextclade executable not found. Make sure it is included in the system $PATH.")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")

def list(full=False):
    """List available Nextclade datasets.

    Parameters
    ----------
    full : bool, default:False
        Determines if only names of datasets will be shown.
    
    Returns
    -------
    None

    """

    if full:
        list_command = "nextclade dataset list"
    else:
        list_command = "nextclade dataset list --only-names"
    run_nextclade(list_command.split(), write=True)

def get(id, dir='', version=''):
    """Download specified Nextclade dataset into target directory.
    
    Parameters
    ----------
    id : str 
        Name of the Nextclade dataset.
    dir : str
        Directory where dataset will be downloaded.
    
    Returns
    -------
    None
    
    """

    # Define target directory where data will be downloaded
    if isinstance(dir, str) and dir.strip(): # Use user-provided directory 
        target_dir = dir 
        if not os.path.isdir(target_dir):
            os.makedirs(target_dir)
    else: # Or use current working directory by default
        target_dir = os.getcwd()
    
    get_command = f"nextclade dataset get --name {id} --output-dir {target_dir}"
    if version != '':
        get_command = get_command + f" --tag {version}"
    run_nextclade(get_command.split())