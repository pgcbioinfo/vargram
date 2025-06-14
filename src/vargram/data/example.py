"""Module to fetch example data."""

from vargram import __version__

import os
import requests
import zipfile
import io

def get_specific_tag(username, reponame, tag):
    """Download a given release of the VARGRAM repo. 

    Parameters
    ----------
    username : str
        The owner of the VARGRAM repository.
    reponame : str
        The name of the VARGRAM repository.
    tag : str
        The specific release of VARGRAM.

    Returns
    -------
    requests.Response
        The Response object returned by a GET request to the URL of the release.
        
    """
    
    url = f'https://github.com/{username}/{reponame}/archive/refs/tags/{tag}.zip'
    response = requests.get(url)
    return response

def get_tag_names(username, reponame):
    """Gets all release tags of VARGRAM repository.

    Parameters
    ----------
    username : str
        The GitHub username or owner of the VARGRAM repository.
    reponame : str
        The name of the VARGRAM repository.

    Returns
    -------
    list
        List of all release tags.

    """

    url = f'https://api.github.com/repos/{username}/{reponame}/tags'
    response = requests.get(url)
    response.raise_for_status()
    return [tag['name'] for tag in response.json()]

def get(dir=''):
    """Downloads test data into target directory.
    
    Parameters
    ----------
    dir : str 
        Directory where test data will be downloaded.

    Returns
    -------
    None

    Raises
    ------
    requests.HTTPError
        If the GitHub request fails.
    """
    username = 'pgcbioinfo'
    reponame = 'vargram'

    # Define target directory where data will be downloaded
    if dir != '' and isinstance(dir, str): # Use user-provided directory 
        target_dir = dir 
        if not os.path.isdir(target_dir):
            os.makedirs(target_dir)
    else: # Or use current working directory by default
        target_dir = os.getcwd()

    # Get zipped data from VARGRAM repo based on installed package version
    downloaded_release = f"v{__version__}"
    all_tags = get_tag_names(username, reponame)
    if downloaded_release in all_tags: # i.e. the current version has not yet been pushed to the public repo
        all_tags.remove(downloaded_release)
    for tag in [downloaded_release] + all_tags:
        try:
            response = get_specific_tag(username, reponame, tag)
            response.raise_for_status()
            break
        except Exception:
            continue
        
    # Unzip repo
    with zipfile.ZipFile(io.BytesIO(response.content)) as zipped_repo:

        # Go through each archive member (i.e. filepaths) in zipfile
        # Skip members that aren't inside the test data directory (source_dir)
        source_dir = os.path.join(f"vargram-{tag[1:]}/","tests/test_data/")
        for path in zipped_repo.namelist():
            if path.startswith(source_dir) and  path != source_dir: 

                # Define target file path of member
                int_path = path[len(source_dir):] 
                target_path = os.path.join(target_dir, int_path)

                # Create member subdirectories
                if path.endswith('/'): 
                    if not os.path.exists(target_path):
                        os.mkdir(target_path)
                    continue
                
                # And extract files to target directory
                with zipped_repo.open(path) as source, open(target_path, "wb") as target:
                    target.write(source.read())
    
    if tag[1:] != __version__:
        print(f"Successfully downloaded vargram-{tag} (installed version: v{__version__}) test data to {target_dir}.")