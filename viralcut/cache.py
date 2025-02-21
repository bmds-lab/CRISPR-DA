'''
cache.py

This file handles access to the local 'cache'.
The cache is simply a folder that contains the files downloaded from the NCBI Datasets v2 REST API (dataset.py)
'''

import shutil
from pathlib import Path

from . import config

def add_entry(id: str):
    '''
    This method creates a new folder in the cache.
    Not it will override an existing entry with the same id.

    Arguments: 
    id (str): The id of the folder to create. It should be unique

    Returns:
    None
    '''
    cache = Path(config.CACHE) / id
    if cache.exists():
        shutil.rmtree(cache)
    cache.mkdir()
    return cache

def remove_entry(id: str):
    '''
    This method remove a folder in the cache.

    Arguments: 
    id (str): The id of the folder to remove.

    Returns:
    None
    '''
    cache = Path(config.CACHE) / id
    shutil.rmtree(cache)

def get_missing_entires(ids: list[str]):
    '''
    This method checks the cache to see if any ids do not have cache folders.

    Arguments: 
    id (list[str]): The ids that will be checked for in the cache.

    Returns:
    missing (list[str]): The ids that do not have folders in the cache
    '''
    cache = Path(config.CACHE)
    cache_misses = []
    for id in ids:
        if not (cache / id).exists():
            cache_misses.append(id)
    return cache_misses

def get_file(id: str, fileExtension: str):
    '''
    This method will return the file ending in the provided 
    file extension from the folder matching the id provided.
    Useful for finding specific file types without knowing the full name.

    Arguments: 
    id (str): The folder id to check.
    fileExtension (str): The suffix of the file to retrieve.

    Returns:
    file (Path): The path object of the file
    '''
    cache = Path(config.CACHE) / id
    files = [x for x in cache.glob(f'*{fileExtension}')]
    if len(files) == 0:
        raise RuntimeError(f'Could not file ending with {fileExtension} in entry {id}')
    elif len(files) > 1:
        raise RuntimeError(f'Multiple files ending with {fileExtension} found in {id}\nFiles: {files}')
    return files[0]

def get_missing_files(ids: list[str], fileExtension: str):
    '''
    This method will check mulitple folders to see if a file 
    ending with the provided file extension exists.

    Arguments: 
    ids (list[str]): The ids of the folders to check.
    fileExtension (str): The suffix of the file to check for.

    Returns:
    missing (list[str]): the ids of the folders that do not contian the file ending with 'fileExtension'.
    '''
    cache_misses = []
    for id in ids:
        try:
            _ = get_file(id, fileExtension)
        except:
            cache_misses.append(id)
    return cache_misses