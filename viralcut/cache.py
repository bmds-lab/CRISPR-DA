import shutil
from pathlib import Path

from . import config

def add_entry(id: str):
    cache = Path(config.CACHE) / id
    if cache.exists():
        shutil.rmtree(cache)
    cache.mkdir()
    return cache

def remove_entry(id: str):
    cache = Path(config.CACHE) / id
    shutil.rmtree(cache)

def get_missing_entires(ids: list[str]):
    cache = Path(config.CACHE)
    cache_misses = []
    for id in ids:
        if not (cache / id).exists():
            cache_misses.append(id)
    return cache_misses

def get_file(id: str, fileExtension: str):
    cache = Path(config.CACHE) / id
    files = [x for x in cache.glob(f'*{fileExtension}')]
    if len(files) == 0:
        raise RuntimeError(f'Could not file ending with {fileExtension} in entry {id}')
    elif len(files) > 1:
        raise RuntimeError(f'Multiple files ending with {fileExtension} found in {id}\nFiles: {files}')
    return files[0]

def get_missing_files(ids: list[str], fileExtension: str):
    cache_misses = []
    for id in ids:
        try:
            _ = get_file(id, fileExtension)
        except:
            cache_misses.append(id)
    return cache_misses

def check_cache_integrity():
    print()