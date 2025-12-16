from importlib.resources import files, as_file
from configparser import ConfigParser
from subprocess import run, DEVNULL
from pathlib import Path
import shutil
import os

configFile = Path.home() / ".crispr-da-config.ini"
config = ConfigParser()
config.read(configFile)

def get_config(option):
    global config
    try:
        return config.get('Main', option)
    except:
        raise RuntimeError(f"Failed to get option '{option}' from config")

def secure_config_opener(path, flags):
    return os.open(path, flags, 0o600)

# TODO: Add more verbose print statements
def run_config(force=False, default=False):
    global config, configFile
    # Build ISSL bin
    with as_file(files('crispr_da').joinpath('resources')) as fp:
        resource_dir = fp
    createIndexBin = resource_dir / 'ISSLCreateIndex'
    getOfftargetsBin = resource_dir / 'ISSLGetOfftargets'
    scoreOfftargetsBin = resource_dir / 'ISSLScoreOfftargets'

    if (not (createIndexBin.exists() and getOfftargetsBin.exists() and scoreOfftargetsBin.exists())) or force:
        print("Bin files missing, running build")

        with as_file(files('crispr_da').joinpath('ISSL')) as fp:
            ISSL_dir = fp
        build_dir = ISSL_dir / 'build'
        build_dir.mkdir(parents=True, exist_ok=True)
        # TODO: Capture output for error logging
        run(f"cmake -B {build_dir} -S {ISSL_dir}", shell=True, check=True, stdout=DEVNULL, stderr=DEVNULL)
        run(f"cmake --build {build_dir}", shell=True, check=True, stdout=DEVNULL, stderr=DEVNULL)
        run(f"cmake --install {build_dir} --prefix {resource_dir}", shell=True, check=True, stdout=DEVNULL, stderr=DEVNULL)
        shutil.rmtree(str(build_dir))
        print("Sucessfully built")

    # Get config.ini
    if 'Main' not in config.sections():
        config.add_section('Main')
    if not default:
        # Set cache location
        print('Enter a new value or hit enter to use the value in the sqaure brackets')
        cacheLocation = input(f'Cache location [{config.get('Main', 'Cache', fallback='')}]: ')
        if len(cacheLocation) > 0:
            config.set('Main','Cache', cacheLocation)
        
        # Set bactch size
        batchSizeNCBI = input(f'NCBI API request batch size (Recommended - 100) [{config.get('Main', 'NCBIBatchSize', fallback='')}]: ')
        if len(batchSizeNCBI) > 0:
            config.set('Main','NCBIBatchSize', batchSizeNCBI)
    else:
        config.set('Main','Cache', str(Path.home() / 'CRISPR-DA-Cache'))
        config.set('Main','NCBIBatchSize', "100")

    # TODO: Move NCBI API key here (save to hidden file or keep as export)

    # Export update config
    oldMask = os.umask(0)
    with open(configFile, 'w', opener=secure_config_opener) as outFile:
        config.write(outFile)
    os.umask(oldMask)
    
    # Reload config file
    config = ConfigParser()
    config.read(configFile)