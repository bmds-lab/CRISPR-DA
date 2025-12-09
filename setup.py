import os
import shutil
import pathlib
import setuptools
from setuptools.command.install import install

class MyInstall(install):
    def run(self):
        cwd = pathlib.Path().absolute()
        ISSL_dir = cwd / 'crisprda/ISSL'
        build_dir = ISSL_dir / 'build'
        build_dir.mkdir(parents=True, exist_ok=True)
        os.chdir(str(build_dir))
        self.spawn(['cmake', str(ISSL_dir)])
        self.spawn(['cmake', '--build', '.'])
        os.chdir(str(cwd))
        super().run()
        os.chdir(str(build_dir))
        self.spawn(['cmake', '--install', '.', '--prefix', str(cwd / 'crisprda/resources')])
        os.chdir(str(cwd))
        shutil.rmtree(str(build_dir))

with open('README.md', 'r', encoding='utf-8') as fh:
    long_description = fh.read()

setuptools.setup(
    name='crispr-da',
    version='0.0.1',
    author='Jake Bradford, Dimitri Perrin',
    author_email='dimitri.perrin qut edu au',
    description='For designing CRISPR-Cas9 sgRNA when many viral genomes are to be considered.',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/bmds-lab/crispr-da',
    project_urls={
        'Bug Tracker': 'https://github.com/bmds-lab/crispr-da/issues',
        'Lab website': 'http://biomedicaldatascience.com/'
    },
    classifiers=[
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)'
        'Operating System :: OS Independent',
        'Programming Language :: Python :: 3',
    ],
    packages=setuptools.find_packages(where='crisprda', exclude=['resources']),
    cmdclass={
        'install': MyInstall,
    },
    python_requires='>=3.12',
    entry_points = {
        'console_scripts': [
            'crisprda=crisprda.utils.cli:main'
        ],
    },
    include_package_data=True,
    package_data={'': ['data/*']},
    install_requires=[
       'ncbi-datasets-pylib>=13.22.0',
       'pandas>=1.4.1',
       'ete3>=3.1.3',
       'requests>=2.32.5'
    ]
)