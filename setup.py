import setuptools

with open('README.md', 'r', encoding='utf-8') as fh:
    long_description = fh.read()

setuptools.setup(
    name='viralcut',
    version='0.0.1',
    author='Jake Bradford, Dimitri Perrin',
    author_email='dimitri.perrin qut edu au',
    description='For designing CRISPR-Cas9 sgRNA when many viral genomes are to be considered.',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/bmds-lab/ViralCut',
    project_urls={
        'Bug Tracker': 'https://github.com/bmds-lab/ViralCut/issues',
        'Lab website': 'http://biomedicaldatascience.com/'
    },
    classifiers=[
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)'
        'Operating System :: OS Independent',
        'Programming Language :: Python :: 3',
    ],
    packages=setuptools.find_packages(where='viralcut'),
    python_requires='>=3.8',
    entry_points = {
        'console_scripts': [
            'ViralCut=viralcut.utils.cli:main'
        ],
    },
    include_package_data=True,
    package_data={'': ['data/*']},
    install_requires=[
       'ncbi-datasets-pylib>=13.22.0',
       'pandas>=1.4.1'
    ]
)