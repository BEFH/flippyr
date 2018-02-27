
"""A simple python script tosearch for allele switches, strand flips,
multiallelic sites, ambiguous sites, and indels. The output is in the form of
a .bim-like table and a log file.

See:
https://github.com/BEFH/flippyr
"""

from setuptools import setup, find_packages

setup(
    name='flippyr',
    version='0.2.0',
    description='Find reference mismatches in PLINK filesets and fix them.',
    long_description=('A simple python script tosearch for allele switches,\n'
                      'strand flips, multiallelic sites, ambiguous sites,\n'
                      'and indels. The output is in the form of a .bim-like\n'
                      'table and a log file.'),
    url='https://github.com/BEFH/flippyr',
    author='Brian Fulton-Howard',
    author_email='brian.fulton-howard@mssm.edu',
    license='MIT',
    keywords='bim PLINK FASTA genomics qc',
    py_modules=['flippyr'],
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Operating System :: POSIX :: Linux',
        'Operating System :: MacOS :: MacOS X',
        'Operating System :: Microsoft',
    ],
    install_requires=['pyfaidx','pandas'],
    entry_points={
        'console_scripts': ['flippyr=flippyr:main'],
    },
)
