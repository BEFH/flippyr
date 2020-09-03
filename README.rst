Flippyr for Python
==================
**Note: Version 0.4.0** Previous versions' log files swapped the number of invalid and ambiguous variants. This did not effect the output genotypes or table, but it has been fixed.

This package is designed to align a PLINK_ fileset with a FASTA_ reference genome. It identifies the following issues:

.. _FASTA: https://en.wikipedia.org/wiki/FASTA_format
.. _PLINK: https://www.cog-genomics.org/plink2

- **Strand flipping:** *Flippyr* identifies variants on the reverse strand in relation to the reference genome.
- **Allele reversal:** PLINK operates on the minor (A1) and major (A2) alleles in the fileset unless the ``--keep-allele-order`` or ``--real-ref-alleles`` switches are used. *Flippyr* identifies variants where A2 does not match the reference allele or its complement, but A1 does.
- **Ambiguous alleles:** *Flippyr* identifies alleles where the which are complements. It is impossible to know whether these are correct, strand flipped, or allele reversed, so they are often removed. These include ``AT``, ``TA``, ``CG`` and ``GC``.
- **Indels:** *Flippyr* identifies insertions and deletions for removal from the fileset. There is no consistant format for representation of indels in plink filesets, so *Flippyr* does not check flipping and reversal on indels and instead marks them for deletion.
- **Unmatched alleles:** *Flippyr* marks loci where the reference or reference complement do not match either of the alleles in the PLINK fileset. This could be due to an incorrectly called variant, a missing allele or a chromosome not in the reference genome. PLINK allows chromosomes above 22 which include mitochondrial and various sex chromosome regions; unless the reference FASTA includes those, they will be marked as missing.
- **Multiallelic sites:** *Flippyr* identifies loci with more than one variant and marks them for deletion. Please use output files if this is not the desired behavior.

Use cases:
----------

There are several common use cases for *Flippyr*:

- Pre-imputation quality control: Imputation software and services like the Michigan Imputation Server require all variants to be aligned to the forward reference allele and disallow indels and multiallelic sites. This software prepares your PLINK files for imputation. You will also need to split your sample into individual chromosomes and create sorted Gzipped Variant Call format (``.vcf.gz``) files.
- Harmonization: If you are combining multiple data sets for joint analysis, you may wish to harmonize them to the same reference.

Dependancies:
-------------
*Flippyr* is compatible with Python version 3 and above, and requires ``pyfaidx`` and ``pandas``. These dependancies are automatically installed when using ``pip`` or ``easy_install``.

In order for *Flippyr* to actually repair the fileset, PLINK_ must be installed and on your ``PATH``.

Usage:
------
Command line:

.. code::

  usage: flippyr.py [-h] [-s] [-p] [-o OUTPUTPREFIX]
                    [--outputSuffix OUTPUTSUFFIX] [-m] [-i]
                    fasta bim

  A simple python script to search for allele switches, strand flips,
  multiallelic sites, ambiguous sites, and indels. The output is in the form of
  a .bim-like table and a log file.

  positional arguments:
    fasta                 Fasta file containing all chromosomes in the plink
                          fileset
    bim                   .bim file from binary plink fileset.

  optional arguments:
    -h, --help            show this help message and exit
    -s, --silent          Supress output to stdout.
    -p, --plink           Run the plink command.
    -o OUTPUTPREFIX, --outputPrefix OUTPUTPREFIX
                          Change output file prefix.
    --outputSuffix OUTPUTSUFFIX
                          Change output file suffix for plink file.
    -m, --keepMultiallelic
                          Do not delete multiallelic sites.
    -i, --keepIndels      Do not delete insertions/deletions.

Python:

.. code:: python

 import flippyr
  '''Run flippyr and get bim-like dataframe of results, as well as log.
  The inputs are the name of the FASTA file, the name of the PLINK bim file,
  and a boolian True if you do not want output to stdout.'''
  bim, log = flippyr.run(fasta,bim,silent=False)
  '''Alternatively, get the same outputs as with flippyr.run, but also save
  the same files as when running on the command line and output the PLINK
  command. There is an extra boolean argument that runs PLINK.'''
  bim, log, runPlink = flippyr.run(fasta,bim,plink=False,silent=False)

Only the ``.bim`` file is required for most operations, but the entire fileset must be in the directory to run PLINK.

The FASTA files must have the chromosome number (1,...,22,X,Y,M) at the beginning of the id lines. The number can be prefixed with "chr" and there can be more information following a space.
