#!/usr/bin/env python3
"""usage: flippyr.py [-h] [-s] [-p] [-o OUTPUTPREFIX]
                  [--outputSuffix OUTPUTSUFFIX] [-m] [-i]
                  [--plinkMem PLINKMEM] fasta bim

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
  --plinkMem            Set the memory limit for plink."""

import sys
import argparse
import os
import re

import pandas as pd
from pyfaidx import Fasta


def get_ref(chrom, pos, fasta):
    '''Read in fasta file'''
    def chr_replace(chrom, orig, new): # map PLINK chrom to fasta chrom
        chrom_nochr = [re.sub("chr", "", x) for x in chrom]
        return [new if x in orig else x for x in chrom_nochr]
    complements = {"A": "T", "T": "A", "C": "G", "G": "C"}
    fa_ = Fasta(fasta, sequence_always_upper=True,
                as_raw=True, read_ahead=900000)
    fa_.records = {re.sub("chr", "", x): y for x, y in
                   zip(fa_.records.keys(), fa_.records.values())}
    # steps to keep sex chromosomes and mitochondria by mapping sequences #
    mitochroms = ['26', 'M', 'MT', '0M']
    if 'M' in fa_.records:
        chrom = chr_replace(chrom, mitochroms, 'M')
    elif 'MT' in fa_.records:
        chrom = chr_replace(chrom, mitochroms, 'MT')
    chrom = chr_replace(chrom, ['23', '25', 'XY'], 'X')
    chrom = chr_replace(chrom, ['24'], 'Y')
    # done mapping sequences #
    infa = set(fa_.records.keys())  # chr in fasta file.
    cond = [chr_ in infa for chr_ in chrom]  # test if chr in fasta
    ref = [fa_[x][int(y) - 1] if z else "N"
           for x, y, z in zip(chrom, pos, cond)]
    comp = [complements.get(x) for x in ref]  # get complements of ref allele
    fa_.close()  # Close fasta file
    return ref, comp, infa


def build_table(fasta, bim):
    '''Parse bim and fasta'''
    df = pd.read_csv(bim, sep="\t", header=None, usecols=[0, 1, 3, 4, 5],
                     names=["chr", "ID", "position", "minor", "major"],
                     dtype={"chr": str, "ID": str, "position": int,
                            "minor": str, "major": str}, engine="c")
    inbim = set(df.chr)  # get set of chr in bim file.
    df['ref'], df['complement'], infa = get_ref(
        df.chr.values, df.position.values, fasta)
    if not inbim.issubset(infa):
        log = "".join(["\033[1;31mWarning:\033[0;31m Fasta file does not ",
                       "contain all chromosomes. Variants\non any chromosome ",
                       "not in the fasta file will be marked as missing.",
                       "\033[0m"])
    else:
        log = ""
    return df, log


def bad_alt(a1a2):
    '''Look for any bad alleles.'''
    valid = ['A', 'T', 'C', 'G']
    if len(a1a2[0]) > 1:
        if len(a1a2[1]) > 1:
            return True
        else:
            return a1a2[1] not in valid
    return any([a not in valid for a in a1a2])


def test_allele(major, minor, ref, complement):
    '''conditions for flipping, reversal, etc.'''
    if bad_alt([major, minor]):
        result = (1, "invalid")
    elif major + minor in ["AT", "TA", "GC", "CG"]:
        result = (2, "ambiguous")
    elif major == ref:
        result = (0, "match")
    elif major == complement:
        result = (3, "strand")
    elif minor == ref:
        result = (4, "allele")
    elif minor == complement:
        result = (5, "strand and allele")
    else:
        result = (6, "no match")
    return result


def test(df):
    '''function to test ref and alt'''
    df["outcome"], df["explanation"] = zip(*[test_allele(w, x, y, z)
        for w, x, y, z in zip(df.major.values, df.minor.values,
        df.ref.values, df.complement.values)])
    df["indel"] = [len(x) != 2 for x in df[['minor', 'major']].sum(axis=1)]
    df["multiallelic"] = df[["chr", "position"]].duplicated(keep=False)
    counts = [0 if v is None else v for v in map(
        df.outcome.value_counts().get, range(7))]
    multi_sum = sum(df.multiallelic)
    indel_sum = sum(df.indel)
    counts[0] -= multi_sum + indel_sum
    counts.append(multi_sum)
    counts.append(indel_sum)
    counts.insert(0, df.shape[0])

    log = ["\033[1mThere are the following sites:",
           "\033[1;32m[{}]\033[0m total",
           "\033[1;32m[{}]\033[0m correct biallelic SNP",
           "\033[1;32m[{}]\033[0m invalid alternate",
           "\033[1;32m[{}]\033[0m ambiguous",
           "\033[1;32m[{}]\033[0m strand flipped",
           "\033[1;32m[{}]\033[0m allele switched",
           "\033[1;32m[{}]\033[0m strand flipped & allele switched",
           "\033[1;32m[{}]\033[0m unmatched",
           "",
           "\033[1;32m[{}]\033[0m multiallelic",
           "\033[1;32m[{}]\033[0m insertion/deletion"]
    log = "\n".join(log).format(*counts)
    return df, log


class output():
    "Print output to the command line."

    def __init__(self, silent=False, sep="\n"):
        self.windows = (os.name == "nt")
        self.silent = silent
        self.sep = sep
        self.log = []

    def __len__(self):
        return len(self.log)

    def __getitem__(self, item):
        return self.log[item]

    def __repr__(self):
        return repr(self.log)

    def new(self, text, stderr=False, postfix="", otr=False):
        "print text and add to log"
        stripped = self.strip(text)
        if not self.silent:
            if self.windows:
                text = stripped
            if stderr:
                print(text, file=sys.stderr)
            else:
                print(text)
        if not otr:
            self.log.append(stripped + postfix)

    @staticmethod
    def strip(text):
        "remove ANSI escapes"
        return re.sub("\\033[\[\;0-9]*m", "", text)

    def write(self, fname):
        '''write log to file'''
        log = self.sep.join(self.log)
        with open(fname, "w") as f:
            f.write(log + "\n")


def run(fasta, bim, silent=False):
    log = output(silent=silent)
    log.new("\033[1mLoading files...\033[0m", otr=True)
    bim, out = build_table(fasta, bim)
    if not out:
        log.new("\n\033[1mFinding misalignments...\033[0m\n", otr=True)
    bim, out = test(bim)
    log.new(out)
    return bim, log


def writeFiles(fasta, bim, outname, plink=False, silent=False,
               p_suff="_flipped", multi=False, indel=False,
               mem="auto"):
    # Initialize plink command
    runPlink = ("plink -bfile {a} --make-bed --out {b} "
                "--real-ref-alleles").format(
                                             a=re.sub("\.bim", "", bim),
                                             b=outname + p_suff)

    bim, log = run(fasta, bim, silent)
    bim.to_csv(outname + ".log.tab", sep="\t", index=False)

    # Write file with ids to delete:
    dels = [state in [1, 2, 6] for state in bim.outcome]
    # [invalid, ambiguous, or no match]
    if not indel:  # not -i or --keepIndels
        dels = [d or indel_ for d, indel_ in zip(dels, bim.indel)]
    if not multi:  # not -m or --keepMultiallelic
        dels = [d or multi_ for d, multi_ in zip(dels, bim.multiallelic)]

    fname = outname + ".delete"
    if any(dels):
        bim[dels]["ID"].to_csv(fname, sep="\t", index=False, header=False)
        runPlink += " --exclude {}".format(fname)
        bim.drop(bim[dels].index, inplace=True)  # drop deleted vars from df
    else:
        open(fname, 'a').close()
    dels = None

    # Write file with ids to flip:
    flips = [i in [3, 5] for i in bim.outcome]
    fname = outname + ".flip"
    if any(flips):
        bim[flips]["ID"].to_csv(fname, sep="\t", index=False, header=False)
        runPlink += " --flip {}".format(fname)
    else:
        open(fname, 'a').close()
    flips = None

    # Write file with ref alleles and IDs to allele correct:
    allele = [i in [4, 5] for i in bim.outcome]
    fname = outname + ".allele"
    if any(allele):
        bim[allele][["ref", "ID"]].to_csv(fname, sep="\t",
                                          index=False, header=False)
        runPlink += " --a2-allele {} 1 2".format(fname)
        if mem != "auto":
            runPlink += " --memory {}".format(mem)
    else:
        open(fname, 'a').close()
    flips = None

    if plink:
        import subprocess
        output_ = subprocess.run(runPlink, shell=True,
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.STDOUT)
        log.new("\n\n" + output_.stdout.decode("utf-8"))
    else:
        log.new("\nRun {}.runPlink to remove ambiguous or\n".format(outname) +
                "unmatched sites, flip reverse strand sites and use " +
                "reference alleles.")
    if os.name != "nt":
        runPlink = "#!/usr/bin/env bash\n\n" + runPlink
    # Write log file.
    log.write("{}.log".format(outname))
    fname = "{}.runPlink".format(outname)
    with open(fname, "w") as f:
        f.write(runPlink)
    os.chmod(fname, 0o755)
    return bim, log, runPlink


def main():
    description = ['A simple python script to search for allele switches,',
                   'strand flips, multiallelic sites, ambiguous sites,',
                   'and indels. The output is in the form of a .bim-like',
                   'table and a log file.']
    parser = argparse.ArgumentParser(description=" ".join(description))
    parser.add_argument("fasta", help="""Fasta file containing all chromosomes
    in the plink fileset""")
    parser.add_argument("bim", help=".bim file from binary plink fileset.")
    parser.add_argument("-s", "--silent", action="store_true",
                        help="Supress output to stdout.")
    parser.add_argument("-p", "--plink", action="store_true",
                        help="Run the plink command.")
    parser.add_argument("--plinkMem", type=int, default=-9,
                        help="Set the memory limit for plink.")
    parser.add_argument("-o", "--outputPrefix", type=str, default="0",
                        help="Change output file prefix.")
    parser.add_argument("--outputSuffix", type=str, default="_flipped",
                        help="Change output file suffix for plink file.")
    parser.add_argument("-m", "--keepMultiallelic", action="store_true",
                        help="Do not delete multiallelic sites.")
    parser.add_argument("-i", "--keepIndels", action="store_true",
                        help="Do not delete insertions/deletions.")
    args = parser.parse_args()
    if args.outputPrefix == "0":
        outname = os.path.splitext(args.bim)[0]
    else:
        outname = args.outputPrefix

    mem = args.plinkMem if args.plinkMem != -9 else "auto"

    writeFiles(args.fasta, args.bim, outname, plink=args.plink,
               silent=args.silent, p_suff=args.outputSuffix,
               multi=args.keepMultiallelic, indel=args.keepIndels,
               mem=mem)


if __name__ == "__main__":
    main()
