#!/usr/bin/env python3
import pysam
from pysam import VariantFile
import sys
import argparse
from collections import Counter
import re
import multiprocessing
import itertools
from datetime import date

parser = argparse.ArgumentParser(
    description='VCF generation script from BAM file')

parser.add_argument('-b',
                    '--bam',
                    help="BAM file to turn into VCF",
                    type=str,
                    required=True)

parser.add_argument('-o',
                    '--out',
                    help='Output path for tab separated alignment result file',
                    type=str,
                    required=True)

parser.add_argument('-c',
                    '--cores',
                    help='Number of cores to use for processing',
                    default=1,
                    type=int,
                    required=False)

parser.add_argument('-d',
                    '--mindepth',
                    help='Minimal depth at which to not consider any alternative alleles',
                    default=10,
                    type=int,
                    required=False)

parser.add_argument('-af',
                    '--minAF',
                    help='Minimal allele frequency to output',
                    default=0.01,
                    type=float,
                    required=False)

parser.add_argument('-r',
                    '--reference',
                    help='Reference genome, must be the same as for BAM',
                    type=str,
                    required=True)

parser.add_argument('--debug',
                    help='Debug flag to print out all positions and variants',
                    action='store_true')

args = parser.parse_args()

insert_finder = re.compile("(.*)\+\d+(.*)")


def print_header(ref, outfile):
    print('##fileformat=VCFv4.0', file=outfile)
    today = date.today()
    print('##fileDate=' + today.strftime("%Y%m%d"), file=outfile)
    print('##source=' + ' '.join(sys.argv), file=outfile)
    print('##reference=' + ref, file=outfile)
    print('##contig=' + ref, file=outfile)
    print("""##INFO=<ID=DP,Number=1,Type=Integer,Description="Raw Depth">
##INFO=<ID=AF,Number=1,Type=Float,Description="Allele Frequency">
##INFO=<ID=DP4,Number=4,Type=Integer,Description="Counts for ref-forward bases, ref-reverse, alt-forward and alt-reverse bases">
##INFO=<ID=INDEL,Number=1,Type=Integer,Description="Indicates that the variant is an indel with 1 or not with 0">
##INFO=<ID=MAJOR,Number=1,Type=Integer,Description="Indicates that the variant is the major variant with 1 or not with 0">""",
          file=outfile)
    print('##FILTER=<ID=minaf' + str(
        args.minAF) + ',Description="Allele frequency below indicated minimum">',
          file=outfile)
    print('##FILTER=<ID=mindp' + str(
        args.mindepth) + ',Description="Total depth below indicated minimum">',
          file=outfile)
    print('#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO',
          file=outfile)


def work(start, stop):
    bamfile = pysam.AlignmentFile(args.bam, "rb")
    ref = bamfile.references[0]
    pileup = bamfile.pileup(contig=ref, start=start, stop=stop,
                            ignore_orphans=False, min_mapping_quality=0,
                            min_base_quality=0, truncate=True)

    return ([
        parse_column(p.reference_pos, p.get_query_sequences(add_indels=True),
                     p.get_num_aligned()) for p in pileup])


def parse_column(ref_pos, allele_list, num_aln):
    if num_aln < args.mindepth:
        return (None)

    FWD_allele = Counter()
    REV_allele = Counter()
    COMB_allele = Counter()

    ref_allele = reference_seq[ref_pos]

    # Initialize the reference allele to 0 in case it is not present in any of the reads
    FWD_allele[(ref_allele.upper(), ref_allele.upper())] = 0
    REV_allele[(ref_allele.upper(), ref_allele.upper())] = 0
    COMB_allele[(ref_allele.upper(), ref_allele.upper())] = 0

    def add_allele(ref, alt):
        if (ref.islower() | alt.islower()):
            REV_allele[(ref.upper(), alt.upper())] += 1
            COMB_allele[(ref.upper(), alt.upper())] += 1
        else:
            FWD_allele[(ref, alt)] += 1
            COMB_allele[(ref, alt)] += 1

    for var in allele_list:
        # Ignore positions that represent deletions
        if '*' in var:
            continue

        # - means next nucleotide is a deletion
        if '-' in var:
            # Get the nucleotides that were deleted
            if var.islower():
                var = reference_seq[ref_pos:(ref_pos + len(var) - 2)].lower()
            else:
                var = reference_seq[ref_pos:(ref_pos + len(var) - 2)]
            add_allele(var, ref_allele)
        # + means next nucleotide is a insertion
        elif '+' in var:
            var = ''.join(insert_finder.match(var).groups())
            add_allele(ref_allele, var)
        else:
            add_allele(ref_allele, var)

    # Sort alleles by counts
    alleles = sorted(COMB_allele, key=COMB_allele.get, reverse=True)
    counts = [str(FWD_allele[k]) + ',' + str(REV_allele[k]) for k in alleles]
    alt_freq = [(FWD_allele[k] + REV_allele[k]) / num_aln for k in alleles]

    # Find and save the index of the reference allele
    for n, allele in enumerate(alleles):
        if allele[0] == allele[1]:
            ref_index = n
            break

    return ([str(num_aln), ref_pos, alleles, counts, alt_freq, ref_index])


if __name__ == '__main__':
    bamfile = pysam.AlignmentFile(args.bam, "rb")
    ref = bamfile.references[0]
    reference_seq = pysam.Fastafile(args.reference).fetch(ref)

    # Make an array of start-stop intervals to parallelize processing
    genome_length = len(reference_seq)
    split = int(genome_length / args.cores)
    split_ranges = [[i * split + 1, (i + 1) * split] for i in range(args.cores)]
    # Adjust the last "stop" to be the genome length
    split_ranges[-1][1] = genome_length

    with multiprocessing.Pool(processes=args.cores) as p:
        resultlist = p.starmap(work, iter(split_ranges))

    with open(args.out, "w") as outfile:
        print_header(ref, outfile)
        for result in itertools.chain.from_iterable(resultlist):
            if result == None:
                continue

            num_aln, ref_pos, alleles, counts, alt_freq, ref_index = result

            for n, alt in enumerate(alleles):
                # Print every position and every variant if debugging is on
                if args.debug != True:  # not true because flag arg is None if not provided and True if provided
                    # Do not print reference allele as a variant
                    if n == ref_index:
                        continue

                    # Do not print minor alleles below minAF
                    if alt_freq[n] < args.minAF:
                        continue

                # Print an INDEL flag if the alt is an indel
                indel = ';INDEL=0' if (len(alt[0]) + len(
                    alt[1])) == 2 else ';INDEL=1'

                # Print a MAJOR flag if this is the major allele (if it is at the top position in the AF sorted alleles list)
                major = ';MAJOR=1' if n == 0 else ';MAJOR=0'

                # Print vcf lines, add 1 to the reference position to make the coordinates 1 based
                print(ref, ref_pos + 1, ".", alt[0], alt[1], ".", "PASS",
                      "DP=" + num_aln + ";AF=" + "{:.6f}".format(
                          alt_freq[n]) + ";DP4=" + counts[ref_index] + "," +
                      counts[n] + indel + major, sep='\t', file=outfile)