#!/usr/bin/env python
from pysam import VariantFile
import argparse

parser = argparse.ArgumentParser(description='Script to filter nanopore derived VCF files')

parser.add_argument('-i',
                    '--infile',
                    help='Input vcffile',
                    type=str,
                    required=True)


parser.add_argument('-o',
                    '--outfile',
                    #metavar='File',
                    help='Filtered vcffile',
                    type=str,
                    required=True)

parser.add_argument('-af',
                    '--minAF',
                    help='Minimal allele frequency to output',
                    default=0.25,
                    type=float,
                    required = False)

args = parser.parse_args()

if __name__ == '__main__':
    vcffile = VariantFile(args.infile)

    with open(args.outfile, 'w') as outfile:
        filtered_variant_dict = {}
        for rec in vcffile:
            #Skip minor allele variants
            if rec.info['MAJOR'] == 0:
                continue
            #Skip positions with an allele frequency below threshold
            if rec.info['AF'] < args.minAF:
                continue
            #Ignore out-of-frame indels
            if not (len(rec.alleles[0])-len(rec.alleles[1])) % 3 == 0:
                continue

            #Add record when it passes all filters
            filtered_variant_dict[rec.pos] = rec

        print(vcffile.header, file=outfile, end='')

        #Iterate over the complete genome (untill the last variant) to be able to skip variants which are within deletion gaps
        pos = 1
        if len(filtered_variant_dict) > 0:
            end = max(filtered_variant_dict)
            while pos <= end:
                if pos in filtered_variant_dict:
                    print(filtered_variant_dict[pos], file=outfile, end='')
                    #If ref is bigger than alt, we have a deletion and increment pos by the deletion size to skip any variants inside the accepted deletion
                    ref_len=len(filtered_variant_dict[pos].alleles[0])
                    alt_len=len(filtered_variant_dict[pos].alleles[1])
                    if ref_len > alt_len:
                        pos += (ref_len-alt_len)
                #Increment pos
                pos += 1