# Script to convert SNP VCF to Zarr database

import argparse
import sys
import os
from itertools import repeat

import pandas as pd
import multiprocessing as mp

import zarr
import numcodecs
import allel

def args():

    parser = argparse.ArgumentParser(
            description = 'Converts a VCF to an on-disk Zarr database',
            usage = 'python 3.7 vcf_to_zarr.py [options]')
    parser.add_argument('-i', '--vcf_in', required=True, type=str,
            help='Path to input VCF on disk')
    parser.add_argument('-t', '--tabix-exec', required=True,type=str,
            help='Full path to tabix executable')
    parser.add_argument('-c', '--chrom_path', required=True, type=str,
            help='Path to text file with chromosome names (one per line)')

    args = parser.parse_args()

    return args.vcf_in, args.tabix_exec, args.chrom_path


def chroms_to_list(chromosome_file):
    """Converts text file with chromosome names to list

    Args:
        chromosomes (str): Path to text file with chromsomes

    Returns: 
        chrom_list (:obj:`list` of :obj:`str`): List with
        chromosomes names as elements.
    """

    chrom_list = pd.read_table(chromosome_file, header=None).iloc[:,0].tolist()
    
    return chrom_list


def vcf_to_zarr(vcf_in, tabix_exec, chrom):
    """Convert on-disk VCF to on-disk Zarr database using
    scikit-allele and zarr modules

    Zarr database written to same directory as input VCF
    
    Args:
        vcf_in (str): Path to input VCF on disk
        tabix_exec (str): Full path to tabix executable
        chrom (str): Chromosome for which Zarr database should be created

    Returns:
        None
    """
    vcf_path = os.path.dirname(vcf_in)

    # allel.vcf_to_zarr returns a directory with Zarr databse
    # Set Zarr database outdir
    zarr_base = os.path.basename(vcf_in).split('.')[0]
    zarr_out = vcf_path + '/' + zarr_base + '.zarr'

    # Rename 'numalt' field. Required by Zarr to distinguish `NUMALT` from `numalt`
    # `numalt` is automatically computed by scikit-allel
    rename_dict = {'variants/numalt':'variants/numalt_sci'}

    # Use vcf_to_zarr function from scikit-allel to create zarr database
    # Currently optimized for biallelic SNP VCF but easy to extend functionality
    allel.vcf_to_zarr(
            input=vcf_in,
            output=zarr_out,
            overwrite=True,
            group=chrom,
            rename_fields=rename_dict,
            fields='*',
            alt_number=1,
            tabix=tabix_exec,
            region=chrom,
            compressor=numcodecs.Blosc(cname='zstd', clevel=1, shuffle=False)
            )


def main():

    vcf_in, tabix_exec, chrom_path = args()
    chrom_list = chroms_to_list(chrom_path)

    iterator = zip(repeat(vcf_in), repeat(tabix_exec), chrom_list)

    with mp.Pool(processes = len(chrom_list)) as pool:
        pool.starmap(vcf_to_zarr, iterator)
   

if __name__ == '__main__':
    main()
    

