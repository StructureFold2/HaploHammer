#!/usr/bin/env python
#David Tack

'''

'''

#Imports
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
import glob
import string
import itertools
import collections
import argparse


#Functions
def read_in_fasta(afasta):
    '''Reads in a fasta file to a dictionary'''
    fasta_dict = {}
    fasta_sequences,fasta_dict = SeqIO.parse(open(afasta),'fasta'),{}
    for fasta in fasta_sequences:
        fasta_dict[fasta.id] = str(fasta.seq)
    return fasta_dict

def generate_names():
    '''Generates arbitary tags'''
    A = list(string.ascii_uppercase)
    B = [''.join(pair) for pair in itertools.product(string.ascii_uppercase, repeat=2)]
    C = [''.join(pair) for pair in itertools.product(string.ascii_uppercase, repeat=3)]
    return A+B+C

def generate_haps(sort_freq=False):
    '''Generate haplotype information'''
    hap_seqs,hap_genos,hap_counts = {},{},{}
    for fasta in sorted(glob.glob('*.fa')):
        fasta_name,fasta_as_dict = fasta.strip('.fa'),read_in_fasta(fasta)
        for seq_name,sequence in fasta_as_dict.items():
            d_key = hash(sequence)
            hap_seqs.setdefault(seq_name, {})[d_key]= sequence
            hap_genos.setdefault(fasta_name,{})[seq_name]=d_key
            hap_counts.setdefault(seq_name,collections.Counter())[d_key]+=1
    
    #Make conversion dictionary
    conversion,generic_names = {},generate_names()
    
    if sort_freq:
        for k, v in hap_counts.items():
            gene,allele_codes = k, v.most_common()
            for i in range(0,len(allele_codes)):
                conversion.setdefault(gene,{})[allele_codes[i][0]]=generic_names[i]
    
    else:
        for k, v in hap_seqs.items():
            gene,allele_codes = k, sorted(v.keys())
            for i in range(0,len(allele_codes)):
                conversion.setdefault(gene,{})[allele_codes[i]]=generic_names[i]

    #Return components
    return hap_seqs,hap_genos,conversion

def write_nested_alleles_as_fasta(info,tl_dict,outfyle,LW=80):
    '''Writes out the <.fasta> file of all alleles, names appropriately'''
    with open(outfyle,'w') as g:
        for gene, alleles in info.items():
            for allele,sequence in alleles.items():
                name = '_'.join([gene,tl_dict[gene][allele]])
                g.write('>' + name + '\n')
                for i in xrange(0,len(sequence),LW):
                    g.write(sequence[i:i+LW] + '\n') 

def write_csv_grid(info,tl_dict,outfyle='out.csv'):
    '''writes as csv matrix'''
    common_keys = sorted(list(set([item for sublist in [q.keys() for q in info.values()] for item in sublist])))
    with open(outfyle,'w') as g:
        header = ','.join(['line']+common_keys)+'\n'
        g.write(header)
        for k, v in sorted(info.items()):
            line = ','.join([k]+[tl_dict[herp][v[herp]] for herp in common_keys])+'\n'
            g.write(line)

def check_extension(astring,extension):
    '''Checks and fixes things to have the proper extension'''
    out_string = astring if astring.endswith(extension) else astring + extension
    return out_string

#Main Function
def main():
    parser = argparse.ArgumentParser(description='Assembles a pan-transcriptome <.fasta> and matrix')
    parser.add_argument('-sort',action='store_true',default=False,help='Alleles are named by frequency')
    parser.add_argument('-fasta',default='out.fasta', help='Specify output <.fasta> name')
    parser.add_argument('-matrix',default='out.csv', help='Specify output <.csv> name')
    args = parser.parse_args()
    
    #Call, collect, organize haplotypes
    seqs, genos, conver = generate_haps(args.sort)
    
    #Write out all alleles
    write_nested_alleles_as_fasta(seqs,conver,check_extension(args.fasta,'.fasta'))
    
    #Write out allele matrix
    write_csv_grid(genos,conver,check_extension(args.matrix,'.csv'))


if __name__ == '__main__': 
    main()
