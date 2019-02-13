#!/usr/bin/env python
#David Tack

'''
Does an OK job of converting CDS based fasta to peptide based fasta.

'''

#Imports
import glob
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
import string
import argparse

#Functions


#######

def quik_trim(a_string):
    '''Trims to Triplet'''
    new_string = a_string
    while len(new_string) %3 != 0:
        new_string = new_string[:-1]
    return new_string

def translate_to_peptides_and_hack(aseq):
    '''Translates'''
    translation_table = {'TTT': 'F', 'TCT': 'S', 'TAT': 'Y', 'TGT': 'C', 'TTC': 'F', 
                        'TCC': 'S', 'TAC': 'Y', 'TGC': 'C', 'TTA': 'L', 'TCA': 'S', 'TAA': '*', 'TGA': '*', \
                        'TTG': 'L', 'TCG': 'S', 'TAG': '*', 'TGG': 'W', 'CTT': 'L', 'CCT': 'P', 'CAT': 'H', \
                        'CGT': 'R', 'CTC': 'L', 'CCC': 'P', 'CAC': 'H', 'CGC': 'R', 'CTA': 'L', 'CCA': 'P', \
                        'CAA': 'Q', 'CGA': 'R', 'CTG': 'L', 'CCG': 'P', 'CAG': 'Q', 'CGG': 'R', 'ATT': 'I', \
                        'ACT': 'T', 'AAT': 'N', 'AGT': 'S', 'ATC': 'I', 'ACC': 'T', 'AAC': 'N', 'AGC': 'S', \
                        'ATA': 'I', 'ACA': 'T', 'AAA': 'K', 'AGA': 'R', 'ATG': 'M', 'ACG': 'T', 'AAG': 'K', \
                        'AGG': 'R', 'GTT': 'V', 'GCT': 'A', 'GAT': 'D', 'GGT': 'G', 'GTC': 'V', 'GCC': 'A', \
                        'GAC': 'D', 'GGC': 'G', 'GTA': 'V', 'GCA': 'A', 'GAA': 'E', 'GGA': 'G', 'GTG': 'V', \
                        'GCG': 'A', 'GAG': 'E', 'GGG': 'G'}
    try:#Easy Mode - no illegal codons
        translated_seq = ''.join([translation_table[aseq[i:i+3]] for i in range(0,len(aseq),3)])
        truncated_seq = translated_seq.split('*')[0]+'*'#got it!
        return truncated_seq
    except KeyError:
        aa = ''
        for i in range(0,len(aseq),3):
            new_codon = translation_table[aseq[i:i+3]] if aseq[i:i+3] in translation_table else '?'
            aa = aa + new_codon
        truncated_seq = aa.split('*')[0]+'*'
        return aa

def regressive_translation(fasta_dict):
    ''''''
    new = {}
    for k, v in fasta_dict.items():
        base_sequence = v.translate(None,string.punctuation)
        true_sequence = base_sequence if len(base_sequence) % 3 == 0 else quik_trim(base_sequence)
        new[k] = translate_to_peptides_and_hack(true_sequence)
    return new

def write_fasta(fastadict,outfile,LW=80):
    '''Writes a <.fasta>'''
    with open(outfile, 'w') as out:
        for k,v in fastadict.items():
            out.write('>' + k + '\n')
            for i in xrange(0,len(v),LW):
                out.write(v[i:i+LW] + '\n')

def read_fasta(genome_fasta):
    '''Reads in a <.fasta>'''
    fasta_sequences,fasta_dict =SeqIO.parse(open(genome_fasta),'fasta'),{}
    for fasta in fasta_sequences:
        fasta_dict[fasta.id] = str(fasta.seq)
    return fasta_dict

def add_suffix(a_file_name,suffix):
    '''Suffixes an outfile simply'''
    components = a_file_name.split('.')
    base_name,extension = components[:-1],components[-1]
    return '.'.join([('_').join([('.').join(base_name),suffix]),extension])
    
def change_extension(a_file_name,extension):
    '''Changes a file extension'''
    components = a_file_name.split('.')
    return '.'.join(['.'.join(components[:-1]),extension])

#Main Function
def main():
    parser = argparse.ArgumentParser(description='Translates CDS <.fasta> files to peptide <.fasta> files')
    parser.add_argument('-no_truncate',action='store_false',default = True,help = 'Do not truncate peptide at premature stop codons')
    parser.add_argument('-single',default = None, help = '[default = current directory] Operate on this single file')
    parser.add_argument('-in_ext',default = 'fa', help = '[default = fa] <.fasta> file extension')
    parser.add_argument('-suffix',default = 'peptide', help = '[default = peptide] Suffix for out files')
    #parser.add_argument('-codon_table',default = None, help = 'non-standard codon table',dest='table')
    args = parser.parse_args()
    
    #Que of files to process
    target_fyles = [args.single] if args.single else glob.glob('*.'+args.in_ext.strip(string.punctuation))
    
    #Process files...
    for fyle in target_fyles:
        data = read_fasta(args.single)
        new_data = regressive_translation(data)
        new_name = add_suffix(args.single,args.suffix)
        write_fasta(new_data,new_name)

if __name__ == '__main__':
    main()
