#!/usr/bin/env python

'''
Converts a letter-based haplotype matrix to classic GWAS format.
'''

#Imports
import collections
import argparse

#Functions
def read_in_haplo_matrix(afile):
    ''''''
    matrix = {}
    with open(afile,'r') as f:
        for line in f:
            bits = line.strip().split(',')
            if bits[0] == 'line':
                locus_names = bits[1:]
            else:
                accession,genotypes = bits[0],bits[1:]
                for n, g in zip(locus_names,genotypes):
                    matrix.setdefault(n, {})[accession] = g
    return matrix

def data_compressor(nested_dictionary):
    ''''''
    spin_to_win = {}
    for name, sub_dict in nested_dictionary.items():
        turnover = {}
        for accession, allele in sub_dict.items():
            if allele in turnover:
                turnover[allele].append(accession)
            else:
                turnover[allele] = [accession]
        spin_to_win[name] = turnover
    return spin_to_win
            
def filter_by_freq(data,percentage):
    ''''''
    END = {}
    for name, sub_dict in data.items():
        
        #Generate Frequencies,master accession list
        sub_freqs,all_accessions = collections.Counter(),[]
        for allele, accessions in sub_dict.items():
            sub_freqs[allele] = len(accessions)
            all_accessions.extend(accessions)
        
        #List of Alleles past frequency threshold
        keep_freqs = []
        for x, b in sub_freqs.items():
            if float(b)/sum(sub_freqs.values()) > percentage:
                keep_freqs.append(x)
        
        for item in keep_freqs:
            ex_name = '_'.join([name,item])
            onez = sub_dict[item]
            zeroz = list(set(all_accessions).difference(set(onez)))
            for line_accession in onez:
                END.setdefault(line_accession, {})[ex_name] = 1
            for line_accession in zeroz:
                END.setdefault(line_accession, {})[ex_name] = 0      
    return END


def write_out_hapwas_matrix(information,outfyle):
    ''''''
    common_keys = sorted(list(set([item for sublist in [q.keys() for q in information.values()] for item in sublist])))
    with open(outfyle,'w') as g:
        header = ','.join(['line']+common_keys)+'\n'
        g.write(header)
        for k, v in sorted(information.items()):
            line = ','.join([k]+[str(v[dkey]) for dkey in common_keys])+'\n'
            g.write(line)

def check_extension(astring,extension):
    '''Checks and fixes things to have the proper extension'''
    out_string = astring if astring.endswith(extension) else astring + extension
    return out_string

def main():
    parser = argparse.ArgumentParser(description='Converts a letter-based haplotype matrix to classic GWAS format.')
    parser.add_argument('csv', type=str, help='Haplotype Matrix file to reformat')
    parser.add_argument('-freq',type=float, default=0.05, help='[default = 0.05] alelle frequency threshold for retention')
    parser.add_argument('-name',type=str,default=None, help='Change the name of the outfile, overrides default')
    args = parser.parse_args()
    
    #Generate out name
    out_name = check_extension(args.name,'.csv') if args.name else check_extension(args.csv,'.csv').replace('.csv','_gwasmatrix.csv')
    
    raw_data = read_in_haplo_matrix(args.csv)
    processed_data = data_compressor(raw_data)
    end_data = filter_by_freq(processed_data,args.freq)
    write_out_hapwas_matrix(end_data,out_name)
    

if __name__ == '__main__': 
    main()
