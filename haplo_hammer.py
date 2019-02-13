#!/usr/bin/env python
#David Tack

'''
Description.
'''

#Imports
import os
import glob
import copy
import bisect
import argparse
import annotation_parser
#
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC


#Classes
class Transcript_Feature(object):
    '''Serves as a general class for transcript features'''
    def __init__(self,transcript_name,chromosome,start,end,strand,feature_type,feature_number,sequence=''):
        self.transcript_name = transcript_name
        self.chromosome = chromosome
        self.start = int(start)
        self.end = int(end)
        self.strand = strand
        self.feature_type = feature_type
        self.feature_number = int(feature_number)
        self.feature_name = '~'.join([transcript_name,feature_type,str(feature_number)])
        self.sequence = sequence
        self.modifications = []
        self.index_string = ((self.start,self.end),self.feature_name)

    def add_sequence(self,sequence):
        self.sequence = sequence

    def add_modification(self,modification):
        self.modifications.append(modification)

    def apply_modifications(self):
        transcript_frame = list(self.sequence)
        for modification in self.modifications:
            position, reference, alternate = modification
            local_position = position - self.start
            #SNP
            if len(reference) ==1 and len(alternate) == 1:
                transcript_frame[local_position] = alternate
            #Insertion
            elif len(reference)==1 and len(alternate) > 1:
                transcript_frame[local_position] = alternate
            #Deletion
            elif len(reference) >1 and len(alternate) == 1:
                deletion_spots = [local_position+numb+1 for numb in range(0,len(reference[1:]))]
                for spot in deletion_spots:
                    try:
                        transcript_frame[spot] = ''
                    #It is possible that a deletion runs into another feature.
                    except IndexError:
                        continue
        self.modified_sequence = ''.join(transcript_frame)

class Transcript(object):
    '''Holds Transcript Features'''
    def __init__(self):
        self.features = []
    
    def add_feature(self,feature):
        self.features.append(feature)
        self.strand = feature.strand
        self.chromosome = feature.chromosome
        self.transcript_name = feature.transcript_name
        
    def order_by_coordinates(self):
        self.features.sort(key=lambda x: x.start)
    
    def order_by_features(self):
        self.features.sort(key=lambda x: (x.feature_number, x.feature_type))
    
    def restrict_types(self,region):
        self.features = [f for f in self.features if f.feature_type == region]
    
    def collect_modification(self,modification,target):
        for feature in self.features:
            if feature.feature_name == target:
                feature.add_modification(modification)
    
    def assemble_transcript(self):
        self.features.sort(key=lambda x: x.start)
        end_seq = ''.join([feature.modified_sequence for feature in self.features])
        out_seq = end_seq if self.strand == '+' else str(Seq(end_seq).reverse_complement())
        return out_seq

#Functions
def read_in_fasta(genome_fasta):
    '''Reads in a fasta to a dictionary'''
    fasta_dict = {}
    fasta_sequences,fasta_dict = SeqIO.parse(open(genome_fasta),'fasta'),{}
    for fasta in fasta_sequences:
        fasta_dict[fasta.id] = str(fasta.seq)
    return fasta_dict

def extract_sequence(fasta_dict,chromosome,start,end):
    '''Pulls a sequence from a chromosome'''
    try:
        return fasta_dict[chromosome][start-1:end]
    except KeyError:
        return ''

def generate_poly_site(Pos,Ref,Alt,Alt_freqs):
    '''distills polymorphic choices'''
    potentials = Alt.split(',')
    simple_subs = dict([(('A','G'),'R'),(('C','T'),'Y'),(('A','T'),'W'),
                        (('C','G'),'S'),(('A','C'),'M'),(('G','T'),'K'),
                        (('C','G','T'),'B'),(('A','C','T'),'H'),
                        (('A','G','T'),'D'),(('A','G','C'),'V'),
                        (('A','C','G','T'),'N')])
    
    if len(potentials) == 1:
        if map(len,[Alt,Ref]) == [1,1]:
            key = tuple(sorted([Ref,Alt]))
            replace = simple_subs[key] if key in simple_subs else 'X'
            return (int(Pos),Ref,replace)
        else:
            return None

    if len(potentials) == 2:
       if map(len,potentials+[Ref]) ==[1,1,1]:
           key = tuple(sorted(potentials))
           replace = simple_subs[key] if key in simple_subs else 'X'
           return (int(Pos),Ref,replace)
       else:
           if all([len(Ref) > len(X) for X in potentials]):
               dd = sorted(potentials,key = lambda d: len(d), reverse=True)[0]
               return (int(Pos),Ref,dd)

           if all([len(Ref) < len(X) for X in potentials]):
               ii = sorted(potentials,key = lambda i: len(i))[0]
               return (int(Pos),Ref,ii)
    else:
        return None

def create_bisectable_index(transcript_dictionary):
    '''Takes a transcript dictionary, transforms to a bisectable index'''
    speedy_dict = {}
    for v in transcript_dictionary.values():
        for feature in v.features:
            if feature.chromosome in speedy_dict:
                speedy_dict[feature.chromosome].append(feature.index_string)
            else:
                speedy_dict[feature.chromosome] = [feature.index_string]
    for chromosome in speedy_dict.values():
        chromosome.sort()
    return speedy_dict

def initialize_features(indexes,index_attribute,genome_fasta):
    '''Build Features'''
    annotation,sequences,sub_index = {},read_in_fasta(genome_fasta), getattr(indexes,index_attribute)
    for name, features_list in sub_index.items():
        strand,chromosome = indexes.strand_map[name],indexes.chromosome_map[name]
        for feat in features_list:
            sequence = extract_sequence(sequences,chromosome,feat.start,feat.end)
            if sequence:
                feature = Transcript_Feature(name,chromosome,feat.start,feat.end,strand,feat.feat_class,feat.feat_numb,sequence)
                if name in annotation:
                    annotation[name].add_feature(feature)
                else:
                    annotation[name] = Transcript()
                    annotation[name].add_feature(feature)
    return annotation

def collect_modifications(transcripts_dict,bisectable_map,vcf_fyle,replace=None,f_range=10):
    '''Walk through a VCF file, passing all alterations through the bisectable map, collect in transcripts'''
    with open(vcf_fyle,'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            else:
                CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO,FORMAT,BO = line.strip().split()
                gff_chromosome = replace[CHROM] if replace else CHROM
                modification = (int(POS),REF,ALT)
                local_index = bisect.bisect(bisectable_map[gff_chromosome], ((int(POS), ''), ''))
                local_features = bisectable_map[gff_chromosome][max(0, local_index-f_range):local_index+f_range]
                hits = [feat for feat in local_features if feat[0][0]<= int(POS) and int(POS) <=feat[0][1]]
                if hits:
                    for hit in hits:
                        feature_name = hit[1]
                        transcript_name,f_type,f_numb = hit[1].split('~')
                        transcripts_dict[transcript_name].collect_modification(modification,feature_name)

def collect_poly_modifications(transcripts_dict,bisectable_map,vcf_fyle,replace=None,f_range=10):
    '''Walk through a VCF file, passing all alterations through the bisectable map, collect in transcripts'''
    with open(vcf_fyle,'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            else:
                CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO,FORMAT,BO = line.strip().split()
                gff_chromosome = replace[CHROM] if replace else CHROM
                alt_freqs = dict([q.split('=') for q in INFO.split(';') if '=' in q])['AF'].split(',')
                chaos_flag  = True if len(alt_freqs) !=1 or int(float(alt_freqs[0])) != 1 else False
                modification = (int(POS),REF,ALT) if not chaos_flag else generate_poly_site(POS,REF,ALT,alt_freqs)
                if modification:
                    local_index = bisect.bisect(bisectable_map[gff_chromosome], ((int(POS), ''), ''))
                    local_features = bisectable_map[gff_chromosome][max(0, local_index-f_range):local_index+f_range]
                    hits = [feat for feat in local_features if feat[0][0]<= int(POS) and int(POS) <=feat[0][1]]
                    if hits:
                        for hit in hits:
                            feature_name = hit[1]
                            transcript_name,f_type,f_numb = hit[1].split('~')
                            transcripts_dict[transcript_name].collect_modification(modification,feature_name)

def transcripts_to_fasta_dict(transcript_dictionary):
    '''Applies all modifications, converts to straight transcript:seq dictionary'''
    fasta_dictionary = {}
    for transcript in transcript_dictionary.values():
        fasta_name = transcript.transcript_name
        #Apply the modifications
        for feature in transcript.features:
            feature.apply_modifications()
        #Generate the modified transcript
        fasta_dictionary[fasta_name]= transcript.assemble_transcript()
    return fasta_dictionary
        
def write_fasta(fastadict,outfile,LW=80):
    '''Derp'''
    with open(outfile, 'w') as out:
        for k,v in sorted(fastadict.items()):
            out.write('>' + k + '\n')
            for i in xrange(0,len(v),LW):
                out.write(v[i:i+LW] + '\n') 

def read_convert_txt(convert_txt):
    '''If for some reason the <.vcf> chromosomes have a different nomenclature, adjust'''
    convert = {}
    with open(convert_txt,'r') as f:
        for line in f:
            if line.strip():
                vcf_chrm,annotation_chrm = line.strip().split(':')
                convert[vcf_chrm] = annotation_chrm
    return convert

def ez_rename(in_fyle,suffix):
    '''Renames an output based on the input name'''
    return '_'.join(['.'.join(in_fyle.split('.')[:-1]),suffix])+'.fa'

#Main Function
def main():
    parser = argparse.ArgumentParser(description='Assembles haplotype specific <.fasta> of sequence features with <.vcf> data')
    parser.add_argument('genome',type=str,help='Genome <.fasta>')
    parser.add_argument('annotation',type=str,help='Annotation <.gff/.gtf/.gff3>')
    parser.add_argument('feature',type=str.upper, default = None,choices = ['RNA','CDS','FP','TP'],help = 'Target Features')
    parser.add_argument('-method',type=str.upper, default='FIXED',choices= ['FIXED','POLY'],help = '[default = FIXED] VCF Filtering Used')
    parser.add_argument('-outdir',type=str,default = 'out_fasta',help='[default = out_fasta] Out Directory')
    parser.add_argument('-single',type=str, default = None,help='[<.vcf>] Run on single <.vcf>, rather than directory')
    parser.add_argument('-chr_names',type=str, default = None,help='[<.txt>] Convert <.vcf> chromosomes')
    args = parser.parse_args()
    
    #Gather Info and extras
    selected_index = {'RNA':'exon_map','CDS':'cds_map','FP':'fp_map','TP':'tp_map'}[args.feature]
    chrom_convert = read_convert_txt(args.chr_names) if args.chr_names else None
    fyle_que = sorted(glob.glob('*.vcf')) if not args.single else [args.single]
    collectors = {'FIXED':collect_modifications,'POLY':collect_poly_modifications}

    #Assemble Indexes
    full_indexes = annotation_parser.OmicsMap(args.annotation)
    full_indexes.build_check()
    transcripts = initialize_features(full_indexes,selected_index,args.genome)
    bisect_map = create_bisectable_index(transcripts)

    #Out Directory Stuff
    if not os.path.isdir(args.outdir):
        os.mkdir(args.outdir)
    
    #Process File Que
    for fyle in fyle_que:
        #Generate Out File Name
        out_file = ez_rename(fyle,args.feature)
        
        #Replicate the transcripts
        copy_transcripts = copy.deepcopy(transcripts)

        #Secret Collect!
        collectors[args.method](copy_transcripts,bisect_map,fyle,chrom_convert)
        
        #Convert, then write!
        fasta_sequences = transcripts_to_fasta_dict(copy_transcripts)
        write_fasta(fasta_sequences,os.path.join(args.outdir,out_file))

if __name__ == '__main__':
    main()