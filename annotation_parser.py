#!/usr/bin/env python

'''
Parses <.gtf>/<.gff>/<.gff3> to easy to use indexes.
Formatting disparity may render this incompatible with some files.
'''

class OmicsFeature(object):
    '''Generic Feature'''
    def __init__(self,start=0,end=0,feat_class=None,transcript=None,feat_numb=0):
        self.start = int(start)
        self.end = int(end)
        self.feat_class = feat_class
        self.transcript = transcript
        self.feat_numb = int(feat_numb)
    
    def assign_number(self,number):
        self.feat_numb = int(number)
    
    def as_tuple(self):
        return (self.start,self.end,self.feat_class,self.feat_numb,self.transcript)
    
    def extend_right(self,extension):
        self.end = self.end + int(extension)
        
    def extend_left(self,extension):
        self.start = self.start - int(extension)

class OmicsMap(object):
    '''Takes a <.gtf> or <.gff>, should produce all relevant information stored Pythonically'''
    def __init__(self,annotation_file = ''):
        self.build_file = annotation_file
        
        #Basic Indexes
        self.feature_map = {}    #Transcript: [(OmicsFeature) ...]
        self.chromosome_map = {} #Transcript,Gene: Chromosome
        self.strand_map = {}     #Transcript,Gene: Strand
        self.gene_map = {}       #Gene: (OmicsFeature)
        
        #Derived Indexes
        self.exon_map = {}       #Transcript: [(Exon) ...]
        self.intron_map = {}     #Transcript: [(Intron)...]
        self.cds_map = {}        #Transcript: [(CDS)...]
        self.fp_map = {}         #Transcript: [(FPTUR)...]
        self.tp_map = {}         #Transcript: [(TPTUR)...]
        
        #Special Maps
        self.ex_map = {}         #Transcript: [(Exon)(Intron)....]
        self.utr_map = {}        #Transcript: [(UTR)...]
        
        if annotation_file:
            with open(annotation_file,'r') as f:
                for line in f:
                    if line.strip() and not line.startswith('#'):
                        chrm,source,item_type,start,end,qual,strand,phase,attribute = line.split('\t')
                        
                        if item_type == 'gene':
                            gene = gene_parser(attribute)
                            feature = OmicsFeature(start,end,item_type,gene)
                            if gene:
                                self.strand_map[gene] = strand
                                self.chromosome_map[gene] = chrm
                                self.gene_map[gene] = feature
                        
                        else:                       
                            transcript = attribute_parser(attribute)
                            feature = OmicsFeature(start,end,item_type,transcript)
                            if transcript:
                                self.strand_map[transcript] = strand
                                self.chromosome_map[transcript] = chrm
                                if transcript in self.feature_map:
                                    self.feature_map[transcript].append(feature)
                                else:
                                    self.feature_map[transcript] = [feature]
                
            #Create Sub-Maps
            for transcript, features in self.feature_map.items():
                exons = [e for e in features if e.feat_class == 'exon']
                coding_seqs = [c for c in features if c.feat_class == 'CDS']
                fp_utrs = [fp for fp in features if fp.feat_class.lower() == 'five_prime_utr']
                tp_utrs = [tp for tp in features if tp.feat_class.lower() == 'three_prime_utr']
                
                #<.gtf>s sometimes don't include the stop codon as 'CDS', Need to buffer +3/-3
                stop_codons = [sc for sc in features if sc.feat_class == 'stop_codon']
                extra_codon = True if stop_codons else False

                if self.strand_map[transcript] == '+':
                    
                    if exons:
                        exons.sort(key=lambda x: x.start)
                        for number, exon in enumerate(exons,1):
                            exon.assign_number(number)
                        introns = [OmicsFeature(exons[i-1].end+1,exons[i].start-1,'intron',transcript,i) for i in range(1,len(exons))]
                        self.exon_map[transcript] = exons
                        self.intron_map[transcript] = introns
                        self.ex_map[transcript] = sorted(exons+introns,key=lambda x: x.start)
                 
                    for generic_features in [coding_seqs,fp_utrs,tp_utrs]:
                        if generic_features:
                            generic_features.sort(key=lambda x: x.start)
                            for number, g_feature in enumerate(generic_features,1):
                                g_feature.assign_number(number)
                    
                    if coding_seqs:
                        #Correct for <.gtf> stop codons
                        if extra_codon:
                            coding_seqs[-1].extend_right(3)
                        self.cds_map[transcript] = coding_seqs
                    if fp_utrs:
                        self.fp_map[transcript] = fp_utrs         
                    if tp_utrs:    
                        self.tp_map[transcript] = tp_utrs        
                    if tp_utrs or fp_utrs:
                        self.utr_map[transcript] = sorted(fp_utrs+tp_utrs,key=lambda x:x.start)

                if self.strand_map[transcript] == '-':
                    
                    if exons:
                        exons.sort(key=lambda x: x.start,reverse=True)
                        for number, exon in enumerate(exons,1):
                            exon.assign_number(number)
                        exons.sort(key=lambda x: x.start)
                        introns = [OmicsFeature(exons[i-1].end+1,exons[i].start-1,'intron',transcript,len(exons)-i) for i in range(1,len(exons))]
                        self.exon_map[transcript] = exons
                        self.intron_map[transcript] = introns
                        self.ex_map[transcript] = sorted(exons+introns,key=lambda x: x.start)

                    for generic_features in [coding_seqs,fp_utrs,tp_utrs]:
                        if generic_features:
                            generic_features.sort(key=lambda x: x.start, reverse=True)
                            for number, g_feature in enumerate(generic_features,1):
                                g_feature.assign_number(number)
                            generic_features.sort(key=lambda x: x.start)

                    if coding_seqs:
                        #Correct for <.gtf> stop codons
                        if extra_codon:
                            coding_seqs[0].extend_left(3)
                        self.cds_map[transcript] = coding_seqs
                    if fp_utrs:
                        self.fp_map[transcript] = fp_utrs         
                    if tp_utrs:    
                        self.tp_map[transcript] = tp_utrs        
                    if tp_utrs or fp_utrs:
                        self.utr_map[transcript] = sorted(fp_utrs+tp_utrs,key=lambda x:x.start)
 
    def build_check(self):
        '''Prints out basic statst to verify the annotation build'''
        print ''
        print 'Used {} as annotation source...'.format(self.build_file)
        print 'Build contains:'
        print '{} genes located on {} chromosomes'.format(str(len(self.gene_map)),str(len(set(self.chromosome_map.values()))))
        print '{} transcripts containing {} exons'.format(str(len(self.exon_map)),str(sum(len(qbert) for qbert in self.exon_map.values())))
        print ''
        

def attribute_parser(attribute_field):
    '''Generic Transcript Field Extractor'''
    fields,transcript = iter([field.strip() for field in attribute_field.split(';')]),None
    while not transcript:
        try:
            field = fields.next()
            if 'transcript_id' in field:
                transcript= field.split()[1].strip('\"')
            elif 'Parent=' in field:
                if ':' in field:
                    transcript = field.split(':')[-1]
                else:
                    transcript = field.split(',')[0].split('=')[1]
        except StopIteration:
            return None
    return transcript  

def gene_parser(attribute_field):
    '''Gene Field Extractor'''
    fields,gene = iter([field.strip() for field in attribute_field.split(';')]),None
    while not gene:
        try:
            field = fields.next()
            if 'gene_id' in field:
                gene = field.split()[1].strip('\"')
            elif 'ID=' in field:
                if ':' in field:
                    gene = field.split(':')[-1]
                else:
                    gene = field.split('=')[1]
        except StopIteration:
            return None
    return gene


def main():
    pass

if __name__ == '__main__':
    main()