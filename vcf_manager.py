#!/usr/bin/env python

'''
Runs vcftools in a batch fashion
'''

#Imports
import argparse
import glob
import subprocess
import os

def drive_vcf_tools(in_vcf,out_vcf,out_dir,quality,filter_type,filter_value='1'):
    '''drives vcftools'''
    command = ' '.join(['vcftools','--vcf',in_vcf,'--minQ',quality,'--recode','--recode-INFO-all',filter_type,filter_value,'--stdout','>',os.path.join(out_dir,out_vcf)])
    subprocess.call(command,shell=True)
    subprocess.call(' '.join(['mv','out.log',os.path.join(out_dir,out_vcf.replace('.vcf','.log'))]),shell=True)

def main():
    parser = argparse.ArgumentParser(description='Drives vcftools, filters for HaploHammer')
    parser.add_argument('-single',type=str,default=None,help='Operate on a single file, rather than current dir')
    parser.add_argument('-quality',type=str,default='30',help='[default = 30] vcftools \'--minQ\'')
    parser.add_argument('-in_ext',type=str,choices = ['gz','vcf'], help = 'Operate on files with this extension')
    parser.add_argument('-method',choices=['fixed','poly'],default='fixed',help='[default = fixed] Changes to keep')
    parser.add_argument('-outdir',type=str,default='filtered_vcf',help='[default = filtered_vcf] Out Directory')
    parser.add_argument('-rm_vcf',action='store_true', help='remove unfiltered vcf after filtering')
    parser.add_argument('-zip_out',action='store_true', help='compress new VCF after filtering')
    args = parser.parse_args()
    
    #Collect files to be operated on, set filter
    vcf_que = [args.single] if args.single else sorted(glob.glob('*.'+args.in_ext))
    run_type = {'fixed':'--non-ref-af','poly':'--non-ref-ac-any'}[args.method]
    
    #Out Directory Stuff
    if not os.path.isdir(args.outdir):
        os.mkdir(args.outdir)
    
    #Run function on files
    for vcf in vcf_que:
        working_file = vcf
        
        #Unzip if appropriate
        if working_file.endswith('.gz'):
            print '\033[94mUnzipping\033[0m {}'.format(vcf)
            subprocess.call(' '.join(['gunzip','-c',working_file,'>',working_file.strip('.gz')]),shell=True)
            working_file = working_file.strip('.gz')
        
        #Filter file
        out_suffix = '_'.join(['',args.method,'minQ'+args.quality])+'.vcf'
        suffixed_name = working_file.replace('.vcf',out_suffix)
        print '\033[92mFiltering\033[0m {} \033[92m -> \033[0m {}'.format(working_file,os.path.join(args.outdir,suffixed_name))
        drive_vcf_tools(working_file,suffixed_name,args.outdir,args.quality,run_type)
        
        #Zip Output
        if args.zip_out:
            print '\033[93mZipping\033[0m {}'.format(suffixed_name)
            subprocess.call(' '.join(['gzip',os.path.join(args.outdir,suffixed_name)]),shell=True)  
        
        #Remove unzipped VCF
        if args.rm_vcf:
            print '\033[91mDeleting\033[0m {}'.format(working_file)
            subprocess.call(' '.join(['rm',working_file]),shell=True)
        
        #Cleanup
        print ''

if __name__ == '__main__':
    main()