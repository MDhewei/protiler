#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 3 13:32:33 2019
@author: Wei He
@email: whe3@mdanderson.org
Call HS region from CRISPR tiling screen data and predict HS region from common protein features
"""

__version__ = "1.0.0"
import os, argparse,logging, sys,pkg_resources
from ProTiler_call import *
from ProTiler_predict import *



def ProTilerMain():
    ## Set logging format
    logging.basicConfig(level=logging.DEBUG,  
                    format='%(levelname)s:%(asctime)s @%(message)s',  
                    datefmt='%a, %d %b %Y %H:%M:%S',
                    filemode='a')
    ## Add aruguments for user input with command lines.
    parser = argparse.ArgumentParser(description='ProTiler: Call functionally essential protein \
             regions from CRISPR tiling screen data.')
    
    ## Add subparser: protiler call 
    subparsers =  parser.add_subparsers(help='commands to run ProTiler',dest='subcmd')
    subp_call = subparsers.add_parser('call',help='call HS region from CRISPR tiling screen data.')
    
    ## Add required group of parameters for call
    req_group = subp_call.add_argument_group(title='Required arguments for ProTiler call.',description='')
    req_group.add_argument('-i','--inputfile',type=str,help='The inputfile contains information of CRISPR tiling   screensincluding symbol of target gene(s),sgRNA information, targeted residue position and CRISPR scores. Accept .txt and .cvs fileformats',required=True)
    req_group.add_argument('-g','--gene_id',type=str,help='the symbol of target gene, for example: CREBBP',required=True)
    req_group.add_argument('-s','--score_columns',type=str,help='colums number(s) of input crispr score, for example:1,2,3',required=True)
    
    ## Add optional group of parameters for call
    opt_group = subp_call.add_argument_group(title='Optional arguments for ProTiler call.',description='')
    opt_group.add_argument('-f','--half_size',type=str,help='The number of neiboring signals from each side selected to filter inefficient sgRNAs',default='5')
    opt_group.add_argument('-t1','--threshold1',type=str,help='Threshold to supress the filtered signals',default='2')
    opt_group.add_argument('-t2','--threshold2',type=str,help='Threshold to detect changing points using TGUH method',default='1.5')
    opt_group.add_argument('-o','--outputdir',type=str,help='Directory to output files,if no directory is given, the result will be generated in current working directory',default='ProTilerOutput')
    
    ## Add sub parser protiler predict
    subp_predict = subparsers.add_parser('predict',help='predict HS region from public protein data.Only gene id is required.')
    
    ## Add required group of parameters for predict
    req_group2 = subp_predict.add_argument_group(title='Required arguments for ProTiler predict.',description='')
    req_group2.add_argument('-l','--gene_list',type=str,help='A list of candidate genes you for whih want to predict their HS regions,for example: CREBBP,EP300',required=True)
    
    ## Add optional group of parameters for predict.
    opt_group2 = subp_predict.add_argument_group(title='Optional arguments for ProTiler predict.',description='')
    opt_group2.add_argument('-b1','--bandwidth1',type=int,help='bandwidth for PTMs KDE',default=4)
    opt_group2.add_argument('-b2','--bandwidth2',type=int,help='bandwidth for SIFT KDE',default=26)
    opt_group2.add_argument('-m','--gamma',type=int,help='gamma paramether for SVM model',default=10)
    opt_group2.add_argument('-c','--penalty',type=int,help='Penalty paramether for SVM model',default=0.01)
    opt_group2.add_argument('-o','--outputdir',type=str,help='Directory to output files',default='ProTilerOutput')
    
    args = parser.parse_args() ## Obtain parameters from user input with commandline
    
    ## To generate the outputdir for saving result files
    try:
        os.mkdir(args.outputdir)
        logging.info('Creat the outputdir {} to place result files'.format(args.outputdir))
    except OSError:
        logging.warning('outputdir {} already exist'.format(args.outputdir))
    
    outputdir = os.path.join(os.getcwd(),args.outputdir) ## Set the path to outputdir
    
    ## To load and read genome-wide annotation files
    #RequiredFilePath = pkg_resources.resource_filename(__name__, 'StaticFiles')
    
    logging.info('Loading and reading genome-wide exons information files...')
    exons_dic_file = os.path.join(RequiredFilePath,'exons_dict')
    exons_dic = json.loads( open(exons_dic_file).read())
    
    logging.info('Loading and reading genome-wide domain annotation...')
    domain_dic_file =  os.path.join(RequiredFilePath,'domain_dic') 
    domain_genome_dic = json.loads( open(domain_dic_file).read())
    
    logging.info('Loading and reading genome-wide secondary structure annotation files...')
    ss_dic_genome = os.path.join(RequiredFilePath,'ss_genome_dic_new')
    ss_genome_dic = json.loads(open(ss_dic_genome).read())
    
    logging.info('Loading and reading genome-wide conservation annotation files...')
    sift_dic_genome = os.path.join(RequiredFilePath,'sift_genome_dic')
    sift_genome_dic = json.loads(open(sift_dic_genome).read())
    
    logging.info('Loading and reading genome-wide post-translational modification file...')
    ac_file = os.path.join(RequiredFilePath,'ac_dic')
    me_file = os.path.join(RequiredFilePath,'me_dic')
    po_file = os.path.join(RequiredFilePath,'po_dic')
    ub_file = os.path.join(RequiredFilePath,'ub_dic') 
    ac_dic = json.loads( open(ac_file).read())
    me_dic = json.loads( open(me_file).read())
    po_dic = json.loads( open(po_file).read())
    ub_dic = json.loads( open(ub_file).read())
    
    
    if args.subcmd == None:
        parser.print_help()
        sys.exit(0)
        
    if args.subcmd == 'call':
        ## Obtain parameters from user input
        Inputfile =  args.inputfile; gene = args.gene_id
        columns = args.score_columns; size = args.half_size
        th1 = args.threshold1; th2 = args.threshold2
        
        with open(Inputfile,'r') as f:
            firstline = f.readline()
            if 'Symbol' not in firstline or 'AA' not in firstline:
                logging.error('Input table should contain columns with name: Symbol \
                              and AA to present gene symbol and residue positions')
                sys.exit(-1)
        
        CallHSRegions(Inputfile,gene,columns,size,th1,th2)
        Visualization(gene,args.outputdir,exons_dic,sift_genome_dic,domain_genome_dic,
                      ss_genome_dic,po_dic,ac_dic,me_dic,ub_dic)
        
    if args.subcmd == 'predict':
        gene_str = args.gene_list; gene_list = gene_str.split(',')
        bd1 = args.bandwidth1; bd2 = args.bandwidth2
        c = args.gamma; g = args.penalty
        TrainingFile = os.path.join(RequiredFilePath,'SVM_Input.csv')
        
        for gene in gene_list:
            logging.info('Predicting HS regions for %s...'%gene)
            
            ## jump out the current loop if gene is not in the current database
            if exons_dic.has_key(gene)==False:
                logging.warning('Annotations of '+ gene +' is not available in the database. Try another alias name')
                continue
            
            df_hs = HSprediction(gene,TrainingFile,bd1,bd2,c,g,exons_dic,domain_genome_dic,ss_genome_dic,
                                 sift_genome_dic,po_dic,ac_dic,me_dic,ub_dic,outputdir)
            Plot(df_hs,gene,outputdir)

if __name__ == '__main__':
    ProTilerMain()   