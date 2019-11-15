#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
SVM model for genome-wide HS-region prediction
@author: Wei He
@E-mail: whe3@mdanderson.org
@Date: 9/18/2018
"""

########import all the pacakges needed###########          
from __future__ import division                             
import pandas as pd                                         
import matplotlib.pyplot as plt                             
from pandas import DataFrame                                
import numpy as np                                          
import seaborn as sns                                                 
import math,os,json,logging,sys,pkg_resources


########Packages for machine learning############
from sklearn.neighbors.kde import KernelDensity    
from sklearn.ensemble import BaggingClassifier
from sklearn.multiclass import OneVsRestClassifier
from sklearn.svm import SVC
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import KFold
from sklearn.preprocessing import scale

#RequiredFilePath = pkg_resources.resource_filename(__name__, 'StaticFiles')
###################### Define functions for HS prediction and visualization #####################

'''This function is to get KDE of PTMs distribution of certain gene
   
   o input: (1) gene id: the symbol of target gene, eg: CREBBP
            (2) PTM dic: A python dict recording genome-wide Post-translational 
                modification sites annotation 
            (3) aa_list: all the residue positions along protein encoded by certain gene.
            (4) bw: bandwidth for kernel density estimation
            
   o output: return kde value of input PTMs sites distribution.
'''
def GetPTMList(gene,PTM_dic,aa_list,bw):
    x_grid = np.array(aa_list).reshape(-1,1)
    kde_value = [0]*len(x_grid)
    
    if gene in PTM_dic:
        ptm_list = []
        for i in PTM_dic[gene]:
            if int(i[1]) >=2:
                ptm_list.append(i[0])
        resNum = max(aa_list)
        x = np.array(ptm_list).reshape(-1,1)
        if len(x) != 0:
            kde_value = kde_sklearn(x,x_grid,bandwidth=bw)
    return kde_value



'''This function is to get KDE of SIFT score at each residue position
   
   o input: (1)con_list: a list of SIFT score at each residue position
            (2)aa_list: residue list of certain gene
            (3)bw: bandwidth for kde training
   
   o output: This function returns a list of kde values for SIFT score list.
'''
def GetSIFTKde(con_list,aa_list,bw):
    con_list = list(pd.cut(con_list,100,labels=False))
    x_grid = np.array(aa_list).reshape(-1,1)
    con_kde_value = [0]*len(con_list)
    ls = []
    for i in range(0,len(con_list)):
        ls += [(i+1)]*int(con_list[i])
    #print ls
    x = np.array(ls).reshape(-1,1)
    
    con_kde_value = kde_sklearn(x, x_grid,bw)
    return con_kde_value



'''This function is to get residue list of certain gene 

   o input: the exon(start,end) list of certain gene 
   
   o output: residue list of certain protein [1,resMax]
    
'''
def GetResList(exon_list):
    resLoc = 0 ;resNumList = [1]
    for i in range(len(exon_list)):
        resLoc += math.ceil((exon_list[i][1] - exon_list[i][0])/3)
        resNumList.append(resLoc)
    resMax = int(resNumList[-1])
    return range(1,resMax+1)


'''This function is to get genomic locatin of a certain residue
   
   o input: (1) exon_list: A list of exon(start,end) of certain gene.
            (2) aa_pos: the residue position which you want to find 
                corresponding genomic location.
            (3) starnd: The strand of gene on the genome. '+' or '-'
   
   o output: the genomic loci correponding to the input residue position.
'''
def GetGenomicLocation(exon_list,aa_pos,strand):
    resLoc = 0; gene_loci = 0
    for i in range(len(exon_list)):
        resLoc += math.ceil((exon_list[i][1] - exon_list[i][0])/3)
        if resLoc >= aa_pos: 
            if strand == '+':
                gene_loci = exon_list[i][1] - (resLoc - aa_pos)*3
            else: 
                gene_loci = exon_list[i][0] + (resLoc - aa_pos)*3
            break
    return gene_loci
            
        
'''Get the kde given a list of one-dimension scores or positions
   The function will return returns the log-likelihood of the samples
'''
def kde_sklearn(x, x_grid, bandwidth):
    '''Kernel Density Estimation with Scikit-learn'''
    kde_skl = KernelDensity(kernel='gaussian',bandwidth=bandwidth)
    kde_skl.fit(x)
    #score_samples''() returns the log-likelihood of the samples
    log_pdf = np.exp(kde_skl.score_samples(x_grid))
    return log_pdf



''' This is the main function here to gain all the information for a certain gene for HS 
    region prediction including PTMs, SIFT score, domain information, Secondary structure,

    o input: (1) gene: The symbol of target gene, eg: CREBBP
             (2) bd1: the bandwidth for PTMs KDE training.
             (3) bd2: the bandwidth for SIFT score training.
             (4) exons_dic: A python dict recording genome-wide exon information annptation.
             (5) doamin_genome_dic: A python dict recording genome-wide domain annotation.
             (6) ss_genome_dic: A python dict recording genome-wide secondary structure annotation.
             (7) sift_genome_dic: A python dict recording genome-wide SIFT score annotation.
             (8) po_dic: A python dict recording genome-wide phosphoruylation sites annotation.
             (9) ac_dic: A python dict recording genome-wide acetylation sites annotation.
             (10) me_dic: A python dict recording genome-wide methylation sites annotation.
             (11) ub_dic: A python dict recording genome-wide ubiquitylation sites annotation.
            
            
     o output:  The function will return a dataframe containing information of all the features for 
                the protein encoded by target gene.
'''
def GetFeatureTable(gene,bd1,bd2,exons_dic,domain_genome_dic,ss_genome_dic,sift_genome_dic,po_dic,
                    ac_dic,me_dic,ub_dic):
        
    logging.info('Start annotation process for %s'%gene)
    logging.info('Get exons information of target gene...')
    if gene in exons_dic:
        strand = exons_dic[gene]['strand']
        if strand == '+':
            exon_list = exons_dic[gene]['exons']
        else:
            exon_list = exons_dic[gene]['exons']
            exon_list.reverse()
    
    aa_list = GetResList(exon_list)
    aa_index = [i/len(aa_list) for i in aa_list]
    
    logging.info('Extract residues located within annotated pfam protein domain...')
    dom_ls = []
    if gene in domain_genome_dic:
        for dom in domain_genome_dic[gene].keys():
            start = domain_genome_dic[gene][dom]['start']
            end = domain_genome_dic[gene][dom]['end']
            dom_ls += range(start,end+1)
    else:
        logging.warning('Domain annotation is not available for %s,try another name.'%gene)
    ## For each features, use a list to present values at each residue position
    gene_list = []; con_list = []; domain_list=[]; po_kde_list=[]; B_list=[]
    ac_kde_list = []; me_kde_list=[]; ub_kde_list=[];iso_list=[]; H_list=[]
    
    logging.info('Get conservation(SIFT) score at each residue position...')
    if gene in sift_genome_dic:
        con_list = sift_genome_dic[gene]
        con_kde_value = GetSIFTKde(con_list,aa_list,bd2)
    else:
        con_kde_value = [0]*len(aa_list)
        logging.warning('SIFT score is not available for %s,try another name.'%gene)
    
    logging.info('Do loop to gain secondary structure and domain information at each residue...')
    for res in aa_list:  
        gene_list.append(gene)
        
        ## To get secondary structure prediction at each residue position
        if gene in ss_genome_dic:
            if res <= len(ss_genome_dic[gene]):
                ss = ss_genome_dic[gene][res-1]
                if ss == 'H': H_list.append(1);B_list.append(0)
                if ss == 'E': H_list.append(0);B_list.append(1)
                if ss == 'C': H_list.append(0);B_list.append(0)
            else: H_list.append(0);B_list.append(0)
        else: H_list.append(0);B_list.append(0)
        
        ## To get domain annotation at each residue position
        if res in dom_ls:
            domain_list.append(1)
        else:domain_list.append(0)
    
    logging.info('Smoothing the Conservation score and PTMs distribution using kernel density estimation')
    po_kde_value = GetPTMList(gene,po_dic,aa_list,bd1)
    ac_kde_value = GetPTMList(gene,ac_dic,aa_list,bd1)
    me_kde_value = GetPTMList(gene,me_dic,aa_list,bd1)
    ub_kde_value = GetPTMList(gene,ub_dic,aa_list,bd1)
    
    
    logging.info('Generating the dataframe to record all the feature information of the protein'+ 
                 'at each residue location...')
    df = DataFrame({'Gene':gene_list,'AA_index':aa_index, 'Helix':H_list, 'Sheet':B_list,'con_kde':con_kde_value,
                          'phosphorylation':po_kde_value,'acetylation':ac_kde_value,'methylation':me_kde_value,
                           'ubiquitylation':ub_kde_value,'In_domain':domain_list},
                           columns=['Gene','AA_index','con_kde','Helix','Sheet','phosphorylation',
                                    'acetylation','methylation','ubiquitylation','In_domain'])
    
    
    ## Use log(value+sudocount) to present kde scores
    s1 = 4.94065645841e-324; s2 = 2.00614554652e-81 
    df['phosphorylation'] = [np.log10(i+s1) for i in df['phosphorylation']]
    df['acetylation'] = [np.log10(i+s1) for i in df['acetylation']]
    df['ubiquitylation'] = [np.log10(i+s1) for i in df['ubiquitylation']]
    df['methylation'] = [np.log10(i+s1) for i in df['methylation']]
    df['con_kde'] = [np.log10(i+s2) for i in df['con_kde']]
    
    return df



'''
## Predict HS region using SVM model: Input a dataframe for training and a dataframe for quering
## Returns: 
   # (1)pred_score: a list of scores at each residue position; 
   # (2)pred_class: a list of 1(hs)or0(non-hs)
'''
def SVMpredictor(df_train,df_query,c,g):
    n_estimators= 100  # Using a bagging SVM to train the model
    clf = OneVsRestClassifier(BaggingClassifier(SVC(C=c,kernel='rbf', probability=True,gamma=g), 
                                                    max_samples=1.0/n_estimators, n_estimators=n_estimators,random_state=10))
    ## Scaling the training set and testing set
    df_train_scale = scale(df_train.iloc[:,2:])
    df_test_scale = scale(df_query.iloc[:,1:])
    
    clf.fit(df_train_scale, df_train['class'])
    pred_score = list(clf.decision_function(df_test_scale))
    pred_class = clf.predict(df_test_scale)
    
    return pred_score,pred_class
    

    
'''This funciton is to predict HS regions of the target protein given the feature table.

   o input:  (1) gene: The symbol of target gene, eg: CREBBP
             (2) bd1: the bandwidth for PTMs KDE training.
             (3) bd2: the bandwidth for SIFT score training.
             (4) c: The penalty paramether for SVM model.
             (5) g: The gamma parameter for SVM model.
             (6) exons_dic: A python dict recording genome-wide exon information annptation.
             (7) doamin_genome_dic: A python dict recording genome-wide domain annotation.
             (8) ss_genome_dic: A python dict recording genome-wide secondary structure annotation.
             (9) sift_genome_dic: A python dict recording genome-wide SIFT score annotation.
             (10) po_dic: A python dict recording genome-wide phosphoruylation sites annotation.
             (11) ac_dic: A python dict recording genome-wide acetylation sites annotation.
             (12) me_dic: A python dict recording genome-wide methylation sites annotation.
             (13) ub_dic: A python dict recording genome-wide ubiquitylation sites annotation.
             (14) outputdir: The directroy to save the result files.
             
    o output: The function will return a dataframe recording the features and SVM predicted HS region
              and corresponding scores. Also write the dataframe into a .csv file.
'''
def HSprediction(gene,TrainingFile,bd1,bd2,c,g,exons_dic,domain_genome_dic,ss_genome_dic,sift_genome_dic,
                 po_dic,ac_dic,me_dic,ub_dic,outputdir):
    
    ## SVM_Input.csv is the table file of training data with features and labels.
    
    df_train = pd.read_csv(TrainingFile).loc[:,['Gene','class','AA_index','con_kde','Helix','Sheet',
              'phosphorylation','acetylation','methylation','ubiquitylation','In_domain']]
    
    df_query = GetFeatureTable(gene,bd1,bd2,exons_dic,domain_genome_dic,ss_genome_dic,
                               sift_genome_dic,po_dic,ac_dic,me_dic,ub_dic)
    
    df_test_scale = scale(df_query.iloc[:,1:])
    pred_score,pred_class = SVMpredictor(df_train,df_query,c,g)
    
    df_query['SVM_Score'] = pred_score
    df_query['SVM_Class'] = pred_class
    df_query.to_csv(os.path.join(outputdir,gene+'_HS_Region_file.csv'))
    return df_query
    

''' This function is to plot bar track to present continous values like SIFT score, PTMs KDE, etc.
    
    o input: (1) resMax: the length of the target protein.
             (2) cont_list: A list of continous values.
             (3) pos: the position to locate the bar track.
             (4) name: the title for the bar track.
             (5) c: Color used for drawing the color bar: 'Reds','Greens','Oranges',etc.
             (6) s: True or False to determine the sorting order.
             
    o output: This function will plot the bar track given a list of continous values.
'''
def PlotContlist(resMax,cont_list,pos,name,c,s):
    plt.text(-(resMax/7),pos+0.2,name,fontsize=10)
    pair_list = [(i+1,cont_list[i]) for i in range(len(cont_list))]
    pair_list = sorted(pair_list,key=lambda x:x[1],reverse=s)
    color_list = sns.color_palette(c,len(cont_list))
    
    zero_list = [i for i in pair_list if i[1]==min(cont_list)]
    non_zero = [i for i in pair_list if i[1]!=min(cont_list)]
    plt.bar([i[0] for i in zero_list],height=[0.6]*len(zero_list),width=1,bottom=pos,color=color_list[0])
    plt.bar([i[0] for i in non_zero],height=[0.6]*len(non_zero),width=1,bottom=pos,color=color_list[len(zero_list):])    



''' This function is to plot bar track to present binary values like domain annotation,secondary structure, etc.
    
    o input: (1) resMax: the length of the target protein.
             (2) bi_list: A list of binary values.
             (3) pos: the position to locate the bar track.
             (4) name: the title for the bar track.
             (5) c: Color used for drawing the color bar: 'Reds','Greens','Oranges',etc.
             
    o output: This function will plot the bar track given a list of binary values.
'''
def PlotBilist(resMax,bi_list,pos,name,c):
    plt.text(-(resMax/7),pos+0.2,name,fontsize=10)
    loc_list = [i for i in range(resMax) if bi_list[i]==1]
    plt.bar(resMax/2,0.6,width=resMax,bottom=pos,color='silver',alpha=0.2)
    plt.bar(loc_list,[0.6]*len(loc_list),width=1.5,bottom=pos,color=c)


''' This function is to plot predicted HS regions together with other protein annotations
    
    o input: (1) df_hs: the dataframe generated by GetFeatureTable function.
             (2) gene: the symbol of target gene, eg: CREBBP.
             (3) outputdir: the directory to save result files
    
    o output: XX.png: the figure for the visualization of HS regions and protein annotations.

'''
def Plot(df_hs,gene,outputdir):
    resMax = df_hs.shape[0]
    svm_class = list(df_hs['SVM_Class'])
    svm_score = list(df_hs['SVM_Score'])
    con_list = list(df_hs['con_kde'])
    helix_list = list(df_hs['Helix'])
    sheet_list = list(df_hs['Sheet'])
    domain_list = list(df_hs['In_domain'])
    pho_list = list(df_hs['phosphorylation'])
    ace_list = list(df_hs['acetylation'])
    met_list = list(df_hs['methylation'])
    ubi_list = list(df_hs['ubiquitylation'])
    
    plt.figure(figsize=(10,8),dpi=300)
    plt.title('Predicted CKHS region and protein annotations for ' + gene)
    logging.info('Plot predicted HS regions...')
    PlotBilist(resMax,svm_class,8.0,'HS_Region','darkred')
    
    logging.info('Plot predicted SVM score...')
    PlotContlist(resMax,svm_score,7.0,'SVM_Score','Reds',False)
    
    logging.info('Plot SIFT score KDE...')
    PlotContlist(resMax,con_list,6.0,'SIFT_KDE','Oranges',True)
    
    logging.info('Plot pfam domain regions...')
    PlotBilist(resMax,domain_list,5.0,'Pfam_Domain','y')
    
    logging.info('Plot secondary structures...')
    PlotBilist(resMax,helix_list,4.0,'','r')
    PlotBilist(resMax,sheet_list,4.0,'','b')
    plt.text(-(resMax/7),4.2,'Helix',fontsize=10,color='r')
    plt.text(-(resMax/14),4.2,'Sheet',fontsize=10,color='b')
    
    logging.info('Plot PTMs KDE...')
    PlotContlist(resMax,pho_list,3.0,'Pho_KDE','Purples',False)
    PlotContlist(resMax,ubi_list,2.0,'Ubi_KDE','Oranges',False)
    PlotContlist(resMax,ace_list,1.0,'Ace_KDE','Reds',False)
    PlotContlist(resMax,met_list,0.0,'Met_KDE','Greens',False)
    
    ax = plt.gca()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    x = np.linspace(0,resMax);y = 0*x
    plt.plot(x,y,color='black')
    plt.yticks([])
    plt.xticks([i-i%10 for i in range(0,int(resMax),int(resMax/100)*10)])
    plt.xlabel('AA_Location')
    logging.info('Saving figures...')
    plt.savefig(outputdir +'/PredictedHSregion_'+gene+'.png')
    logging.info('Finished!')
    plt.show()
    
