#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Functions for HS region calling and visulization
@author: Wei He
@E-mail: whe3@mdanderson.org
@Date: 10/18/2018
"""

#################### Import all the packages ##################
from __future__ import division                
import pandas as pd
import math,json,os,pkg_resources  
import matplotlib.pyplot as plt
import matplotlib; matplotlib.pyplot.switch_backend('agg') 
from matplotlib.gridspec import GridSpec
import numpy as np 
import seaborn as sns
import os,sys,logging,subprocess
import pandas as pd
from sklearn.neighbors.kde import KernelDensity 

RequiredFilePath = pkg_resources.resource_filename(__name__, 'StaticFiles')
############################# Functions defined to call and visualize HS regions #################################

'''This function is to call HS region using TGUH method

   o input:  (1) inputfile: The table file recording CRISPR tiling screen data:.txt and .csv are accepted
             (2) gene: the symbol of target gene, eg. CREBBP
             (3) columns: the column number in the table file representing CRISPR scores
             (4) size: the number of neighboring sgRNAs selected to filter inefficient guides
             (5) th1: the threshold used to supress outliers
             (6) th2: the threshold for TGUH method to detect changing points
             
   o Output: Tis function mainly run the Rscript 'ProTiler_TGUH.R', will result in three main files:
             (1) 'XX_Score.csv': Record the functional scores at each residue after filter
             (2) 'XX_Est.csv': Record all the mean score for all the segments along the protein
             (3) 'XX_Segments.csv': the segments call by TGUH together with other information for each segment

'''
def CallHSRegions(inputfile,gene,columns,size,th1,th2):
    logging.info('Segmentation using TGUH method...')
    if inputfile.upper().endswith('TXT')==False and inputfile.upper().endswith('CSV')==False:
        logging.error('Input file format is not supported, .csv or .txt is required!')
        sys.exit(-1)
    col = 'c('+columns+')'
    tguh_script = os.path.join(RequiredFilePath,'ProTiler_TGUH.R')
    logging.info('Runing R script:%s'%tguh_script)
    cmd = ['Rscript',tguh_script,inputfile,gene,col,size,th1,th2]
    ret = subprocess.call(cmd, shell=False)
    if ret:
        logging.error('{} exited with code {}.'.format(cmd[1],ret))
        sys.exit(ret)



'''This function is plot the scatter of scores for all the sgRNAs and the segments call by TGUH method

   o input:  (1) ResMax: The length of protein encoded by target gene. 
             (2) gene: the symbol of target gene, eg. CREBBP
             (3) Two table files recording functional scores for each residue and segment generated 
                 by CallHSRegion function: gene+'_Score.csv' and gene+'XX_Est.csv'
             
   o Output: Tis function will generate a scatter plot presenting all the fucntional scores for each 
             sgRNA at certain residue position. With inefficient sgRNAs plotted grey and efficient sgRNAs
             plotted black. Also plot several lines presenting segments of the signals.
             
'''
def PlotScatterOfScore(resMax,gene):
    ## Read the score information in table files
    df_est = pd.read_csv(gene+'_Est.csv')
    df_score = pd.read_csv(gene+'_Score.csv')
    df_score_1 = df_score[df_score['is.efficient']==1]
    df_score_0 = df_score[df_score['is.efficient']==0]
    
    est_list = df_est['aa.est']
    score_list_1 = df_score_1['v.filtered']
    aa_list_1 = df_score_1['aa.pos']
    score_list_0 = df_score_0['aver.score']
    aa_list_0 = df_score_0['aa.pos']
    #aa_list = list(df_score['aa.pos'])
    #score_list = df_score['aver.score']
    
    ## Plot scatter plot
    #plt.ylabel('Tiling CRISPR Z-Score',fontsize=8)
    plt.scatter(aa_list_1,score_list_1,s=4,marker='o',label='Efficient')
    plt.scatter(aa_list_0,score_list_0,s=4,marker='o',c='silver',label='Inefficient')
    x = np.linspace(0,resMax);y = 0*x
    
    ## Plot segment lines 
    plt.plot(range(0,len(est_list)),est_list,c='r',label='Seg_pattern')
    plt.plot(x,y,color='black',linewidth=0.8)
    plt.text(-(resMax/6.5),(max(score_list_1)+min(score_list_1))/2,'Z-score:',fontsize=14)
    
    ## Remove the intermediate files
    os.system('rm '+gene+'_Est.csv')
    os.system('rm '+gene+'_Score.csv')
    
def MergeHSregion(hs_list,start_index=0):
    for i in range(start_index,len(hs_list)-1):
        if int(hs_list[i+1][1]) - int(hs_list[i][2]) < 5:
            hs_new = [hs_list[i][0],hs_list[i][1],hs_list[i+1][2],hs_list[i][3]+hs_list[i+1][3],
                      (hs_list[i][4]+hs_list[i+1][4])/2, hs_list[i][5]+hs_list[i+1][5]]
            hs_list[i] = hs_new
            del hs_list[i+1]
            #print hs_list
            return MergeHSregion(hs_list,start_index=i)
            
    return hs_list
        
            
'''This function is plot the HS region called by our method in a bar track

   o input:  (1) resMax: The length of protein encoded by target gene. 
             (2) gene: the symbol of target gene, eg. CREBBP
             (3) outputdir: the directory to save the result files
             (4) Table files recording all the informations for HS regions: 'XX_Segments.csv'
             
   o Output: Tis function will generate a track recording HS regions(dark color) 
             and non-HS region(light color)
             
'''
def PlotHSregion(resMax,gene,outputdir):
    ## Read the HS region information from table file.
    df_hs = pd.read_csv(gene+'_Segments.csv')
    #print df_hs
    hs_list = []
    for i in range(0,df_hs.shape[0]):
        hs_list.append([df_hs['Gene'][i],df_hs['AA.start'][i],df_hs['AA.end'][i],df_hs['n'][i],
                        df_hs['m'][i],df_hs['lenght'][i]])
    
    hs_list = MergeHSregion(hs_list)
    df_hs =  pd.DataFrame(hs_list,columns=['Gene','AA.start','AA.end','n','m','length'])
    plt.text(-(resMax/6.5),7.2,'CKHS_Regions:',fontsize=14)
    
    ## Generate a blank track with light color
    plt.bar(resMax/2,0.6,width=resMax,bottom=7.0,color='silver',alpha=0.2)
    #plt.bar(start + (end-start)/2, 0.4,width=end-start,bottom = 2.0,facecolor='darkred',alpha=0.3)
    for i in range(df_hs.shape[0]):
        start = list(df_hs['AA.start'])[i]
        end = list(df_hs['AA.end'])[i]
        score = list(df_hs['m'])[i]
        plt.bar(start + (end-start)/2, 0.6 ,width=end-start,bottom = 7.0,facecolor='darkred',alpha=0.3)
        #x = np.linspace(start+5,end-5);y = 0*x+3.5
    
    df_hs.to_csv(outputdir+'/'+gene+'_CKHS_Regions.csv')
    ## move the file recording HS regions into outputfir
    os.system('rm '+gene+'_Segments.csv')

      
'''This function is plot the pfam domain annotation in a bar track.

   o input:  (1) resMax: The length of protein encoded by target gene. 
             (2) gene: the symbol of target gene, eg. CREBBP
             (3) domain_dic: A python dict recording domain annotations 
                 for all the gene in human genome.The structure is: 
                {gene_id:{domain1:{'start':N1,'end':N2},domain2:{'start':N3,'end':N4},...}
             
   o Output: Tis function will generate a track recording domain annotations(different colors)
             
'''
def PlotDomain(resMax,gene,domain_genome_dic): 
    ## Generate a blank track with light color
    plt.bar(resMax/2,0.6,width=resMax,bottom=4.0,color='silver',alpha=0.2)
    plt.text(-(resMax/6.5),4.2,'Pfam_Domains:',fontsize=14)
    
    ## Make sure the domain_dic contains the target genes
    if gene in domain_genome_dic:
       domain_dic = domain_genome_dic[gene]
       color_list = ['red','blue','purple','yellowgreen','orange','green','brown','pink',
                     'darkcyan','coral','darkorange','darkred','gray','silver']
       
       #color_list[i%len(color_list)]
       dom_list = [(dom,domain_dic[dom]['start'],domain_dic[dom]['end']) for dom in domain_dic.keys()]
       dom_list = sorted(dom_list,key=lambda x:x[1])
       for i in range(len(dom_list)):
           dom = dom_list[i][0]
           dom_loc = dom_list[i][1]
           dom_wid = dom_list[i][2]-dom_list[i][1]
           ## Plot the domain region with different color and text annotation
           plt.bar(dom_loc + dom_wid/2, 0.6 ,width=dom_wid,bottom = 4.0,
                   facecolor=color_list[0],alpha=0.3)
           ## Add text to annotate domain names 
           if len(dom) < 5:
               plt.text(dom_loc+dom_wid/2.3,3.5,dom,fontsize=12,horizontalalignment='center')
           else:
               plt.text(dom_loc+dom_wid/2.3,3.5,dom[0:8],fontsize=12,horizontalalignment='center')
                
    else:
         logging.info('Domain annotation is not avilable for %s'%gene)



'''This function is to get the residue list corresponding to each exons of target 
   gene and length of protein encoded by target gene.
   
   o input: (1) gene: the symbol(name) of taget gene, eg: CREBBP
            (2) exons_dic: a python dict recording genome-wide information about
                exons(start,end),strand,exon numbers,etc.
   
   O ouput: This function will return:
            (1) resNumList: a list of [residue.start,residue.end] corresponding 
                           to all the exons of the gene.
            (2) resMax: the length of protein encoded by target gene.
   
'''   
def GetResList(gene,exons_dic):
    ## Make sure that the exons_dic contains infotmation of target gene
    if gene in exons_dic: 
        #print gene
        strand = exons_dic[gene]['strand']
        if strand == '+':
            exon_list = exons_dic[gene]['exons']
        else:
            exon_list = exons_dic[gene]['exons']
            exon_list.reverse()
            
        resLoc = 0 ;resNumList = [1]
        for i in range(len(exon_list)):
            resLoc += math.ceil((exon_list[i][1] - exon_list[i][0])/3)
            resNumList.append(resLoc)
        resMax = resNumList[-1]
    return resNumList,resMax     


'''This function is to plot all the exons of the target gene in a bar track
   
   o input: two outputs of the function GetResList resNumList and resMax are 
            inputs of this function
   
   O ouput: This function will generate a track recording exons annotations
   
'''
def PlotExons(resNumList,resMax):
    plt.text(-(resMax/6.5),1.2,'Exons:',fontsize=14)
    ## Use dark and color to discriminate adjacent exons.
    for i in range(len(resNumList)-1):
        exon_wid = resNumList[i+1]-resNumList[i]
        if (i+1)%2 == 0: color = 'silver'
        else: color = 'black'
        plt.bar(resNumList[i]+exon_wid/2,0.6,width=exon_wid,bottom = 1.0, 
                facecolor = color, alpha=0.2)
        #plt.annotate(str(i+1),xy=(resNumList[i]+exon_wid/2,0),
        #xytext=(resNumList[i]+exon_wid/4,0),fontsize=8)
        
        
   
'''This function is to get KDE of PTMs distribution of certain gene
    
    o input: (1)aa_list: all the sgrna cutting sites along certain gene
             (2)bw: bandwidth for kde
    o output: return kde value at each sites
'''
def GetPtmKDE(ptm_list,resMax,bw):
    x_grid = np.array(range(1,int(resMax)+1)).reshape(-1,1)
    kde_value = [0]*len(x_grid)
    x = np.array(ptm_list).reshape(-1,1)
    if len(x) != 0:
        kde_value = kde_sklearn(x,x_grid,bandwidth=bw)
    return kde_value        



'''This function is to plot bar tracks recording different PTMs distributions.
   
   o input: (1) resMax: The length of protein encoded by target gene.
            (2) gene: Symbol of target gene, eg: CREBBP
            (3) po_dic: a python dic recording all the phosphorylation sites for whole genome.
            (4) ac_dic: a python dic recording all the acetylation sites for whole genome.
            (5) me_dic: a python dic recording all the methylation sites for whole genome.
            (6) ub_dic: a python dic recording all the ubiquitylation sites for whole genome.
            
   o output: This function will plot four bar tracks recording distribution of different PTMs
'''
def PlotPMS(resMax,gene,po_dic,ac_dic,me_dic,ub_dic):
    plt.text(-(resMax/6.5),-1.8,'Acetylation:',fontsize=14)
    plt.bar(resMax/2,0.6,width=resMax,bottom=-2.0,color='silver',alpha=0.2)
    if gene in ac_dic:
       ac_loc_list = [i[0] for i in ac_dic[gene] if i[1]>=2]
       ac_ref_list = [i[1] for i in ac_dic[gene] if i[1]>=2]
       #print ac_loc_list,ac_ref_list
       #ac_kde_list = GetPtmKDE(ac_loc_list,resMax,4)
       plt.bar(ac_loc_list,[0.6]*len(ac_loc_list),width=2,bottom=-2.0,color='r',label='Acety_site')
       #plt.bar(ac_loc_list,ac_ref_list,width=resMax/800,color='black')
    else:
        logging.warning('No acetylation sites in the protein!')
    
    plt.text(-(resMax/6.5),-4.8,'Methylation:',fontsize=14)
    plt.bar(resMax/2,0.6,width=resMax,bottom=-5.0,color='silver',alpha=0.2)
    if gene in me_dic:
       me_loc_list = [i[0] for i in me_dic[gene] if i[1]>=2]
       me_ref_list = [i[1] for i in me_dic[gene] if i[1]>=2]
       plt.bar(me_loc_list,[0.6]*len(me_loc_list),width=2,color='g',bottom=-5.0,label='Methy_site')
       #plt.bar(me_loc_list,me_ref_list,width=resMax/800,color='black')
    else:   
        logging.warning('No Methylation sites in the protein!')

    plt.text(-(resMax/6.5),-0.3,'Phosphorylation:',fontsize=14)
    plt.bar(resMax/2,0.6,width=resMax,bottom=-0.5,color='silver',alpha=0.2)
    if gene in po_dic:
       po_loc_list = [i[0] for i in po_dic[gene] if i[1]>=2]
       po_ref_list = [i[1] for i in po_dic[gene] if i[1]>=2]
       plt.bar(po_loc_list,[0.6]*len(po_loc_list),width=2,color='b',bottom=-0.5,label='Phos_site')
       #plt.bar(po_loc_list,po_ref_list,width=resMax/800,color='black')
    else:   
        logging.warning('No Phosphorylation sites in the protein!')

    plt.text(-(resMax/6.5),-3.3,'Ubiquitination:',fontsize=14)
    plt.bar(resMax/2,0.6,width=resMax,bottom=-3.5,color='silver',alpha=0.2)
    if gene in ub_dic:
       ub_loc_list = [i[0] for i in ub_dic[gene] if i[1]>=2]
       ub_ref_list = [i[1] for i in ub_dic[gene] if i[1]>=2]
       plt.bar(ub_loc_list,[0.6]*len(ub_loc_list),width=2,color='y',bottom=-3.5,label='Ubiq_site')
    else:   
        logging.warning('No ubiquitylation sites in the protein!')

    

        
'''This function is to plot bar track recording predicted secondary structures.
   
   o input: (1) resMax: The length of protein encoded by target gene.
            (2) gene: Symbol of target gene, eg: CREBBP
            (3) ss_genome_dic: a python dic recording all predicted secondary 
                structures for whole genome.
            
   o output: This function will plot bar track recording predicted secondary 
             structures: alpha helix(red), beta sheet(blue)
'''
def PlotSecondaryStructure(resMax,gene,ss_genome_dic):
    plt.text(-(resMax/6.5),2.7,'SS Predict:',fontsize=14)
    plt.bar(resMax/2,0.6,width=resMax,bottom=2.5,color='silver',alpha=0.2)
    plt.annotate('Helix',xy=(-(resMax/6.5),2.0),color='r',fontsize=12)
    plt.annotate('Sheet',xy=(-(resMax/12.5),2.0),color='b',fontsize=12)
    if gene in ss_genome_dic:
       ss_list = ss_genome_dic[gene]
       h_pos_list = []; e_pos_list=[]; c_pos_list=[]
       for res,ss in enumerate(ss_list,1):
           if ss == 'H':
               h_pos_list.append(res)
           elif ss == 'E':
               e_pos_list.append(res)
           elif ss == 'C':
               c_pos_list.append(res)
       ## red bar represent alpha helix and blue bar represent beta sheet
       plt.bar(h_pos_list,[0.6]*len(h_pos_list),width=1,color='r',bottom=2.5,label='Alpha helix')
       plt.bar(e_pos_list,[0.6]*len(e_pos_list),width=1,color='b',bottom=2.5,label='Beta sheet')
       #plt.scatter(c_pos_list,[1.1]*len(c_pos_list),marker='^',s=8,c='r')
       #plt.legend(loc='lower left',fontsize='small')
    else:
        logging.warning('Secondary structure is not available for %s.'%gene)
    

'''Get the kde for a given list of scores'''
def kde_sklearn(x, x_grid, bandwidth):
    ## Kernel Density Estimation with Scikit-learn
    kde_skl = KernelDensity(kernel='gaussian',bandwidth=bandwidth)
    kde_skl.fit(x)
    ##score_samples''() returns the log-likelihood of the samples
    log_pdf = np.exp(kde_skl.score_samples(x_grid))
    return log_pdf


'''This function is to get KDE of SIFT score at each residue position

    o input: (1) con_list: a list of SIFT score at each residue position
             (2) resMax: the length of the protein.
             (3) bw: bandwidth for kde training
             
    o output  This function returns a list of kde values
'''
def GetSIFTKde(con_list,resMax,bw):
    con_list = list(pd.cut(con_list,100,labels=False))
    x_grid = np.array(range(1,int(resMax)+1)).reshape(-1,1)
    #con_kde_value = [0]*len(con_list)
    ls = []
    for i in range(0,len(con_list)):
        ls += [(i+1)]*int(con_list[i])
    #print ls
    x = np.array(ls).reshape(-1,1)
    
    con_kde_value = kde_sklearn(x, x_grid,bw)
    return con_kde_value


'''This function is to plot bar track recording SIFT score at each residue position
   
   o input: (1) resMax: The length of protein encoded by target gene.
            (2) gene: Symbol of target gene, eg: CREBBP
            (3) con_genome_dic: a python dic recording all SIFT score at each
                residue position for whole genome.
            
   o output: This function will plot bar track recording SIFT score at each residue
             position. More conserved darker color.
'''
def PlotConScore(resMax,gene,sift_genome_dic):
    plt.text(-(resMax/6.5),5.7,'SIFT Score:',fontsize=14)
    start = int(-(resMax/10));end = int(-(resMax/17.5));Len=end-start
    
    #plt.bar(range(start,end),height=[0.35]*(end-start),width=1,bottom=5.0,color=sns.color_palette('Reds',Len))
    #plt.annotate('High',xy=(-(resMax/7.5),5.0),fontsize=8)
    #plt.annotate('Low',xy=(-(resMax/19),5.0),fontsize=8)
    
    if  gene in sift_genome_dic:
        score_list = sift_genome_dic[gene]
        if int(resMax) <= len(score_list):
            score_list = score_list[0:int(resMax)]
        else:
            score_list += [np.mean(score_list)]*(int(resMax)-len(score_list))
        #kde_list = GetSIFTKde(score_list,resMax,26)
        #PlotKDE(resMax,gene,kde_list,4.75,'SIFT_KDE:','Oranges',True)
        pair_list = [(i+1,score_list[i]) for i in range(len(score_list))]
        pair_list = sorted(pair_list,key=lambda x:x[1],reverse=True)
        pos_new = [x[0] for x in pair_list]
        plt.bar(pos_new,height=[0.6]*len(pair_list),width=1,bottom=5.5,color=sns.color_palette('Reds',len(score_list)))
    else:
        logging.warning('SIFT score is not avilable for %s.'%gene)
    


'''This function is to plot bar track recording KDE for SIFT score and PTMs distribution
   
   o input: (1) resMax: The length of protein encoded by target gene.
            (2) gene: Symbol of target gene, eg: CREBBP
            (3) kde_list: the kernel density estimation of certain list of value
            (3) pos: the position to draw the bar track
            (4) name: the name of the bar track, eg: 'SIFT_KDE:'
            (5) color: the color used for drawing the KDE score. ('Reds','Oranges','Greens',etc.)
            (6) a: False or True to judge the sorting order.
            
   o output: This function will plot bar track for the given kde value list.
'''
def PlotKDE(resMax,gene,kde_list,pos,name,color,a):
    logging.info('Annotate the title of track..')
    plt.text(-(resMax/8),pos+0.2,name,fontsize=8) 
    
    logging.info('Sort the kde list by values...')
    pair_list = [(i,kde_list[i-1]) for i in range(1,len(kde_list)+1)]
    pair_list = sorted(pair_list,key=lambda x:x[1],reverse=a)
        
    color_list = sns.color_palette(color,len(kde_list))
    logging.info('Figure out 0 values in the list...')
    ## In case that many position got same value equal to 0, but color still change as they sort randomaly
    zero_list = [i for i in pair_list if i[1]==min(kde_list)]
    non_zero = [i for i in pair_list if i[1]!=min(kde_list)]
    logging.info('Plot the color bars recording kde value...')  
    plt.bar([i[0] for i in zero_list],height=[0.6]*len(zero_list),width=1,bottom=pos,color=color_list[0])
    plt.bar([i[0] for i in non_zero],height=[0.6]*len(non_zero),width=1,bottom=pos,color=color_list[len(zero_list):])    
    

'''This function is the intergation of all the functions above to plot the HS region together with 
   other protein annotations.
   
   o input: (1) gene: Symbol of target gene, eg: CREBBP
            (2) outputdir: The directory to save the result files.
            (3) exons_dic: a python dict recording genome-wide information about
                exons(start,end),strand,exon numbers,etc.
            (4) con_genome_dic: a python dict recording genome-wide SIFT scores.
            (5) domain_dic: a python dict recording genome-wide domain annotations.
            (6) ss_genome_dic: a python dict recording genome-wide secondary structure annotations.
            (7) po_dic: a python dict recording genome-wide phosphorylation sites. 
            (8) ac_dic: a python dict recording genome-wide acetylation sites.
            (9) me_dic: a python dict recording genome-wide methylation sites.
           (10) po_dic: a python dict recording genome-wide ubiquitylation sites.
    
    o output: XX.png file: figure for the visualization of CRISPR screen scores, HS regions, protein annotations. 
   
'''
def Visualization(gene,outputdir,exons_dic,sift_genome_dic,domain_genome_dic,ss_genome_dic,po_dic,ac_dic,me_dic,ub_dic):
    if gene not in exons_dic:
        logging.error('The information of ' +Gene +' is not available in the database.' 
                      +'Visualization stop... Try another gene alias name')
        sys.exit(-1)
    resList,resMax = GetResList(gene,exons_dic)
    plt.figure(figsize=(12,8),dpi=300)
    gs = GridSpec(2, 1, height_ratios=[1,2.5])
    logging.info('Plotting scatters of z-scores and segments...')
    plt.rcParams["font.size"] = 14
    plt.rcParams["font.family"] = "Arial"
    plt.subplot(gs[0])
    plt.title('CKHS profile and annotations for '+ gene,fontsize=16)
    plt.xlim(-resMax/6,resMax*51/50)
    plt.xticks([])
    plt.yticks(fontsize=12)
    PlotScatterOfScore(resMax,gene)
    ax = plt.gca()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_position(('data',-0.1))
    
    plt.subplot(gs[1])
    plt.yticks([])
    plt.ylim(-6.0,8.0)
    plt.xlim(-resMax/6,resMax*51/50)
    plt.xticks([i-i%10 for i in range(0,int(resMax),int(resMax/100)*10)])
    logging.info('Plotting HS regions...')
    PlotHSregion(resMax,gene,outputdir)
    logging.info('Plotting Exon regions...')
    PlotExons(resList,resMax)
    logging.info('Plotting conservation score...')
    PlotConScore(resMax,gene,sift_genome_dic)
    logging.info('Plotting domain information...')
    PlotDomain(resMax,gene,domain_genome_dic)
    logging.info('Plotting secondary structures...')
    PlotSecondaryStructure(resMax,gene,ss_genome_dic)
    logging.info('Plotting post-translaitonal modifications...')
    PlotPMS(resMax,gene,po_dic,ac_dic,me_dic,ub_dic)
    ax = plt.gca()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    
    ax.spines['left'].set_visible(False)
    x = np.linspace(0,resMax);y = 0*x-6
    plt.plot(x,y,color='black')
    
    plt.text(-(resMax/6.5),-6,'AA Location:',fontsize=14)
    #plt.xlabel('Residue location',fontsize=10)
    plt.tight_layout(h_pad=0)
    logging.info('Saving figures...')
    plt.savefig(outputdir +'/Segmentfigure_'+gene+'.png')
    logging.info('Finished!')
    plt.show()


