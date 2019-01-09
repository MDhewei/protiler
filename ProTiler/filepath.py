#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
import file paths to static files and outputdir
Created on Tue May 29 16:51:23 2018
@author: Wei He
"""
import os,json,logging
import pkg_resources

   
RequiredFilePath = pkg_resources.resource_filename(__name__, 'StaticFiles')
#print RequiredFilePath
#ResultFiles = os.getcwd()

#RequiredFilePath = os.path.join(os.getcwd(),'StaticFiles')
#ResultFiles = os.path.join(os.getcwd(),'ResultFiles')
#domain_file = os.path.join(RequiredFilePath,'pfam_genomic_position.txt')
#exons_dic_file = os.path.join(RequiredFilePath,'exons_dict')
#ac_file = os.path.join(RequiredFilePath,'ac_dic')
#me_file = os.path.join(RequiredFilePath,'me_dic')
#po_file = os.path.join(RequiredFilePath,'po_dic')
#ub_file = os.path.join(RequiredFilePath,'ub_dic') 
#su_file = os.path.join(RequiredFilePath,'su_dic') 
#con_file =  os.path.join(RequiredFilePath,'Conservation_dic') 
#ss_dic_file = os.path.join(RequiredFilePath,'ss_dic')
#con_kde_file = os.path.join(RequiredFilePath,'con_smooth_dic') 

#ac_kde_file = os.path.join(RequiredFilePath,'ac_kde_dic')
#me_kde_file = os.path.join(RequiredFilePath,'me_kde_dic')
#po_kde_file = os.path.join(RequiredFilePath,'po_kde_dic')
#ub_kde_file = os.path.join(RequiredFilePath,'ub_kde_dic') 

#exons_dic = json.loads( open(exons_dic_file).read())
#print 'Loading and reading genome-wide domain annotation...'
#domain_dic_file =  os.path.join(RequiredFilePath,'domain_dic') 
#print 'Done!'
#domain_all = json.loads( open(domain_dic_file).read())
#ac_dic = json.loads( open(ac_file).read())
#me_dic = json.loads( open(me_file).read())
#po_dic = json.loads( open(po_file).read())
#ub_dic = json.loads( open(ub_file).read())
#con_dic = json.loads( open(con_file).read())
#ss_dic = json.loads( open(ss_dic_file).read())

#ac_kde_dic = json.loads(open(ac_kde_file).read())
#me_kde_dic = json.loads(open(me_kde_file).read())
#po_kde_dic = json.loads(open(po_kde_file).read())
#ub_kde_dic = json.loads(open(ub_kde_file).read())
#con_kde_dic = json.loads(open(con_kde_file).read())



