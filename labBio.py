# -*- coding: utf-8 -*-
"""
Created on Tue Jan 26 15:33:18 2016

@author: Duarte
"""

from Bio import SeqIO
from Bio import Entrez
from Bio.Blast import NCBIWWW
import urllib

import shutil
import os.path

def get_genome_zone(start,stop,filename):
    Entrez.email = "duartelpmacau@hotmail.com"
    handle = Entrez.efetch(db="nucleotide", rettype="gb", retmode="text", seq_start=start, seq_stop=stop, id="15638995")
    file=open(filename,"w")#creating a GenBank file
    file.write(handle.read())
    file.close()
    handle.close()
    #moving the file to another directory
    path=os.getcwd()
    src = path+os.sep+filename #source folder
    dst = path+os.sep+"res"+os.sep#destination folder
    if not (os.path.exists(dst)):
        os.makedirs(dst)
    shutil.move(src, dst)#moving file to another directory
    record = SeqIO.read(dst+filename, "genbank") 
    return record
    
    
def removerGenBank():
    path=os.getcwd()+os.sep+"res"+os.sep
    shutil.rmtree(path)
  
def blast(GI_numb,filename):
    result_handle = NCBIWWW.qblast("blastp", "nr", GI_numb)
    save_file = open(filename, "w")
    save_file.write(result_handle.read())
    save_file.close()
    result_handle.close()
    #moving the file to another directory
    path=os.getcwd()
    src = path+os.sep+filename #source folder #source folder
    dst = path+os.sep+"res"+os.sep#destination folder
    shutil.move(src, dst)
    
def myblast(GI_numb,filename):
    url = 'http://www.uniprot.org/mapping/'
    params = {
    'from':'P_GI',
    'to':'ID',
    'format':'tab',
    'query':GI_numb
    }
    
    
    data = urllib.parse.urlencode(params)
    binario = data.encode("ASCII")
    request = urllib.request.Request(url,binario)
    contact = "duartelpmacau@hotmail.com" # Please set your email address here to help us debug in case of problems.
    request.add_header('User-Agent', 'Python %s' % contact)
    response = urllib.request.urlopen(request)
    html = "".join(map(chr,response.read()))
    
    save_file = open(filename, "w")
    save_file.write(html)
    save_file.close()
    #moving the file to another directory
    path=os.getcwd()
    src = path+os.sep+filename #source folder #source folder
    dst = path+os.sep+"res"+os.sep#destination folder
    shutil.move(src, dst)    
    
    
if __name__ == "__main__":
    removerGenBank()
    record=get_genome_zone(1015256,1137818,"genbank.gb")
    #print(record)
    myblast(499184939,"proteina.txt")
    for i in range(len(record.features)):
        feature_i = record.features[i]
#==============================================================================
#         if feature_i.type == "gene":
#             print(feature_i.qualifiers.get("db_xref"))
#             geneID =
#==============================================================================
    