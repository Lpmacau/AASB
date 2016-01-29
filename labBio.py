# -*- coding: utf-8 -*-
"""
Created on Tue Jan 26 15:33:18 2016

@author: Duarte
"""

from Bio import SeqIO
from Bio import Entrez
from Bio.Blast import NCBIWWW
import urllib
import xml.etree.ElementTree as ET
import time

import shutil
import os.path

def get_genome_zone(GI_Number,start,stop,filename,tagRes,path):
    Entrez.email = "duartelpmacau@hotmail.com"
    handle = Entrez.efetch(db="nucleotide", rettype="gb", retmode="text", seq_start=start, seq_stop=stop, id=str(GI_Number))
    file=open(path+filename,"w")#creating a GenBank file
    file.write(handle.read())
    file.close()
    handle.close()
    
    record = SeqIO.read(path+filename, "genbank") 
    return record
    
    
def removerGenBank():
    path=os.getcwd()+os.sep+"res"+os.sep
    if os.path.exists(path):
        shutil.rmtree(path)
    
def myblast(GI_numb,filename,tagRes,path):
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
    try:
        response = urllib.request.urlopen(request)
        html = "".join(map(chr,response.read()))
       
        save_file = open(path+filename, "w")
        save_file.write(html)
        save_file.close()
        return 1
    except:
        return -1
      
    
def getXML(UniProtID,tagRes,path):
    filename = UniProtID+".xml"
    url = "http://www.uniprot.org/uniprot/"+filename
    
    try:
        request = urllib.request.Request(url)    
        response = urllib.request.urlopen(request)
        data = response.read()
        
        data_str = "".join(map(chr,data))
        save_file = open(path+filename, "w")
        save_file.write(data_str)
        save_file.close()
        return 1
    except: 
        return -1
   
def recursiveXMLParse(XMLChild):
    # Aceder as informações (escolher as que estao na tabela)
    #print(XMLChild.tag,XMLChild.attrib)
    for child in XMLChild:
        recursiveXMLParse(child)
    
def parseXML(UniProtID,tagRes,path):
    filename = path+UniProtID+".xml"
    tree = ET.parse(filename)
    root = tree.getroot()
    for child in root:
        recursiveXMLParse(child)
    
    
def getUniProtID(GI_Number,tagRes,path):
    filename = path+GI_Number+".txt"
    file = open(filename,"r")
    
    for line in file.readlines():
        #procurar por TREPA e fazer return do ID, se nao houvre, retornar -1
        strs = line.split()
        protID = strs[1]
        if "TREPA" in protID:
            uniProtID = protID.split("_")
            file.close()
            print(uniProtID[0])
            return uniProtID[0]
    file.close()
    return -1

if __name__ == "__main__":
    tagRes = time.strftime("%H")+ time.strftime("%M") + time.strftime("%S")
    path = os.getcwd()
    os.mkdir("res-"+tagRes)
    if tagRes not in path:
        os.chdir(path+os.sep+"res-"+tagRes)
        path = path+os.sep+"res-"+tagRes+os.sep
    else:
        path = path+os.sep
    
    record=get_genome_zone(15638995,1015256,1137818,"genbank.gb",tagRes,path)
      
    file = open("res.csv","w")
    file.write("Start,Stop,Strand,GeneID,\"Locus Tag\",\"Protein Product\",\"#Aminoacids\",\"Product (Protein Name)\",UniProtID,\"Subcellular Location\",\"GeneOntology Terms\",\"GeneOntology Identifiers\"\n")
    try:
        for i in range(len(record.features)):
            feature_i = record.features[i]
            if feature_i.type == "CDS":
                start = feature_i.location.start
                end = feature_i.location.end
                strandNum = feature_i.location.strand
                if strandNum == 1:
                    strand = "+"
                else:
                    strand = "-"
                #print(feature_i)
                
                
                
                tag = feature_i.qualifiers.get("db_xref")[0]
                try:
                    tag2 =feature_i.qualifiers.get("db_xref")[1]
                except:
                    tag2 = ""
                if "GeneID" in tag:
                    GeneID = tag[7:]
                else:
                    if tag2 != "":
                        GeneID = tag2[7:]
                try:
                    locus_tag = feature_i.qualifiers.get("locus_tag")[0]            
                    protein_id = feature_i.qualifiers.get("protein_id")[0]          
                    product = feature_i.qualifiers.get("product")[0]
                    nAmino = len(feature_i.qualifiers.get("translation")[0])
                    file.write("\""+str(start)+"\","+"\""+str(end)+"\","+"\""+strand+"\","+"\""+str(GeneID)+"\","+"\""+locus_tag+"\","+"\""+protein_id+"\","+"\""+str(nAmino)+"\","+"\""+product+"\n")
                except:
                    print("Não encontrei")
            
                if "GI" in tag :
                    GI_Number = tag[3:]
                    error = myblast(GI_Number,str(GI_Number)+".txt",tagRes,path)
                    if error != -1:
                        print("GI_Number "+GI_Number)                
                        uniProtID = getUniProtID(GI_Number,tagRes,path)
                        print("Uniprot "+str(uniProtID))
                        if uniProtID != -1:
                            if getXML(uniProtID,tagRes,path) != -1:
                                parseXML(uniProtID,tagRes,path)
                else:
                    if "GI" in tag2:
                        GI_Number = tag2[3:]
                        error = myblast(GI_Number,str(GI_Number)+".txt",tagRes,path)
                        if error != -1:
                            print("GI_Number "+GI_Number)                
                            uniProtID = getUniProtID(GI_Number,tagRes,path)
                            print("Uniprot "+str(uniProtID))
                            if uniProtID != -1:
                                if getXML(uniProtID,tagRes,path) != -1:
                                    parseXML(uniProtID,tagRes,path)
    except:
        print("-")
    file.close()       
    