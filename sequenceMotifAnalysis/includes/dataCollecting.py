'''
Created on 07.11.2021

@author: grustentier
'''

from sequenceMotifAnalysis.includes.stringUtils import cleanString


def collectFastaData(fastaFilePaths):
    fastaData = []

    for fastaFilePath in fastaFilePaths:
        data = getFastaData(fastaFilePath)
        if data is not None:
            fastaData.extend(data)
          
    return fastaData


def collectTmhmmData(tmhmmFilePaths):
    tmhmmData = []
    
    for tmhmmFilePath in tmhmmFilePaths:
        data = getTmhmmData(tmhmmFilePath)
        if data is not None:
            tmhmmData.extend(data)
            
    return tmhmmData 


def findCorrespondingTmhmmData(fasta_id, tmhmmData):    
    for data in tmhmmData: 
        tmhmm_id = data["id"]  
        if str(fasta_id).lower() == str(tmhmm_id).lower():
            return data
        
    return None


def findCorrespondingFastaData(tmhmm_id, fastaData):    
    for data in fastaData: 
        fasta_id = data["id"]  
        if str(fasta_id).lower() == str(tmhmm_id).lower():
            return data
        
    return None


def getFastaData(filePath):
    file = open(filePath, 'r')
    lines = file.readlines()    
    fastaData = []

    for line in lines:
        line = str(line).replace("\n", "")
        if str(line).startswith(">"):
            fastaData.append({"id":line[1:], "sequence":""}) 
        else:
            fastaData[-1]["sequence"] += line
    
    file.close() 
    return fastaData


def getTmhmmData(filePath):
    file = open(filePath, 'r')
    lines = file.readlines()    
    tmhmmData = []
    TMHMM_PROTEIN = None
    
    for line in lines:
        line = str(line).replace("\n", "")
        if str(line).startswith("<pre>"):            
            TMHMM_PROTEIN = {"id":str(line).split(" ")[1], "areas":[]}
            
        if str(line).startswith("</pre>"):
            tmhmmData.append(TMHMM_PROTEIN)
            TMHMM_PROTEIN = None
            
        if TMHMM_PROTEIN is not None and str(line).startswith(TMHMM_PROTEIN["id"]):
            line = cleanString(line)
            split = str(line).strip(" ").split("-")
            TMHMM_PROTEIN["areas"].append({"topology":str(split[-3]).lower(), "from":split[-2], "to":split[-1]})
    
    file.close() 
    return tmhmmData
