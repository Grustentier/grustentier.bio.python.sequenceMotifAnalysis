'''
Created on 10.10.2021

@author: grustentier
'''

print("""


  ___                                                 __  __         _     _    __       _                   _               _               
 / __|  ___   __ _   _  _   ___   _ _    __   ___    |  \/  |  ___  | |_  (_)  / _|     /_\    _ _    __ _  | |  _  _   ___ (_)  _ _    __ _ 
 \__ \ / -_) / _` | | || | / -_) | ' \  / _| / -_)   | |\/| | / _ \ |  _| | | |  _|    / _ \  | ' \  / _` | | | | || | |_ / | | | ' \  / _` |
 |___/ \___| \__, |  \_,_| \___| |_||_| \__| \___|   |_|  |_| \___/  \__| |_| |_|     /_/ \_\ |_||_| \__,_| |_|  \_, | /__| |_| |_||_| \__, |
                |_|                                                                                              |__/                  |___/ 
  _                ___                     _                  _     _                                                                        
 | |__   _  _     / __|  _ _   _  _   ___ | |_   ___   _ _   | |_  (_)  ___   _ _                                                            
 | '_ \ | || |   | (_ | | '_| | || | (_-< |  _| / -_) | ' \  |  _| | | / -_) | '_|                                                           
 |_.__/  \_, |    \___| |_|    \_,_| /__/  \__| \___| |_||_|  \__| |_| \___| |_|                                                             
         |__/                                                                                                                                

""")

'''
Analyzing the variable positions of given sequence motifs with a regEx XYn.
X represents the starting amino acid and Y the ending by n-1 variable positions.
For example, a GG4 motif with n-1 is represented by three variable positions (GxxxG).
These variable positions have to be considered statistically within the topologies transmembrane (tm), non-transmembrane (ntm) and transition (trans).
'''

import os
import re
import math
import seaborn
import argparse
from cv2 import data 
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

parser = argparse.ArgumentParser(description='Code for analysis of variable sequence motif positions  for different topologies.')
parser.add_argument('--fasta_input_dir', default='./testdata/fasta', help='Path to the input dir including fasta files.')
parser.add_argument('--tmhmm_input_dir', default='./testdata/tmhmm', help='Path to the input dir including tmhmm files, generated extensive and  with no graphics.')
#parser.add_argument('--export_dir', default='', type=str, help='The export/output directory for exporting heatmap')
parser.add_argument('--max_variable_positions', default=9, type=int, help='The export/output directory')
parser.add_argument('--sort_by', default=None, type=str, help='Sorting amin acid presentation by properties like: alphabetical or hydrophob.')
parser.add_argument('--regex', default='PG10,LF10,PG9,LF9,VF8,LF8,GY8,GA7,AG7,AA7,GG7,LY6,VG6,SA6,PG6,AL6,PG5,GS5,LG5,AG5,GN4,IV4,IL4,GS4,GG4,SG4,VL4,AS4,GA4,AG4,SA3,AA3,GL3', type=str, help='Comma separted REGEXES like XXn representing a starting and a ending aminoacid by n - 1 variable position between both X.')
 
arguments = parser.parse_args() 

AMINO_ACIDS_ONE_LETTER_CODE_ALPHABETICAL = ["A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y"]
AMINO_ACIDS_ONE_LETTER_CODE_HYDROPHOB =    ["A","C","F","I","L","M","P","T","V","W","Y","D","E","G","H","K","N","Q","R","S"]
AMINO_ACIDS_ONE_LETTER_CODE_SMALL =        ["A","C","D","G","N","P","V","S","T","E","F","H","I","K","L","M","Q","R","W","Y"]
AMINO_ACIDS_ONE_LETTER_CODE_PLOAR =        ["C","D","E","H","K","N","Q","R","S","T","W","Y","A","F","G","I","L","M","P","V"]
 
def cleanString(string,replacement='-'):
    a =  re.sub('[^a-zA-Z0-9.?]',replacement,string) 
    return re.sub(replacement+'+', replacement, a)

def collectFilePaths(directory):
    filePaths = []
    for root, _, files in os.walk(directory):  
        if root.endswith(os.sep) is False:root+=os.sep
        filePaths.extend([root + file for file in files])
    return filePaths

def getFastaData(filePath):
    file = open(filePath, 'r')
    lines = file.readlines()
    
    fastaData = []
    for line in lines:
        line = str(line).replace("\n","")
        if str(line).startswith(">"):
            fastaData.append({"id":line[1:],"sequence":""} ) 
        else:
            fastaData[-1]["sequence"]+=line
    
    file.close() 
    return fastaData

def getTmhmmData(filePath):
    file = open(filePath, 'r')
    lines = file.readlines()
    
    tmhmmData = []
    TMHMM_PROTEIN = None
    for line in lines:        
        line = str(line).replace("\n","")
        if str(line).startswith("<pre>"):            
            TMHMM_PROTEIN = {"id":str(line).split(" ")[1],"areas":[]}
        if str(line).startswith("</pre>"):
            tmhmmData.append(TMHMM_PROTEIN)
            TMHMM_PROTEIN = None
            
        if TMHMM_PROTEIN is not None and str(line).startswith(TMHMM_PROTEIN["id"]):
            line = cleanString(line)
            split = str(line).strip(" ").split("-")
            TMHMM_PROTEIN["areas"].append({"topology":str(split[-3]).lower(),"from":split[-2],"to":split[-1]})
    
    file.close() 
    return tmhmmData

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

def findCorrespondingTmhmmData(fasta_id,tmhmmData):    
    for data in tmhmmData: 
        tmhmm_id = data["id"]  
        if str(fasta_id).lower() == str(tmhmm_id).lower():
            return data
    return None

def getTopology(start,end,tmhmmData):  
    for area in tmhmmData["areas"]:
        start_tmhmm = int(area["from"]) - 1
        end_tmhmm = int(area["to"]) - 1 
        topology = area["topology"].strip("")          
        
        if start >= start_tmhmm and start <= end_tmhmm and end >= start_tmhmm and end <= end_tmhmm:
            if topology == "outside" or topology == "inside":
                return "ntm"
            elif topology == "tmhelix": 
                return "tm"
            else: return "new_unknow_topology"
        if start >= start_tmhmm and start <= end_tmhmm and end > end_tmhmm:
            return "trans"
        if start < start_tmhmm and end >= start_tmhmm and end <= end_tmhmm:
            return "trans"
        
    return "unknown"
    
def getPossibleMotifs(fastaData,tmhmmData):
    possibleMotifs = {"tm":{},"ntm":{},"trans":{}}
    for data in fastaData: 
        sequence = data["sequence"] 
        fasta_id = data["id"]
        current_tmhmmData = findCorrespondingTmhmmData(fasta_id,tmhmmData)
        if current_tmhmmData is None:print("No tmhmm data for to fasta id:",fasta_id)  
        
        for i in range(0,len(sequence)):
            for j in range(3,arguments.max_variable_positions+1):
                if (i + j) >= len(sequence):continue
                regEx = sequence[i]+sequence[i+j]+str(j)                 
                if regEx not in REGEXES: continue                
                motifSeq = sequence[i:i+j+1]                                
                if len(motifSeq)>3:                         
                    topology2Assign = getTopology(i,i+j,current_tmhmmData)
                    if topology2Assign not in possibleMotifs.keys():possibleMotifs[topology2Assign]={}
                    if regEx not in possibleMotifs[topology2Assign].keys():possibleMotifs[topology2Assign][regEx] = [] 
                    possibleMotifs[topology2Assign][regEx].append(motifSeq) 
         
    return possibleMotifs 

def createClusterMap(dataFrame,x_axis_labels,y_axis_labels): 
    if len(dataFrame) > 0:
        custom_color_map = LinearSegmentedColormap.from_list(name='custom_navy',colors=[(0/255, 0/255, 255/255),(255/255, 0/255, 0/255)])
        seaborn.set(font_scale=0.1) 
        heatmap = seaborn.heatmap(dataFrame, cmap=custom_color_map, cbar=True, annot=False, annot_kws={"size": 2},xticklabels=x_axis_labels,yticklabels=y_axis_labels)
        heatmap.set_yticklabels(heatmap.get_yticklabels(), rotation=0)  
        
        plt.show() 
        #plt.savefig("TODO/heatmap.png",dpi=600)
        
def getAminoAcidLetters():
    aminoAcidLetters = None
    if arguments.sort_by == "hydrophob":
        aminoAcidLetters = AMINO_ACIDS_ONE_LETTER_CODE_HYDROPHOB
    elif arguments.sort_by == "small":
        aminoAcidLetters = AMINO_ACIDS_ONE_LETTER_CODE_SMALL
    elif arguments.sort_by == "polar":
        aminoAcidLetters = AMINO_ACIDS_ONE_LETTER_CODE_PLOAR
    else:aminoAcidLetters = AMINO_ACIDS_ONE_LETTER_CODE_ALPHABETICAL        
    
    return aminoAcidLetters

def printProgress(steps,maximum):
    output = ""
    maxSteps2Console = 20
    for _ in range(0,int((steps/maximum)*maxSteps2Console)):
        output +="."
    print("["+output+"]", str(int(round((steps/maximum)*100,0)))+"%") 

if __name__ == "__main__":  
    fastaFilePaths = collectFilePaths(arguments.fasta_input_dir)
    assert len(fastaFilePaths)>0,"No fasta files have been found!!!"
    tmhmmFilePaths = collectFilePaths(arguments.tmhmm_input_dir)
    assert len(tmhmmFilePaths)>0,"No tmhmm files have been found!!!"
    assert len(arguments.regex)>0,"No regex found have been found!!!"
    
    REGEXES = []
    if str(arguments.regex).find(",") != -1:
        REGEXES = str(arguments.regex).upper().strip(" ").split(",")
    else:
        REGEXES.append(str(arguments.regex).upper().strip(" "))
        
    AMINO_ACIDS_ONE_LETTER_CODE = getAminoAcidLetters()
     
    maxSteps = 4   
    printProgress(0,maxSteps)
    print("Parsing and collecing fast files...")
    fastaData = collectFastaData(fastaFilePaths)   
    printProgress(1,maxSteps) 
    print("Parsing and collecing tmhmm files...")
    tmhmmData = collectTmhmmData(tmhmmFilePaths)
    printProgress(2,maxSteps)    
    print("Gathering motifs from sequences...")      
    possibleMotifs = getPossibleMotifs(fastaData,tmhmmData)
    printProgress(3,maxSteps)   
    
    print("Generating heatmap...")  
    dataFrame = []
    y_axis_labels = []  
    for topology in possibleMotifs.keys():
        for regEx in possibleMotifs[topology]: 
            registrations = []
            
            for varPosIndex in range(0,int(regEx[2:])-1):
                registrations.append({})
                for letter in AMINO_ACIDS_ONE_LETTER_CODE:                    
                    registrations[-1][letter] = 0                    
                #y_axis_labels.append(regEx+"-"+str(varPosIndex+1)+"-"+topology)
                rowLabel = regEx+"-"+topology
                if rowLabel not in y_axis_labels:
                    y_axis_labels.append(rowLabel) 
                else:y_axis_labels.append("")
             
            for motif in possibleMotifs[topology][regEx]:
                varPosOfMotif = motif[1:-1] 
                for pos in range(0,len(varPosOfMotif)):
                    if varPosOfMotif[pos] in AMINO_ACIDS_ONE_LETTER_CODE:
                        registrations[pos][varPosOfMotif[pos]]+=1
                    
            if len(possibleMotifs[topology][regEx]) > 0:
                for reg in registrations: 
                    dataFrameContent = []
                    for aminoAcid in AMINO_ACIDS_ONE_LETTER_CODE:
                        z = reg[aminoAcid]
                        n = len(possibleMotifs[topology][regEx])
                        quotient = z/n
                        #if quotient != 0.0:quotient = math.log1p(quotient)
                        dataFrameContent.append(quotient)
                    dataFrame.append(dataFrameContent)
    
    printProgress(4,maxSteps)   
    createClusterMap(dataFrame,AMINO_ACIDS_ONE_LETTER_CODE,y_axis_labels) 
    