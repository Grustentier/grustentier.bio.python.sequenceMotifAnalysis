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
 __   __               _          _      _           __  __         _     _    __     ___              _   _     _                           
 \ \ / /  __ _   _ _  (_)  __ _  | |__  | |  ___    |  \/  |  ___  | |_  (_)  / _|   | _ \  ___   ___ (_) | |_  (_)  ___   _ _    ___        
  \ V /  / _` | | '_| | | / _` | | '_ \ | | / -_)   | |\/| | / _ \ |  _| | | |  _|   |  _/ / _ \ (_-< | | |  _| | | / _ \ | ' \  (_-<        
   \_/   \__,_| |_|   |_| \__,_| |_.__/ |_| \___|   |_|  |_| \___/  \__| |_| |_|     |_|   \___/ /__/ |_|  \__| |_| \___/ |_||_| /__/        
                                                                                                                                             
  _                ___                     _                  _     _                                                                        
 | |__   _  _     / __|  _ _   _  _   ___ | |_   ___   _ _   | |_  (_)  ___   _ _                                                            
 | '_ \ | || |   | (_ | | '_| | || | (_-< |  _| / -_) | ' \  |  _| | | / -_) | '_|                                                           
 |_.__/  \_, |    \___| |_|    \_,_| /__/  \__| \___| |_||_|  \__| |_| \___| |_|                                                             
         |__/                                                                                                                                


""")

'''
Analyzing the variable positions of given sequence motifs with regEx like XYn.
X represents the starting amino acid and Y the ending by n-1 variable positions.
For example, a GG4 motif with n-1 is represented by three variable x positions (GxxxG).
These variable positions have to be considered statistically within the topologies transmembrane (tm), non-transmembrane (ntm) and transition (trans).

Steffen Grunert, Florian Heinke, Dirk Labudde, 
"Structure Topology Prediction of Discriminative Sequence Motifs in Membrane Proteins with Domains of Unknown Functions", 
Structural Biology, vol. 2013, Article ID 249234, 10 pages, 2013. https://doi.org/10.1155/2013/249234
'''

import os
import re
import math
import numpy
import pandas
import seaborn
import argparse
from cv2 import data
from scipy import stats
import matplotlib.pyplot as plt 
from sklearn.manifold import MDS
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from matplotlib.colors import LinearSegmentedColormap

parser = argparse.ArgumentParser(description='Code for analysis of variable sequence motif positions  for different topologies.')
parser.add_argument('--fasta_input_dir', default='testdata'+os.sep+'fasta'+os.sep+'rhodopsins', help='Path to the input dir including fasta files.')
parser.add_argument('--tmhmm_input_dir', default='testdata'+os.sep+'tmhmm'+os.sep+'rhodopsins', help='Path to the input dir including tmhmm files, generated extensive and  with no graphics.')
#parser.add_argument('--export_dir', default='', type=str, help='The export/output directory for exporting heatmap')
parser.add_argument('--max_variable_positions', default=9, type=int, help='The export/output directory')
parser.add_argument('--sort_by', default=None, type=str, help='Sorting amin acid presentation by properties like: alphabetical or hydrophob.')
parser.add_argument('--regex', default='PG10,LF10,PG9,LF9,VF8,LF8,GY8,GA7,AG7,AA7,GG7,LY6,VG6,SA6,PG6,AL6,PG5,GS5,LG5,AG5,GN4,IV4,IL4,GS4,GG4,SG4,VL4,AS4,GA4,AG4,SA3,AA3,GL3', type=str, help='Comma separted REGEXES like XXn representing a starting and a ending aminoacid by n - 1 variable position between both X.')
#parser.add_argument('--regex', default='GG4,SG4,AA3', type=str, help='Comma separted REGEXES like XXn representing a starting and a ending aminoacid by n - 1 variable position between both X.')
 
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
     
def createHeatMaps(dataFrames): 
    fig, a = plt.subplots(nrows=len(dataFrames.keys()))
    fig.subplots_adjust(wspace=0.01)
    seaborn.set(font_scale=0.5) 
    index = 0
    
    for topology in dataFrames.keys(): 
        if len(dataFrames[topology]["dataFrame"]) > 0:
            custom_color_map = LinearSegmentedColormap.from_list(name='custom_navy',colors=[(0/255, 0/255, 255/255),(255/255, 0/255, 0/255)])
            heatmap = seaborn.heatmap(dataFrames[topology]["dataFrame"],ax=a[index], cmap=custom_color_map, cbar=True, annot=False, annot_kws={"size": 1},xticklabels=AMINO_ACIDS_ONE_LETTER_CODE,yticklabels=dataFrames[topology]["yLabels"])
            heatmap.set_yticklabels(heatmap.get_yticklabels(), rotation=0, fontsize = 6)  
        index = index + 1     
    
    plt.show()
        
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

lastProgressValue = None 
def printProgress(steps,maximum):
    output = ""
    maxSteps2Console = 20
    for _ in range(0,int((steps/maximum)*maxSteps2Console)):output +="."
    value = int(round((steps/maximum)*100,0))
    global lastProgressValue
    if lastProgressValue != value:
        print("["+output+"]", str(value)+"%")
        lastProgressValue = value 
    
def getHeatmapDataFrames(possibleMotifs):
    dataFrames = {}
    
    for topology in possibleMotifs.keys():
        dataFrames[topology]={"yLabels":[],"dataFrame":[]}  
        for regEx in possibleMotifs[topology]: 
            registrations = []
            
            for _ in range(0,int(regEx[2:])-1):
                registrations.append({})
                for letter in AMINO_ACIDS_ONE_LETTER_CODE:                    
                    registrations[-1][letter] = 0                    
                #y_axis_labels.append(regEx+"-"+str(varPosIndex+1)+"-"+topology)
                rowLabel = regEx+"-"+topology
                if rowLabel not in dataFrames[topology]["yLabels"]:
                    dataFrames[topology]["yLabels"].append(rowLabel) 
                else:dataFrames[topology]["yLabels"].append("")
             
            for motifSequence in possibleMotifs[topology][regEx]:
                varPosOfMotif = motifSequence[1:-1] 
                for pos in range(0,len(varPosOfMotif)):
                    if varPosOfMotif[pos] in AMINO_ACIDS_ONE_LETTER_CODE:
                        registrations[pos][varPosOfMotif[pos]]+=1
                    
            if len(possibleMotifs[topology][regEx]) > 0:
                for pos in range(0,len(registrations)): 
                    dataFrameContent = []
                    for aminoAcid in AMINO_ACIDS_ONE_LETTER_CODE: 
                        #quotient = registrations[pos][aminoAcid] / len(possibleMotifs[topology][regEx])
                        quotient = registrations[pos][aminoAcid] / AMINO_ACID_OCCURRENCES[pos]
                        if quotient != 0.0:quotient = math.log(quotient)
                        dataFrameContent.append(quotient)
                    dataFrames[topology]["dataFrame"].append(dataFrameContent)
                    
    return dataFrames

def getPositionSpecificStatistics(possibleMotifs):
    positionSpecificStatistics= [] 
    
    for topology in possibleMotifs.keys():         
        for regEx in possibleMotifs[topology]:             
            regExsStatistic = []             
            
            for pos in range(0,int(regEx[2:])-1): 
                regExsStatistic.append({"regEx":regEx,"position":pos,"topology":topology,"occurences":[0 for _ in AMINO_ACIDS_ONE_LETTER_CODE],"occurence-ratios":[0.0 for _ in AMINO_ACIDS_ONE_LETTER_CODE]}) 
             
            for motifSequence in possibleMotifs[topology][regEx]:
                varPosOfMotif = motifSequence[1:-1] 
                for pos in range(0,len(varPosOfMotif)):
                    if varPosOfMotif[pos] in AMINO_ACIDS_ONE_LETTER_CODE:
                        indexInArray = AMINO_ACIDS_ONE_LETTER_CODE.index(varPosOfMotif[pos])
                        regExsStatistic[pos]["occurences"][indexInArray]+=1
                        
            for currentRegExsStatistic in regExsStatistic:
                for pos in range(0, len(currentRegExsStatistic["occurence-ratios"])):                    
                    '''
                        Within the corresponding research article log(P(a|pos|topology)/P(a|nature))...
                        Here a simpler information base approach. 
                    '''                        
                    #quotient = currentRegExsStatistic["occurences"][pos] / len(possibleMotifs[topology][regEx])
                    quotient = currentRegExsStatistic["occurences"][pos] / AMINO_ACID_OCCURRENCES[pos]                    
                    if quotient != 0.0:quotient = math.log(quotient)
                    currentRegExsStatistic["occurence-ratios"][pos] = quotient
             
            positionSpecificStatistics.extend(regExsStatistic)     
            
    return positionSpecificStatistics

def getSpearmanDistanceMatrix(positionSpecificStatistics):
    data = {}
    
    for i in range(0,len(positionSpecificStatistics)):
        comatrix = [] 
        for j in range(0,len(positionSpecificStatistics)):
            correlation,_ = stats.spearmanr(positionSpecificStatistics[i]["occurence-ratios"], positionSpecificStatistics[j]["occurence-ratios"])
            comatrix.append(correlation)   
        dataKey = positionSpecificStatistics[i]["regEx"]+"-"+str(positionSpecificStatistics[i]["position"])
        data[dataKey] = comatrix   
        
    data["topology"] = [positionSpecificStatistics[i]["topology"] for i in range(0,len(positionSpecificStatistics))]
    return pandas.DataFrame(data)

def cluster(dataFrames,cluster = 3):
    fig, a = plt.subplots(nrows=4,ncols=len(dataFrames))
    fig.subplots_adjust(wspace=0.2,hspace=0.4)
    seaborn.set(font_scale=0.5)     
    col = 0
    row = 0
    
    targets = ['tm', 'ntm', 'trans']
    colors = ['r', 'g', 'b', 'c', 'm']
        
    for dataFrame in dataFrames: 
        features = dataFrame.columns[0:-1]        
        x = dataFrame.loc[:, features].values
        #y = dataFrame.loc[:,['topology']].values
        
        ''' PCA scatter plot '''
        pca = PCA(2)         
        pcaComponents = pca.fit_transform(x)
        pcaComponentsDataFrame = pandas.DataFrame(data = pcaComponents, columns = ['pc1', 'pc2'])
        finalPCADataFrame = pandas.concat([pcaComponentsDataFrame, dataFrame[['topology']]], axis = 1)
        for target, color in zip(targets,colors):
            indicesToKeep = finalPCADataFrame['topology'] == target
            a[row,col].scatter(finalPCADataFrame.loc[indicesToKeep, 'pc1'], finalPCADataFrame.loc[indicesToKeep, 'pc2'], c = color, s = 50)
        a[row,col].set_xlabel('Principal Component 1', fontsize = 6)
        a[row,col].set_ylabel('Principal Component 2', fontsize = 6)
        a[row,col].set_title('PCA')
        a[row,col].legend(targets)
        row = row + 1
        
        ''' KMeans scatter plot based on PCA components '''
        kmeans = KMeans(n_clusters= cluster)        
        labels = kmeans.fit_predict(pcaComponents)
        centroids = kmeans.cluster_centers_
        classes = numpy.unique(labels) 
                
        for clazz in labels:         
            filtered_label = pcaComponents[labels == clazz]
            a[row,col].scatter(filtered_label[:,0], filtered_label[:,1],color=colors[list(classes).index(clazz)])
        
        a[row,col].scatter(centroids[:,0] , centroids[:,1] , s = 80, color = 'k')
        a[row,col].set_title('kmeans based on PCA components')
        a[row,col].legend(targets)
        row = row + 1
         
        ''' MDS scatter plot '''
        mds = MDS(n_components=2,metric=True,n_init=4,max_iter=300,verbose=0,eps=0.001,n_jobs=None,random_state=42,dissimilarity='euclidean')
        mdsComponents = mds.fit_transform(x)   
        mdsComponentsDataFrame = pandas.DataFrame(data = mdsComponents, columns = ['pc1', 'pc2'])
        finalMDSDataFrame = pandas.concat([mdsComponentsDataFrame, dataFrame[['topology']]], axis = 1)
        for target, color in zip(targets,colors):
            indicesToKeep = finalMDSDataFrame['topology'] == target
            a[row,col].scatter(finalMDSDataFrame.loc[indicesToKeep, 'pc1'], finalMDSDataFrame.loc[indicesToKeep, 'pc2'], c = color, s = 50)        
        a[row,col].set_xlabel('Principal Component 1', fontsize = 6)
        a[row,col].set_ylabel('Principal Component 2', fontsize = 6)
        a[row,col].set_title('MDS')
        a[row,col].legend(targets)
        row = row + 1
        
        ''' KMeans scatter plot based on PCA components '''
        kmeans = KMeans(n_clusters= cluster)        
        labels = kmeans.fit_predict(mdsComponents)
        centroids = kmeans.cluster_centers_
        classes = numpy.unique(labels) 
                
        for clazz in labels:         
            filtered_label = mdsComponents[labels == clazz]
            a[row,col].scatter(filtered_label[:,0], filtered_label[:,1],color=colors[list(classes).index(clazz)])
        
        a[row,col].scatter(centroids[:,0] , centroids[:,1] , s = 80, color = 'k')
        a[row,col].set_title('kmeans based on MDS components')
        a[row,col].legend(targets) 
        
        row = 0
        col = 1
        
    plt.show()
    
def getAminoAcidOccurrencesinNature(fastaData):
    occurrences = [0 for _ in AMINO_ACIDS_ONE_LETTER_CODE]
    
    for data in fastaData: 
        for letter in data["sequence"]: 
            if str(letter).upper() in AMINO_ACIDS_ONE_LETTER_CODE: 
                occurrences[AMINO_ACIDS_ONE_LETTER_CODE.index(str(letter).upper())]+=1
                
    return occurrences        

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
     
    maxSteps = 7   
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
    
    print("Determining amino acid occurences in nature")
    AMINO_ACID_OCCURRENCES = getAminoAcidOccurrencesinNature(fastaData)
    
    print("Generating some statistics...")
    positionSpecificStatistics = getPositionSpecificStatistics(possibleMotifs) 
    printProgress(4,maxSteps) 
      
    print("Generating some distance matrices...")    
   
    df1 = {AMINO_ACIDS_ONE_LETTER_CODE[aaIndex]:[positionSpecificStatistics[i]["occurence-ratios"][aaIndex] for i in range(0,len(positionSpecificStatistics))] for aaIndex in range(0,len(AMINO_ACIDS_ONE_LETTER_CODE))}
    df1["topology"] = [positionSpecificStatistics[i]["topology"] for i in range(0,len(positionSpecificStatistics))]
    printProgress(5,maxSteps)
    dataFrames = [pandas.DataFrame(df1),getSpearmanDistanceMatrix(positionSpecificStatistics)]
    
    print("Clustering...")    
    cluster(dataFrames)
    printProgress(6,maxSteps) 
    
    print("Generating heatmap...")  
    dataFrames = getHeatmapDataFrames(possibleMotifs)         
    printProgress(7,maxSteps) 
      
    print("FINISHED...")
    createHeatMaps(dataFrames)