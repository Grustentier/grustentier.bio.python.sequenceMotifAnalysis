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
import sys
import math
import numpy
import pandas
import seaborn
import argparse
import xml.dom.minidom
import matplotlib.pyplot as plt
import xml.etree.ElementTree as ET
from cv2 import data
from scipy import stats
from sklearn.manifold import MDS
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA 
from matplotlib.colors import LinearSegmentedColormap

parser = argparse.ArgumentParser(description='Code for analysis of variable sequence motif positions  for different topologies.')
parser.add_argument('--fasta_input_dir', default='testdata'+os.sep+'fasta'+os.sep+'rhodopsins', help='Path to the input dir including fasta files.')
parser.add_argument('--tmhmm_input_dir', default='testdata'+os.sep+'tmhmm'+os.sep+'rhodopsins', help='Path to the input dir including tmhmm files, generated extensive and  with no graphics.')
parser.add_argument('--sort_by', default=None, type=str, help='Selection: alphabetical, hydrophob or None; Sorting amino acid presentation by properties.')
'''Gerstein pattern (regExes) '''
#parser.add_argument('--regexes', default='PG10,LF10,PG9,LF9,VF8,LF8,GY8,GA7,AG7,AA7,GG7,LY6,VG6,SA6,PG6,AL6,PG5,GS5,LG5,AG5,GN4,IV4,IL4,GS4,GG4,SG4,VL4,AS4,GA4,AG4,SA3,AA3,GL3', type=str, help='Comma separted REGEXES like XXn representing a starting and a ending aminoacid by n - 1 variable position between both X.')
''' Other regExes after running analyseConsecutiveMotifs.py ''' 
parser.add_argument('--regexes', default='WP5, LL7, YL3, DF7, LL5, YP8, DA3, KF3, TL3, LL4, LL6, WL7, LL3, DK4, DP5, WY3, PW3, TL4, YA7, RP9, LA4, WL6, KG4, KL7, WT4, PL6, PL3, DV6, RD3, WT3, IL5, PL8, YT7, YL9, GL3, DG8, DP6, AF4, FL4, GF4, GA3, YL6, LF8, DM3, IG3, WF5, MG4, LA3, TL7, LY4, RT8, KL8, YT6, TL9, MT3, GA4, LK5, LP3, DT4, DG7, DL8, AL8, DT5, YK8, GL4, LE4, LV7, WP4, WL8, DL7, RT7', type=str, help='Comma separted REGEXES like XXn representing a starting and a ending aminoacid by n - 1 variable position between both X.')
parser.add_argument('--amino_acid_abundance', default=None, type=str, help='Selection: pdbtm, fasta or None; pdbtm: amino acid abundance from currently solved transmembrane proteins from pdbtm: 6420 (alpha: 5899 , beta: 473, version: 2021-10-08, from Tusnady (http://pdbtm.enzim.hu/)); fasta: amino acid abundance from current fasta dataset; None: amino acid abundance from on https://en.wikipedia.org/wiki/Amino_acid;')

 
arguments = parser.parse_args() 

AMINO_ACIDS_ONE_LETTER_CODE_ALPHABETICAL = ["A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y"]
AMINO_ACIDS_ONE_LETTER_CODE_HYDROPHOB =    ["A","C","F","I","L","M","P","T","V","W","Y","D","E","G","H","K","N","Q","R","S"]
AMINO_ACIDS_ONE_LETTER_CODE_SMALL =        ["A","C","D","G","N","P","V","S","T","E","F","H","I","K","L","M","Q","R","W","Y"]
AMINO_ACIDS_ONE_LETTER_CODE_PLOAR =        ["C","D","E","H","K","N","Q","R","S","T","W","Y","A","F","G","I","L","M","P","V"]

def cleanString(string,replacement='-'):
    a =  re.sub('[^a-zA-Z0-9.?]',replacement,string) 
    return re.sub(replacement+'+', replacement, a)

def cluster(dataFrames,cluster = 3):
    fig, a = plt.subplots(nrows=4,ncols=len(dataFrames))
    fig.subplots_adjust(wspace=0.2,hspace=0.4)
    seaborn.set(font_scale=0.5)     
    col = 0
    row = 0
    
    targets = ['tm', 'ntm', 'trans']
    colors = ['#cc4c48', '#4ab3c9', '#8cb555', 'c', 'm']
        
    for index in range(0,len(dataFrames)):         
        dataFrame = dataFrames[index]
        if dataFrame.empty is True:continue
        
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
        
        printProgress(index,len(dataFrames)-1,"Creating cluster plots...")  
        
    plt.show()

def collectFilePaths(directory):
    filePaths = []
    for root, _, files in os.walk(directory):  
        if root.endswith(os.sep) is False:root+=os.sep
        filePaths.extend([root + file for file in files])
    return filePaths

def collectFastaData(fastaFilePaths):
    fastaData = []

    for index in range(0,len(fastaFilePaths)): 
        fastaFilePath = fastaFilePaths[index]
        data = getFastaData(fastaFilePath)
        if data is not None:
            fastaData.extend(data)  
                  
        printProgress(index,len(fastaFilePaths)-1,"Parsing and collecting fasta files...")
            
    return fastaData

def collectTmhmmData(tmhmmFilePaths):
    tmhmmData = []
    
    for index in range(0,len(tmhmmFilePaths)): 
        tmhmmFilePath = tmhmmFilePaths[index]
        data = getTmhmmData(tmhmmFilePath)
        if data is not None:
            tmhmmData.extend(data)
            
        printProgress(index,len(tmhmmFilePaths)-1,"Parsing and collecting tmhmm files...")
            
    return tmhmmData 

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
            
        printProgress(index,len(dataFrames.keys())-1,"Creating heatmap...") 
        index = index + 1     
    
    plt.show()
    
def exportWinners(motifWinners):
    motifTopologiesFilePath = '.'+os.sep+'inputdata'+os.sep+'motif-topologies.xml'
    if os.path.exists(motifTopologiesFilePath) is False:
        root = ET.Element('motifs')
        with open(motifTopologiesFilePath, "wb") as f:
            f.write(ET.tostring(root))
    
    tree=ET.parse(motifTopologiesFilePath)
    root=tree.getroot()
     
    index = 0
    for regEx in motifWinners.keys():
        winnerTopology = max(motifWinners[regEx].items(), key = lambda k : k[1])[0] 
         
        motifElement = root.find(".//*[@regEx='"+regEx+"']")
        if motifElement is None:
            motifElement = ET.Element('motif')
                 
        motifElement.set("regEx", regEx)
        motifElement.set("topology", winnerTopology)   
        root.append(motifElement)       
        
        printProgress(index,len(motifWinners.keys())-1,"Exporting topology winners...") 
        index = index + 1       
        
    with open(motifTopologiesFilePath, "wb") as f:
        f.write(ET.tostring(root))
    
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

def findCorrespondingTmhmmData(fasta_id,tmhmmData):    
    for data in tmhmmData: 
        tmhmm_id = data["id"]  
        if str(fasta_id).lower() == str(tmhmm_id).lower():
            return data
        
    return None

def getFastaData(filePath):
    file = open(filePath, 'r')
    lines = file.readlines()    
    fastaData = []

    for index in range(0,len(lines)):
        line = lines[index]
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
    
    for index in range(0,len(lines)):
        line = lines[index]
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
    minVarPos, maxVarPos = getMinMaxVarPos()

    for index in range(0,len(fastaData)):
        data = fastaData[index]
        sequence = data["sequence"] 
        fasta_id = data["id"]
        current_tmhmmData = findCorrespondingTmhmmData(fasta_id,tmhmmData)
        if current_tmhmmData is None:print("No tmhmm data for to fasta id:",fasta_id)  
        
        for i in range(0,len(sequence)):
            for j in range(minVarPos+1,maxVarPos+2):
                if (i + j) >= len(sequence):continue
                regEx = sequence[i]+sequence[i+j]+str(j)  
                if regEx not in REGEXES: continue    
                motifSeq = sequence[i:i+j+1]         
                if len(motifSeq)>3:                         
                    topology2Assign = getTopology(i,i+j,current_tmhmmData)
                    if topology2Assign not in possibleMotifs.keys():possibleMotifs[topology2Assign]={}
                    if regEx not in possibleMotifs[topology2Assign].keys():possibleMotifs[topology2Assign][regEx] = [] 
                    possibleMotifs[topology2Assign][regEx].append(motifSeq) 
                    
        printProgress(index,len(fastaData)-1,"Gathering motifs from sequences...")
         
    return possibleMotifs 
    
def getHeatmapDataFrames(possibleMotifs):
    dataFrames = {}
    
    index = 0
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
                    
        printProgress(index,len(possibleMotifs.keys())-1,"Collecting heatmap data...") 
        index = index + 1
                    
    return dataFrames

def getPositionSpecificStatistics(possibleMotifs):
    positionSpecificStatistics= [] 
    
    index = 0
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
            
        printProgress(index,len(possibleMotifs.keys())-1,"Determining position specific statistics...")  
        index = index + 1 
            
    return positionSpecificStatistics

def getDataFrameFromOccurrences(positionSpecificStatistics):
    data = {AMINO_ACIDS_ONE_LETTER_CODE[aaIndex]:[positionSpecificStatistics[i]["occurence-ratios"][aaIndex] for i in range(0,len(positionSpecificStatistics))] for aaIndex in range(0,len(AMINO_ACIDS_ONE_LETTER_CODE))}
    data["topology"] = [positionSpecificStatistics[i]["topology"] for i in range(0,len(positionSpecificStatistics))]
    return pandas.DataFrame(data)

def getDataFrameFromSpearman(positionSpecificStatistics):
    data = {}

    for i in range(0,len(positionSpecificStatistics)):
        comatrix = [] 
        for j in range(0,len(positionSpecificStatistics)):
            correlation,_ = stats.spearmanr(positionSpecificStatistics[i]["occurence-ratios"], positionSpecificStatistics[j]["occurence-ratios"])
            comatrix.append(correlation)   
        dataKey = positionSpecificStatistics[i]["regEx"]+"-"+str(positionSpecificStatistics[i]["position"])
        data[dataKey] = comatrix   
        
        printProgress(i,len(positionSpecificStatistics)-1,"Calculating spearman matrix")
        
    data["topology"] = [positionSpecificStatistics[i]["topology"] for i in range(0,len(positionSpecificStatistics))]
    return pandas.DataFrame(data)    
    
def getAminoAcidOccurrencesinNatureFromFasta(fastaData):
    print("... from current fastaData") 
    occ = [0 for _ in AMINO_ACIDS_ONE_LETTER_CODE]

    for index in range(0,len(fastaData)):
        data = fastaData[index]
        for letter in data["sequence"]: 
            if str(letter).upper() in AMINO_ACIDS_ONE_LETTER_CODE: 
                occ[AMINO_ACIDS_ONE_LETTER_CODE.index(str(letter).upper())]+=1
        
        printProgress(index,len(fastaData)-1,"Determining amino acid occurrences in nature from fasta data")
        
    return occ

def getAminoAcidOccurrencesinNatureFromWiki(): 
    print('Determining amino acid occurrences in nature based on on https://en.wikipedia.org/wiki/Amino_acid')
    wiki = {}
    wiki["W"] = 0.0125
    wiki["C"] = 0.0138
    wiki["H"] = 0.0226
    wiki["M"] = 0.0232
    wiki["Y"] = 0.0291
    wiki["F"] = 0.0387
    wiki["Q"] = 0.039
    wiki["N"] = 0.0393
    wiki["P"] = 0.0502
    wiki["K"] = 0.0519
    wiki["D"] = 0.0549
    wiki["I"] = 0.0549
    wiki["T"] = 0.0553
    wiki["R"] = 0.0578
    wiki["E"] = 0.0632
    wiki["V"] = 0.0703
    wiki["G"] = 0.0125
    wiki["S"] = 0.0714
    wiki["A"] = 0.0876
    wiki["L"] = 0.0968
    
    return [wiki[aa] for aa in AMINO_ACIDS_ONE_LETTER_CODE]

def getAminoAcidOccurrencesinNatureFromTusnady():  
    print("Determining amino acid occurrences in nature based on Tusnady (http://pdbtm.enzim.hu/)")
    try:       
        occurrences = [0 for _ in AMINO_ACIDS_ONE_LETTER_CODE]
        doc = xml.dom.minidom.parse('.'+os.sep+"inputdata"+os.sep+"pdbtmalpha.xml")
        sequences2One = ""
        
        for seq in doc.getElementsByTagName('SEQ'):
            unstripped = seq.firstChild.nodeValue
            stripped = re.sub('[^a-zA-Z]',"",unstripped)
            sequences2One+=stripped
            
        for letter in sequences2One: 
            if str(letter).upper() in AMINO_ACIDS_ONE_LETTER_CODE: 
                occurrences[AMINO_ACIDS_ONE_LETTER_CODE.index(str(letter).upper())]+=1
                
        for i in range(0,len(occurrences)):
            occurrences[i] = occurrences[i] / len(sequences2One)                
                
        return occurrences
    except:
        return getAminoAcidOccurrencesinNatureFromWiki()    
    
def getAminoAcidOccurrencesinNature(fastaData): 
    if str(arguments.amino_acid_abundance).lower() == "fasta":        
        return getAminoAcidOccurrencesinNatureFromFasta(fastaData)
    elif str(arguments.amino_acid_abundance).lower() == "pdbtm":        
        return getAminoAcidOccurrencesinNatureFromTusnady() 
    else:
        return getAminoAcidOccurrencesinNatureFromWiki()
    
def getMinMaxVarPos():
    minVarPos = 10000000000
    maxVarPos = 0
    
    for regEx in REGEXES:
        varPos = int(regEx[2:])-1
        if varPos < minVarPos: minVarPos = varPos  
            
        if varPos > maxVarPos: maxVarPos = varPos 
    return minVarPos,maxVarPos    

def getAminoAcidOccurrences2MotifPosition(regEx,position,positionSpecificStatistics):
    occurrences = []
    for stat in positionSpecificStatistics:
        if stat["regEx"] == regEx and stat["position"] == position: 
            return stat["occurence-ratios"]
    return occurrences

def determineWinners(possibleMotifs,positionSpecificStatistics):
    motifWinners = {}
    
    index = 0
    for t in possibleMotifs.keys():
        for regEx in possibleMotifs[t].keys():
            if regEx not in motifWinners.keys():motifWinners[regEx] = {}
            if t not in motifWinners[regEx].keys():motifWinners[regEx][t] = 0
            
            for motif in possibleMotifs[t][regEx]:
                varPosSeq = motif[1:-1]
                for pos in range(0,len(varPosSeq)):
                    letterAtPos = varPosSeq[pos] 
                    if letterAtPos not in AMINO_ACIDS_ONE_LETTER_CODE:continue
                    oc = getAminoAcidOccurrences2MotifPosition(regEx,pos,positionSpecificStatistics)
                    assert len(oc) > 0, "wrong assingment"
                    motifWinners[regEx][t]+=oc[AMINO_ACIDS_ONE_LETTER_CODE.index(letterAtPos)]     
                        
        printProgress(index,len(possibleMotifs.keys())-1,"Determining topology winners...") 
        index = index + 1                
                
    return motifWinners
    
def printProgress(steps,maximum,name="todo", bar_length = 20, width = 20):  
    percent = float(steps) / maximum
    arrow = '-' * int(round(percent*bar_length) - 1) + '>'
    spaces = ' ' * (bar_length - len(arrow))
    sys.stdout.write("\r{0: <{1}} : [{2}]{3}%".format(\
                     name, width, arrow + spaces, int(round(percent*100))))
    sys.stdout.flush()    
    
    if steps >= maximum:     
        sys.stdout.write('\n\n')  
    
if __name__ == "__main__":  
    fastaFilePaths = collectFilePaths(arguments.fasta_input_dir)
    assert len(fastaFilePaths)>0,"No fasta files have been found!!!"
    tmhmmFilePaths = collectFilePaths(arguments.tmhmm_input_dir)
    assert len(tmhmmFilePaths)>0,"No tmhmm files have been found!!!"
    assert len(arguments.regexes)>0,"No regex found have been found!!!"
    
    REGEXES = []
    if str(arguments.regexes).find(",") != -1:
        REGEXES = str(arguments.regexes).upper().replace(" ","").split(",")
    else:
        REGEXES.append(str(arguments.regexes).upper().replace(" ",""))        
        
    AMINO_ACIDS_ONE_LETTER_CODE = getAminoAcidLetters() 
    
    ''' Collecting fasta data '''
    fastaData = collectFastaData(fastaFilePaths)
    
    ''' Collecting tmhmm data  '''    
    tmhmmData = collectTmhmmData(tmhmmFilePaths)
    
    ''' Gathering motifs from sequences ''' 
    possibleMotifs = getPossibleMotifs(fastaData,tmhmmData)
    
    ''' Determining amino acid occurrences in nature from different sources'''  
    AMINO_ACID_OCCURRENCES = getAminoAcidOccurrencesinNature(fastaData)
    
    ''' Determining position specific statistics ''' 
    positionSpecificStatistics = getPositionSpecificStatistics(possibleMotifs) 
    
    ''' Determining topology winners ''' 
    motifWinners = determineWinners(possibleMotifs,positionSpecificStatistics)
    
    ''' Exporting topology winners '''  
    exportWinners(motifWinners)
    
    ''' Generating panda.DataFrames from position specific statistics '''      
    dataFrames = [getDataFrameFromOccurrences(positionSpecificStatistics),getDataFrameFromSpearman(positionSpecificStatistics)]
    
    ''' Clustering data and presenting within scatter plots '''
    cluster(dataFrames)
    
    ''' Generating dataFrame for heatmap''' 
    dataFrames = getHeatmapDataFrames(possibleMotifs)
    
    ''' Creating heatmap '''   
    createHeatMaps(dataFrames)