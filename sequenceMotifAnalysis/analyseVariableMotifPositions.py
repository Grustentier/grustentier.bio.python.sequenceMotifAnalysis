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
Based on position specific amino acid occurrences, the topology state (tm, ntm or trans) will be predicted using a simple winner takes it all formular.
These topology information (export to ./inputdata/motif-topologies.xml) can be used to label high occurrent motifs within the resulting graph generated 
by script analyseConsecutiveMotifs.py for highlighting topology unspecific motifs.

Steffen Grunert, Florian Heinke, Dirk Labudde, 
"Structure Topology Prediction of Discriminative Sequence Motifs in Membrane Proteins with Domains of Unknown Functions", 
Structural Biology, vol. 2013, Article ID 249234, 10 pages, 2013. https://doi.org/10.1155/2013/249234
'''

import os
import sys
import math
import argparse
import matplotlib.pyplot as plt
import xml.etree.ElementTree as ET
sys.path.append(os.path.join(os.path.dirname(__file__), "..")) 

from sequenceMotifAnalysis.includes import letters
from sequenceMotifAnalysis.includes import fileUtils as fU
from sequenceMotifAnalysis.includes import dataFrames as dF
from sequenceMotifAnalysis.includes import stringUtils as sU
from sequenceMotifAnalysis.includes import dataCollecting as dC
from sequenceMotifAnalysis.includes import dataVisualization as dV

parser = argparse.ArgumentParser(description='Code for analysis of variable sequence motif positions  for different topologies.')
parser.add_argument('--fasta_input_dir', default=os.path.dirname(__file__) + os.sep + 'inputdata' + os.sep + 'fasta' + os.sep + 'DUF', help='Path to the input dir including fasta files.')
parser.add_argument('--tmhmm_input_dir', default=os.path.dirname(__file__) + os.sep + 'inputdata' + os.sep + 'tmhmm' + os.sep + 'DUF', help='Path to the input dir including tmhmm files, generated extensive and  with no graphics.')
'''Gerstein pattern (regExes) '''
parser.add_argument('--regexes', default='PG10,LF10,PG9,LF9,VF8,LF8,GY8,GA7,AG7,AA7,GG7,LY6,VG6,SA6,PG6,AL6,PG5,GS5,LG5,AG5,GN4,IV4,IL4,GS4,GG4,SG4,VL4,AS4,GA4,AG4,SA3,AA3,GL3', type=str, help='Comma separted REGEXES like XXn representing a starting and a ending aminoacid by n - 1 variable position between both X.')
''' Other regExes maybe coming from console after running analyseConsecutiveMotifs.py to specify/predict the structure topology of given motifs ''' 
# parser.add_argument('--regexes', default='WP5, LL7, YL3, DF7, LL5, YP8, DA3, KF3, TL3, LL4, LL6, WL7, LL3, DK4, DP5, WY3, PW3, TL4, YA7, RP9, LA4, WL6, KG4, KL7, WT4, PL6, PL3, DV6, RD3, WT3, IL5, PL8, YT7, YL9, GL3, DG8, DP6, AF4, FL4, GF4, GA3, YL6, LF8, DM3, IG3, WF5, MG4, LA3, TL7, LY4, RT8, KL8, YT6, TL9, MT3, GA4, LK5, LP3, DT4, DG7, DL8, AL8, DT5, YK8, GL4, LE4, LV7, WP4, WL8, DL7, RT7', type=str, help='Comma separted REGEXES like XXn representing a starting and a ending aminoacid by n - 1 variable position between both X.')
parser.add_argument('--display_heatmap', default=True, type=str, help='True or 1 for displaying heatmap representing amino acid occurrences at variable position of motifs defined at --regexes.')
parser.add_argument('--display_clustermap', default=True, type=str, help='True or 1 for displaying clustermaps representing clustered variable positions of motifs defined at --regexes.')
parser.add_argument('--plotting_3d', default=False, type=str, help='True or 1 for displaying cluster results within 3d plots')
parser.add_argument('--sort_by', default=None, type=str, help='Selection: alphabetical, hydrophob or None; Sorting amino acid presentation by properties.')
arguments = parser.parse_args() 

    
def exportWinners(motifWinners, filePath=os.path.dirname(__file__) + os.sep + 'inputdata' + os.sep + 'motif-topologies.xml'):    
    if os.path.exists(filePath) is False:
        root = ET.Element('motifs')
        with open(filePath, "wb") as f:
            f.write(ET.tostring(root))
    
    tree = ET.parse(filePath)
    root = tree.getroot()
     
    for regEx in motifWinners.keys():
        winnerTopology = max(motifWinners[regEx].items(), key=lambda k : k[1])[0] 
         
        motifElement = root.find(".//*[@regEx='" + regEx + "']")
        if motifElement is None:
            motifElement = ET.Element('motif')
                 
        motifElement.set("regEx", regEx)
        motifElement.set("topology", winnerTopology)   
        root.append(motifElement)   
        
    with open(filePath, "wb") as f:
        f.write(ET.tostring(root))


def getTopology(start, end, tmhmmData):  
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

    
def getPossibleMotifs(fastaData, tmhmmData):
    possibleMotifs = {"tm":{}, "ntm":{}, "trans":{}}
    minVarPos, maxVarPos = getMinMaxVarPos()

    for data in fastaData:
        sequence = data["sequence"] 
        fasta_id = data["id"]
        current_tmhmmData = dC.findCorrespondingTmhmmData(fasta_id, tmhmmData)
        if current_tmhmmData is None:print("No tmhmm data for to fasta id:", fasta_id)  
        
        for i in range(0, len(sequence)):
            for j in range(minVarPos + 1, maxVarPos + 2):
                if (i + j) >= len(sequence):continue
                regEx = sequence[i] + sequence[i + j] + str(j)  
                if regEx not in REGEXES: continue    
                motifSeq = sequence[i:i + j + 1]       
                  
                if len(motifSeq) > 3:                         
                    topology2Assign = getTopology(i, i + j, current_tmhmmData)
                    if topology2Assign not in possibleMotifs.keys():possibleMotifs[topology2Assign] = {}
                    if regEx not in possibleMotifs[topology2Assign].keys():possibleMotifs[topology2Assign][regEx] = [] 
                    possibleMotifs[topology2Assign][regEx].append(motifSeq) 
         
    return possibleMotifs 


def getPositionSpecificStatistics(possibleMotifs):
    positionSpecificStatistics = [] 
    
    for topology in possibleMotifs.keys():         
        for regEx in possibleMotifs[topology]:             
            regExsStatistic = []             
            
            for pos in range(0, int(regEx[2:]) - 1): 
                regExsStatistic.append({"regEx":regEx, "position":pos, "topology":topology, "occurrences":[0 for _ in AMINO_ACID_LETTERS], "occurrence-ratios":[0.0 for _ in AMINO_ACID_LETTERS]}) 
             
            for motifSequence in possibleMotifs[topology][regEx]:
                varPosOfMotif = motifSequence[1:-1] 
                for pos in range(0, len(varPosOfMotif)):
                    if varPosOfMotif[pos] in AMINO_ACID_LETTERS:
                        indexInArray = AMINO_ACID_LETTERS.index(varPosOfMotif[pos])
                        regExsStatistic[pos]["occurrences"][indexInArray] += 1
            positionSpecificStatistics.extend(regExsStatistic)
        
    for i in range(0, len(positionSpecificStatistics)):
        for ocPos in range(0, len(positionSpecificStatistics[i]["occurrences"])):
            quotient = 0.0                      
            if AMINO_ACID_OCCURRENCES[ocPos] != 0.0:
                quotient = (positionSpecificStatistics[i]["occurrences"][ocPos] / len(possibleMotifs[positionSpecificStatistics[i]["topology"]][positionSpecificStatistics[i]["regEx"]])) / AMINO_ACID_OCCURRENCES[ocPos]                    
            if quotient != 0.0:quotient = math.log2(quotient)
            positionSpecificStatistics[i]["occurrence-ratios"][ocPos] = quotient
            
    return positionSpecificStatistics


def getAminoAcidOccurrencesinNature(): 
    print('\n\nDetermining amino acid occurrences in nature based on on https://en.wikipedia.org/wiki/Amino_acid\n')
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
    
    return [wiki[aa] for aa in AMINO_ACID_LETTERS]

    
def getMinMaxVarPos():
    minVarPos = 10000000000
    maxVarPos = 0
    
    for regEx in REGEXES:
        varPos = int(regEx[2:]) - 1
        if varPos < minVarPos: minVarPos = varPos  
            
        if varPos > maxVarPos: maxVarPos = varPos 
    return minVarPos, maxVarPos    


def getOccurrencesRatios(regEx, position, topology, positionSpecificStatistics):
    occurrences = []
    for stat in positionSpecificStatistics:
        if stat["regEx"] == regEx and stat["position"] == position and stat["topology"] == topology: 
            return stat["occurrence-ratios"]
    return occurrences


def determineWinners(possibleMotifs, positionSpecificStatistics):
    motifWinners = {}
    
    for topology in possibleMotifs.keys():
        for regEx in possibleMotifs[topology].keys():
            if regEx not in motifWinners.keys():motifWinners[regEx] = {}
            if topology not in motifWinners[regEx].keys():motifWinners[regEx][topology] = 0
            
            for motif in possibleMotifs[topology][regEx]:
                varPosSeq = motif[1:-1]
                for pos in range(0, len(varPosSeq)):
                    letterAtPos = varPosSeq[pos] 
                    if letterAtPos not in AMINO_ACID_LETTERS:continue
                    oc = getOccurrencesRatios(regEx, pos, topology, positionSpecificStatistics)
                    assert len(oc) > 0, "wrong assingment"
                    motifWinners[regEx][topology] += oc[AMINO_ACID_LETTERS.index(letterAtPos)]   
                
    return motifWinners

    
def printProgress(steps, maximum, name="progress", bar_length=20, width=20):  
    if maximum == 0:
        percent = 1.0
    else:
        percent = float(steps) / maximum
    arrow = '-' * int(round(percent * bar_length) - 1) + '>'
    spaces = ' ' * (bar_length - len(arrow))
    # sys.stdout.write("\r{0: <{1}} : [{2}]{3}%".format(name, width, arrow + spaces, int(round(percent*100))))
    sys.stdout.write("\r{0: <{1}}[{2}]{3}%".format("", 0, arrow + spaces, int(round(percent * 100))))
    sys.stdout.flush()    
    
    if steps >= maximum:     
        sys.stdout.write('\n\n')  

    
if __name__ == "__main__":
    fastaFilePaths = fU.collectFilePaths(arguments.fasta_input_dir)
    assert len(fastaFilePaths) > 0, "No fasta files have been found!!!"
    tmhmmFilePaths = fU.collectFilePaths(arguments.tmhmm_input_dir)
    assert len(tmhmmFilePaths) > 0, "No tmhmm files have been found!!!"
    assert len(arguments.regexes) > 0, "No regex found have been found!!!"
    
    REGEXES = []
    if str(arguments.regexes).find(",") != -1:
        REGEXES = str(arguments.regexes).upper().replace(" ", "").split(",")
    else:
        REGEXES.append(str(arguments.regexes).upper().replace(" ", ""))        
        
    AMINO_ACID_LETTERS = letters.getAminoAcidLetters(arguments.sort_by) 
    
    steps = 7
    step = 0
    if sU.asBoolean(arguments.display_heatmap) is True:steps += 2
    if sU.asBoolean(arguments.display_clustermap) is True:steps += 3
    
    printProgress(step, steps)
    step = step + 1
    
    ''' Collecting fasta data '''
    fastaData = dC.collectFastaData(fastaFilePaths)
    printProgress(step, steps)
    step = step + 1
    
    ''' Collecting tmhmm data  '''    
    tmhmmData = dC.collectTmhmmData(tmhmmFilePaths)
    printProgress(step, steps)
    step = step + 1
    
    ''' Gathering motifs from sequences ''' 
    possibleMotifs = getPossibleMotifs(fastaData, tmhmmData)
    printProgress(step, steps)
    step = step + 1
    
    ''' Determining amino acid occurrences in nature from different sources'''  
    AMINO_ACID_OCCURRENCES = getAminoAcidOccurrencesinNature()
    printProgress(step, steps)
    step = step + 1
    
    ''' Determining position specific statistics ''' 
    positionSpecificStatistics = getPositionSpecificStatistics(possibleMotifs) 
    printProgress(step, steps)
    step = step + 1
    
    ''' Determining topology winners ''' 
    motifWinners = determineWinners(possibleMotifs, positionSpecificStatistics)
    printProgress(step, steps)
    step = step + 1
    
    ''' Exporting topology winners '''
    exportWinners(motifWinners)
    printProgress(step, steps)
    step = step + 1
    
    if sU.asBoolean(arguments.display_heatmap) is True:
        ''' Generating dataFrame for heatmap''' 
        dataFrames = dF.getDataFrames4Heatmap(possibleMotifs, AMINO_ACID_LETTERS, AMINO_ACID_OCCURRENCES)
        printProgress(step, steps)
        step = step + 1
        
        ''' Creating heatmap '''   
        dV.createHeatMaps(dataFrames, AMINO_ACID_LETTERS)
        printProgress(step, steps)
        step = step + 1
    
    if sU.asBoolean(arguments.display_clustermap) is True:
        ''' Generating panda.DataFrames from position specific statistics '''      
        ''' Clustering data and presenting within scatter plots '''
        df1 = [{"title":"Clustering based on raw data (position specific amino acid occurrences)", "dataFrame":dF.getDataFrameFromOccurrences(positionSpecificStatistics, AMINO_ACID_LETTERS)}]
        
        if sU.asBoolean(arguments.plotting_3d) is True:
            dV.cluster3d(df1)
        else:
            dV.cluster(df1)
            
        printProgress(step, steps)
        step = step + 1
        
        df2 = [{"title":"Clustering based on spearman rank correlation", "dataFrame":dF.getDataFrameFromSpearman(positionSpecificStatistics)}]
       
        if sU.asBoolean(arguments.plotting_3d) is True:
            dV.cluster3d(df2)
        else:
            dV.cluster(df2)
            
        printProgress(step, steps)
        step = step + 1
       
        df3 = [{"title":"Clustering based on pearson rank correlation", "dataFrame":dF.getDataFrameFromPearson(positionSpecificStatistics)}]
        
        if sU.asBoolean(arguments.plotting_3d) is True:
            dV.cluster3d(df3)
        else:
            dV.cluster(df3)
            
        printProgress(step, steps)
        step = step + 1
    
    if sU.asBoolean(arguments.display_heatmap) is True or sU.asBoolean(arguments.display_clustermap) is True:
        plt.show()
        
