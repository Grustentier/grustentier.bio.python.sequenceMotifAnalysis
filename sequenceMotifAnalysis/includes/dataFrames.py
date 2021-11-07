'''
Created on 07.11.2021

@author: grustentier
'''
import math
import pandas
from scipy import stats


def getDataFrames4Heatmap(possibleMotifs, AMINO_ACID_LETTERS, AMINO_ACID_OCCURRENCES):
    dataFrames = {}    
   
    for topology in possibleMotifs.keys():
        dataFrames[topology] = {"yLabels":[], "dataFrame":[]}  
        for regEx in possibleMotifs[topology]: 
            registrations = []
            
            for _ in range(0, int(regEx[2:]) - 1):
                registrations.append({})
                for letter in AMINO_ACID_LETTERS:                    
                    registrations[-1][letter] = 0                    
                # y_axis_labels.append(regEx+"-"+str(varPosIndex+1)+"-"+topology)
                rowLabel = regEx + "-" + topology
                if rowLabel not in dataFrames[topology]["yLabels"]:
                    dataFrames[topology]["yLabels"].append(rowLabel) 
                else:dataFrames[topology]["yLabels"].append("")
             
            for motifSequence in possibleMotifs[topology][regEx]:
                varPosOfMotif = motifSequence[1:-1] 
                for pos in range(0, len(varPosOfMotif)):
                    if varPosOfMotif[pos] in AMINO_ACID_LETTERS:
                        registrations[pos][varPosOfMotif[pos]] += 1
                    
            if len(possibleMotifs[topology][regEx]) > 0:
                for pos in range(0, len(registrations)): 
                    dataFrameContent = []
                    for aminoAcid in AMINO_ACID_LETTERS: 
                        # quotient = registrations[pos][aminoAcid] / len(possibleMotifs[topology][regEx])
                        quotient = (registrations[pos][aminoAcid] / len(possibleMotifs[topology][regEx])) / AMINO_ACID_OCCURRENCES[AMINO_ACID_LETTERS.index(aminoAcid)]
                        if quotient != 0.0:quotient = math.log(quotient)
                        dataFrameContent.append(quotient)
                    dataFrames[topology]["dataFrame"].append(dataFrameContent)
                    
    return dataFrames


def getDataFrameFromOccurrences(positionSpecificStatistics, AMINO_ACID_LETTERS):
    data = {AMINO_ACID_LETTERS[aaIndex]:[positionSpecificStatistics[i]["occurrence-ratios"][aaIndex] for i in range(0, len(positionSpecificStatistics))] for aaIndex in range(0, len(AMINO_ACID_LETTERS))}
    data["topology"] = [positionSpecificStatistics[i]["topology"] for i in range(0, len(positionSpecificStatistics))]
    return pandas.DataFrame(data)


def getDataFrameFromSpearman(positionSpecificStatistics):
    data = {}

    for i in range(0, len(positionSpecificStatistics)):
        comatrix = [] 
        for j in range(0, len(positionSpecificStatistics)):
            correlation, _ = stats.spearmanr(positionSpecificStatistics[i]["occurrence-ratios"], positionSpecificStatistics[j]["occurrence-ratios"])
            comatrix.append(correlation)   
        dataKey = positionSpecificStatistics[i]["regEx"] + "-" + str(positionSpecificStatistics[i]["position"]) + "-" + positionSpecificStatistics[i]["topology"]
        data[dataKey] = comatrix  
        
    data["topology"] = [positionSpecificStatistics[i]["topology"] for i in range(0, len(positionSpecificStatistics))]
    return pandas.DataFrame(data) 


def getDataFrameFromPearson(positionSpecificStatistics):
    data = {}

    for i in range(0, len(positionSpecificStatistics)):
        comatrix = [] 
        for j in range(0, len(positionSpecificStatistics)):
            correlation, _ = stats.pearsonr(positionSpecificStatistics[i]["occurrence-ratios"], positionSpecificStatistics[j]["occurrence-ratios"])
            comatrix.append(correlation)   
        dataKey = positionSpecificStatistics[i]["regEx"] + "-" + str(positionSpecificStatistics[i]["position"]) + "-" + positionSpecificStatistics[i]["topology"]
        data[dataKey] = comatrix 
        
    data["topology"] = [positionSpecificStatistics[i]["topology"] for i in range(0, len(positionSpecificStatistics))]
    return pandas.DataFrame(data)
