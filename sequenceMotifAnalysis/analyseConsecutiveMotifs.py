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
  _  _          _            __  __         _     _    __     ___          _                _     _                                          
 | || |  _  _  | |__   ___  |  \/  |  ___  | |_  (_)  / _|   |   \   ___  | |_   ___   __  | |_  (_)  ___   _ _                              
 | __ | | || | | '_ \ |___| | |\/| | / _ \ |  _| | | |  _|   | |) | / -_) |  _| / -_) / _| |  _| | | / _ \ | ' \                             
 |_||_|  \_,_| |_.__/       |_|  |_| \___/  \__| |_| |_|     |___/  \___|  \__| \___| \__|  \__| |_| \___/ |_||_|                            
                                                                                                                                             
  _                ___                     _                  _     _                                                                        
 | |__   _  _     / __|  _ _   _  _   ___ | |_   ___   _ _   | |_  (_)  ___   _ _                                                            
 | '_ \ | || |   | (_ | | '_| | || | (_-< |  _| / -_) | ' \  |  _| | | / -_) | '_|                                                           
 |_.__/  \_, |    \___| |_|    \_,_| /__/  \__| \___| |_||_|  \__| |_| \___| |_|                                                             
         |__/                                                                                                                                


""")

'''
Analysis of the statistical occurrence of consecutive sequence motifs in membrane proteins.
It is a simple alternative to the HMM logo for visualizing frequently occurring motifs.
This script is not yet finished. In further steps, I will use topology-specific motif information from the analyseVariableMotifPositions.py script to label topology-unspecific motifs.
For example, a common motif in transmembrane has non-transmembrane amino acid properties like shown and described in paper:

Grunert, S., Labudde, D. 
Graph representation of high-dimensional alpha-helical membrane protein data. 
BioData Mining 6, 21 (2013). https://doi.org/10.1186/1756-0381-6-21
'''

import os
import re
import argparse
from cv2 import data
import networkx as nx
import matplotlib.pyplot as plt 

parser = argparse.ArgumentParser(description='Code for analysis of variable sequence motif positions  for different topologies.')
parser.add_argument('--fasta_input_dir', default='.'+os.sep+'testdata'+os.sep+'fasta'+os.sep+'rhodopsins', help='Path to the input dir including fasta files.')
parser.add_argument('--tmhmm_input_dir', default='.'+os.sep+'testdata'+os.sep+'tmhmm'+os.sep+'rhodopsins', help='Path to the input dir including tmhmm files, generated extensive and  with no graphics.') 
parser.add_argument('--min_variable_positions', default=3, type=int, help='The min number of min variable x- position of a given sequence motif XYn by n-1 variable postions like GG4 = GxxxG')
parser.add_argument('--max_variable_positions', default=9, type=int, help='The max number of min variable x- position of a given sequence motif XYn by n-1 variable postions like GG4 = GxxxG')
parser.add_argument('--max_nodes', default=30, type=int, help='The max number of graph nodes to visualize') 
parser.add_argument('--topology', default="tm", type=str, help='The topology to analyse. (tm for transmembrane or ntm for none-transmembrane')
parser.add_argument('--as_tree', default=False, type=str, help='True or 1 for displaying graph as tree.')
 
arguments = parser.parse_args() 

AMINO_ACIDS_ONE_LETTER_CODE_ALPHABETICAL = ["A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y"]

def boolean_string(s):
    if str(s).lower() not in ['false', 'true', '1', '0']:
        raise ValueError('Not a valid boolean string')
    return str(s).lower() == 'true' or str(s).lower() == '1'
        
def cleanString(string,replacement='-'):
    a =  re.sub('[^a-zA-Z0-9.?]',replacement,string) 
    return re.sub(replacement+'+', replacement, a)

def collectFilePaths(directory):
    filePaths = []
    for root, _, files in os.walk(directory):  
        if root.endswith(os.sep) is False:
            root+=os.sep
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
    
    printProgress(0,len(fastaFilePaths))
    for index in range(0,len(fastaFilePaths)): 
        fastaFilePath = fastaFilePaths[index]
        data = getFastaData(fastaFilePath)
        if data is not None:
            fastaData.extend(data)        
        printProgress(index,len(fastaFilePaths))
            
    return fastaData

def collectTmhmmData(tmhmmFilePaths):
    tmhmmData = []
    
    printProgress(0,len(tmhmmFilePaths))
    for index in range(0,len(tmhmmFilePaths)): 
        tmhmmFilePath = tmhmmFilePaths[index]
        data = getTmhmmData(tmhmmFilePath)
        if data is not None:
            tmhmmData.extend(data)
        printProgress(index,len(tmhmmFilePaths))
            
    return tmhmmData 

def findCorrespondingFastaData(tmhmm_id,fastaData):    
    for data in fastaData: 
        fasta_id = data["id"]  
        if str(fasta_id).lower() == str(tmhmm_id).lower():
            return data
        
    return None
    
def collectTopologySpecficSubSequences(fastaData,tmhmmData): 
    subsequences = []
    
    printProgress(0,len(tmhmmData))
    for index in range(0,len(tmhmmData)):
        data_tmhmm = tmhmmData[index]
        tmhmm_id = data_tmhmm["id"]
        current_fastaData = findCorrespondingFastaData(tmhmm_id,fastaData)
        if current_fastaData is None:
            continue

        for area in data_tmhmm["areas"]:
            append = False
            if str(arguments.topology).lower() == "tm" and area["topology"] == "tmhelix":
                append = True                
            elif str(arguments.topology).lower() == "ntm" and (area["topology"] == "inside" or area["topology"] == "outside"):
                append = True
            else:append = False
            
            if append is True:       
                start_tmhmm = int(area["from"]) - 1
                end_tmhmm = int(area["to"]) - 1 
                subsequences.append(current_fastaData["sequence"][start_tmhmm:end_tmhmm])          
        printProgress(index,len(tmhmmData))          
                            
    return subsequences     

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

def createConsecutiveMotifsStatistics(topologySpecificSubSequences): 
    statistics = {}
    
    printProgress(0,len(topologySpecificSubSequences))
    for sequenceIndex in range(0,len(topologySpecificSubSequences)):
        sequence = topologySpecificSubSequences[sequenceIndex]        
        for i in range(0,len(sequence)):            
            for j in range(arguments.min_variable_positions,arguments.max_variable_positions+1): 
                if (i + j) >= len(sequence):
                    continue    
                regEx1 = sequence[i]+sequence[i+j]+str(j)  
                motifSeq1 = sequence[i:i+j+1]                 
                for l in range(arguments.min_variable_positions,arguments.max_variable_positions+1):  
                    if (i+j+1 + l) >= len(sequence):
                        continue
                    regEx2 = sequence[i+j+1]+sequence[i+j+1+l]+str(l) 
                    motifSeq2 = sequence[i+j+1:i+j+1+l+1]                        
                    cmKey = regEx1 + "-" + regEx2                     
                    if cmKey not in statistics:
                        statistics[cmKey] = {"sequence":sequence,"regEx1":regEx1,"regEx2":regEx2,"motifSequence1":motifSeq1,"motifSequence2":motifSeq2,"occurrence":1}
                    else:
                        statistics[cmKey]["occurrence"] += 1                             
        printProgress(sequenceIndex,len(topologySpecificSubSequences))
        
    return statistics
        
def getTop(g,node):
    nodes = g.predecessors(node) 
    if len(list(nodes)) == 0: 
        return node
    for n in nodes:
        return getTop(g,n)        
    return None

def getDiGraph(statistics):
    graph =  nx.DiGraph()
         
    if arguments.max_nodes is not None:
        statistics = statistics[0:arguments.max_nodes]
        
    for listEntry in statistics: 
        try:
            regEx1 = listEntry[1]["regEx1"]
            regEx2 = listEntry[1]["regEx2"] 
            graph.add_node(regEx1)
            graph.add_node(regEx2)
            #weight = listEntry[1]["occurrence"] * int(listEntry[1]["regEx1"][2:]) * int(listEntry[1]["regEx2"][2:])
            graph.add_edge(regEx1,regEx2,weight=listEntry[1]["occurrence"]) 
        except: 
            pass 
        
    return graph 

def createGraph(statistics):
    statistics = sorted(statistics.items(), key=lambda x: x[1]["occurrence"], reverse=True)
    graph = getDiGraph(statistics)
    
    if arguments.as_tree is True:
        layout = nx.nx_pydot.graphviz_layout(graph,prog='dot')
    else:  
        layout = nx.nx_pydot.graphviz_layout(graph)
        
    labels = nx.get_edge_attributes(graph,'weight')
    nx.draw_networkx_edge_labels(graph,layout,edge_labels=labels)
    
    nodeColor = "black"
    if str(arguments.topology).lower() == "tm":
        nodeColor = "#cd4d48"
    if str(arguments.topology).lower() == "ntm":
        nodeColor = "#8cb555"
    
    cent = nx.centrality.betweenness_centrality(graph,weight=None,normalized=False,endpoints=True) 
    nx.draw(graph,layout,width=1,linewidths=1,node_size=[v*500 for v in cent.values()],node_color=nodeColor, edge_color='silver',alpha=0.9,labels={node:node for node in graph.nodes()})
    plt.axis('off')    

if __name__ == "__main__":  
    fastaFilePaths = collectFilePaths(arguments.fasta_input_dir)
    assert len(fastaFilePaths)>0,"No fasta files have been found!!!"
    tmhmmFilePaths = collectFilePaths(arguments.tmhmm_input_dir)
    assert len(tmhmmFilePaths)>0,"No tmhmm files have been found!!!"       
        
    AMINO_ACIDS_ONE_LETTER_CODE = AMINO_ACIDS_ONE_LETTER_CODE_ALPHABETICAL
    
    print("Parsing and collecting fasta files...")
    fastaData = collectFastaData(fastaFilePaths)     
    
    print("Parsing and collecting tmhmm files...")
    tmhmmData = collectTmhmmData(tmhmmFilePaths)    
    
    print("Gathering motifs from sequences...")      
    topologySpecificSubSequences = collectTopologySpecficSubSequences(fastaData,tmhmmData)
    
    print("Collecting consecutive Motifs...")
    statistics = createConsecutiveMotifsStatistics(topologySpecificSubSequences) 
    
    print("Creating graph...")    
    createGraph(statistics)  
    plt.show()    
    