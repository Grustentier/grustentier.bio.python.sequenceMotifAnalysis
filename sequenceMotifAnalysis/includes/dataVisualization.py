'''
Created on 07.11.2021

@author: grustentier
'''
import pandas
import seaborn
import matplotlib.pyplot as plt
from sklearn.manifold import MDS
# from sklearn.cluster import KMeans
from sklearn.decomposition import PCA 
from matplotlib.colors import LinearSegmentedColormap

TOPOLOGIES = ['tm', 'ntm', 'trans']
COLORS = ['#cc4c48', '#4ab3c9', '#8cb555', 'c', 'm']


# Visualization of variable x-positions of selected sequnece motifs based on topology specific amino acid distribution using PCA and MDS within 2D plots
# Each dot within a plot represent a x-position of one of the selected sequence motifs
def cluster(dataFrames, cluster=3):    
        
    for index in range(0, len(dataFrames)):   
        row = 0  
        fig, a = plt.subplots(nrows=2, ncols=1)
        fig.subplots_adjust(wspace=0.2, hspace=0.4)
        seaborn.set(font_scale=0.5)
        fig.canvas.set_window_title(dataFrames[index]["title"])   
       
        dataFrame = dataFrames[index]["dataFrame"]
        if dataFrame.empty is True:continue
        
        features = dataFrame.columns[0:-1]        
        x = dataFrame.loc[:, features].values
        # y = dataFrame.loc[:,['topology']].values
        ''' PCA scatter plot '''
        pca = PCA(2)         
        pcaComponents = pca.fit_transform(x)
        pcaComponentsDataFrame = pandas.DataFrame(data=pcaComponents, columns=['pc1', 'pc2'])
        finalPCADataFrame = pandas.concat([pcaComponentsDataFrame, dataFrame[['topology']]], axis=1)
        for target, color in zip(TOPOLOGIES, COLORS):
            indicesToKeep = finalPCADataFrame['topology'] == target
            a[row].scatter(finalPCADataFrame.loc[indicesToKeep, 'pc1'], finalPCADataFrame.loc[indicesToKeep, 'pc2'], c=color, s=50)
        a[row].set_xlabel('Principal Component 1', fontsize=6)
        a[row].set_ylabel('Principal Component 2', fontsize=6)
        a[row].set_title('PCA')
        a[row].legend(TOPOLOGIES)
        row = row + 1
        
        ''' KMeans scatter plot based on PCA components 
        kmeans = KMeans(n_clusters= cluster)        
        labels = kmeans.fit_predict(pcaComponents)
        centroids = kmeans.cluster_centers_
        classes = numpy.unique(labels) 
                
        for clazz in labels:         
            filtered_label = pcaComponents[labels == clazz]
            a[row].scatter(filtered_label[:,0], filtered_label[:,1],color=colors[list(classes).index(clazz)])
        
        a[row].scatter(centroids[:,0] , centroids[:,1] , s = 80, color = 'k')
        a[row].set_title('kmeans based on PCA components')
        a[row].legend(targets)
        row = row + 1
        '''
         
        ''' MDS scatter plot '''
        mds = MDS(n_components=2, metric=True, n_init=4, max_iter=300, verbose=0, eps=0.001, n_jobs=None, random_state=42, dissimilarity='euclidean')
        mdsComponents = mds.fit_transform(x)   
        mdsComponentsDataFrame = pandas.DataFrame(data=mdsComponents, columns=['pc1', 'pc2'])
        finalMDSDataFrame = pandas.concat([mdsComponentsDataFrame, dataFrame[['topology']]], axis=1)
        for target, color in zip(TOPOLOGIES, COLORS):
            indicesToKeep = finalMDSDataFrame['topology'] == target
            a[row].scatter(finalMDSDataFrame.loc[indicesToKeep, 'pc1'], finalMDSDataFrame.loc[indicesToKeep, 'pc2'], c=color, s=50)        
        a[row].set_xlabel('Principal Component 1', fontsize=6)
        a[row].set_ylabel('Principal Component 2', fontsize=6)
        a[row].set_title('MDS')
        a[row].legend(TOPOLOGIES)
        row = row + 1
        
        ''' KMeans scatter plot based on PCA components 
        kmeans = KMeans(n_clusters= cluster)        
        labels = kmeans.fit_predict(mdsComponents)
        centroids = kmeans.cluster_centers_
        classes = numpy.unique(labels) 
                
        for clazz in labels:         
            filtered_label = mdsComponents[labels == clazz]
            a[row].scatter(filtered_label[:,0], filtered_label[:,1],color=colors[list(classes).index(clazz)])
        
        a[row].scatter(centroids[:,0] , centroids[:,1] , s = 80, color = 'k')
        a[row].set_title('kmeans based on MDS components')
        a[row].legend(targets) 
        '''
        
    # plt.show()


# Visualization of variable x-positions of selected sequnece motifs based on topology specific amino acid distribution using PCA and MDS within 3D plots
# Each dot within a plot represent a x-position of one of the selected sequence motifs   
def cluster3d(dataFrames):
    
    fig = plt.figure()
    for index in range(0, len(dataFrames)): 
        a = fig.add_subplot(211, projection='3d')
        fig.subplots_adjust(wspace=0.2, hspace=0.4)
        seaborn.set(font_scale=0.5)
        fig.canvas.set_window_title(dataFrames[index]["title"])   
       
        dataFrame = dataFrames[index]["dataFrame"]
        if dataFrame.empty is True:continue
        
        features = dataFrame.columns[0:-1]        
        x = dataFrame.loc[:, features].values
        # y = dataFrame.loc[:,['topology']].values
        ''' PCA scatter plot '''
        pca = PCA(3)         
        pcaComponents = pca.fit_transform(x)
        pcaComponentsDataFrame = pandas.DataFrame(data=pcaComponents, columns=['pc1', 'pc2', 'pc3'])
        finalPCADataFrame = pandas.concat([pcaComponentsDataFrame, dataFrame[['topology']]], axis=1)
        for target, color in zip(TOPOLOGIES, COLORS):
            indicesToKeep = finalPCADataFrame['topology'] == target
            a.scatter(finalPCADataFrame.loc[indicesToKeep, 'pc1'], finalPCADataFrame.loc[indicesToKeep, 'pc2'], finalPCADataFrame.loc[indicesToKeep, 'pc3'], c=color, s=50)
        a.set_xlabel('Principal Component 1', fontsize=6)
        a.set_ylabel('Principal Component 2', fontsize=6)
        a.set_zlabel('Principal Component 3', fontsize=6)
        a.set_title('PCA')
        a.legend(TOPOLOGIES)
    
        a = fig.add_subplot(212, projection='3d')
        ''' MDS scatter plot '''
        mds = MDS(n_components=3, metric=True, n_init=4, max_iter=300, verbose=0, eps=0.001, n_jobs=None, random_state=42, dissimilarity='euclidean')
        mdsComponents = mds.fit_transform(x)   
        mdsComponentsDataFrame = pandas.DataFrame(data=mdsComponents, columns=['pc1', 'pc2', 'pc3'])
        finalMDSDataFrame = pandas.concat([mdsComponentsDataFrame, dataFrame[['topology']]], axis=1)
        for target, color in zip(TOPOLOGIES, COLORS):
            indicesToKeep = finalMDSDataFrame['topology'] == target
            a.scatter(finalMDSDataFrame.loc[indicesToKeep, 'pc1'], finalMDSDataFrame.loc[indicesToKeep, 'pc2'], finalMDSDataFrame.loc[indicesToKeep, 'pc3'], c=color, s=50)        
        a.set_xlabel('Principal Component 1', fontsize=6)
        a.set_ylabel('Principal Component 2', fontsize=6)
        a.set_zlabel('Principal Component 3', fontsize=6)
        a.set_title('MDS')
        a.legend(TOPOLOGIES)
        
    # plt.show()


# Visualization of position specific amino acid distribution data of selected sequence motif within customized heatmap    
def createHeatMaps(dataFrames, AMINO_ACID_LETTERS): 
    fig, a = plt.subplots(nrows=len(dataFrames.keys()))
    fig.subplots_adjust(wspace=0.01)
    seaborn.set(font_scale=0.5) 
    
    index = 0    
    for topology in dataFrames.keys(): 
        if len(dataFrames[topology]["dataFrame"]) > 0:
            custom_color_map = LinearSegmentedColormap.from_list(name='custom_navy', colors=[(0 / 255, 0 / 255, 255 / 255), (255 / 255, 0 / 255, 0 / 255)])
            heatmap = seaborn.heatmap(dataFrames[topology]["dataFrame"], ax=a[index], cmap=custom_color_map, cbar=True, annot=False, annot_kws={"size": 1}, xticklabels=AMINO_ACID_LETTERS, yticklabels=dataFrames[topology]["yLabels"])
            heatmap.set_yticklabels(heatmap.get_yticklabels(), rotation=0, fontsize=6)            
        
        index = index + 1     
    
    # plt.show()
