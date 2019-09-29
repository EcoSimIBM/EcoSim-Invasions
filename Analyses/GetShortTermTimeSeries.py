from __future__ import division
import sys
import numpy as np
import pandas as pd
import sklearn
from sklearn import metrics
#from sklearn.tree import DecisionTreeClassifier
#from sklearn.linear_model import LogisticRegressionCV
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import MinMaxScaler
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from sklearn.utils import resample
from rfpimp import *
from sklearn.ensemble.forest import _generate_unsampled_indices
import eli5
from eli5.sklearn import PermutationImportance
#anova
#import statsmodels.api as sm
#from statsmodels.formula.api import ols
import scipy.stats as stats
import scikit_posthocs as sp
from scipy.stats import rankdata

np.set_printoptions(threshold=11000)
np.set_printoptions(suppress=True)

#skip the header line
numRuns = 10
numLevels = 5
treatments = ["InvasionsDefaultAcceptors", "InvasionsDefaultAcceptorsNiches", "InvasionsNichesAcceptors", "InvasionsNichesAcceptorsNiches"]
typs = ["NGDT"]
fileN = "Results_Inv_Prey_Female"
resultsColumns = ["total population", "nbPerCell"]
startRow = 15100
maxRow = 20100

future = 60
window = 5

numberOfComparisonsp1=0
numberOfComparisonsp05=0
numberOfComparisonsp01=0

lowbounds = {}
upperbounds = {}

tempColumns = resultsColumns[:]  
tempColumns.append("class")

for columnName in tempColumns:
    lowbounds[columnName] = 1000000
    upperbounds[columnName] = -1000000
 
for typ in typs:
    framesListAllMeans = []
    framesListAllSDevs = []
    framesListAllMedians = []
    framesListAll2p5 = []
    framesListAll97p5 = []
    for t in treatments:
        print("Processing treatment " + t)
        importances=[]
        sys.stdout.flush()
        for n in range(1, 1 + numLevels):
            sys.stdout.flush()
            print("\tProcessing level " + str(n))
            arr = []
            for r in range(1, 1 + numRuns):
                print("\t\tProcessing run number " + str(r))
                df = pd.read_csv(t + str(typ) + str(n) + "/run" + str(r) + "/" + fileN + ".csv", sep=', ')
                col = df[resultsColumns]
                actualMax = min(df.shape[0] - (future + 1), maxRow-1)
                rowsProcessed = 0
                
                for timestep in range(startRow-1, actualMax, 100):
                    sys.stdout.flush()
                    if(col["total population"].values[timestep-1] == 0):
                        toGet = col.values[timestep:timestep+future]
                        
                        toGet[:, 0] = np.log10(toGet[:, 0]+1)
                        
                        toCompareTo = col.values[timestep + future]
                        numberIndividualsTimestep = toGet[0, 0]
                        numberIndividualsTimestepPlus99 = toCompareTo[0]
                        deltaNumberIndivs = numberIndividualsTimestepPlus99 - numberIndividualsTimestep
                        didPopulationIncrease = 1 if deltaNumberIndivs > 0 else 0
                        didPopulationPersist = 1 if numberIndividualsTimestepPlus99 > 0 else 0
                        percentageIncrease = numberIndividualsTimestepPlus99 / numberIndividualsTimestep
                        
                        if (didPopulationPersist):
                            arr.append(np.array(toGet))
                        rowsProcessed = rowsProcessed + 1
            sys.stdout.flush()
            tempColumns = resultsColumns[:]  
            
            arr = np.array(arr)
            
            theMeans = np.nanmean(arr, axis=0)
            
            theSDevs = np.nanstd(arr, axis=0)
           
            frMeans = pd.DataFrame(theMeans, columns=tempColumns)
            frSDev = pd.DataFrame(theSDevs, columns=tempColumns)
            
            for columnName in tempColumns:
                lowbounds[columnName] = min(np.nanmin(frMeans[[columnName]].values), lowbounds[columnName])
                upperbounds[columnName] = max(np.nanmax(frMeans[[columnName]].values), upperbounds[columnName])
            
            framesListAllMeans.append(frMeans)
            framesListAllSDevs.append(frSDev)
          
    tempColumns = resultsColumns[:]  
    ynames = tempColumns[:]
    ynames[0] = "Abundance"
    ynames[1] = "Compactness"
    
    byGeneticDiversityDDa = []
    byGeneticDiversityDDa.append(framesListAllMeans[0])
    byGeneticDiversityDDa.append(framesListAllMeans[1])
    byGeneticDiversityDDa.append(framesListAllMeans[2])
    byGeneticDiversityDDa.append(framesListAllMeans[3])
    byGeneticDiversityDDa.append(framesListAllMeans[4])
    
    byGeneticDiversityDNa = []
    byGeneticDiversityDNa.append(framesListAllMeans[5])
    byGeneticDiversityDNa.append(framesListAllMeans[6])
    byGeneticDiversityDNa.append(framesListAllMeans[7])
    byGeneticDiversityDNa.append(framesListAllMeans[8])
    byGeneticDiversityDNa.append(framesListAllMeans[9])
    
    byGeneticDiversityNDa = []
    byGeneticDiversityNDa.append(framesListAllMeans[10])
    byGeneticDiversityNDa.append(framesListAllMeans[11])
    byGeneticDiversityNDa.append(framesListAllMeans[12])
    byGeneticDiversityNDa.append(framesListAllMeans[13])
    byGeneticDiversityNDa.append(framesListAllMeans[14])
    
    byGeneticDiversityNNa = []
    byGeneticDiversityNNa.append(framesListAllMeans[15])
    byGeneticDiversityNNa.append(framesListAllMeans[16])
    byGeneticDiversityNNa.append(framesListAllMeans[17])
    byGeneticDiversityNNa.append(framesListAllMeans[18])
    byGeneticDiversityNNa.append(framesListAllMeans[19])
    
    byGeneticDiversityDDas = []
    byGeneticDiversityDDas.append(framesListAllSDevs[0])
    byGeneticDiversityDDas.append(framesListAllSDevs[1])
    byGeneticDiversityDDas.append(framesListAllSDevs[2])
    byGeneticDiversityDDas.append(framesListAllSDevs[3])
    byGeneticDiversityDDas.append(framesListAllSDevs[4])
    
    byGeneticDiversityDNas = []
    byGeneticDiversityDNas.append(framesListAllSDevs[5])
    byGeneticDiversityDNas.append(framesListAllSDevs[6])
    byGeneticDiversityDNas.append(framesListAllSDevs[7])
    byGeneticDiversityDNas.append(framesListAllSDevs[8])
    byGeneticDiversityDNas.append(framesListAllSDevs[9])
    
    byGeneticDiversityNDas = []
    byGeneticDiversityNDas.append(framesListAllSDevs[10])
    byGeneticDiversityNDas.append(framesListAllSDevs[11])
    byGeneticDiversityNDas.append(framesListAllSDevs[12])
    byGeneticDiversityNDas.append(framesListAllSDevs[13])
    byGeneticDiversityNDas.append(framesListAllSDevs[14])
    
    byGeneticDiversityNNas = []
    byGeneticDiversityNNas.append(framesListAllSDevs[15])
    byGeneticDiversityNNas.append(framesListAllSDevs[16])
    byGeneticDiversityNNas.append(framesListAllSDevs[17])
    byGeneticDiversityNNas.append(framesListAllSDevs[18])
    byGeneticDiversityNNas.append(framesListAllSDevs[19])
    
    
    #these are both 2D lists that contain the data, across run types, by genetic diversity level
    byGeneticDiversitya = []
    byGeneticDiversitya.append(byGeneticDiversityDDa)
    byGeneticDiversitya.append(byGeneticDiversityDNa)
    byGeneticDiversitya.append(byGeneticDiversityNDa)
    byGeneticDiversitya.append(byGeneticDiversityNNa)
    
    byGeneticDiversityas = []
    byGeneticDiversityas.append(byGeneticDiversityDDas)
    byGeneticDiversityas.append(byGeneticDiversityDNas)
    byGeneticDiversityas.append(byGeneticDiversityNDas)
    byGeneticDiversityas.append(byGeneticDiversityNNas)
    
    #need to iterate over all columns and create corresponding charts
    
    for i in range(4):
        byGDa = byGeneticDiversitya[i]
        byGDas = byGeneticDiversityas[i]
        colCompleted = 0
        for columnName in tempColumns:
            plt.clf()
            #get the required column by genetic diversity
            
            label = ynames[colCompleted]
            byGeneticDiversityaFR0 = byGDa[0]
            byGeneticDiversityaFR1 = byGDa[1]
            byGeneticDiversityaFR2 = byGDa[2]
            byGeneticDiversityaFR3 = byGDa[3]
            byGeneticDiversityaFR4 = byGDa[4]
            
            byGeneticDiversityasFR0 = byGDas[0]
            byGeneticDiversityasFR1 = byGDas[1]
            byGeneticDiversityasFR2 = byGDas[2]
            byGeneticDiversityasFR3 = byGDas[3]
            byGeneticDiversityasFR4 = byGDas[4]
            
            #these are now dataframes... get the column I want
            
            byGeneticDiversityaFR0 = byGeneticDiversityaFR0[[columnName]].values
            byGeneticDiversityaFR1 = byGeneticDiversityaFR1[[columnName]].values
            byGeneticDiversityaFR2 = byGeneticDiversityaFR2[[columnName]].values
            byGeneticDiversityaFR3 = byGeneticDiversityaFR3[[columnName]].values
            byGeneticDiversityaFR4 = byGeneticDiversityaFR4[[columnName]].values
            
            byGeneticDiversityasFR0 = byGeneticDiversityasFR0[[columnName]].values
            byGeneticDiversityasFR1 = byGeneticDiversityasFR1[[columnName]].values
            byGeneticDiversityasFR2 = byGeneticDiversityasFR2[[columnName]].values
            byGeneticDiversityasFR3 = byGeneticDiversityasFR3[[columnName]].values
            byGeneticDiversityasFR4 = byGeneticDiversityasFR4[[columnName]].values
            
            ylim=(lowbounds[columnName], upperbounds[columnName])
            fig, ax = plt.subplots()
        
            box = ax.get_position()
            ax.set_position([box.x0 + 0.05, box.y0, box.width * 0.8, box.height])
            
            if colCompleted == 0:
                plt.ylim(0, 6)
            else:
                plt.ylim(ylim)
            plt.xticks(fontsize='large')
            plt.yticks(fontsize='large')
            plt.xlabel('Time Step', fontsize='x-large', fontweight='bold')
            plt.ylabel(label, fontsize='x-large', fontweight='bold')
            
            p1 = ax.plot(range(60), byGeneticDiversityaFR0, color='r')
            p2 = ax.plot(range(60), byGeneticDiversityaFR1, color='y')
            p3 = ax.plot(range(60), byGeneticDiversityaFR2, color='b')
            p4 = ax.plot(range(60), byGeneticDiversityaFR3, color='g')
            p5 = ax.plot(range(60), byGeneticDiversityaFR4, color='m')
            f1 = ax.fill_between(range(60), np.clip((byGeneticDiversityaFR0 - byGeneticDiversityasFR0).flatten(), 0, None), np.clip((byGeneticDiversityaFR0 + byGeneticDiversityasFR0).flatten(), 0, None), color='r', alpha=0.1)
            f2 = ax.fill_between(range(60), np.clip((byGeneticDiversityaFR1 - byGeneticDiversityasFR1).flatten(), 0, None), np.clip((byGeneticDiversityaFR1 + byGeneticDiversityasFR1).flatten(), 0, None), color='y', alpha=0.1)
            f3 = ax.fill_between(range(60), np.clip((byGeneticDiversityaFR2 - byGeneticDiversityasFR2).flatten(), 0, None), np.clip((byGeneticDiversityaFR2 + byGeneticDiversityasFR2).flatten(), 0, None), color='b', alpha=0.1)
            f4 = ax.fill_between(range(60), np.clip((byGeneticDiversityaFR3 - byGeneticDiversityasFR3).flatten(), 0, None), np.clip((byGeneticDiversityaFR3 + byGeneticDiversityasFR3).flatten(), 0, None), color='g', alpha=0.1)
            f5 = ax.fill_between(range(60), np.clip((byGeneticDiversityaFR4 - byGeneticDiversityasFR4).flatten(), 0, None), np.clip((byGeneticDiversityaFR4 + byGeneticDiversityasFR4).flatten(), 0, None), color='m', alpha=0.1)
            ax.legend((p1[0], p2[0], p3[0], p4[0], p5[0]), ('L1', 'L2', 'L3', 'L4', 'L5'), bbox_to_anchor=(1, 0.5), loc='center left')
            plt.savefig('ByGDSTTimeSeries' + columnName + str(i) + '.png')
            colCompleted += 1