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
typerinos = ["NGDT"]
fileN = "Results_Inv_Prey_Female"
resultsColumns = ["total population", "nbPerCell"]                    
startRow = 15100
maxRow = 20100

startOffset = 5
window = 5
future = 60

numberOfComparisonsp1=0
numberOfComparisonsp05=0
numberOfComparisonsp01=0

lowbounds = {}
upperbounds = {}

tempColumns = resultsColumns[:]  
tempColumns.append("class")
tempColumns.append("percentageIncrease")

for columnName in tempColumns:
    lowbounds[columnName] = 1000000
    upperbounds[columnName] = -1000000
 
def custom_legend(colors, labels, linestyles=None):
        """ Creates a list of matplotlib Patch objects that can be passed to the legend(...) function to create a custom
            legend.

        :param colors: A list of colors, one for each entry in the legend. You can also include a linestyle, for example: 'k--'
        :param labels:  A list of labels, one for each entry in the legend.
        """

        if linestyles is not None:
            assert len(linestyles) == len(colors), "Length of linestyles must match length of colors."

        h = list()
        for k,(c,l) in enumerate(zip(colors, labels)):
            clr = c
            ls = 'solid'
            if linestyles is not None:
                ls = linestyles[k]
            patch = matplotlib.patches.Patch(color=clr, label=l, linestyle=ls)
            h.append(patch)
        return h
        
def grouped_boxplot(data, group_names=None, subgroup_names=None, ax=None, subgroup_colors=None, box_width=0.6, box_spacing=1.0, ylims=(0, 1)):
    """ Draws a grouped boxplot. The data should be organized in a hierarchy, where there are multiple
        subgroups for each main group.

    :param data: A dictionary of length equal to the number of the groups. The key should be the
                group name, the value should be a list of arrays. The length of the list should be
                equal to the number of subgroups.
    :param group_names: (Optional) The group names, should be the same as data.keys(), but can be ordered.
    :param subgroup_names: (Optional) Names of the subgroups.
    :param subgroup_colors: A list specifying the plot color for each subgroup.
    :param ax: (Optional) The axis to plot on.
    """
    
    #use the keys from the dict if none are explicity provided
    if group_names is None:
        group_names = data.keys()

    #get current axis if none is provided
    if ax is None:
        ax = plt.gca()
    plt.sca(ax)
    
    plt.ylim(ylims)
    box = ax.get_position()
    ax.set_position([box.x0 + 0.05, box.y0, box.width * 0.8, box.height])

    #number of subgroups
    nsubgroups = np.array([len(v) for v in data.values()])
    #ensure that the number of subgroups matches across all data
    assert len(np.unique(nsubgroups)) == 1, "Number of subgroups for each property differ!"
    nsubgroups = nsubgroups[0]

    #if a subgroup_colors list is provided, use it.. otherwise, initialize with random colours
    if subgroup_colors is None:
        subgroup_colors = list()
        for k in range(nsubgroups):
            subgroup_colors.append(np.random.rand(3))
    else:
        assert len(subgroup_colors) == nsubgroups, "subgroup_colors length must match number of subgroups (%d)" % nsubgroups

    def _decorate_box(_bp, _d):
        plt.setp(_bp['boxes'], lw=0, color='k')
        plt.setp(_bp['whiskers'], lw=3.0, color='k')

        # fill in each box with a color
        assert len(_bp['boxes']) == nsubgroups
        for _k,_box in enumerate(_bp['boxes']):
            _boxX = list()
            _boxY = list()
            for _j in range(5):
                _boxX.append(_box.get_xdata()[_j])
                _boxY.append(_box.get_ydata()[_j])
            _boxCoords = list(zip(_boxX, _boxY))
            _boxPolygon = plt.Polygon(_boxCoords, facecolor=subgroup_colors[_k])
            ax.add_patch(_boxPolygon)

        # draw a black line for the median
        for _k,_med in enumerate(_bp['medians']):
            _medianX = list()
            _medianY = list()
            for _j in range(2):
                _medianX.append(_med.get_xdata()[_j])
                _medianY.append(_med.get_ydata()[_j])
                plt.plot(_medianX, _medianY, 'k', linewidth=3.0)

            # draw a black dot for the mean
            plt.plot([np.mean(_med.get_xdata())], [np.mean(_d[_k])], color='w', marker='.',
                      markeredgecolor='k', markersize=16)
        
    cpos = 1
    label_pos = list()
    for k in group_names:
        d = data[k]
        nsubgroups = len(d)
        pos = np.arange(nsubgroups) + cpos
        label_pos.append(pos.mean())
        bp = plt.boxplot(d, positions=pos, widths=box_width)
        _decorate_box(bp, d)
        cpos += nsubgroups + box_spacing

    plt.xlim(0, cpos-1)
    plt.xticks(label_pos, group_names)

    if subgroup_names is not None:
        leg = custom_legend(subgroup_colors, subgroup_names)
        plt.legend(handles=leg, bbox_to_anchor=(1, 0.5), loc='center left')

 
for typerino in typerinos:
    #iterate over all treatments
    framesList0 = []
    framesList1 = []
    framesListAll = []
    for t in treatments:
        print("Processing treatment " + t)
        #importances list
        importances=[]
        sys.stdout.flush()
        for n in range(1, 1 + numLevels):
            sys.stdout.flush()
            print("\tProcessing level " + str(n))
            #make a numpy array
            arr = []
            #make an output file
            #iterate over all runs
            for r in range(1, 1 + numRuns):
                print("\t\tProcessing run number " + str(r))
                df = pd.read_csv(t + str(typerino) + str(n) + "/run" + str(r) + "/" + fileN + ".csv", sep=', ')
                col = df[resultsColumns]
                actualMax = min(df.shape[0] - (future + 1), maxRow-1)
                rowsProcessed = 0
                for timestep in range(startRow-1, actualMax, 100):
                    sys.stdout.flush()
                    if(col["total population"].values[timestep-1] == 0):
                        toGet = col.values[timestep+startOffset:timestep+startOffset+window]
                        averages = np.nanmean(toGet, axis=0)
                        toCompareTo = col.values[timestep + future]
                        numberIndividualsTimestep = toGet[0, 0]
                        numberIndividualsTimestepPlus99 = toCompareTo[0]
                        deltaNumberIndivs = numberIndividualsTimestepPlus99 - numberIndividualsTimestep
                        entropyTimestep = toGet[0, 1]
                        entropyTimestepPlus99 = toCompareTo[1]
                        didPopulationIncrease = 1 if deltaNumberIndivs > 0 else 0
                        didPopulationPersist = 1 if numberIndividualsTimestepPlus99 > 0 else 0
                        percentageIncrease = numberIndividualsTimestepPlus99 / numberIndividualsTimestep
                        toGet = np.ndarray.tolist(averages)
                        toGet.append(didPopulationPersist)
                        toGet.append(percentageIncrease)
                        arr.append(np.array(toGet))
                        rowsProcessed = rowsProcessed + 1
            sys.stdout.flush()
            tempColumns = resultsColumns[:]
            tempColumns.append("class")
            tempColumns.append("percentageIncrease")
            fr = pd.DataFrame(arr)
            fr.columns = tempColumns
            fr = fr.dropna()
            
            
            for columnName in tempColumns:
                lowbounds[columnName] = min(fr[[columnName]].values.min(), lowbounds[columnName])
                upperbounds[columnName] = max(fr[[columnName]].values.max(), upperbounds[columnName])
            
            need0 = fr[fr["class"] == 0]
            framesList0.append(need0)
            need1 = fr[fr["class"] == 1]
            
            framesList1.append(need1)
            framesListAll.append(fr)
    
    
    tempColumns = resultsColumns[:]  
    tempColumns.append("class")
    tempColumns.append("percentageIncrease")
    
    labelList = ["Abundance", "Compactness"]                    

    labelList.append("class")
    labelList.append("percentage increase")
    
    byGeneticDiversityDD0 = []
    byGeneticDiversityDD0.append(framesList0[0])
    byGeneticDiversityDD0.append(framesList0[1])
    byGeneticDiversityDD0.append(framesList0[2])
    byGeneticDiversityDD0.append(framesList0[3])
    byGeneticDiversityDD0.append(framesList0[4])
    
    byGeneticDiversityDN0 = []
    byGeneticDiversityDN0.append(framesList0[5])
    byGeneticDiversityDN0.append(framesList0[6])
    byGeneticDiversityDN0.append(framesList0[7])
    byGeneticDiversityDN0.append(framesList0[8])
    byGeneticDiversityDN0.append(framesList0[9])
    
    byGeneticDiversityND0 = []
    byGeneticDiversityND0.append(framesList0[10])
    byGeneticDiversityND0.append(framesList0[11])
    byGeneticDiversityND0.append(framesList0[12])
    byGeneticDiversityND0.append(framesList0[13])
    byGeneticDiversityND0.append(framesList0[14])
    
    byGeneticDiversityNN0 = []
    byGeneticDiversityNN0.append(framesList0[15])
    byGeneticDiversityNN0.append(framesList0[16])
    byGeneticDiversityNN0.append(framesList0[17])
    byGeneticDiversityNN0.append(framesList0[18])
    byGeneticDiversityNN0.append(framesList0[19])
    
    byGeneticDiversityDD1 = []
    byGeneticDiversityDD1.append(framesList1[0])
    byGeneticDiversityDD1.append(framesList1[1])
    byGeneticDiversityDD1.append(framesList1[2])
    byGeneticDiversityDD1.append(framesList1[3])
    byGeneticDiversityDD1.append(framesList1[4])
    
    byGeneticDiversityDN1 = []
    byGeneticDiversityDN1.append(framesList1[5])
    byGeneticDiversityDN1.append(framesList1[6])
    byGeneticDiversityDN1.append(framesList1[7])
    byGeneticDiversityDN1.append(framesList1[8])
    byGeneticDiversityDN1.append(framesList1[9])
    
    byGeneticDiversityND1 = []
    byGeneticDiversityND1.append(framesList1[10])
    byGeneticDiversityND1.append(framesList1[11])
    byGeneticDiversityND1.append(framesList1[12])
    byGeneticDiversityND1.append(framesList1[13])
    byGeneticDiversityND1.append(framesList1[14])
    
    byGeneticDiversityNN1 = []
    byGeneticDiversityNN1.append(framesList1[15])
    byGeneticDiversityNN1.append(framesList1[16])
    byGeneticDiversityNN1.append(framesList1[17])
    byGeneticDiversityNN1.append(framesList1[18])
    byGeneticDiversityNN1.append(framesList1[19])
    
    byGeneticDiversityDDa = []
    byGeneticDiversityDDa.append(framesListAll[0])
    byGeneticDiversityDDa.append(framesListAll[1])
    byGeneticDiversityDDa.append(framesListAll[2])
    byGeneticDiversityDDa.append(framesListAll[3])
    byGeneticDiversityDDa.append(framesListAll[4])
    
    byGeneticDiversityDNa = []
    byGeneticDiversityDNa.append(framesListAll[5])
    byGeneticDiversityDNa.append(framesListAll[6])
    byGeneticDiversityDNa.append(framesListAll[7])
    byGeneticDiversityDNa.append(framesListAll[8])
    byGeneticDiversityDNa.append(framesListAll[9])
    
    byGeneticDiversityNDa = []
    byGeneticDiversityNDa.append(framesListAll[10])
    byGeneticDiversityNDa.append(framesListAll[11])
    byGeneticDiversityNDa.append(framesListAll[12])
    byGeneticDiversityNDa.append(framesListAll[13])
    byGeneticDiversityNDa.append(framesListAll[14])
    
    byGeneticDiversityNNa = []
    byGeneticDiversityNNa.append(framesListAll[15])
    byGeneticDiversityNNa.append(framesListAll[16])
    byGeneticDiversityNNa.append(framesListAll[17])
    byGeneticDiversityNNa.append(framesListAll[18])
    byGeneticDiversityNNa.append(framesListAll[19])
    
    
    AllSSframes0 = framesList0[0]
    AllSSframes0 = pd.concat([AllSSframes0, framesList0[1]], sort=False)
    AllSSframes0 = pd.concat([AllSSframes0, framesList0[2]], sort=False)
    AllSSframes0 = pd.concat([AllSSframes0, framesList0[3]], sort=False)
    AllSSframes0 = pd.concat([AllSSframes0, framesList0[4]], sort=False)
    AllNSframes0 = framesList0[5]
    AllNSframes0 = pd.concat([AllNSframes0, framesList0[6]], sort=False)
    AllNSframes0 = pd.concat([AllNSframes0, framesList0[7]], sort=False)
    AllNSframes0 = pd.concat([AllNSframes0, framesList0[8]], sort=False)
    AllNSframes0 = pd.concat([AllNSframes0, framesList0[9]], sort=False)
    AllSNframes0 = framesList0[10]
    AllSNframes0 = pd.concat([AllSNframes0, framesList0[11]], sort=False)
    AllSNframes0 = pd.concat([AllSNframes0, framesList0[12]], sort=False)
    AllSNframes0 = pd.concat([AllSNframes0, framesList0[13]], sort=False)
    AllSNframes0 = pd.concat([AllSNframes0, framesList0[14]], sort=False)
    AllNNframes0 = framesList0[15]
    AllNNframes0 = pd.concat([AllNNframes0, framesList0[16]], sort=False)
    AllNNframes0 = pd.concat([AllNNframes0, framesList0[17]], sort=False)
    AllNNframes0 = pd.concat([AllNNframes0, framesList0[18]], sort=False)
    AllNNframes0 = pd.concat([AllNNframes0, framesList0[19]], sort=False)
    
    AllSSframes1 = framesList1[0]
    AllSSframes1 = pd.concat([AllSSframes1, framesList1[1]], sort=False)
    AllSSframes1 = pd.concat([AllSSframes1, framesList1[2]], sort=False)
    AllSSframes1 = pd.concat([AllSSframes1, framesList1[3]], sort=False)
    AllSSframes1 = pd.concat([AllSSframes1, framesList1[4]], sort=False)
    AllNSframes1 = framesList1[5]
    AllNSframes1 = pd.concat([AllNSframes1, framesList1[6]], sort=False)
    AllNSframes1 = pd.concat([AllNSframes1, framesList1[7]], sort=False)
    AllNSframes1 = pd.concat([AllNSframes1, framesList1[8]], sort=False)
    AllNSframes1 = pd.concat([AllNSframes1, framesList1[9]], sort=False)
    AllSNframes1 = framesList1[10]
    AllSNframes1 = pd.concat([AllSNframes1, framesList1[11]], sort=False)
    AllSNframes1 = pd.concat([AllSNframes1, framesList1[12]], sort=False)
    AllSNframes1 = pd.concat([AllSNframes1, framesList1[13]], sort=False)
    AllSNframes1 = pd.concat([AllSNframes1, framesList1[14]], sort=False)
    AllNNframes1 = framesList1[15]
    AllNNframes1 = pd.concat([AllNNframes1, framesList1[16]], sort=False)
    AllNNframes1 = pd.concat([AllNNframes1, framesList1[17]], sort=False)
    AllNNframes1 = pd.concat([AllNNframes1, framesList1[18]], sort=False)
    AllNNframes1 = pd.concat([AllNNframes1, framesList1[19]], sort=False)

    byGeneticDiversity0 = []
    byGeneticDiversity1 = []
    byGeneticDiversitya = []
    byGeneticDiversity0.append(byGeneticDiversityDD0)
    byGeneticDiversity0.append(byGeneticDiversityDN0)
    byGeneticDiversity0.append(byGeneticDiversityND0)
    byGeneticDiversity0.append(byGeneticDiversityNN0)
    byGeneticDiversity1.append(byGeneticDiversityDD1)
    byGeneticDiversity1.append(byGeneticDiversityDN1)
    byGeneticDiversity1.append(byGeneticDiversityND1)
    byGeneticDiversity1.append(byGeneticDiversityNN1)
    byGeneticDiversitya.append(byGeneticDiversityDDa)
    byGeneticDiversitya.append(byGeneticDiversityDNa)
    byGeneticDiversitya.append(byGeneticDiversityNDa)
    byGeneticDiversitya.append(byGeneticDiversityNNa)
    
    for i in range(4):
        byGD0 = byGeneticDiversity0[i]
        byGD1 = byGeneticDiversity1[i]
        byGDa = byGeneticDiversitya[i]
        
        colnum=0
        for columnName in tempColumns:
            label = labelList[colnum]
            colnum += 1
            plt.clf()
            byGeneticDiversity0FR0 = byGD0[0]
            byGeneticDiversity0FR1 = byGD0[1]
            byGeneticDiversity0FR2 = byGD0[2]
            byGeneticDiversity0FR3 = byGD0[3]
            byGeneticDiversity0FR4 = byGD0[4]
            
            byGeneticDiversity1FR0 = byGD1[0]
            byGeneticDiversity1FR1 = byGD1[1]
            byGeneticDiversity1FR2 = byGD1[2]
            byGeneticDiversity1FR3 = byGD1[3]
            byGeneticDiversity1FR4 = byGD1[4]
            
            byGeneticDiversityaFR0 = byGDa[0]
            byGeneticDiversityaFR1 = byGDa[1]
            byGeneticDiversityaFR2 = byGDa[2]
            byGeneticDiversityaFR3 = byGDa[3]
            byGeneticDiversityaFR4 = byGDa[4]
            
            byGeneticDiversity0FR0 = byGeneticDiversity0FR0[[columnName]].values
            byGeneticDiversity0FR1 = byGeneticDiversity0FR1[[columnName]].values
            byGeneticDiversity0FR2 = byGeneticDiversity0FR2[[columnName]].values
            byGeneticDiversity0FR3 = byGeneticDiversity0FR3[[columnName]].values
            byGeneticDiversity0FR4 = byGeneticDiversity0FR4[[columnName]].values
            
            byGeneticDiversity1FR0 = byGeneticDiversity1FR0[[columnName]].values
            byGeneticDiversity1FR1 = byGeneticDiversity1FR1[[columnName]].values
            byGeneticDiversity1FR2 = byGeneticDiversity1FR2[[columnName]].values
            byGeneticDiversity1FR3 = byGeneticDiversity1FR3[[columnName]].values
            byGeneticDiversity1FR4 = byGeneticDiversity1FR4[[columnName]].values
            
            byGeneticDiversityaFR0 = byGeneticDiversityaFR0[[columnName]].values
            byGeneticDiversityaFR1 = byGeneticDiversityaFR1[[columnName]].values
            byGeneticDiversityaFR2 = byGeneticDiversityaFR2[[columnName]].values
            byGeneticDiversityaFR3 = byGeneticDiversityaFR3[[columnName]].values
            byGeneticDiversityaFR4 = byGeneticDiversityaFR4[[columnName]].values
            
            ph0 = [np.ravel(byGeneticDiversity0FR0).tolist(), np.ravel(byGeneticDiversity0FR1).tolist(), np.ravel(byGeneticDiversity0FR2).tolist(), np.ravel(byGeneticDiversity0FR3).tolist(), np.ravel(byGeneticDiversity0FR4).tolist()]
            ph1 = [np.ravel(byGeneticDiversity1FR0).tolist(), np.ravel(byGeneticDiversity1FR1).tolist(), np.ravel(byGeneticDiversity1FR2).tolist(), np.ravel(byGeneticDiversity1FR3).tolist(), np.ravel(byGeneticDiversity1FR4).tolist()]
            pha = [np.ravel(byGeneticDiversityaFR0).tolist(), np.ravel(byGeneticDiversityaFR1).tolist(), np.ravel(byGeneticDiversityaFR2).tolist(), np.ravel(byGeneticDiversityaFR3).tolist(), np.ravel(byGeneticDiversityaFR4).tolist()]
            
            try:
                f_value, p_value = stats.kruskal(byGeneticDiversity0FR0, byGeneticDiversity0FR1, byGeneticDiversity0FR2, byGeneticDiversity0FR3, byGeneticDiversity0FR4)
                if (p_value <=0.1):
                    numberOfComparisonsp1 += 1
                    if (p_value <=0.05):
                        numberOfComparisonsp05 += 1
                    if (p_value <=0.01):
                        numberOfComparisonsp01 += 1
                    print("KW for failures in column " + str(columnName) + " RT " + str(i) + "*** " + str(p_value))
                    try:
                        pvals = sp.posthoc_conover(ph0, p_adjust='holm')
                        truth = np.logical_and(pvals <= 0.1, pvals >= 0)
                        if (np.any(truth)):
                            print("significant comparison:")
                            print(pvals)
                    except Exception as e:
                        print('Could not compute posthoc conover: ' + str(e))
            except:
                print("")
            try:
                f_value, p_value = stats.kruskal(byGeneticDiversity1FR0, byGeneticDiversity1FR1, byGeneticDiversity1FR2, byGeneticDiversity1FR3, byGeneticDiversity1FR4)
                if (p_value <=0.1):
                    numberOfComparisonsp1 += 1
                    if (p_value <=0.05):
                        numberOfComparisonsp05 += 1
                    if (p_value <=0.01):
                        numberOfComparisonsp01 += 1
                    print("KW for successes in column " + str(columnName) + " RT " + str(i) + "*** " + str(p_value))
                    try:
                        pvals = sp.posthoc_conover(ph1, p_adjust='holm')
                        truth = np.logical_and(pvals <= 0.1, pvals >= 0)
                        if (np.any(truth)):
                            print("significant comparison:")
                            print(pvals)
                    except Exception as e:
                        print('Could not compute posthoc conover: ' + str(e))
            except:
                print("")
            try:
                f_value, p_value = stats.kruskal(byGeneticDiversityaFR0, byGeneticDiversityaFR1, byGeneticDiversityaFR2, byGeneticDiversityaFR3, byGeneticDiversityaFR4)
                if (p_value <=0.1):
                    numberOfComparisonsp1 += 1
                    if (p_value <=0.05):
                        numberOfComparisonsp05 += 1
                    if (p_value <=0.01):
                        numberOfComparisonsp01 += 1
                    print("KW for all in column " + str(columnName) + " RT " + str(i) + "*** " + str(p_value))
                    try:
                        pvals = sp.posthoc_conover(pha, p_adjust='holm')
                        truth = np.logical_and(pvals <= 0.1, pvals >= 0)
                        if (np.any(truth)):
                            print("significant comparison:")
                            print(pvals)
                    except Exception as e:
                        print('Could not compute posthoc conover: ' + str(e))
            except:
                print("")
            dataByGD = {
                'L1': [byGeneticDiversity0FR0, byGeneticDiversity1FR0, byGeneticDiversityaFR0],
                'L2': [byGeneticDiversity0FR1, byGeneticDiversity1FR1, byGeneticDiversityaFR1],
                'L3': [byGeneticDiversity0FR2, byGeneticDiversity1FR2, byGeneticDiversityaFR2],
                'L4': [byGeneticDiversity0FR3, byGeneticDiversity1FR3, byGeneticDiversityaFR3],
                'L5': [byGeneticDiversity0FR4, byGeneticDiversity1FR4, byGeneticDiversityaFR4]
            } 
            
            ylim=(lowbounds[columnName], upperbounds[columnName])
            grouped_boxplot(dataByGD, group_names=['L1', 'L2', 'L3', 'L4', 'L5'], subgroup_names=['Fail', 'Success', 'All'], subgroup_colors=['red', 'blue', 'magenta'], ylims=ylim)
            
            plt.xticks(fontsize='large')
            plt.yticks(fontsize='large')
            plt.xlabel('Genetic Diversity Level', fontsize='x-large', fontweight='bold')
            plt.ylabel(label, fontsize='x-large', fontweight='bold')
            
            plt.savefig('ByGD' + columnName + str(i) + '.png')
    
    colnum=0
    for columnName in tempColumns:
    
        label = labelList[colnum]
        colnum += 1
        plt.clf()
        byRunType0FR0 = AllSSframes0
        byRunType0FR1 = AllNSframes0
        byRunType0FR2 = AllSNframes0
        byRunType0FR3 = AllNNframes0
        
        byRunType1FR0 = AllSSframes1
        byRunType1FR1 = AllNSframes1
        byRunType1FR2 = AllSNframes1
        byRunType1FR3 = AllNNframes1
        
        byRunType0FR0 = byRunType0FR0[[columnName]].values
        byRunType0FR1 = byRunType0FR1[[columnName]].values
        byRunType0FR2 = byRunType0FR2[[columnName]].values
        byRunType0FR3 = byRunType0FR3[[columnName]].values
        
        byRunType1FR0 = byRunType1FR0[[columnName]].values
        byRunType1FR1 = byRunType1FR1[[columnName]].values
        byRunType1FR2 = byRunType1FR2[[columnName]].values
        byRunType1FR3 = byRunType1FR3[[columnName]].values
        
        ph0 = [np.ravel(byRunType0FR0).tolist(), np.ravel(byRunType0FR1).tolist(), np.ravel(byRunType0FR2).tolist(), np.ravel(byRunType0FR3).tolist()]
        ph1 = [np.ravel(byRunType1FR0).tolist(), np.ravel(byRunType1FR1).tolist(), np.ravel(byRunType1FR2).tolist(), np.ravel(byRunType1FR3).tolist()]
            
        
        try:
            f_value, p_value = stats.kruskal(byRunType0FR0, byRunType0FR1, byRunType0FR2, byRunType0FR3)
            if (p_value <=0.1):
                numberOfComparisonsp1 += 1
                if (p_value <=0.05):
                    numberOfComparisonsp05 += 1
                if (p_value <=0.01):
                    numberOfComparisonsp01 += 1
                print("KW for failures in column " + str(columnName) + " RT " + str(i) + "*** " + str(p_value))
                
                try:
                    pvals = sp.posthoc_conover(ph0, p_adjust='holm')
                    truth = np.logical_and(pvals <= 0.1, pvals >= 0)
                    if (np.any(truth)):
                        print("significant comparison:")
                        print(pvals)
                except Exception as e:
                    print('Could not compute posthoc conover: ' + str(e))
                
        except:
            print("")
        try:
            f_value, p_value = stats.kruskal(byRunType1FR0, byRunType1FR1, byRunType1FR2, byRunType1FR3)
            if (p_value <=0.1):
                numberOfComparisonsp1 += 1
                if (p_value <=0.05):
                    numberOfComparisonsp05 += 1
                if (p_value <=0.01):
                    numberOfComparisonsp01 += 1
                print("KW for successes in column " + str(columnName) + " RT " + str(i) + "*** " + str(p_value))
                
                try:
                    pvals = sp.posthoc_conover(ph1, p_adjust='holm')
                    truth = np.logical_and(pvals <= 0.1, pvals >= 0)
                    if (np.any(truth)):
                        print("significant comparison:")
                        print(pvals)
                except Exception as e:
                    print('Could not compute posthoc conover: ' + str(e))
               
        except:
            print("")
        
        dataByRunType = {
            'S->S': [byRunType0FR0, byRunType1FR0],
            'N->S': [byRunType0FR1, byRunType1FR1],
            'S->N': [byRunType0FR2, byRunType1FR2],
            'N->N': [byRunType0FR3, byRunType1FR3]
        
        }
        
        diff = upperbounds[columnName] - lowbounds[columnName]
        ylim=(lowbounds[columnName] - diff * 0.1, upperbounds[columnName] + diff * 0.1)
        
        grouped_boxplot(dataByRunType, group_names=['S->S', 'N->S', 'S->N', 'N->N'], subgroup_names=['Fail', 'Success'], subgroup_colors=['red', 'blue'], ylims=ylim)
        
        plt.xticks(fontsize='large')
        plt.yticks(fontsize='large')
        plt.xlabel('Run Type', fontsize='x-large', fontweight='bold')
        plt.ylabel(label, fontsize='x-large', fontweight='bold')
        
        plt.savefig('ByRT' + columnName + '.png')

print('numberOfComparisonsp1')
print(numberOfComparisonsp1)
print('numberOfComparisonsp05')
print(numberOfComparisonsp05)
print('numberOfComparisonsp01')
print(numberOfComparisonsp01)            