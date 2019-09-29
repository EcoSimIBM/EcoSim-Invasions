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
import statsmodels.api as sm
from statsmodels.formula.api import ols
import scipy.stats as stats

np.set_printoptions(threshold=11000)
np.set_printoptions(suppress=True)

#skip the header line
numRuns = 10
numLevels = 5
treatments = ["InvasionsDefaultAcceptors", "InvasionsDefaultAcceptorsNiches", "InvasionsNichesAcceptors", "InvasionsNichesAcceptorsNiches"]
typs = ["NGDT"]
fileN = "Results_Inv_Prey_Female"
resultsColumns = ["total population", "act_EatRatio", "act_EatFailedRatio", "act_ReproduceRatio", "act_ReproduceFailedRatio", "act_EscapeRatio", "act_SearchFoodRatio", "act_SearchFoodFailedRatio", "act_SocializeRatio", "act_SocializeFailedRatio", "act_WaitRatio", "act_Move2StrongestPreyRatio", "act_Move2StrongestPreyCellRatio", "act_Move2StrongestPreyCellFailedRatio", "act_Move2WeakestPreyCellRatio", "act_Move2WeakestPreyCellFailedRatio", "strength", "nbPerCell", "arcs", "MaxSpeed", "EnergyTransferred", "ageDeath"]
#cut columns: act_ExplorationRatio, dist_evolved, act_Move2StrongestPreyFailedRatio, energy, vision, MaxEnergy, speed, eat attempts, reproduction efficiency
startRow = 15100
maxRow = 20100

startOffset = 5
window = 5
future = 60
  
for typ in typs:
    arr = []
    for t in treatments:
        print("Processing treatment " + t)
        #importances list
        importances=[]
        sys.stdout.flush()
        for n in range(1, 1 + numLevels):
            sys.stdout.flush()
            print("\tProcessing level " + str(n))
            #make an output file
            #iterate over all runs
            for r in range(1, 1 + numRuns):
                print("\t\tProcessing run number " + str(r))
                #load up the file
                df = pd.read_csv(t + str(typ) + str(n) + "/run" + str(r) + "/" + fileN + ".csv", sep=', ')
                col = df[resultsColumns]
                actualMax = min(df.shape[0] - (future + 1), maxRow-1)
                rowsProcessed = 0
                
                for timestep in range(startRow-1, actualMax, 100):
                    sys.stdout.flush()
                    if(col["total population"].values[timestep-1] == 0):
                        #get the row I need with all the data I want (above)
                        toGet = col.values[timestep+startOffset:timestep+startOffset+window]
                        averages = np.nanmean(toGet, axis=0)
                        
                        eatAttempts = averages[1] + averages[2]
                        eatEfficiency = averages[1] / eatAttempts
                        reprodAttempts = averages[3] + averages[4]
                        reprodEfficiency = averages[3] / reprodAttempts
                        averages[1] = eatEfficiency
                        averages[2] = eatAttempts
                        averages[3] = reprodEfficiency
                        averages[4] = reprodAttempts
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
    tempColumns[1] = "eatEfficiency"
    tempColumns[2] = "eatAttempts"
    tempColumns[3] = "reprodEfficiency"
    tempColumns[4] = "reprodAttempts"
    tempColumns.append("class")
    tempColumns.append("percentageIncrease")
    fr = pd.DataFrame(arr)
    fr.columns = tempColumns
    fr = fr.dropna()
    mlX = fr.values[:,1:-2]
    mly = fr.values[:,-2]
    X_train, X_test, y_train, y_test = train_test_split(mlX, mly, test_size=0.3, random_state=42)
    
    frColumns = resultsColumns[1:]
    frColumns[0] = "eatEfficiency"
    frColumns[1] = "eatAttempts"
    frColumns[2] = "reprodEfficiency"
    frColumns[3] = "reprodAttempts"
    frColumns.append('class')
    
    toFrame = np.concatenate((X_train, np.expand_dims(y_train, axis=1)), axis=1)
    frTr = pd.DataFrame(toFrame, columns=frColumns)
    
    #first get the samples of class 0
    c1 = frTr[frTr['class']== 0]
    #then get the samples of class 1
    c2 = frTr[frTr['class'] == 1]
    XTr1 = c1
    XTr2 = c2
    XY_train = pd.concat([XTr1, XTr2], axis=0)
    X_train = XY_train.values[:, :-1]
    y_train = XY_train.values[:,-1]
    
    print(X_train)
    
    rho, pvals = stats.spearmanr(X_train)
    
    rhoDF = pd.DataFrame(rho, columns=frColumns[:-1], index=frColumns[:-1])
    pvals = pd.DataFrame(pvals, columns=frColumns[:-1], index=frColumns[:-1])
    
    rhoDF.to_csv("SpearmansRhoDF.csv")
    pvals.to_csv("SpearmansPvals.csv")
    
    print("rho:")
    print(rho)
    print("pvals:")
    print(pvals)
    
    print("X_train.shape:")
    print(X_train.shape)
    print("y_train.shape:")
    print(y_train.shape)
    
    tc = frColumns[:-1]
    
    X_train = pd.DataFrame(X_train, columns=tc)
    y_train = pd.DataFrame(y_train)
    
    X_train = X_train.drop(['reprodEfficiency', 'eatAttempts'], axis=1)
    
    X_test = pd.DataFrame(X_test, columns=tc)
    X_test = X_test.drop(['reprodEfficiency', 'eatAttempts'], axis=1)
    
    base_rf = RandomForestClassifier(n_estimators=1000, n_jobs=-1, oob_score=True) #max_depth=4, 
    clf = clone(base_rf)
    clf.fit(X_train, y_train.values.ravel())
    
    print("ELI5 FEATURE IMPORTANCE:")
    perm = PermutationImportance(clf, n_iter=10, scoring=metrics.make_scorer(metrics.roc_auc_score, needs_proba=True)).fit(X_test, y_test)
    
    print("perm.results")
    res = np.array(perm.results_)
    print(eli5.format_as_text(eli5.explain_weights(perm, feature_names = X_train.columns.tolist())))
    importances.append(res)
    print("Training score:")
    print(clf.score(X_train, y_train))
    print("Testing score:")
    print(clf.score(X_test, y_test))
    
    predictions = clf.predict(X_test)
    cm = metrics.confusion_matrix(y_test, predictions)
    print("Confusion Matrix:")
    print(cm)
    print("Accuracy:",metrics.accuracy_score(y_test, predictions))
    print("Precision:",metrics.precision_score(y_test, predictions))
    print("Recall:",metrics.recall_score(y_test, predictions))
    
    y_pred_proba = clf.predict_proba(X_test)[::,1]
    fpr, tpr, _ = metrics.roc_curve(y_test,  y_pred_proba)
    auc = metrics.roc_auc_score(y_test, y_pred_proba)
    f1Score = metrics.f1_score(y_test, predictions) 
    plt.plot(fpr,tpr,label="data 1, auc="+str(auc))
    plt.legend(loc=4)
    print("auc:")
    print(auc)
    
    print("f1Score:")
    print(f1Score)
    plt.savefig("ROC" + str(t) + "-" + str(n) + ".png", format="png")