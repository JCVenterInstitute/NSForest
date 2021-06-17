def NS_Forest(adata, clusterLabelcolumnHeader = "louvain", rfTrees = 1000, Median_Expression_Level = 0, Genes_to_testing = 6, betaValue = 0.5):
    #adata = scanpy object
    #rfTrees = Number of trees
    #Median_Expression_Level = median expression level for removing negative markers
    #Genes_to_testing = How many top genes ranked by binary score will be evaluated in permutations by fbeta-score (as the number increases the number of permutation rises exponentially!)
    #betaValue = Set values for fbeta weighting. 1 is default f-measure. close to zero is Precision, greater than 1 weights toward Recall

    #libraries
    import numpy as np
    import pandas as pd
    import numexpr
    import itertools
    from subprocess import call
    import scanpy as sc
    from sklearn.ensemble import RandomForestClassifier
    from sklearn.tree import DecisionTreeClassifier
    from sklearn import tree
    from sklearn.metrics import fbeta_score
    from sklearn.metrics import accuracy_score
    from sklearn.metrics import confusion_matrix
    import graphviz
    import time
    
    # Functions
    def randomForest(adata,dataDummy,column,rfTrees,threads): #Runs Random forest on the binary dummy variables; outputs all genes ranked Gini Index
        x_train = adata.X
        names = adata.var_names
        y_train = dataDummy[column]
        rf = RandomForestClassifier(n_estimators=rfTrees, n_jobs=threads, random_state=123456)
        rf.fit(x_train, y_train)
        Ranked_Features = sorted(zip([round(x, 8) for x in rf.feature_importances_], names),reverse=True)
        return Ranked_Features    
        
    def rankInformative(Ranked_Features,column,rankedDict,howManyInformativeGenes2test): #subsets list according to howManyInformativeGenes2test parameter
        RankedList = []
        midcounter = 0
        for x in Ranked_Features:
            midcounter +=1
            RankedList.append(x[1])
            if midcounter==howManyInformativeGenes2test:
                break
        rankedDict[column] = RankedList
        return RankedList             

    def negativeOut(x, column,medianValues,Median_Expression_Level): # Removes genes with median expression < Median_Expression_Level parameter
        Positive_RankedList_Complete = []
        for i in x:
            if medianValues.loc[column, i] > Median_Expression_Level:
                print(i)
                print(medianValues.loc[column, i])
                Positive_RankedList_Complete.append(i)
            else:
                print(i)
                print(medianValues.loc[column, i])
                print("Is Right Out!")
        return Positive_RankedList_Complete

    def binaryScore(Positive_RankedList_Complete, InformativeGenes, medianValues, column, clusters2Loop, Ranked_Features, Genes_to_testing,Binary_store_DF): # Takes top ranked positive genes (number according to Genes_to_testing) and computes Binary score for each gene
        Positive_RankedList = list(Positive_RankedList_Complete[0:InformativeGenes])
        Median_RF_Subset = medianValues.loc[:, Positive_RankedList]
        Rescaled_Matrix = pd.DataFrame()

        for i in Positive_RankedList:
            Target_value = medianValues.loc[column, i]
            Rescaled_values = Median_RF_Subset[[i]].divide(Target_value)
            Rescaled_Matrix = pd.concat([Rescaled_Matrix,Rescaled_values],axis=1)
        difference_matrix = Rescaled_Matrix.apply(lambda x: 1-x, axis=1)
        difference_matrix_clean1 = difference_matrix.where(difference_matrix >= 0,other=0)
        difference_matrix_clean = difference_matrix_clean1.where(difference_matrix > 0, 0)
        ColumnSums = difference_matrix_clean.sum(0)
        rescaled = ColumnSums/clusters2Loop
       
        # Double sort so that for ties, the RF ranking prevails!     
        Ranked_Features_df = pd.DataFrame(Ranked_Features)
        Ranked_Features_df.rename(columns={1: 'Symbol'}, inplace=True)
        Ranked_Features_df_indexed = Ranked_Features_df.set_index("Symbol")
        rescaled_df = pd.DataFrame(rescaled)
        binaryAndinformation_Ranks = rescaled_df.join(Ranked_Features_df_indexed,lsuffix='_scaled', rsuffix='_informationGain')
        binaryAndinformation_Ranks.sort_values(by=['0_scaled','0_informationGain'],ascending= [False, False], inplace = True)
        Binary_ranked_Genes = binaryAndinformation_Ranks.index.tolist()
        Binary_RankedList = list(Binary_ranked_Genes[0:Genes_to_testing])
        Binary_scores = rescaled.to_dict()
        Binary_store_DF = Binary_store_DF.append(binaryAndinformation_Ranks)
        
        return Binary_RankedList,Binary_store_DF


    def DT_cutOffs(x, column, dataDummy): # For each gene in the top binary gene, function finds optimal decision tree cutoff for F-beta testing
        cut_dict = {}
        for i in x:
            filename = str(i)
            y_train = dataDummy[column]
            x_train = adata[:,i].X
            X = x_train[:, None]
            clf = tree.DecisionTreeClassifier(max_leaf_nodes=2)
            clf = clf.fit(x_train, y_train) 
            threshold = clf.tree_.threshold
            cut_dict[i] = threshold[0]
        return cut_dict

    def queryGenerator(Binary_RankedList, cut_dict): # Builds dict to create queries for F-beta testing
        queryList = []
        for i in Binary_RankedList:
            str1 = i
            current_value = cut_dict.get(str1)
            queryString1 = str(str1.replace("-", "_").replace(".", "_"))+'>='+ str(current_value)
            queryList.append(queryString1)
        return queryList

    def permutor(x): # creates all combinations of queries built above
        binarylist2 = x
        combs = []
        for i in range(1, len(x)+1):
            els = [list(x) for x in itertools.combinations(binarylist2, i)]
            combs.extend(els)
        return combs

    def fbetaTest(x, column, adata, Binary_RankedList,testArray, betaValue): # uses queries to perform F-beta testing at the betaValue set in parameters
        fbeta_dict = {}
        subset_adata = adata[:,Binary_RankedList]
        Subset_dataframe = pd.DataFrame(data = subset_adata.X, index = subset_adata.obs_names, columns = subset_adata.var_names)
        Subset_dataframe.columns = Subset_dataframe.columns.str.replace("-", "_").str.replace(".", "_")

        for list in x:
            testArray['y_pred'] = 0
            betaQuery = '&'.join(list)            
            Ineq1 = Subset_dataframe.query(betaQuery)
            testList = Ineq1.index.tolist()
            testArray.loc[testList, 'y_pred'] = 1
            f1 = fbeta_score(testArray['y_true'], testArray['y_pred'], average= 'binary', beta=betaValue)
            tn, fp, fn, tp = confusion_matrix(testArray['y_true'], testArray['y_pred']).ravel()
            ### strip betaQuery and normalize        
            dictName = column+"&"+betaQuery.replace("_", "-")
            fbeta_dict[dictName] = f1, tn, fp, fn, tp
        return fbeta_dict
    
    def ReportReturn(max_grouped_df):  # Cleaning up results to return as dataframe
        for column in max_grouped_df.columns[8:14]:
            max_grouped_df[column] = max_grouped_df[column].str.replace('nan', '')
        max_grouped_df["NSForest_Markers"] = max_grouped_df[max_grouped_df.columns[8:14]].values.tolist()
        max_grouped_df = max_grouped_df[['clusterName',"f-measure",'markerCount','NSForest_Markers','True Positive','True Negative','False Positive','False Negative',1,2,3,4,5,6,"index"]]
        
        for i in max_grouped_df.index:
            cleanList = [string for string in max_grouped_df.loc[i,'NSForest_Markers'] if string != ""]
            max_grouped_df.at[i, 'NSForest_Markers'] = cleanList
        Results = max_grouped_df
        return Results
  
  
    #Parameters of interest
    
    #Random Forest parameters
    threads = -1   #Number of threads to use, -1 is the greedy option where it will take all available CPUs/RAM

    #Filtering and ranking of genes from random forest parameters
    howManyInformativeGenes2test = 15 #How many genes from the GiniRanking move on for further testing...
    
    #How many top genes from the Random Forest ranked features will be evaluated for binariness 
    InformativeGenes = 10 
    
    #Main function# 
    
    #Creates dummy columns for one vs all Random Forest modeling
    dataDummy = pd.get_dummies(adata.obs[clusterLabelcolumnHeader], columns=[clusterLabelcolumnHeader], prefix = "", prefix_sep = "")
    
    #Creates matrix of cluster median expression values
    medianValues = pd.DataFrame(columns=adata.var_names, index=adata.obs[clusterLabelcolumnHeader].cat.categories)                                                                                                 
    
    ClusterList = adata.obs[clusterLabelcolumnHeader].unique()    
    for clust in ClusterList: #adata.obs.Clusters.cat.categories:    
        subset_adata = adata[adata.obs[clusterLabelcolumnHeader].isin([clust]),:]
        Subset_dataframe = pd.DataFrame(data = subset_adata.X, index = subset_adata.obs, columns = subset_adata.var_names)
        medianValues.loc[clust] = Subset_dataframe.median()
    
    medianValues.to_csv('NSForest3_medianValues.csv')

    ##Use Mean
    #for clust in adata.obs.Clusters.cat.categories:
        #medianValues.loc[clust] = adata[adata.obs[clusterLabelcolumnHeader].isin([clust]),:].X.mean(0)
    
    clusters2Loop = len(dataDummy.columns)-1 
    print (clusters2Loop)
    #gives us the top ten features from RF
    rankedDict = {}
    f1_store_1D = {}
    Binary_score_store_DF = pd.DataFrame()
    DT_cutoffs_store = {}
    
    for column in dataDummy.columns:
        print(column)
        Binary_store_DF = pd.DataFrame()
        
        #Run Random Forest and get a ranked list
        Ranked_Features = randomForest(adata, dataDummy, column, rfTrees, threads)
        RankedList = rankInformative(Ranked_Features,column,rankedDict,howManyInformativeGenes2test)
        
        #Setup testArray for f-beta evaluation
        testArray = dataDummy[[column]]
        testArray.columns = ['y_true']
        
        #Rerank according to expression level and binary score
        Positive_RankedList_Complete = negativeOut(RankedList, column, medianValues, Median_Expression_Level)
        print(Positive_RankedList_Complete)
        
        outputlist = binaryScore(Positive_RankedList_Complete, InformativeGenes, medianValues, column, clusters2Loop, Ranked_Features, Genes_to_testing,Binary_store_DF)
        Binary_RankedList = outputlist[0]
        Binary_score_store_DF_extra = outputlist[1].assign(clusterName = column)         
        Binary_score_store_DF = Binary_score_store_DF.append(Binary_score_store_DF_extra)
        
        #Get expression cutoffs for f-beta testing
        cut_dict = DT_cutOffs(Binary_RankedList, column, dataDummy)
        DT_cutoffs_store[column] = cut_dict
        
        #Generate expression queries and run those queries using fscore() function
        queryInequalities = queryGenerator(Binary_RankedList, cut_dict)
        FullpermutationList = permutor(queryInequalities)
        print(len(FullpermutationList))
        f1_store = fbetaTest(FullpermutationList, column, adata, Binary_RankedList, testArray, betaValue)
        f1_store_1D.update(f1_store)
        
    #Report generation and cleanup for file writeouts 
    f1_store_1D_df = pd.DataFrame() #F1 store gives all results.
    f1_store_1D_df = pd.DataFrame.from_dict(f1_store_1D)
    Results_df = f1_store_1D_df.transpose()
    Results_df.columns = ["f-measure", "True Negative", "False Positive", "False Negative", "True Positive"]
    Results_df['markerCount'] = Results_df.index.str.count('&')
    Results_df.reset_index(level=Results_df.index.names, inplace=True)
    Results_df_done= Results_df['index'].apply(lambda x: pd.Series(x.split('&')))
    NSForest_Results_Table=Results_df.join(Results_df_done)
    NSForest_Results_Table_Fin = pd.DataFrame()
    NSForest_Results_Table_Fin = NSForest_Results_Table[NSForest_Results_Table.columns[0:8]]
    
    for i, col in enumerate(NSForest_Results_Table.columns[8:15]):
        splitResults = NSForest_Results_Table[col].astype(str).apply(lambda x: pd.Series(x.split('>='))) 
        firstOnly = splitResults[0]
        Ascolumn = firstOnly.to_frame()
        Ascolumn.columns = [col]
        NSForest_Results_Table_Fin = NSForest_Results_Table_Fin.join(Ascolumn)
    
    NSForest_Results_Table_Fin.rename(columns={0:'clusterName'},inplace=True) #rename columns by position
    NSForest_Results_Table_Fin.sort_values(by=['clusterName','f-measure','markerCount'],ascending= [True, False, True], inplace = True)
    print (NSForest_Results_Table_Fin)
    time.perf_counter()
    
    #Write outs
    Binary_score_store_DF.to_csv('NS-Forest_v3_Extended_Binary_Markers_Supplmental.csv')
    NSForest_Results_Table_Fin.to_csv('NS-Forest_v3_Full_Results.csv')
    
    #Subset of full results 
    max_grouped = NSForest_Results_Table_Fin[NSForest_Results_Table_Fin.groupby("clusterName")["f-measure"].transform('max') == NSForest_Results_Table_Fin['f-measure']]
    max_grouped_df = pd.DataFrame(max_grouped)
    
    ##Move binary genes to Results dataframe
    clusters2Genes = pd.DataFrame(columns = ['Gene', 'clusterName'])
    clusters2Genes["clusterName"] = Binary_score_store_DF["clusterName"]
    clusters2Genes["Gene"] = Binary_score_store_DF.index
    GroupedBinarylist = clusters2Genes.groupby('clusterName').apply(lambda x: x['Gene'].unique()) 
    BinaryFinal = pd.DataFrame(columns = ['clusterName','Binary_Genes'])
    BinaryFinal['clusterName'] = GroupedBinarylist.index
    BinaryFinal['Binary_Genes'] = GroupedBinarylist.values
        
    Results = ReportReturn(max_grouped_df)
    #Results["NSForest_Markers"] = Results["NSForest_Markers"].apply(clean_alt_list)
    
    Result = pd.merge(Results, BinaryFinal, on='clusterName')
    Result.to_csv('NSForest_v3_Final_Result.csv')
    ResultUnique = Result.drop_duplicates(subset=["clusterName"]) 
    
    time.perf_counter()
    
    return ResultUnique
