
import pandas as pd
from sklearn.tree import DecisionTreeClassifier
from sklearn.metrics import fbeta_score
from sklearn.metrics import precision_score
from sklearn.metrics import recall_score
from sklearn.metrics import confusion_matrix
import itertools

## construct decision tree for each gene and evaluate the fbeta score in all combinations ==> outputs markers with max fbeta, and all scores
def myDecisionTreeEvaluation(adata, df_dummies, cl, genes_eval, beta = 0.5, combinations = True):
    """\
    Calculating performance metrics for `genes_eval`. 

    Parameters
    ----------
        adata: AnnData
            Annotated data matrix.
        df_dummies: pd.DataFrame
            Dummy dataframe for one vs all model. 
        cl: str
            Specified `cluster_header` value. 
        genes_eval: list
            List of genes to find best combination in sklearn.tree's DecisionTreeClassifier. 
        beta: float (default: 0.5)
            `beta` parameter in sklearn.metrics's fbeta_score. 
        combinations: bool (default: True)
            Whether to find the combination of `genes_eval` with the highest fbeta_score. 
    
    Returns
    -------
    markers: list
        list of markers with returned scores. 
    scores: list
        fbeta, ppv, recall, tn, fp, fn, tp of `markers`. 
    """
    # Training decision tree based on single gene split
    dict_pred = {}
    for i in genes_eval:
        x_train = adata[:,i].to_df() # OG anndata
        y_train = df_dummies[cl]
        tree_clf = DecisionTreeClassifier(max_leaf_nodes=2)
        tree_clf = tree_clf.fit(x_train, y_train) 
        dict_pred[i] = tree_clf.apply(x_train)-1
    df_pred = pd.DataFrame(dict_pred) #cell-by-gene
    
    combs = []
    # Getting every possible subset of genes_eval
    if combinations: 
        for L in range(1, len(genes_eval)+1):
            els = [list(x) for x in itertools.combinations(genes_eval, L)]
            combs.extend(els)
    # Using all available marker genes
    else: 
        combs = [genes_eval]
    # print("COMBINATIONS:", combs)
    
    # Checking if every possible subset of genes_eval leads to correct prediction
    dict_scores = {} 
    for ii in combs:
        y_true = df_dummies[cl]
        # if at least 1 gene is incorrect, 0
        y_pred = df_pred[ii].product(axis=1)
        fbeta = fbeta_score(y_true, y_pred, average='binary', beta=beta)
        ppv = precision_score(y_true, y_pred, average='binary', zero_division=0)
        recall = recall_score(y_true, y_pred, average='binary', zero_division=0)
        tn, fp, fn, tp = confusion_matrix(y_true, y_pred).ravel()
        dict_scores['&'.join(ii)] = fbeta, ppv, recall, tn, fp, fn, tp
    df_scores = pd.DataFrame(dict_scores) #cell-by-genecombo
        
    ## find which combination has the max fbeta
    idx_max = df_scores.idxmax(axis=1)[0] # axis = 1 means axis = columns, [0] is fbeta (first column)
    markers = idx_max.split('&')
    scores = df_scores[idx_max]
    
    return markers, scores
