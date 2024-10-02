
import pandas as pd
import scanpy as sc
from sklearn.ensemble import RandomForestClassifier

def myRandomForest(adata, df_dummies, cl, n_trees, n_jobs, n_top_genes, binary_dummies):
    """\
    Running sklearn.ensemble's RandomForestClassifier on the binary dummy variables. 

    Parameters
    ----------
        adata: AnnData
            Annotated data matrix.
        df_dummies: pd.DataFrame
            Dummy dataframe for one vs all model. 
        cl: str
            Specified `cluster_header` value. 
        n_trees: int (default: 1000)
            `n_estimators` parameter in sklearn.ensemble's RandomForestClassifier. 
        n_jobs: int (default: -1)
            `n_jobs` parameter in sklearn.ensemble's RandomForestClassifier. 
        beta: float (default: 0.5)
            `beta` parameter in sklearn.metrics's fbeta_score. 
        n_top_genes: int (default: 15)
            Taking the top `n_top_genes` genes ranked by sklearn.ensemble's RandomForestClassifier as input for sklearn.tree's DecisionTreeClassifier. 
        binary_dummies: pd.DataFrame
            Dataframe of binary scores filtered by `gene_selection`. 

    Returns
    -------
    top_rf_genes: list 
        The top `n_top_genes` genes ranked by Gini Impurity. 
    """
    x_train = adata.to_df()
    y_train = df_dummies[cl]
    
    ## pre-select genes based on gene_selection criterium
    ind_genes_selected = binary_dummies.loc[cl] == 1
    x_train = x_train.loc[:,ind_genes_selected]   # subset x_train
    print(f"\tPre-selected {x_train.shape[1]} genes to feed into Random Forest.")

    rf_clf = RandomForestClassifier(n_estimators=n_trees, n_jobs=n_jobs, random_state=123456) #<===== criterion="gini", by default
    rf_clf.fit(x_train, y_train)
    ## get feature importance and rank/subset top genes
    top_rf_genes = pd.Series(rf_clf.feature_importances_, index=x_train.columns).sort_values(ascending=False)[:n_top_genes]
    
    return top_rf_genes  
