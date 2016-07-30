
def glm(loom, group_attr, method="sample"):
    """
    Calculate bayesian generalized linear regression with 'group_attr' as the grouping of cells.

    Args:
        group_attr (string):	Name of the col_attr that will be used to group cells.

    Returns:
        Nothing, but (re)creates the following datasets in the file:
            
            'regression_means' 		Matrix of floats, shaped as (genes, labels).
            'regression_stdevs' 	Matrix of floats, shaped as (genes, labels).
            'regression_labels'		Labels in lexical order.
            'binarized'				Matrix of 64-bit integers, shaped as (genes, n).
                                    n is sufficiently large to hold one bit per label
    """
    if not self.col_attrs.__contains__("_TotalRNA") or not self.row_attrs.__contains__("_Excluded"):
        raise RuntimeError("Regression requires feature selection")
        
    if self.file.__contains__("regression_means"):
        del self.file["regression_means"]
    if self.file.__contains__("regression_stdevs"):
        del self.file["regression_stdevs"]
    if self.file.__contains__("regression_labels"):
        del self.file["regression_labels"]

    labels = list(set(self.col_attrs[group_attr]))
    labels.sort()
    self.file["regression_labels"] = labels

    # Create the design matrix
    x = []
    for lbl in labels:
        indicators = (self.col_attrs[group_attr] == lbl).astype('int')
        x.append(indicators)
    x = np.array(x)
    ix = 0

    # Calculate the relative transcriptome sizes
    kappa = self.col_attrs["_TotalRNA"]/self.col_attrs["_TotalRNA"].mean()
    
    # Get the data for relevant genes
    selection = (1-self.row_attrs["_Excluded"]).astype('bool')
    y = self[selection,:].astype('int')

    c = y.shape[1]
    g = y.shape[0]
    k = x.shape[0]
    data = {
        "C": c,
        "G": g,
        "K": k,
        "kappa": kappa,
        "x": x,
        "y": y
    }
    init = {
        "beta": np.zeros((g,k)) + 0.1,
        "r": np.zeros((g,)) + 0.1,
        "basal": np.zeros((g,)) + 0.1
    }
    stan = CmdStan("/Users/Sten/Dropbox/Code/cmdstan/")
    result = stan.fit("bregression", data, init, method=method, debug=True)
    means = np.array([result["beta.%d.%d" % (row,col)].mean() for row in xrange(1,g+1) for col in xrange(1,k+1)]).reshape((g,k))
    stdevs = np.array([result["beta.%d.%d" % (row,col)].std() for row in xrange(1,g+1) for col in xrange(1,k+1)]).reshape((g,k))

    self.file["regression_means"] = means
    self.file["regression_stdevs"] = stdevs
    
    return result
    # TODO: binarize