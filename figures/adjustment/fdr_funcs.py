# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import numpy as np

def fdr(pval, q=0.05, cV=1):
    
    if len(pval.shape) > 1:
        raise Exception('pval must be a vector')
    
    # Sort p-values
    pval  = np.sort(pval)
    oidxR = np.argsort(np.argsort(pval))

    # Number of observations
    V = len(pval);
    
    # Order (indices), in the same size as the pvalues
    idx = np.arange(V)+1
    
    # Line to be used as cutoff
    thrline = idx*q/V/cV;
    
    # Find the largest pval, still under the line
    idxthr = pval <= thrline
    if np.any(idxthr):
        
        # This is the FDR threshold
        thr = np.max(pval[idxthr])
    
        # Deal with the case when all the points under the line
        # are equal to zero, and other points are above the line
        if thr == 0:
            thr = np.max(thrline[idxthr])
    else:
        # Case when it does not cross
        thr = 0
    
    # p-corrected
    pcor = pval*V*cV/idx
    
    # p-adjusted
    padj = np.flip(np.minimum.accumulate(np.flip(pcor, axis=0)), axis=0)
    padj = padj[oidxR]
    
    # Outputs
    return thr, padj

def bky7(pval, q=0.05, cV=1):
    
    if len(pval.shape) > 1:
        raise Exception('pval must be a vector')
    
    # Sort p-values
    pval  = np.sort(pval)
    oidxR = np.argsort(np.argsort(pval))
    
    # Number of observations
    V = len(pval);
    
    # Order (indices), in the same size as the pvalues
    idx = np.arange(V)+1
    
    idxthr = np.full(idx.shape, False)
    for v in range(V):
        # Line to be used as cutoff
        thrline = idx[v:]*q/(V+1-(v+1)*(1-q))/cV
        
        # P-vals that survive the cutoff
        if np.any(pval[v:]<=thrline):
            idxthr[v] = True
        else:
            break
    
    if np.any(idxthr):
        
        # This is the FDR threshold
        thr = np.max(pval[idxthr])
    
        # Deal with the case when all the points under the line
        # are equal to zero, and other points are above the line
        if thr == 0:
            v = idx[idxthr][-1]
            thr = v*q/(V+1-v*(1-q))/cV # note here it's v, not v+1 as above, because here the indices start at 1, whereas above they start at 0
    else:
        # Case when it does not cross
        thr = 0
    
    # p-corrected
    pcor = np.ones(pval.shape);
    for v in range(V):
        # For every p-value, this is the minimum q that will eventually
        # satisfy Definition #7 of Benjamini, Krieger and Yekutieli (2006).
        pcor[v] = np.min(pval[v:]*(V+1-(v+1))*cV/(idx[v:]-(v+1)*pval[v:]));

    # p-adjusted
    #padj = np.flip(np.minimum.accumulate(np.flip(pcor, axis=0)), axis=0)
    padj = np.maximum.accumulate(pcor)
    padj = padj[oidxR]
    
    # Outputs
    return thr, padj

punc = np.array([0.005, 0.010, 0.014, 0.025, 0.042, 0.066,  0.1,   0.12,  0.17,   0.28,  0.36,  0.524, 0.61, 0.68,  0.78,  0.9, 0.96])

fdrthr, fdradj = fdr (punc, q=.2)
bkythr, bkyadj = bky7(punc, q=.2)