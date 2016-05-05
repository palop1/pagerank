import os
import sys
import csv
import time
import numpy as np
import pandas as pd

from sklearn.preprocessing import normalize
from scipy.sparse import csr_matrix, lil_matrix
from time import gmtime, strftime

def powermethod(A,x0,iter=130,d=0.85):
    """ 
    Calculates pagerank using Power Method
    
    
    power method modified to compute
     the maximal real eigenvector 
     of the matrix M built on top of the input matrix A 
     
    """
    #Power Method Algorithm
    n = A.shape[1]
    delta = (1-d)*(np.array([1]*n,dtype='float')/n)
    for i in range(iter):
        x0 = A.dot(x0).dot(d) + delta
    end = time.time()
    with open("PR_2015403258.txt", "w") as f:
        for i in range(len(x0)):
            f.write(str(i) + "," + "{:.8f}".format(x0[i]))
            f.write('\n')

def process_matrix(matrix):
    """
    Returns scaled connectivity matrix version of
    the network matrix in input parameter
    
    """
    for i in range(matrix.shape[1]):
        s =  matrix.getrow(i).sum()
        if not int(matrix[i,:].sum()):
            matrix[i] = np.ones(matrix.shape[1])
    matrix = normalize(matrix, norm='l1', copy=False, axis=1)
    matrix = matrix.transpose()
    return matrix

def get_number_of_sites(name):
    """
    Calculates how many sites are listed in the file
    
    returns integer number of sites
    """
    with open(name) as tsv:
        no_sites = 0
        for line in csv.reader(tsv, dialect="excel-tab"):
            for value in line:
                value = int(value)
                if value > no_sites:
                    no_sites = value
    return no_sites + 1

def input_file_to_matrix(name):
    """
    Changes input file into numpy matrix that maps the link relationships 
    into numpy matrix
    
    Example:
    0	1
    0	2
    1	3
    2	3
    3	0
    
    Returns:
    [[ 0.  1.  1.  0.]
     [ 0.  0.  0.  1.]
     [ 0.  0.  0.  1.]
     [ 1.  0.  0.  0.]]
    
    """
    #no_sites = get_number_of_sites(name)
    no_sites = 4
    matrix = lil_matrix((no_sites,no_sites))
    data = np.loadtxt(name, delimiter="\t")
    for line in data:
        matrix[int(line[0]),int(line[1])] += 1.0
    return matrix

if __name__=='__main__':
    """
    PageRank calculator for dynamic size link files.
    Usage:
    
    python pagerank.py <path to source file>
    """
    print "Timer starts!"
    start = time.time()
    matrix = input_file_to_matrix(str(sys.argv[1]))
    print "File read into matrix"
    end = time.time()
    print(end - start)
    
    print matrix.todense()
    scm = process_matrix(matrix)
    print scm.todense()
    print "Matrix processed"
    x0 = np.ones(scm.shape[1])
    print "x0 matrix created"
    powermethod(scm,x0)
    print "Pagerank calculated"
    
