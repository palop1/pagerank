import sys
import csv
import time
import numpy as np

from scipy.sparse import csc_matrix, lil_matrix, eye
from time import gmtime, strftime


def linearEquations(A,m=0.15):
    """ solving linear equations 
        of the system (I-(1-m)*A)*x = m*s """
    n = A.shape[1]
    C = eye(n,n)-A.dot(1-m)
    b = m*(np.array([1]*n,dtype='float32')/n)
    ranks = np.linalg.solve(C.todense(),b)
    
    with open("PR_2015403258.txt", "w") as f:
        for i in range(len(ranks)):
            f.write(str(i) + "," + "{:.8f}".format(ranks[i]))
            f.write('\n')

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
    
    with open("PR_2015403258.txt", "w") as f:
        for i in range(len(x0)):
            f.write(str(i) + "," + "{:.8f}".format(x0[i]))
            f.write('\n')
    
    
def process_matrix(matrix):
    """
    Returns scaled connectivity matrix version of
    the network matrix in input parameter
    
    """
    
    idxs = [idx for idx in range(matrix.shape[1]) if idx not in list(set(matrix.nonzero()[0]))]
    row_norm = np.bincount(matrix.indices, weights=matrix.data * matrix.data)
    matrix.data /= np.take(row_norm, matrix.indices)
    for idx in idxs:
        matrix[idx] = np.ones(matrix.shape[1])/matrix.shape[1]
    
    """
    for i in range(matrix.shape[1]):
        
        s =  matrix.getrow(i).sum()
        if not int(matrix[i,:].sum()):
            matrix[i] = np.ones(matrix.shape[1])
            s = matrix.shape[1]
        matrix[i] /= s
    """
    matrix = matrix.transpose()
    return matrix

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
    no_sites = 100
    
    with open(name) as tsv:
        matrix = lil_matrix((no_sites,no_sites),dtype=float) 
        for line in csv.reader(tsv, dialect="excel-tab"):
            #line = [float(value) for value in line]
            matrix[float(line[0]),float(line[1])] += 1.0
    """
    matrix = lil_matrix((no_sites,no_sites)) 
    data = np.loadtxt(name, delimiter="\t")
    for line in data:
        matrix[int(line[0]),int(line[1])] += 1.0
    """
    matrix = matrix.tocsc()
    return matrix

if __name__=='__main__':
    """
    PageRank calculator for dynamic size link files.
    Usage:
    
    python pagerank.py <path to source file>
    """
    print strftime("%Y-%m-%d %H:%M:%S", gmtime())
    
    print "File I/O starts "
    start = time.time()
    matrix = input_file_to_matrix(str(sys.argv[1]))
    end = time.time()
    print "File I/O done in:", str(end-start)
    
    print "Matrix processing starts"
    start = time.time()
    srm = process_matrix(matrix)
    x0 = np.ones(srm.shape[1])
    end = time.time()
    print "Matrix processed in: ", str(end-start)
    
    print "Powermethod calculation starts!"
    start = time.time()
    powermethod(srm,x0)
    end = time.time()
    print "Powermethod with 130 iterations done: ", str(end-start)
    
    print "Linear Equations calculation starts!"
    start = time.time()
    linearEquations(srm)
    end = time.time()
    print "Linear Equations done: ", str(end-start)
    
    print strftime("%Y-%m-%d %H:%M:%S", gmtime())
