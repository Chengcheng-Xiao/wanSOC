from __future__ import print_function
import numpy as np
from numpy.linalg import inv

def creat_basis_lm(orb):
    '''
    Creat |lm,s> stated for orbital.
    '''
    if orb=='p':
        basis = []
        for m in [-1,0,1]:
            basis.append([1,m])
    if orb=='d':
        basis = []
        for m in [-2,-1,0,1,2]:
            basis.append([2,m])
    if orb=='s':
        print('we do not consider the soc in s orbital')
    if orb=='f':
        print('for now, soc for f orbital is not added.')
        exit()
    return basis
