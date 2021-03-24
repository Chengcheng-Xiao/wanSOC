from __future__ import print_function
import numpy as np

def creat_basis_lm(orb):
    '''
    Creat |lm,s> stated for orbital.
    '''
    if orb=='s':
        print('we do not consider the soc in s orbital')
        exit()
    elif orb=='p':
        basis = []
        for m in [-1,0,1]:
            basis.append([1,m])
    elif orb=='d':
        basis = []
        for m in [-2,-1,0,1,2]:
            basis.append([2,m])
    elif orb=='f':
        basis = []
        for m in [-3,-2,-1,0,1,2,3]:
            basis.append([3,m])
    else:
        print('wrong l quantum number.')
        exit()
    return basis
