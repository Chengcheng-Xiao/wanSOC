from __future__ import print_function
import numpy as np
from numpy.linalg import inv
from wanSOC.basis import *

# L matrices
def MatLp(basis):
    '''
    <lm|L+|lm>
    order l=1: -2 -1 0 1 2 -> 0 1 2 3 4
    order l=2: -1 0 1      -> 0 1 2
    '''
    MatL = np.zeros([len(basis),len(basis)],dtype=np.float64)
    # assert len(basis)==basis[0][0]
    for i in range(len(basis)):
        l, m = basis[i][0], basis[i][1]
        if m+1 > l:
            continue
        else:
            MatL[m+l+1,m+l] = np.sqrt((l-m)*(l+m+1))
    return MatL

def MatLm(basis):
    '''
    <lm|L-|lm>
    order l=1: -2 -1 0 1 2 -> 0 1 2 3 4
    order l=2: -1 0 1      -> 0 1 2
    '''
    MatL = np.zeros([len(basis),len(basis)],dtype=np.float64)
    # assert len(basis)==basis[0][0]
    for i in range(len(basis)):
        l, m = basis[i][0], basis[i][1]
        if m-1 < -l:
            continue
        else:
            MatL[m+l-1,m+l] = np.sqrt((l+m)*(l-m+1))
    return MatL

def MatLz(basis):
    '''
    <lm|Lz|lm>
    order l=1: -2 -1 0 1 2 -> 0 1 2 3 4
    order l=2: -1 0 1      -> 0 1 2
    '''
    MatL = np.zeros([len(basis),len(basis)],dtype=np.float64)
    # assert len(basis)==basis[0][0]
    for i in range(len(basis)):
        l, m = basis[i][0], basis[i][1]
        MatL[m+l,m+l] = m
    return MatL

# transform matrix
def trans_L_mat(orb):
    '''
    calculate the real spherical harmonics transformation matrix
    '''
    if orb == 'p':
        # old basis: -1 0 1
        # new basis: pz,px,py
        trans=np.mat([[0,1,0],
                      [np.sqrt(1/2.),0,-np.sqrt(1/2.)],
                      [np.complex(0,np.sqrt(1/2.)),0,np.complex(0,np.sqrt(1/2.))]])
    elif orb == 'd':
        # old basis: -2 -1 0 1 2
        # new basis: dz2,dxz,dyz,dx2-dz2,dxy
        trans=np.mat([[0,0,1,0,0],
                      [0,np.sqrt(1/2.),0,-np.sqrt(1/2.),0],
                      [0,np.complex(0,np.sqrt(1/2.)),0,np.complex(0,np.sqrt(1/2.)),0],
                      [np.sqrt(1/2.),0,0,0,np.sqrt(1/2.)],
                      [np.complex(0,np.sqrt(1/2.)),0,0,0,np.complex(0,-np.sqrt(1/2.))]])
    return trans.T


def trans_S_mat(q_axis):
    '''
    calculate the spin quantization axis transformation matrix
    '''
    # normalize dir
    q_axis /= np.linalg.norm(q_axis)
    # up vector
    sp = 1./np.sqrt(2*1.*(q_axis[2]+1.)) * \
         np.array([q_axis[2]+1.,np.complex(q_axis[0], q_axis[1])])
    # down vector
    sm = 1./np.sqrt(2*1.*(q_axis[2]+1.)) * \
         np.array([np.complex(-q_axis[0],q_axis[1]), q_axis[2]+1.])

    trans=np.mat([sp,sm])
    return trans.T

# hamiltonian
def gen_Hsoc(orb, q_axis):
    '''
    generate H_soc.
    '''
    # create lm basis
    basis = creat_basis_lm(orb)
    # create L matrices
    Lz = MatLz(basis)
    Lp = MatLp(basis)
    Lm = MatLm(basis)
    # generate transform matrix
    trans_L = trans_L_mat(orb)
    # transform them into real spherical harmonics
    Lz = np.dot(np.dot(trans_L.H, Lz),trans_L)
    Lp = np.dot(np.dot(trans_L.H, Lp),trans_L)
    Lm = np.dot(np.dot(trans_L.H, Lm),trans_L)


    # set up pauli matrices
    Sp = np.array([[0,1],[0,0]])
    Sm = np.array([[0,0],[1,0]])
    Sz = np.array([[1,0],[0,-1]])
    # set up quantization axis
    trans_S = trans_S_mat(q_axis)
    # transform them into real spherical harmonics
    Sz = np.dot(np.dot(trans_S.H, Sz),trans_S)
    Sp = np.dot(np.dot(trans_S.H, Sp),trans_S)
    Sm = np.dot(np.dot(trans_S.H, Sm),trans_S)

    Hsoc = np.kron(Sp, Lm) + np.kron(Sm, Lp) + np.kron(Sz, Lz)
    return Hsoc

def get_full_Hsoc(orb_chg, hop_spinor):
    '''
    generate full H_soc.
    '''
    # # set up two atoms with d orbtials, first one's star at orbital number 3, second one at 8.
    # orb_chg = [['d',[1,0,0],0.1,3],['d',[1,0,0],0.1,8]]

    full_Hsoc = np.zeros([hop_spinor.shape[1],hop_spinor.shape[2]],dtype=np.complex64)
    for i in range(len(orb_chg)):
        Hsoc = gen_Hsoc(orb_chg[i][0], orb_chg[i][1])*orb_chg[i][2]
        # split Hsoc into four parts
        slice_index = int(Hsoc.shape[0]/2.)
        Hsoc_uu = Hsoc[0:slice_index,0:slice_index]
        Hsoc_ud = Hsoc[0:slice_index,slice_index:]
        Hsoc_du = Hsoc[slice_index:,0:slice_index]
        Hsoc_dd = Hsoc[slice_index:,slice_index:]
        # add them to full_Hsoc
        # locate index
        uu_s = slice(int(orb_chg[i][3]-1), int(orb_chg[i][3]-1+Hsoc.shape[0]/2.),1)
        uu_e = slice(int(orb_chg[i][3]-1), int(orb_chg[i][3]-1+Hsoc.shape[0]/2.),1)
        ud_s = slice(int(orb_chg[i][3]-1), int(orb_chg[i][3]-1+Hsoc.shape[0]/2.),1)
        ud_e = slice(int(orb_chg[i][3]-1+hop_spinor.shape[1]/2.),int(orb_chg[i][3]-1+Hsoc.shape[0]/2.+hop_spinor.shape[1]/2.),1)
        du_s = slice(int(orb_chg[i][3]-1+hop_spinor.shape[1]/2.),int(orb_chg[i][3]-1+Hsoc.shape[0]/2.+hop_spinor.shape[1]/2.),1)
        du_e = slice(int(orb_chg[i][3]-1),int(orb_chg[i][3]-1+Hsoc.shape[0]/2.),1)
        dd_s = slice(int(orb_chg[i][3]-1+hop_spinor.shape[1]/2.),int(orb_chg[i][3]-1+Hsoc.shape[0]/2.+hop_spinor.shape[1]/2.),1)
        dd_e = slice(int(orb_chg[i][3]-1+hop_spinor.shape[1]/2.),int(orb_chg[i][3]-1+Hsoc.shape[0]/2.+hop_spinor.shape[1]/2.),1)
        # print(uu_s)
        full_Hsoc[uu_s, uu_e] = Hsoc_uu
        full_Hsoc[ud_s, ud_e] = Hsoc_ud
        full_Hsoc[du_s, du_e] = Hsoc_du
        full_Hsoc[dd_s, dd_e] = Hsoc_dd

    return full_Hsoc
