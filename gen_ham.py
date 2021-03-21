from __future__ import print_function
import numpy as np
from wanSOC.basis import *
from wanSOC.io import *
from wanSOC.helper import *
from wanSOC.hamiltonian import *

#%% generatr full hamiltonian form non-collinear hamiltonians
hop_spinor_dn, Rlatt_dn, indR0_dn, deg_dn = read_hr('wannier90.dn_hr.dat', spin='dn')
hop_spinor_up, Rlatt_up, indR0_up, deg_up = read_hr('wannier90.up_hr.dat', spin='up')

# Sanity check
assert Rlatt_dn.all()==Rlatt_up.all(), "Real-lattice hopping vecotr mismatch."
assert indR0_dn==indR0_up, "home cell index mismatch."
assert deg_dn.all()==deg_dn.all(), "degeneracy mismatch"

# add hop_spinor to get full spinor hamiltonian
hop_spinor = hop_spinor_dn + hop_spinor_up


#%% we have two atoms each with one set of d orbitals
# [orbital_type,quantization_axis,strength/2,orbital_index(start from 1)]
orb_chg=[['d',[0,0,1],0.1,3],['d',[0,0,1],0.1,8]]
full_Hsoc=get_full_Hsoc(orb_chg, hop_spinor)
# print(full_Hsoc)
# add H_soc to full hamiltonian
hop_spinor[indR0_dn] += full_Hsoc

#%%
write_hr(deg_dn, Rlatt_dn, hop_spinor, Filename='wannier90_hr.dat')
