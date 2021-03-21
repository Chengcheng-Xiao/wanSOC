# wanSOC

Add on-site SOC to collinear Wannier Hamiltonians (`hr.dat` files) from [Wannier90](http://www.wannier.org/)

Support:
- spin-polarized/unpolarized systems.
- change on-site SOC strength.
- change spin quantization axis.

limitations:
- only support full-range m quantum numbers (e.g. l=1, m=-1, 0, 1)
- only support `p`, `d` orbitals (as of now).

## tutorial
Take a look at `gen_ham.py`.
