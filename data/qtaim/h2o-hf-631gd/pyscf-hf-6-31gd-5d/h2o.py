"""h2o HF/6-31G(d), spherical d (5D): Avogadro QTAIM reference data.

Run with: python h2o.py > h2o.out   (written for PySCF 2.14)
Same fixed geometry as geometry.xyz; PySCF does not reorient without symmetry.
"""
import pyscf
from pyscf import gto, scf
from pyscf.tools import molden

mol = gto.M(
    atom="""
O          0.99902       -0.06151        0.07481
H          1.96108       -0.09190        0.05648
H          0.71180       -0.84832       -0.39981
""",
    unit="Angstrom",
    basis="6-31g*",
    cart=False,
    symmetry=False,
    charge=0,
    spin=0,
)
mf = scf.RHF(mol)
mf.conv_tol = 1e-12
mf.kernel()
print("PySCF", pyscf.__version__)
print("Total energy", repr(mf.e_tot))
molden.from_scf(mf, "h2o.molden")
