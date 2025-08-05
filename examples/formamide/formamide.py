import psi4
import resp
from resp.utils import grid_to_mol2

mol = psi4.geometry(
    """
O   -0.03083811  -6.40367569   3.34506494
N    0.82970231 -10.25305705   5.02929344
H    1.84614178  -9.09807179   1.66343556
H    0.09834096  -9.88400850   6.49779514
H    1.48375759 -11.74633681   4.75638605
C    0.91810299  -8.54279041   3.22290639
units bohr
"""
)

mol.update_geometry()
options = {
    "VDW_SCALE_FACTORS": [1.4, 1.6, 1.8, 2.0],
    "VDW_POINT_DENSITY": 1.0,
    "RESP_A": 0.0005,
    "RESP_B": 0.1,
    "METHOD_ESP": "hf",
    "BASIS_ESP": "aug-cc-pvdz",
}

# Call for first stage fit
charges1 = resp.resp([mol], options)
print("Electrostatic Potential Charges")
print(charges1[0])
print("Restrained Electrostatic Potential Charges")
print(charges1[1])

# Change the value of the RESP parameter A
options["RESP_A"] = 0.001

# Add constraint for atoms fixed in second stage fit
# For formamide, we constrain the two amino hydrogens to be equal
options["constraint_group"] = [[4, 5]]  # Atoms 4 and 5 (1-based indexing)
options["grid"] = ["1_%s_grid.dat" % mol.name()]
options["esp"] = ["1_%s_grid_esp.dat" % mol.name()]

# Call for second stage fit
charges2 = resp.resp([mol], options)

# Get RESP charges
print("\nStage Two:\n")
print("RESP Charges")
print(charges2[1])

# Convert to MOL2 format for visualization
print("\nConverting RESP ESP to MOL2 format...")
grid_to_mol2(
    grid_file="1_default_grid.dat",
    esp_file="1_default_grid_esp.dat",
    output_mol2="resp_esp.mol2",
    normalize=False,
)

# Also create a PDB file of the molecule
from resp.utils import molecule_to_pdb

molecule_to_pdb(mol, filename="formamide.pdb", res_name="FMD")
