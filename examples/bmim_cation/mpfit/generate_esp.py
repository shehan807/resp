import psi4
import resp
from resp.utils import compare_grid_esp, create_difference_esp, grid_to_mol2

mol = psi4.geometry(
""" 
 N   36.025      77.846        41.841  
 C   37.174      78.034        40.967  
 C   36.071      77.613        43.178  
 C   34.784      77.564        43.624  
 N   33.914      77.844        42.526  
 C   34.718      77.752        41.432  
 C   32.418      77.791        42.527  
 C   36.942      77.575        39.507  
 C   38.169      78.008        38.735  
 C   39.418      77.195        39.150  
 H   40.207      77.351        38.414  
 H   39.247      76.122        39.238  
 H   39.731      77.519        40.142  
 H   38.064      77.902        37.655  
 H   38.346      79.064        38.942  
 H   36.792      76.496        39.507  
 H   36.077      78.018        39.014  
 H   37.977      77.465        41.436  
 H   37.417      79.096        40.998  
 H   36.948      77.420        43.777  
 H   34.331      77.379        44.587  
 H   34.393      77.825        40.405  
 H   32.064      78.763        42.181  
 H   31.997      77.738        43.531  
 H   31.952      77.043        41.886  
units angstrom
""")

mol.update_geometry()
#mol.set_molecular_charge(1)  # Set charge to +1 for cation
#mol.set_multiplicity(1)      # Set multiplicity (2S+1, where S is total spin)
psi4.set_options({'reference': 'uks'})  # Unrestricted Kohn-Sham for DFT

options = {'VDW_SCALE_FACTORS'  : [1.4, 1.6, 1.8, 2.0],
           'VDW_POINT_DENSITY'  : 1.0,
           'RESP_A'             : 0.0005,
           'RESP_B'             : 0.1,
	   'METHOD_ESP'         : 'pbe0', 
	   'BASIS_ESP'          : 'def2-SVP', 
	  }

# MPFIT-generated charges for BMIM+
# Based on the charge table:
# N: -0.1915, C1: 0.3158, C2: -0.0192, C3: 0.2316, C4: 0.2341, C5: 0.0616
# C6: -0.0251 (4 carbons in butyl chain)
# H1: 0.0919, H2: 0.1644, H3: 0.0328, H4: 0.0157, H5: -0.0181, H6: 0.0213

charges = [-0.1915,    # N (atom 1)
           0.3158,     # C1 (atom 2) 
           -0.0192,    # C2 (atom 3)
           0.2316,     # C3 (atom 4)
           0.2341,     # C4 (atom 5)
           0.0616,     # C5 (atom 6)
           -0.0251,    # C6 (atom 7)
           -0.0251,    # C6 (atom 8)
           -0.0251,    # C6 (atom 9)
           -0.0251,    # C6 (atom 10)
           -0.1915,    # N (atom 11)
           0.0616,     # C5 (atom 12)
           -0.0251,    # C6 (atom 13)
           0.0919,     # H1 (atom 14)
           0.0919,     # H1 (atom 15)
           0.0919,     # H1 (atom 16)
           0.0213,     # H6 (atom 17)
           0.0213,     # H6 (atom 18)
           0.0213,     # H6 (atom 19)
           0.0213,     # H6 (atom 20)
           0.0213,     # H6 (atom 21)
           0.0213,     # H6 (atom 22)
           0.0213,     # H6 (atom 23)
           0.1644,     # H2 (atom 24)
           0.1644,     # H2 (atom 25)
           0.0328,     # H3 (atom 26)
           0.0157,     # H4 (atom 27)
           -0.0181,    # H5 (atom 28)
           -0.0181,    # H5 (atom 29)
           -0.0181,    # H5 (atom 30)
           -0.0181]    # H5 (atom 31)

# (as a test) ESP-generated charges
#charges = [-0.51241923, -0.01963482, 0.11270785, 0.19172582, 0.16623251, 0.29778778, -0.09928506, -0.35601807, 0.21890322]

# Generate ESP grid files from charges
resp.charges_to_esp(mol, charges, options)

print(f"Generated grid files:")
print(f"  1_{mol.name()}_grid.dat")
print(f"  1_{mol.name()}_grid_esp.dat")

# Compare with RESP-generated ESP
print("\nComparing MPFIT-generated ESP with RESP-generated ESP...")
resp_esp = "../resp/1_default_grid_esp.dat"
mpfit_esp = "1_default_grid_esp.dat"

try:
    metrics = compare_grid_esp(resp_esp, mpfit_esp, verbose=True)
    
    # Create difference ESP file
    print("\nCreating difference ESP file...")
    create_difference_esp(resp_esp, mpfit_esp, output_file='difference_grid_esp.dat')
    
    # Convert to MOL2 format for visualization
    print("\nConverting to MOL2 format...")
    grid_to_mol2(grid_file='1_default_grid.dat', 
                 esp_file='1_default_grid_esp.dat',
                 output_mol2='mpfit_esp.mol2',
                 normalize=False)
    
    # Also create MOL2 for the difference
    print("\nCreating MOL2 for difference ESP...")
    grid_to_mol2(grid_file='1_default_grid.dat',
                 esp_file='difference_grid_esp.dat', 
                 output_mol2='difference_esp.mol2',
                 normalize=False)
    
except Exception as e:
    print(f"Error during comparison: {e}")
