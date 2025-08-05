import psi4
import resp
from resp.utils import grid_to_mol2

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
mol.set_molecular_charge(1)  # Set charge to +1 for cation
mol.set_multiplicity(1)      # Set multiplicity (2S+1, where S is total spin)
psi4.set_options({'reference': 'uks'})  # Unrestricted Kohn-Sham for DFT

options = {'VDW_SCALE_FACTORS'  : [1.4, 1.6, 1.8, 2.0],
           'VDW_POINT_DENSITY'  : 1.0,
           'RESP_A'             : 0.0005,
           'RESP_B'             : 0.1,
	   'METHOD_ESP'         : 'pbe0', 
	   'BASIS_ESP'          : 'def2-SVP', 
	  }

# Call for first stage fit
charges1 = resp.resp([mol], options)
print('Electrostatic Potential Charges')
print(charges1[0])
print('Restrained Electrostatic Potential Charges')
print(charges1[1])

# Change the value of the RESP parameter A
options['RESP_A'] = 0.001

# Add constraint for atoms fixed in second stage fit
constraint_charge = []
for i in range(4, 8):
    constraint_charge.append([charges1[1][i], [i+1]])
options['constraint_charge'] = constraint_charge
options['constraint_group'] = [[2, 3, 4]]
options['grid'] = ['1_%s_grid.dat' %mol.name()]
options['esp'] = ['1_%s_grid_esp.dat' %mol.name()]

# Call for second stage fit
charges2 = resp.resp([mol], options)

# Get RESP charges
print("\nStage Two:\n")
print('RESP Charges')
print(charges2[1])

# Convert to MOL2 format for visualization
print("\nConverting RESP ESP to MOL2 format...")
grid_to_mol2(grid_file='1_default_grid.dat',
             esp_file='1_default_grid_esp.dat',
             output_mol2='resp_esp.mol2',
             normalize=False)

# Convert ESP to kJ/mol and create MOL2
print("\nConverting ESP to kJ/mol...")
import numpy as np
esp_data = np.loadtxt('1_default_grid_esp.dat')
hartree_to_kjmol = 2625.5
esp_kjmol = esp_data * hartree_to_kjmol
np.savetxt('1_default_grid_esp_kjmol.dat', esp_kjmol, fmt='%15.10f')
print(f'ESP range in kJ/mol: [{esp_kjmol.min():.2f}, {esp_kjmol.max():.2f}]')

# Create MOL2 for kJ/mol data
print("\nConverting kJ/mol ESP to MOL2 format...")
grid_to_mol2(grid_file='1_default_grid.dat',
             esp_file='1_default_grid_esp_kjmol.dat',
             output_mol2='resp_esp_kjmol.mol2',
             normalize=False)

# Also create a PDB file of the molecule
from resp.utils import molecule_to_pdb
molecule_to_pdb(mol, filename='pyrazole.pdb', res_name='PYZ')
