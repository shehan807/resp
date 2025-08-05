import psi4
import resp
from resp.utils import grid_to_mol2

mol = psi4.geometry(
""" 
 B   32.538      79.590        46.732
 F   33.202      78.374        46.554
 F   33.459      80.605        47.178
 F   31.821      79.897        45.526
 F   31.715      79.417        47.777
units angstrom
""")

mol.update_geometry()
mol.set_molecular_charge(-1)  # Set charge to +1 for cation
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
# For BF4, group all fluorine atoms to have equal charges (don't use individual constraints)
options['constraint_group'] = [[2, 3, 4, 5]]  # Group all 4 fluorine atoms (atoms 2-5)
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
