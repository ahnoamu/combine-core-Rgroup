#!/usr/bin/env python

#  Core - Rgroup convolution                                   #
#  Convolute / Combine molecules given the core and RGroup     #
#  by Arnold Amusengeri, 2022                                  #

#===========================#
# Example use:
# Method 1
# >> python3 combineCoreRgroup.py  -c  "[*:1]N1C(=CC2=CC=CC=[N+]2CC)SC(=C3N(C4=CC=CC=C4S3)C)C1=O"  -r  "[*:1]CC"  -m   m1
# Method 2
# >> python3 combineCoreRgroup.py  -c  "[*:1]N1C(=CC2=CC=CC=[N+]2CC)SC(=C3N(C4=CC=CC=C4S3)C)C1=O"  -r  "[*:1]CC"  -m   m2
#===========================#

import sys 
import getopt

import numpy as np
import pandas as pd

import rdkit.Chem as Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from rdkit.Chem import DataStructs
from rdkit.Chem.MolStandardize import rdMolStandardize
     
#============Parser=========#

# Parse in core, Rgroup, and method strings
     
def core_Rgroup_method(): 
	core = None
	Rgroup = None
	method = None

	argv = sys.argv[1:] 

	try: 
		opts, args = getopt.getopt(argv, "c:r:m:") 

	except: 
		print("Error") 

	for opt, arg in opts: 
		if opt in ['-c']:
			core = arg 
		elif opt in ['-r']: 
			Rgroup = arg 
		elif opt in ['-m']: 
			method = arg 

	print("core = ", core ," \t\t Rgroup = ", Rgroup, " \t\t method = ", method) 

	return core, Rgroup, method


#=========== Methods =========#

# Method 1

def combineCoreRgroupM1(new_core,new_Rgroup):

	try:
		### Weld Rgroup onto the core

		c = Chem.MolFromSmiles(new_core)     
		r = Chem.MolFromSmiles(new_Rgroup)

		combined = Chem.CombineMols(c,r)

		result = Chem.molzip(combined)

		### Result

		[a.SetAtomMapNum(0) for a in result.GetAtoms()] 

		result_SMI = Chem.MolToSmiles(Chem.RemoveHs(result, sanitize=False), kekuleSmiles=True)

		print('\nResult = ', result_SMI, '\n')

		return result_SMI
		
	except:
		print("Error")


# Method 2


def combineCoreRgroupM2(new_core,new_Rgroup):

	try:
		#### Merge core and Rgroup
		
		c = Chem.MolFromSmiles(new_core, sanitize=False)
		r = Chem.MolFromSmiles(new_Rgroup, sanitize=False)

		combi = Chem.CombineMols(c,r)

		### Obtain indices of dummy atoms
		
		dc,dr = np.where(np.array([i.GetAtomicNum() for i in combi.GetAtoms()]) == 0)[0]
				
		### Obtain indices of neighboring atoms 

		nc = combi.GetAtomWithIdx(int(dc)).GetNeighbors()[0].GetIdx()
		nr = combi.GetAtomWithIdx(int(dr)).GetNeighbors()[0].GetIdx()

		### Remove the dummy atoms 
		
		editable = Chem.EditableMol(combi)
		for d in sorted([dc,dr], reverse=True):
			editable.RemoveAtom(int(d))
		
		### adjust atom indices of the neighbouring atoms
		
		nc -= int(dc < nc) + int(dr < nc) # core
		nr -= int(dc < nr) + int(dr < nr) # Rgroup 
		
		### Weld core and Rgroup by adding a single bond between neighbouring atoms
		
		editable.AddBond(nc,nr,Chem.rdchem.BondType.SINGLE)
		
		result = editable.GetMol()
		
		### Result

		[a.SetAtomMapNum(0) for a in result.GetAtoms()]  
		
		result_SMI = Chem.MolToSmiles(Chem.RemoveHs(result, sanitize=False), kekuleSmiles=True)
		
		print('\nResult = ', result_SMI, '\n')
		
		return result_SMI

	except:
		print("Error")


#======= Application =======#

core, Rgroup, method = core_Rgroup_method()

new_core = str(core) 
new_Rgroup = str(Rgroup)

if str(method) == "m1":
	output = combineCoreRgroupM1(new_core,new_Rgroup)
elif str(method) == "m2":
	output = combineCoreRgroupM1(new_core,new_Rgroup)
else:
	print("\nWrong input method. Try m1 or m2 \n")


### The output can be written to file

# with open("out_file.txt", "a") as file1:
# 	file1.writelines(output+"\n")

#========= END =============#

