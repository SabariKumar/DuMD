import openmm as mm
from openmm import Vec3
import openmm.unit as unit
from openmm.app import PDBFile
from pdbfixer.pdbfixer import PDBFixer, proteinResidues
from multiprocessing import Process, Pipe
import os
from pathlib import Path
import sys
from time import perf_counter
from tqdm import tqdm

def prep_pdb(pdb_file, stage: str, param_dict: dict) -> None: #pH = 7.0, add_H = True, water= 'TIP3P-FB', unit_cell_offset = 5, ion_st = 0 )
    start = perf_counter()
    name = pdb_file.split('.')[0] 
    #Read in the pdb, remove unwanted chains. For now, we'll only keep chains with AA residues or water.
    fixer = PDBFixer(pdb_file)
    chains = []
    hasHydrogens = False
    # Grab chains with just protein residues
    for chain in fixer.topology.chains():
        residues = list(r.name for r in chain.residues())
        if any(r in proteinResidues for r in residues):
            chains.append(chain)
        elif any(r == 'HOH' for r in chain.residues()):
            chains.append(chain)
        else:
            content = ', '.join(set(residues))
            hasHeterogen = True
            
    chains = set(chains)
    delete_indices = [i for i, chain in enumerate(fixer.topology.chains()) if chain not in chains]
    fixer.removeChains(delete_indices)
    
    #Add in missing residues/heavy atoms, replace nonstandard residues
    fixer.findMissingResidues()
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    if(param_dict['add_H']):
        fixer.addMissingHydrogens(param_dict['pH'])
    
    #Create cubic solvent box that's 5A bigger than the biggest dimension  
    cell_dims = [max((pos[i] for pos in fixer.positions))-min((pos[i] for pos in fixer.positions)) for i in range(3)]
    offset = [unit.Quantity(param_dict['unit_cell_offset'], unit.nanometer), unit.Quantity(param_dict['unit_cell_offset'], unit.nanometer),
              unit.Quantity(param_dict['unit_cell_offset'], unit.nanometer)]
    dims = [x+y for x, y in zip(cell_dims, offset)]
    fixer.addSolvent([x.value_in_unit(unit.nanometer) for x in dims], positiveIon='Na+', negativeIon='Cl-', ionicStrength=param_dict['ion_st'] * unit.molar)
    PDBFile.writeFile(fixer.topology, fixer.positions, open(name + '_processed.pdb', 'w'))
    end = perf_counter()
    print(f'Processed pdb file {name} in {end - start} seconds')