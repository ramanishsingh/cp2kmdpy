#!/usr/bin/env python
# coding: utf-8

# In[2]:


import mbuild as mb
from cssi_cp2k.classes import SIM as sim
# read this: https://www.cp2k.org/_media/events:2015_cecam_tutorial:watkins_optimization.pdf for more knowledge


# In[3]:


def info_molecule(molecule):
    num_atoms=len(molecule.atoms)
    name_atoms=[];
    mass_atoms=[];
    for i in range(num_atoms):
        name_atoms.append(molecule.atoms[i].element_name)
        mass_atoms.append(molecule.atoms[i].mass)
    molecular_mass=sum(mass_atoms)
    return name_atoms,mass_atoms,molecular_mass


# In[4]:


def remove_duplicate(x):
    return list(dict.fromkeys(x))


# In[5]:


def basis_set(element_symbol):
    #HERE I should define all the elements and all the basis set
    if element_symbol=='H':
        return "TZV2PX-MOLOPT-GTH"
    elif element_symbol=='F':
        return "TZV2PX-MOLOPT-GTH"
    elif element_symbol=='Cl':
        return "TZV2PX-MOLOPT-GTH"
    elif element_symbol=='I':
        return "DZVP-MOLOPT-SR-GTH"
        


# In[6]:


def potential(element_symbol,functional):
    return "GTH-"+functional
        


# In[8]:


def optimize_molecule(molecule,functional,length,dir):
    molecule.save('mol_unopt_coord.xyz',overwrite='True')
    with open('mol_unopt_coord.xyz', 'r') as fin:
        data = fin.read().splitlines(True)
    with open('mol_unopt_coord.xyz', 'w') as fout:
        fout.writelines(data[2:])
    molecule=molecule.to_parmed()
    atom_list,mass_list,total_mass=info_molecule(molecule)
    unique_atom_list=remove_duplicate(atom_list)
    num_atoms=len(atom_list)
    num_unique_atoms=len(unique_atom_list)

    mySim = sim.SIM()

    mySim.GLOBAL.RUN_TYPE = "GEO_OPT"
    mySim.GLOBAL.PROJECT  = "molecule_opt"
    mySim.GLOBAL.PRINT_LEVEL = "LOW"
    #FORCE EVAL SECTION
    mySim.FORCE_EVAL.METHOD='QUICKSTEP'
    mySim.FORCE_EVAL.SUBSYS.CELL.ABC='{L} {L} {L}'.format(L=length)
    mySim.FORCE_EVAL.SUBSYS.COORD.DEFAULT_KEYWORD='mol_unopt_coord.xyz'
    mySim.FORCE_EVAL.SUBSYS.init_atoms(num_atoms);
    for i in range(num_unique_atoms):
        mySim.FORCE_EVAL.SUBSYS.KIND[i+1].SECTION_PARAMETERS=unique_atom_list[i]
        mySim.FORCE_EVAL.SUBSYS.KIND[i+1].BASIS_SET=basis_set(unique_atom_list[i])
        mySim.FORCE_EVAL.SUBSYS.KIND[i+1].POTENTIAL=potential(unique_atom_list[i],functional)

    mySim.FORCE_EVAL.DFT.BASIS_SET_FILE_NAME=dir+'BASIS_MOLOPT'
    mySim.FORCE_EVAL.DFT.POTENTIAL_FILE_NAME=dir+'GTH_POTENTIALS'
    mySim.FORCE_EVAL.DFT.QS.EPS_DEFAULT=1E-10

    mySim.FORCE_EVAL.DFT.MGRID.CUTOFF=400
    mySim.FORCE_EVAL.DFT.MGRID.REL_CUTOFF=40
    mySim.FORCE_EVAL.DFT.MGRID.NGRIDS=5

    mySim.FORCE_EVAL.DFT.XC.XC_FUNCTIONAL.SECTION_PARAMETERS=functional
    mySim.FORCE_EVAL.DFT.XC.VDW_POTENTIAL.POTENTIAL_TYPE='PAIR_POTENTIAL'
    mySim.FORCE_EVAL.DFT.XC.VDW_POTENTIAL.PAIR_POTENTIAL.TYPE='DFTD3'
    mySim.FORCE_EVAL.DFT.XC.VDW_POTENTIAL.PAIR_POTENTIAL.PARAMETER_FILE_NAME='dftd3.dat'
    mySim.FORCE_EVAL.DFT.XC.VDW_POTENTIAL.PAIR_POTENTIAL.REFERENCE_FUNCTIONAL=functional
    mySim.FORCE_EVAL.DFT.XC.VDW_POTENTIAL.PAIR_POTENTIAL.R_CUTOFF=11

    mySim.FORCE_EVAL.DFT.SCF.SCF_GUESS='ATOMIC'
    mySim.FORCE_EVAL.DFT.SCF.MAX_SCF=200
    mySim.FORCE_EVAL.DFT.SCF.EPS_SCF=1E-6

    mySim.MOTION.GEO_OPT.TYPE='MINIMIZATION'
    mySim.MOTION.GEO_OPT.OPTIMIZER='BFGS'
    mySim.MOTION.GEO_OPT.MAX_ITER=200
    mySim.MOTION.GEO_OPT.MAX_DR=1e-3

    mySim.MOTION.CONSTRAINT.FIXED_ATOMS.LIST ='1'
    mySim.write_changeLog(fn="mol_opt-changeLog.out")
    mySim.write_errorLog()
    mySim.write_inputFile(fn='mol_opt.inp')


# In[ ]:




