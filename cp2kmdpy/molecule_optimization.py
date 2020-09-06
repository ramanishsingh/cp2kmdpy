import numpy as np
import mdtraj as md
import mbuild as mb
import setter
from cp2kmdpy import runners
import os
class Molecule_optimization():
    def __init__(self,molecule=None, functional=None,box=None,cutoff=None, scf_tolerance=None,
                 basis_set=[None],basis_set_filename=None, potential_filename=None,fixed_list=None, periodicity=None,n_iter=None):
        self.molecule=molecule;
        self.functional=functional;
        self.box=box;
        self.cutoff=cutoff;
        self.scf_tolerance=scf_tolerance;
        self.basis_set=basis_set;
        self.basis_set_filename=basis_set_filename;
        self.potential_filename=potential_filename;
        self.periodicity=periodicity
        self.fixed_list=fixed_list;
        self.n_iter=n_iter;
    
    def optimize_files(self):
        
        
        molecule=self.molecule;
        functional=self.functional;
        box=self.box;
        cutoff=self.cutoff
        scf_tolerance=self.scf_tolerance
        basis_set=self.basis_set
        basis_set_filename=self.basis_set_filename
        potential_filename=self.potential_filename;
        fixed_list=self.fixed_list;
        periodicity=self.periodicity
        n_iter=self.n_iter
        print('You can change default settings in setter.single_molecule_opt_file')
        setter.single_molecule_opt_files(molecule,functional,box,cutoff,scf_tolerance, basis_set, basis_set_filename,potential_filename, fixed_list, periodicity,n_iter)

    def run_optimization(self):
        runners.run_single_molecule_optimization(self)

    def return_optimized_molecule(self):
        output_pos_file=self.molecule.name+'_opt-pos-1.xyz'
        self.molecule.save('delete.pdb',overwrite=True)
        trj=md.load(output_pos_file, top='delete.pdb')
        self.molecule.xyz=trj.xyz[-1];
        os.remove('delete.pdb')
        return self.molecule


