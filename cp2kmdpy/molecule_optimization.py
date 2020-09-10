import numpy as np
import mdtraj as md
import mbuild as mb
import setter
from cp2kmdpy import runners
import os
class Molecule_optimization():
    r"""Base class for running molecule optimization.
    :param molecule: Molecule whose structure needs to be optimized
    :type molecule: mBuild molecule
    :param functional: DFT XC functional to be used
    :type functional: string
    :param box: Box in which the molecule should be placed
    :type box: mBuild box
    :param cutoff: Plane wave cutoff (Ry) for DFT calculation
    :type cutoff: float, optional
    :param scf_tolerance: Tolerance for each SCF cycle
    :type scf_tolerance: float, optional
    :param basis_set: Basis set for each atomic kind
    :type basis_set: dictionary with key being the atomic symbol and value being the basis set, both strings
    :param basis_set_filename: Filename for the basis set
    :type basis_set_filename: string, optional, defaults to BASIS_MOLOPT
    :param potential_filename: Filename for the pseudopotential to be used
    :type potential_filename: string, optional, defaults to GTH_POTENTIALS
    :param fixed_list: list of elements to be frozen
    :type fixed_list: string, optional, for example: if atom # 1 to 10 shall be fixed => fixed_list='1..10'
    :param periodicity: Periodicity of the box
    :type periodiicity: string, optional, defaults to 'XYZ'
    :param n_iter: Number of iterations for geometru optimization
    :type n_iter: positive integer, optional

    """
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

        if cutoff==None:
            cutoff=900;
            print('cutoff not specified, set as 900')
        if scf_tolerance==None:
            scf_tolerance=1e-6
            print('scf_tolerance not specified, set as 1e-6')
        if basis_set_filename==None:
            basis_set_filename='BASIS_MOLOPT'
            print('basis_set_filename not defined, set as BASIS_MOLOPT')
        if potential_filename==None:
            potential_filename='GTH_POTENTIALS'
            print('potential_filename not specified, set as GTH_POTENTIALS')
        if periodicity==None:
            periodicity='XYZ';
            print('periodicity not specified, set as XYZ')
        if fixed_list==None:
            fixed_list=1
            print('fixed_list not specified, set as 1')


        if n_iter==None:
            n_iter=10
            print('n_iter not specified, set as 10')


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


