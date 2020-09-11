import numpy as np
import mdtraj as md
import mbuild as mb

from cp2kmdpy import runners
import os
import unyt as u




def info_molecule(molecule):
    num_atoms=len(molecule.atoms)
    name_atoms=[];
    mass_atoms=[];
    for i in range(num_atoms):
        name_atoms.append(molecule.atoms[i].element_name)
        mass_atoms.append(molecule.atoms[i].mass)

    return name_atoms,mass_atoms


def remove_duplicate(x):
    return list(dict.fromkeys(x))

def is_cubic(box):
    angles=box.angles
    for angle in angles:
        if math.isclose(angle, 90.0):
            continue
        else:
            return False
            break
    lengths=box.lengths
    a=lengths[0]
    for length in lengths:
        if math.isclose(length, a):
            continue
        else:
            return False
            break
    return True


def basis_set_setter(element_symbol):
    #HERE I should define all the elements and all the basis set
    if element_symbol=='H':
        return "TZV2PX-MOLOPT-GTH"
    elif element_symbol=='F':
        return "TZV2PX-MOLOPT-GTH"
    elif element_symbol=='Cl':
        return "TZV2PX-MOLOPT-GTH"
    elif element_symbol=='I':
        return "DZVP-MOLOPT-SR-GTH"


def potential(element_symbol,functional):
    return "GTH-"+functional

def pressure_ensemble(val):
    if val=='NPE_F' or val=='NPE_I' or val=='NPT_F' or val=='NPT_I':
        return True
    else:
        return False

def temperature_ensemble(val):
    if val=='MSST' or val=='MSST_DAMPED' or val=='NPT_F' or val=='NPT_I' or val=='NVT' or val=='NVT_ADIABATIC':
        return True
    else:
        return False


class MD():
    r"""Base class for running molecule optimization.
    :param molecules: Molecule(s) in the simulation box
    :type molecule: list, each element is an mBuild molecule
    :param functional: DFT XC functional to be used
    :type functional: string
    :param box: Box in which the molecules should be placed
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
    :param simulation_time: Total time of the simulation
    :type simulation_time: unyt quantity, optional
    :param time_step: Time step for integrating the trajectory
    :type time_step: unyt qunatity, example: time_step=0.1*u.fs # when time step is 0.1 femtosecond
    :param ensemble: MD simulation ensemble
    :type ensemble: string, optional
    :param project_name: Name of the project
    :type project_name: string, optional
    :param temperature: Simulation temperature
    :type temperature: unyt quantity
    :param pressure: Simulation pressure
    :type pressure: unyt quantity
    :param n_molecules: Number of molecules in the box of each kind
    :type n_molecules: list of number of molecules of each kind
    :param thermostat: type of thermostat
    :type thermostat: string, optional
    :param traj_type: Output trajectory format
    :type traj_type: string, optional
    :param traj_freq: Trjectory snapshot capture frequency
    :type traj_freq: integer, optional
    :param seed: Random number seed for MD simulation
    :type seed: integer, optional
    """

    def __init__(self,molecules=None, functional=None,box=None,cutoff=None, scf_tolerance=None,
                 basis_set=[None],basis_set_filename=None, potential_filename=None,fixed_list=None, periodicity=None,simulation_time=None,time_step=None,ensemble=None,project_name=None,temperature=None,pressure=None,n_molecules=None,thermostat=None,traj_type=None,traj_freq=None,seed=None,input_filename=None,output_filename=None):
        self.molecules=molecules;
        self.functional=functional;

        self.box=box;
        self.cutoff=cutoff;
        self.scf_tolerance=scf_tolerance;
        self.basis_set=basis_set;
        self.basis_set_filename=basis_set_filename;
        self.potential_filename=potential_filename;
        self.periodicity=periodicity
        self.fixed_list=fixed_list;
        self.simulation_time=simulation_time
        self.time_step=time_step
        self.ensemble=ensemble
        self.project_name=project_name
        self.temperature=temperature
        self.pressure=pressure
        self.n_molecules=n_molecules
        self.thermostat=thermostat
        self.traj_type=traj_type
        self.traj_freq=traj_freq
        self.seed=seed
        self.input_filename=input_filename
        self.output_filename=output_filename
    
    def md_initialization(self):
        molecules=self.molecules;
        functional=self.functional;
        box=self.box;
        cutoff=self.cutoff
        scf_tolerance=self.scf_tolerance
        basis_set=self.basis_set
        basis_set_filename=self.basis_set_filename
        potential_filename=self.potential_filename;
        fixed_list=self.fixed_list;
        periodicity=self.periodicity
        simulation_time=self.simulation_time
        time_step=self.time_step
        ensemble=self.ensemble
        project_name=self.project_name
        temperature=self.temperature
        pressure=self.pressure
        n_molecules=self.n_molecules
        thermostat=self.thermostat
        traj_type=self.traj_type
        traj_freq=self.traj_freq
        seed=self.seed
        
        
        atom_list=[];
        mass_list=[];
        for i in range(len(molecules)):
            current_molecule=mb.clone(molecules[i])
            current_molecule_pmd=current_molecule.to_parmed()
            x,y=info_molecule(current_molecule_pmd);
            atom_list.extend(x)
            mass_list.extend(y)
        unique_atom_list=remove_duplicate(atom_list)
        num_atoms=len(atom_list)
        num_unique_atoms=len(unique_atom_list)
        unique_atom_list.sort()

        if project_name==None:
            self.project_name='sample_project'
            print('project_name not specified, set as sample_project')
        if cutoff==None:
            self.cutoff=600;
            print('cutoff not specified, set as 600')
        if scf_tolerance==None:
            self.scf_tolerance=1e-6
            print('scf_tolerance not specified, set as 1e-6')
        if basis_set_filename==None:
            self.basis_set_filename='BASIS_MOLOPT'
            print('basis_set_filename not defined, set as BASIS_MOLOPT')
        if potential_filename==None:
            self.potential_filename='GTH_POTENTIALS'
            print('potential_filename not specified, set as GTH_POTENTIALS')
        if periodicity==None:
            self.periodicity='XYZ';
            print('periodicity not specified, set as XYZ')
        if simulation_time==None:
            self.simulation_time=1*u.ps #ps
            print('simulation_time not specified, set as 1 ps')

        if time_step==None:
            lightest=min(mass_list);

            if lightest <1.5:
                self.time_step=0.5*u.fs
                print('time_step not specified, time_step set as 0.5 fs as the lighest element has mass {} au'.format(lightest))
            elif (lightest>=1.5) and (lightest<40):
                self.time_step=1*u.fs
                print('time_step not specified, time_step set as 1 fs as the lighest element has mass {} au'.format(lightest))
            if lightest>=40:

                self.time_step=1.5*u.fs
                print('time_step not specified, time_step set as 1.5 fs as the lighest element has mass {} au'.format(lightest))

        if ensemble==None:
            self.ensemble='NVE'
            print('ensemble not specified, set as NVE')
        if traj_type==None:
            self.traj_type='XYZ'
            print('output trajectory format set as XYZ')
        if traj_freq==None:
            self.traj_freq=10
        if seed == None:
            self.seed=0


        if self.input_filename==None:
            self.input_filename=self.project_name+'_md_input.inp'
            print('input_filename not specified, set as {}'.format(self.input_filename))
        if self.output_filename==None:
            self.output_filename=self.project_name+'_md_output.out'
            print('output_filename not specified, set as {}'.format(self.output_filename))
        
        output_pos_filename=self.project_name+"-pos-1.xyz"

        print('Output position filename is {}'.format(output_pos_filename))

        if self.temperature is not None:

            self.temperature=(temperature.to('K')).value
        if self.pressure is not None:
            self.pressure=(pressure.to('bar')).value
        
        self.simulation_time=(self.simulation_time.to('fs')).value
        self.time_step=(self.time_step.to('fs')).value

        print('You can change default settings in setter.md_files')
        
