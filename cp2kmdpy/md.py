import numpy as np
import mdtraj as md
import mbuild as mb
import setter
from cp2kmdpy import runners
import os
import unyt as u

class MD():
    def __init__(self,molecules=None, functional=None,box=None,cutoff=None, scf_tolerance=None,
                 basis_set=[None],basis_set_filename=None, potential_filename=None,fixed_list=None, periodicity=None,simulation_time=None,time_step=None,ensemble=None,project_name=None,temperature=None,pressure=None,n_molecules=None,thermostat=None,traj_type=None,traj_freq=None,seed=None):
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
    
    def md_files(self):
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

        if project_name==None:
            project_name='sample_project'
        if cutoff==None:
            cutoff=600;
        if scf_tolerance==None:
            scf_tolerance=1e-6
        if basis_set_filename==None:
            basis_set_filename='BASIS_MOLOPT'
        if potential_filename==None:
            potential_filename='GTH_POTENTIALS'
        if periodicity==None:
            periodicity='XYZ';
        if simulation_time==None:
            simulation_time=1*u.ps #ps
        if time_step==None:
            time_step=0.1*u.fs #fs

        if ensemble==None:
            ensemble='NVT'
        if traj_type==None:
            traj_type='XYZ'
        if traj_freq==None:
            traj_freq=10

        if temperature is not None:

            temperature=(temperature.to('K')).value
        if pressure is not None:
            pressure=(pressure.to('bar')).value
        simulation_time=(simulation_time.to('ps')).value
        time_step=(time_step.to('fs')).value

        print('You can change default settings in setter.md_files')
        setter.md_files(molecules,functional,box,cutoff,scf_tolerance, basis_set, basis_set_filename,potential_filename, fixed_list, periodicity,simulation_time,time_step,ensemble,project_name,temperature,pressure,n_molecules,thermostat,traj_type,traj_freq,seed)
