import openmm
from multiprocessing import Process, Pipe
from pathlib import Path
import yaml
from collections import Generator
from typing import Union
from time import perf_counter

def run_openmm(pdb_path: Union[str, Path], stage: str, openmm_args: dict, *args, **kwargs) -> str:
    try:
        pdb = openmm.PDBFile(pdb_path)
        base = Path(pdb_path).stem
        forcefield = openmm.ForceField(openmm_args['forcefield_file'], openmm_args['water_file'])
        system = forcefield.createSystem(pdb.topology,
                                         nonbondedMethod = openmm_args['nonbondedMethod'],
                                         nonbondedCutoff = openmm_args['nonbondedCutoff'] * openmm.unit.nanometer,
                                         constraints = openmm_args['constraints'])
        integrator = integrator(openmm_args['temp'] * openmm.unit.kelvin, openmm_args['friction_coeff']/openmm.unit.picosecond, openmm_args['step_size'] * openmm.unit.picoseconds)
        simulation = openmm.Simulation(pdb.topology, openmm_args['system'], openmm_args['integrator'])
        simulation.context.setPositions(pdb.positions)
        if openmm_args.minimizeEnergy:
            simulation.minimizeEnergy()
        simulation.context.setVelocitiesToTemperature(openmm_args['temp'] * openmm.unit.kelvin)
        if openmm_args.n_steps > 10000:
            simulation.reporters.append(openmm.PDBReporter(f'{base}_{stage}.pdb', 5000))

        simulation.reporters.append(openmm.StateDataReporter('{base}_{stage}.log',
                                                             10000, totalSteps = openmm_args['n_steps'],
                                                             step=True, time=True, 
                                                             speed=True, progress=True,
                                                             elapsedTime=True, potentialEnergy=True,
                                                             kineticEnergy=True, totalEnergy=True,
                                                             temperature=True, volume=True,
                                                             density=True, separator = "\t"))
        simulation.step(openmm_args['n_steps'])
        positions = simulation.context.getState(getPositions=True).getPositions()
        openmm.PDBFile.writeFile(simulation.topology, positions,             
                          open(f'{base}_{stage}_final_pos.pdb', 'w'))
        return '{base}_{stage}_final_pos.pdb'

    except ValueError as e:
        print(f"Encountered ValueError of pdb {openmm_args['pdb_path']}! \n")
        print(e)

    except KeyError as e:
        print(f"Stage {stage} is missing required a required attribute. Please correct your job definition yaml. \n")
        print(e)

class OpenMMArgs(Generator):
    # Generator wrapper class for run_openmm arguments.
    # Create instances from "stage" key-level input yaml dicts to set params for 
    # each stage of the production MD pipeline.

    def __init__(self, yaml_path: Union[Path, str]):
        # The yaml definition file should consist of "stages" as the top keys,
        # with params as key/value pairs following
        with open(yaml_path, 'r') as yaml_file:
            params_dict = yaml.safe_load(yaml_file)
        self.params = [(k, v) for k, v in params_dict.items()]
        self.current_stage_index = 0
        self.max_stage_index = len(self.params)

    def send(self, _):
        if self.current_stage_index < self.max_stage_index:
            return_params = self.params[self.current_stage_index]
            self.current_stage_index += 1
            return return_params
        raise StopIteration
    
    def throw(self, type, value=None, traceback=None):
        super().throw(type, value, traceback)