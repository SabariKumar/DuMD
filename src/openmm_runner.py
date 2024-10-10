import openmm
from multiprocessing import Process, Pipe
from pathlib import Path
from time import perf_counter

def run_openmm(stage, openmm_args: OpenMMArgs):
    try:
        pdb = openmm.PDBFile(pdb_path)
        base = Path(pdb_path).stem
        forcefield = openmm.ForceField(forcefield_file, water_file)
        system = forcefield.createSystem(pdb.topology,
                                         nonbondedMethod = nonbondedMethod,
                                         nonbondedCutoff = nonbondedCutoff * openmm.unit.nanometer,
                                         constraints = constraints)
        integrator = integrator(temp*openmm.unit.kelvin, friction_coeff/openmm.unit.picosecond, step_size*openmm.unit.picoseconds)
        simulation = openmm.Simulation(pdb.topology, system, integrator)
        simulation.context.setPositions(pdb.positions)
        if minimizeEnergy:
            simulation.minimizeEnergy()
        simulation.context.setVelocitiesToTemperature(temp * openmm.unit.kelvin)
        simulation.step(n_steps)
        if n_steps > 10000:
            simulation.reporters.append(openmm.PDBReporter(f'{base}_{stage}.pdb', 5000))

        simulation.reporters.append(openmm.StateDataReporter('{base}_{stage}.log',
                                                             10000, totalSteps = n_steps,
                                                             step=True, time=True, 
                                                             speed=True, progress=True,
                                                             elapsedTime=True, potentialEnergy=True,
                                                             kineticEnergy=True, totalEnergy=True,
                                                             temperature=True, volume=True,
                                                             density=True, separator = "\t"))
        simulation.step(n_steps)
        positions = simulation.context.getState(getPositions=True).getPositions()
        openmm.PDBFile.writeFile(simulation.topology, positions,             
                          open(f'{base}_{stage}_final_pos.pdb', 'w'))

    except ValueError as e:
        print(f"Encountered ValueError of pdb {pdb_path}! \n")
        print(e)

    except AttributeError as e:
        print(f"Stage {stage} is missing required a required attribute")
        print(e)

class OpenMMArgs:
    # Wrapper class for run_openmm arguments.
    # Create instances from "stage" key-level input yaml dicts to set params for 
    # each stage of the production MD pipeline.  
    def __init__(self, stage: str, param_dict: dict):
        self.stage = stage
        for key, val in param_dict.items():
            setattr(self, key, val)   