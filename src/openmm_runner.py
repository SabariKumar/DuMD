import openmm
from multiprocessing import Process, Pipe
from pathlib import Path
from time import perf_counter

def run_openmm(pdb_path: str,
               n_steps: int,
               stage: str,
               forcefield_file: str = 'amber14-all.xml',
               water_file: str = 'amber14/tip3p.xml',
               nonbondedMethod = openmm.app.PME,
               nonbondedCutoff: int = 1,
               constraints = openmm.HBonds,
               integrator = openmm.LangevinMiddleIntegrator,
               temp: int = 300,
               friction_coeff: int = 1,
               step_size: float = 0.002,
               minimizeEnergy: bool = True,
               ):
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

    except:
        print(f"Encountered other error of pdb {pdb_path}! \n")
        print(e)