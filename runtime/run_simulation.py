import argparse
from pathlib import *
import argparse
from src.openmm_runner import OpenMMArgs, run_openmm
from src.preprocess_pdb import prep_pdb

def run_md(arg_parse_dict) -> None:
    stage_args = OpenMMArgs(yaml_path=arg_parse_dict.defn_file)
    pdb_path = ''
    for ind, stage_defn in enumerate(stage_args):
        if ind == 0:
            pdb_path = prep_pdb(arg_parse_dict.pdb_file, *stage_defn)
            print(pdb_path)
        else:
            print(pdb_path)
            print(f"Stage defn:\n{stage_defn[0]}")
            stage_name, stage_args = stage_defn
            pdb_path = run_openmm(pdb_path = pdb_path, stage = stage_name, openmm_args = stage_args)

def parse_args() -> dict:
    parser=argparse.ArgumentParser()
    parser.add_argument('pdb_file', type = str)
    parser.add_argument('defn_file', type = str)
    args=parser.parse_args()
    return args

if __name__ == '__main__':
    inputs = parse_args()
    run_md(arg_parse_dict = inputs)
