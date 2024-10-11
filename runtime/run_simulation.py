import argparse
from pathlib import *
import argparse
from src.openmm_runner import OpenMMArgs, run_openmm
from src.preprocess_pdb import prep_pdb

def run_md(arg_parse_dict: dict) -> None:
    stage_args = OpenMMArgs(yaml_path=arg_parse_dict['yaml_path'])
    for ind, stage_defn in enumerate(stage_args):
        if ind == 0:
            prep_pdb(arg_parse_dict['pdb_path'], *stage_defn)
        else:
            run_openmm(*stage_defn)

def parse_args() -> dict:
    parser=argparse.ArgumentParser()
    parser.add_argument('pdb_file', type = str)
    parser.add_argument('defn_file', type = str)
    args=parser.parse_args()
    return args

if __name__ == '__main__':
    inputs = parse_args()
    run_md(arg_parse_dict = inputs)
