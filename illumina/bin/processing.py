#!/usr/bin/env python3
import subprocess
import sentry_sdk
from pathlib import Path
import os
print(os.environ['HOME'])
import argparse

parser = argparse.ArgumentParser(description='Start processing script')

parser.add_argument('-r',
                    '--run_accession',
                    help="run_accession",
                    type=str,
                    required=True)
parser.add_argument('-p',
                    '--projects_accounts_csv',
                    help="projects_accounts_csv",
                    type=str,
                    required=True)
parser.add_argument('-i1',
                    '--input_file_1',
                    help="input_file_1",
                    type=str,
                    required=True)
parser.add_argument('-i2',
                    '--input_file_2',
                    help="input_file_2",
                    type=str,
                    required=True)
parser.add_argument('-f',
                    '--sars2_fasta',
                    help="sars2_fasta",
                    type=str,
                    required=True)
parser.add_argument('-s',
                    '--study_accession',
                    help="study_accession",
                    type=str,
                    required=False)

def run_process(run_accession, projects_accounts_csv, input_file_1, input_file_2, sars2_fasta, study_accession:str='PRJEB45555'):
    # Define the Bash command to run
    bash_command=f"bash {Path.cwd()}/processing.sh {run_accession} {projects_accounts_csv} {input_file_1} {input_file_2} {sars2_fasta} {study_accession}"
    # Run the Bash command and capture the output
    output = subprocess.check_output(bash_command, shell=True)
    # Print the output
    print(output)

if __name__ == '__main__':
    print(f"SENTRY URL: {os.environ.get('SENTRY_URL')}")
    sentry_sdk.init(f"{os.environ['SENTRY_URL']}")
    try:
        run_process(args.run_accession, args.projects_accounts_csv, args.input_file_1, args.input_file_2, args.sars2_fasta, args.study_accession)
    except Exception as e:
        sentry_sdk.capture_exception(e)