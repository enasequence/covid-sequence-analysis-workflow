#!/usr/bin/env python3
import subprocess
import sentry_sdk
from sentry_sdk import push_scope, capture_exception
from sentry_sdk.integrations.asyncio import AsyncioIntegration
from sentry_sdk.integrations.excepthook import ExcepthookIntegration
import os
import argparse

SENTRY_MSG_LIMIT_LEN=1000
print(f"SENTRY URL: {os.environ.get('SENTRY_URL')}")
sentry_sdk.init(f"{os.environ['SENTRY_URL']}", 
                integrations=[
                    AsyncioIntegration(),
                    ExcepthookIntegration(always_run=True),
                ],
                attach_stacktrace=True,)
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

parser.add_argument('-c',
                    '--task_cpus',
                    help="task cpu",
                    type=str,
                    required=True)

parser.add_argument('-s',
                    '--study_accession',
                    help="study_accession",
                    type=str,
                    default='PRJEB45555',
                    required=False)

args = parser.parse_args()

class FailedDownloadException(Exception):
    pass
class FailedProcessException(Exception):
    pass
class FailedConsensusException(Exception):
    pass

def download(run_accession, projects_accounts_csv, input_file_1, input_file_2, study_accession:str='PRJEB45555'):
    # Define the Bash command to run
    bash_command=f"bash download.sh {run_accession} {projects_accounts_csv} {input_file_1} {input_file_2} {study_accession}"
    # Run the Bash command and capture the output
    output = subprocess.run(bash_command, capture_output=True, text=True, shell=True)
    print(output)
    if output.returncode!=0:
        raise FailedDownloadException(output.stderr[-SENTRY_MSG_LIMIT_LEN:])
    return output

def run_process(run_accession, sars2_fasta, task_cpus):
    # Define the Bash command to run
    bash_command=f"bash map_to_ref.sh {run_accession} {sars2_fasta} {task_cpus}"
    # Run the Bash command and capture the output
    output = subprocess.run(bash_command, capture_output=True, text=True, shell=True)
    print(output)
    if output.returncode!=0:
        raise FailedProcessException(output.stderr[-SENTRY_MSG_LIMIT_LEN:])
    return output

def generate_zip(run_accession):
    for file_type in ["_consensus.fasta", ".coverage", "annot.vcf"]:
        bash_command=f"bgzip {run_accession}{file_type}"
        output = subprocess.run(bash_command, capture_output=True, text=True, shell=True)
        print(output)
        print(f"subprocess output stderr limit: {output.stderr[-SENTRY_MSG_LIMIT_LEN:]}")
        if output.returncode!=0:
            raise FailedConsensusException(output.stderr[-SENTRY_MSG_LIMIT_LEN:])
    return 

def generate_consensus(run_accession, sars2_fasta):
    output = subprocess.run(f"vcf_to_consensus.py -dp 10 -af 0.25 -v {run_accession}.vcf.gz -d {run_accession}.coverage -o headless_consensus.fasta -n {run_accession} -r {sars2_fasta}", capture_output=True, text=True, shell=True)
    print(output)
    print(f"subprocess output stderr limit: {output.stderr[-SENTRY_MSG_LIMIT_LEN:]}")
    if output.returncode!=0:
        raise Exception(output.stderr[-SENTRY_MSG_LIMIT_LEN:])
    
    output = subprocess.run(f"fix_consensus_header.py headless_consensus.fasta > {run_accession}_consensus.fasta", capture_output=True, text=True, shell=True)
    print(output)
    print(f"subprocess output stderr limit: {output.stderr[-SENTRY_MSG_LIMIT_LEN:]}")
    if output.returncode!=0:
        raise Exception(output.stderr[-SENTRY_MSG_LIMIT_LEN:])
    

if __name__ == '__main__':
    print(args)
    try:
        download(args.run_accession, args.projects_accounts_csv, args.input_file_1, args.input_file_2, args.study_accession)
        run_process(args.run_accession, args.sars2_fasta, args.task_cpus)
        generate_consensus(args.run_accession, args.sars2_fasta)
        generate_zip(args.run_accession)
    except Exception as e:
        capture_exception(e)
    finally:
        sentry_sdk.flush()