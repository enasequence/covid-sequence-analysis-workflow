#!/usr/bin/env python3
import os
import subprocess
import sys
import time

# DIR where the current script resides
DIR = os.path.dirname(os.path.abspath(__file__))

# set default values for arguments if not provided
pipeline = sys.argv[1] if len(sys.argv) > 1 else 'nanopore'
profile = sys.argv[2] if len(sys.argv) > 2 else 'gls'
root_dir = sys.argv[3] if len(sys.argv) > 3 else 'gs://prj-int-dev-covid19-nf-gls'
snapshot_date = sys.argv[4] if len(sys.argv) > 4 else '2023-03-15'
concurrency = sys.argv[5] if len(sys.argv) > 5 else '100'
batch_size = sys.argv[6] if len(sys.argv) > 6 else '15000'
dataset_name = sys.argv[7] if len(sys.argv) > 7 else 'sarscov2_metadata'
project_id = sys.argv[8] if len(sys.argv) > 8 else 'prj-int-dev-covid19-nf-gls'
test_submission = 'false'

# Row count and batches
table_name = f"{pipeline}_to_be_processed"
sql = f"SELECT count(*) AS total FROM {project_id}.{dataset_name}.{table_name}"
row_count = subprocess.check_output(f"bq --project_id={project_id} --format=csv query --use_legacy_sql=false \"{sql}\" | grep -v total", shell=True).decode('utf-8')

queue_size = 100  # as defined as queueSize in nextflow.config
batches = (int(row_count) + int(batch_size) - 1) // int(batch_size)
num_of_jobs = (int(concurrency) + queue_size - 1) // queue_size

input_dir = os.path.join(DIR, 'data', snapshot_date)
os.makedirs(input_dir, exist_ok=True)

print(f"Pipeline: {pipeline}")
print(f"batches: {batches}")
print(f"num_of_jobs: {num_of_jobs}")
print(f"input_dir: {input_dir}")
time.sleep(5)

for j in range(num_of_jobs):
    if j >= batches:
        break
    # bsub -n 2 -M 8192 -q production
    job_dir = os.path.join(root_dir, snapshot_date, f"{pipeline}_{j}")
    os.makedirs(job_dir, exist_ok=True)
    os.chdir(job_dir)
    offset = j * int(batch_size)
    subprocess.call([f"{DIR}/manage_nf.sh", str(j), batch_size, project_id, dataset_name, table_name,
                     input_dir, pipeline, profile, root_dir, snapshot_date, test_submission])
    
# Move set.archived.sh and update.receipt.sh outside the loop
# Comment out update.receipt.sh
# subprocess.call([f"{DIR}/update.receipt.sh", batch_index, snapshot_date, pipeline, profile, root_dir, dataset_name, project_id])
subprocess.call([f"{DIR}/set.archived.sh", dataset_name, project_id])

sql = f"CREATE OR REPLACE TABLE {dataset_name}.sra_processing AS SELECT DISTINCT * FROM {dataset_name}.sra_processing"
subprocess.call([f"bq --project_id={project_id} --format=csv query --use_legacy_sql=false \"{sql}\""], shell=True)

# max_mem avg_mem swap stat exit_code exec_cwd exec_host
# bjobs -u all -d -o "jobid job_name user submit_time start_time finish_time run_time cpu
