import os
import subprocess

class RunNextflow:
    def __init__(self):
        self.dir = os.path.dirname(os.path.abspath(__file__))
        self.batch_input = ''
        self.pipeline = 'nanopore'
        self.profile = 'codon'
        self.root_dir = 'gs://prj-int-dev-covid19-nf-gls'
        self.batch_index = '0'
        self.snapshot_date = '2023-03-13'
        self.test_submission = 'true'
        self.study_accession = 'PRJEB45555'
        self.dataset_name = 'sarscov2_metadata'
        self.project_id = 'prj-int-dev-covid19-nf-gls'

    def process_batch(self):
        print("")
        print(f"** Processing samples with {self.dir}/{self.pipeline}/{self.pipeline}.nf. **")

        if self.profile == "awsbatch":
            print(f"** Retrieve secrets to {self.dir}/{self.project_id}-sa-credential.json **")
            subprocess.run(['aws', 'secretsmanager', 'get-secret-value', '--secret-id', GOOGLE_APPLICATION_CREDENTIALS_SECRET_ID, '--query', 'SecretString', '--output', 'text', '>', f"{self.dir}/{self.project_id}-sa-credential.json"])
            subprocess.run(['gcloud', 'auth', 'activate-service-account', '--key-file', f"{self.dir}/{self.project_id}-sa-credential.json"])
            subprocess.run(['gcloud', 'config', 'set', 'project', self.project_id])
            project_bucket = "prj-int-dev-ait-eosc-aws-eval"
            subprocess.run(['aws', 's3', 'cp', self.batch_input, f"{self.dir}/data/"])  # download sample index file from s3 to local dir
            subprocess.run(['aws', 's3', 'cp', f"s3://{project_bucket}/{self.dataset_name}/", f"{self.dir}/data/", '--recursive', '--exclude', '"*/*"'])  # download projects_accounts and .fa files
            self.batch_input = f"{self.dir}/data/{os.path.basename(self.batch_input)}"  # local path to sample index file

        pipeline_dir = f"{self.root_dir}/{self.snapshot_date}/{self.pipeline}_{self.batch_index}"
        print(f"** pipeline_dir: {pipeline_dir} **")

        subprocess.run(['nextflow', '-C', f"{self.dir}/nextflow-lib/nextflow.config", 'run', f"{self.dir}/{self.pipeline}/{self.pipeline}.nf", '-profile', self.profile,
                        '--TEST_SUBMISSION', self.test_submission, '--STUDY', self.study_accession,
                        '--CONFIG_YAML', f"{self.dir}/{self.pipeline}/config.yaml",
                        '--SECRETS', f"{self.dir}/data/projects_accounts.csv",
                        '--SARS2_FA', f"{self.dir}/data/NC_045512.2.fa",
                        '--SARS2_FA_FAI', f"{self.dir}/data/NC_045512.2.fa.fai",
                        '--INDEX', self.batch_input,
                        '--OUTDIR', f"{pipeline_dir}/publishDir",
                        '--STOREDIR', f"{pipeline_dir}/storeDir",
                        '-w', f"{pipeline_dir}/workDir",
                        '-with-tower'])

        subprocess.run([f"{self.dir}/update.receipt.sh", self.batch_index, self.snapshot_date, self.pipeline, self.profile, self.root_dir, self.dataset_name, self.project_id])
        subprocess.run([f"{self.dir}/set.archived.sh", self.dataset_name, self.project_id])

        if self.profile not in ['gls', 'awsbatch']:
            subprocess.run(['rm', '-R', f"{pipeline_dir}/workDir"])
            subprocess.run(['rm', '-R', f"{pipeline_dir}/storeDir"])
            subprocess.run(['rm', '-R', f"{pipeline_dir}/publishDir"])

        if self.profile == "awsbatch":
            subprocess.run(['aws', 's3', 'rm', '--recursive', f"{pipeline_dir}/workDir", '--quiet'])
            subprocess.run(['aws', 's3', 'rm', '--recursive', f"{pipeline_dir}/storeDir", '--quiet'])
            subprocess.run(['aws', 's3', 'rm', '--recursive', f"{pipeline_dir}/publishDir", '--quiet'])
            subprocess.run(['aws', 's3', 'rm', '--recursive', pipeline_dir])

if __name__ == "__main__":
    processor = RunNextflow()
    processor.batch_input = sys.argv[1] if len(sys.argv) > 1 else ''
    processor.pipeline = sys.argv[2] if len(sys.argv) > 2 else 'nanopore'
    processor.profile = sys.argv[3] if len(sys.argv) > 3 else 'codon'
    processor.root_dir = sys.argv[4] if len(sys.argv) > 4 else 'gs://prj-int-dev-covid19-nf-gls'
    processor.batch_index = sys.argv[5] if len(sys.argv) > 5 else '0'
    processor.snapshot_date = sys.argv[6] if len(sys.argv) > 6 else '2023-03-13'
    processor.test_submission = sys.argv[7] if len(sys.argv) > 7 else 'true'
    processor.study_accession = sys.argv[8] if len(sys.argv) > 8 else 'PRJEB45555'
    processor.dataset_name = sys.argv[9] if len(sys.argv) > 9 else 'sarscov2_metadata'
    processor.project_id = sys.argv[10] if len(sys.argv) > 10 else 'prj-int-dev-covid19-nf-gls'
    processor.process_batch()
