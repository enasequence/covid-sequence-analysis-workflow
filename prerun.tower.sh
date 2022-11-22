
pipeline='illumina'
profile='gls'
root_dir='gs://prj-int-dev-covid19-nf-gls-poc'
snapshot_date='2022-06-26' 
batch_size='5'
dataset_name='sarscov2_metadata'
project_id='prj-int-dev-covid19-nf-gls' 
table_name="${pipeline}_to_be_processed"
input_dir="./data/${snapshot_date}"; mkdir -p "${input_dir}"
offset=0
sql="SELECT count(*) AS total FROM ${project_id}.${dataset_name}.${table_name}"
row_count=$(bq --project_id="${project_id}" --format=csv query --use_legacy_sql=false "${sql}" | grep -v total)
echo "row_count is: $row_count"

sql="SELECT * FROM ${project_id}.${dataset_name}.${table_name} LIMIT ${batch_size} OFFSET ${offset}"
bq --project_id="${project_id}" --format=csv query --use_legacy_sql=false --max_rows="${batch_size}" "${sql}" \
| awk 'BEGIN{ FS=","; OFS="\t" }{$1=$1; print $0 }' > "${input_dir}/${table_name}.tsv"
gsutil -m cp "${input_dir}/${table_name}.tsv" "gs://${dataset_name}_poc/${table_name}.tsv"