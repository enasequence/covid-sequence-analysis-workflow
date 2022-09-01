#!/bin/bash

pipeline_dir=${1}

echo $pipeline_dir

echo "${pipeline_dir}/workDir"
rm -R "${pipeline_dir}/workDir" &

echo "${pipeline_dir}/storeDir"
rm -R "${pipeline_dir}/storeDir" &

echo "${pipeline_dir}/publishDir"
rm -R "${pipeline_dir}/publishDir" &
