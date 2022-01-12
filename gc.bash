#!/bin/bash

ILLUMINA_N=$1

if [[ $ILLUMINA_N == -1 ]]
then
    FOLDERS=$(gsutil ls gs://prj-int-dev-covid19-nf-gls/2021-12-13 | grep illumina_)
else
    FOLDERS=$(gsutil ls gs://prj-int-dev-covid19-nf-gls/2021-12-13 | grep illumina_$ILLUMINA_N/)
fi

echo -e "Launch for: \n$FOLDERS\n"

for folder in $FOLDERS
do
if [[ $folder == *"illumina"* ]]; then
    echo -e "Current bucket is: $folder\n"
    SUBFOLDERS=$(gsutil ls ${folder::-1}/results)
    for subfolder in $SUBFOLDERS:
    do
        if [[ $subfolder != *".tar.gz"* ]]; then

            RUN_ID_FOLDER=$(echo "${subfolder::-1}" | grep -oP "SRR.*")
            ILLUMINA_ID=$(echo "$subfolder" | grep -oP "(illumina_[[:digit:]]+)")
            OUTPUT_LOCAL_FOLDER="/home/mansurova/GC_FILES/$RUN_ID_FOLDER"

            FILES_1=$(gsutil ls $subfolder | grep -E ".*(annot|bam|coverage|vcf).*" | grep -vE ".*_filtered.*")

            mkdir -p $OUTPUT_LOCAL_FOLDER
            gsutil cp $FILES_1 $OUTPUT_LOCAL_FOLDER
            cd /home/mansurova/GC_FILES/ && tar -czvf $RUN_ID_FOLDER.tar.gz $RUN_ID_FOLDER
            rm -rf $OUTPUT_LOCAL_FOLDER

            OUTPUT_VM_FILE="gs://2022-01-25-upd-bucket/$ILLUMINA_ID/results/$RUN_ID_FOLDER.tar.gz"
            gsutil -m -h "Content-Type:application/octet-stream" cp /home/mansurova/GC_FILES/$RUN_ID_FOLDER.tar.gz $OUTPUT_VM_FILE
            rm -f /home/mansurova/GC_FILES/$RUN_ID_FOLDER.tar.gz

            FILES_2=$(gsutil ls $subfolder | grep -E ".*(_consensus|_filtered).*")

            mkdir -p $OUTPUT_LOCAL_FOLDER
            gsutil cp $FILES_2 $OUTPUT_LOCAL_FOLDER

            OUTPUT_VM_FOLDER="gs://2022-01-25-upd-bucket/$ILLUMINA_ID/results/"
            gsutil -m -h "Content-Type:application/octet-stream" cp -r $OUTPUT_LOCAL_FOLDER $OUTPUT_VM_FOLDER
            rm -rf $OUTPUT_LOCAL_FOLDER
        fi
    done
fi
done
