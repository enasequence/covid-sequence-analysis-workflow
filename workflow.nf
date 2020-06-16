#!/usr/bin/env nextflow

params.RAW_DIR = ""
params.ANALYSIS_DIR = ""
params.HUMAN_IDX = ""
params.SARS2_IDX = ""
params.SARS2_FA = ""

process quality_control_pre {
    cpus 4
    memory '1 GB'

    input:
    path input_file from params.RAW_DIR
    path analysis_dir from params.ANALYSIS_DIR

    output:


    script:
    """
    fastqc -t ${task.cpus} -o ${analysis_dir}fastqc_pre/" -q ${input_file}
    """
}

process trim_reads {
    cpus 19
    memory '1 GB'

    input:
    path input_file from params.RAW_DIR

    output:

    script:
    """
    trimmomatic PE ${input_file}_1 ${input_file}_2 \
    ${input_file}_trim ${input_file}_trim \
    ${input_file}_trim_un ${input_file}_trim_un -summary \
    ${input_file}_trim_summary -threads ${task.cpus} SLIDINGWINDOW:5:30 MINLEN:50
    """
}

process quality_control_post {
    cpus 4
    memory '1 GB'

    input:
    path input_file from params.RAW_DIR
    path analysis_dir from params.ANALYSIS_DIR

    output:


    script:
    """
    fastqc -t ${task.cpus} -o ${analysis_dir}fastqc_post/" -q ${input_file}
    """
}

process align_reads_to_human_genome {
    cpus 19
    memory '20 GB'

    input:
    path human_idx from params.HUMAN_IDX

    output:

    script:
    """
    bowtie2 --very-sensitive-local -p ${task.cpus} -x ${human_idx} \
    --met-file $s -1 $f -2 $r -U $uf','$ur | samtools view -Sb -f 4 > $bam
    """
}

process convert_bam_to_fastq {
    cpus 4
    memory '4 GB'

    input:

    output:

    script:
    """
    samtools bam2fq -N -1 $f -2 $r -s $s $bam > $u
    """
}

process align_reads_to_sars2_genome {
    cpus 19
    memory '20 GB'

    input:

    output:

    script:
    """
    bowtie2 -p $v_threads --no-mixed --no-discordant --met-file $s \
    -x $v_sars2_idx -1 $f -2 $r | samtools view -bST $v_sars2_fa | \
    samtools sort | samtools view -h -F 4 -b > $bam
    samtools index $bam
    """
}

process remove_duplicates {
    cpus 19
    memory '20 GB'

    input:

    output:

    script:
    """
    picard MarkDuplicates I=$ibam O=$obam REMOVE_DUPLICATES=true \
    M=$ibam'_marked_dup_metrics.txt'
    """
}

process check_coverage {
    cpus 19
    memory '1 GB'

    input:

    output:

    script:
    """
    samtools mpileup -a -A -Q 30 -d 1000000 -f $v_sars2_fa $f > $o
    """
}

process make_small_file_with_coverage {
    cpus: 19
    memory: '1 GB'

    input:

    output:

    script:
    """
    cat $f | awk '{print $2,",",$3,",",$4}' > $o
    """
}

process generate_vcf_files {
    cpus: 19
    memory: '1 GB'

    input:

    output:

    script:
    """
    samtools index $f
    lofreq call-parallel --pp-threads ${task.cpus} -f $v_sars2_fa -o $o $f
    bgzip $f
    tabix $g
    bcftools stats $g > $v
    """
}

process create_consensus_sequence {
    cpus: 19
    memory: '1 GB'

    input:

    output:

    script:
    """
    bcftools filter -i "DP>50" $f   -o $o
    bgzip $o
    tabix $g
    bcftools filter -i "AF>0.5" $g > $go
    fa=${go/'.cfiltered_freq.vcf'/'.cons.fa'}
    gx=${go/'.vcf'/'.vcf.gz'}
    n='>'${go/'.cfiltered_freq.vcf'/''}
    bgzip -c $go > $gx
    bcftools index $gx
    bcftools consensus -f $v_sars2_fa $gx > $fa
    sed -i "1s/.*/$n/" $fa
    rm $g
    rm ${g/'.cfiltered.vcf.gz'/'.cfiltered.vcf.gz.tbi'}
    rm $go
    rm ${go/'cfiltered_freq.vcf'/'cfiltered_freq.vcf.gz.csi'}
    rm ${go/'cfiltered_freq.vcf'/'cfiltered_freq.vcf.gz'}
    """
}

process filter_snv {
    cpus: 19
    memory: '1 GB'

    input:

    output:

    script:
    """
    bcftools filter -i "DP>50" $f   -o $o
    bgzip $o
    tabix $g
    bcftools filter -i "AF>0.1" $g > $go
    """
}

process annotate_snps {
    cpus: 19
    memory: '1 GB'

    input:

    output:

    script:
    """
    cat $f | sed "s/^NC_045512.2/NC_045512/" > $m
    java -Xmx4g -jar ~/tools/snpEff/snpEff.jar -v -s $s sars.cov.2 $m > $o
    """
}