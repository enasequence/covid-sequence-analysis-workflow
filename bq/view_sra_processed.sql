CREATE VIEW prj-int-dev-covid19-nf-gls.sarscov2_metadata.sra_processed
AS SELECT
  DISTINCT T1.run_accession,
  T1.instrument_platform,
  T1.first_public,
  T2.analysis_date
FROM
  `prj-int-dev-covid19-nf-gls.sarscov2_metadata.sra_index` T1,
  `prj-int-dev-covid19-nf-gls.sarscov2_metadata.analysis_archived` T2
WHERE
  -- T1.fastq_ftp IS NOT NULL AND
  -- (REGEXP_CONTAINS(T1.fastq_ftp, r'^ftp.sra.ebi.ac.uk.*.fastq.gz;ftp.sra.ebi.ac.uk.*_1.fastq.gz;ftp.sra.ebi.ac.uk.*_2.fastq.gz$') OR
    --   REGEXP_CONTAINS(T1.fastq_ftp, r'^ftp.sra.ebi.ac.uk.*_1.fastq.gz;ftp.sra.ebi.ac.uk.*_2.fastq.gz$') OR
    --   REGEXP_CONTAINS(T1.fastq_ftp, r'^ftp.sra.ebi.ac.uk.*_1.fastq.gz$')) AND
  T1.run_accession = T2.run_accession
ORDER BY
  T1.run_accession DESC