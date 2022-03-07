CREATE VIEW prj-int-dev-covid19-nf-gls.sarscov2_metadata.nanopore_to_be_processed
AS SELECT
  *
FROM
  `prj-int-dev-covid19-nf-gls.sarscov2_metadata.sra_index` T1
WHERE
  T1.instrument_platform = 'OXFORD_NANOPORE'
  AND T1.fastq_ftp IS NOT NULL
  AND REGEXP_CONTAINS(T1.fastq_ftp, r'^ftp.sra.ebi.ac.uk.*_1.fastq.gz$')
  AND T1.run_accession NOT IN (
  SELECT
    T2.run_ref
  FROM
    `prj-int-dev-covid19-nf-gls.sarscov2_metadata.analysis_archived` T2)
  AND T1.run_accession NOT IN (
  SELECT
    T3.run_accession
  FROM
    `prj-int-dev-covid19-nf-gls.sarscov2_metadata.sra_processing` T3)
  ORDER BY run_accession DESC
