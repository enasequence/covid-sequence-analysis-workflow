CREATE VIEW prj-int-dev-covid19-nf-gls.sarscov2_metadata.sra_processed
AS SELECT
     *
   FROM
     `prj-int-dev-covid19-nf-gls.sarscov2_metadata.sra_index` T1
   WHERE
     T1.run_accession IN (
     SELECT
       T2.run_ref
     FROM
       `prj-int-dev-covid19-nf-gls.sarscov2_metadata.analysis_archived` T2)
     OR T1.run_accession IN (
     SELECT
       T4.run_id
     FROM
       `prj-int-dev-covid19-nf-gls.sarscov2_metadata.submission_metadata` T4)
   ORDER BY
     T1.run_accession DESC