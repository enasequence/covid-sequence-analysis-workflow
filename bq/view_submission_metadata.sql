CREATE VIEW prj-int-dev-covid19-nf-gls.sarscov2_metadata.submission_metadata
AS SELECT
  DISTINCT T2.run_accession AS run_id,
  T2.instrument_platform AS platform,
  T2.instrument_model AS model,
  T2.first_public,
  T2.first_created,
  T2.country,
  T2.collection_date,
  T1.time_submitted,
  T1.snapshot_date,
  -- IFNULL(T1.snapshot_date, '2022-06-27') AS snapshot_date,
  T3.analysis_accession,
  T3.sample_accession,
  T3.submitted_bytes,
  T3.analysis_date
FROM (`prj-int-dev-covid19-nf-gls.sarscov2_metadata.submission_receipts` T1
  FULL JOIN
    `prj-int-dev-covid19-nf-gls.sarscov2_metadata.analysis_archived` T3
  ON
    REGEXP_EXTRACT(T1.file_submitted, r'^([A-Z0-9]+)_[a-z]+.[a-z]+.gz$') = T3.run_ref)
INNER JOIN
  `prj-int-dev-covid19-nf-gls.sarscov2_metadata.sra_index` T2
ON
  T2.run_accession = T3.run_ref

  -- AS SELECT
  --   DISTINCT T2.run_accession AS run_id,
  --   T2.instrument_platform AS platform,
  --   T2.instrument_model AS model,
  --   T2.first_public,
  --   T2.first_created,
  --   T2.country,
  --   T2.collection_date,
  --   T1.time_submitted,
  --   T1.snapshot_date
  -- FROM
  --   `prj-int-dev-covid19-nf-gls.sarscov2_metadata.submission_receipts` T1,
  --   `prj-int-dev-covid19-nf-gls.sarscov2_metadata.sra_index` T2
  -- WHERE
  --   T2.run_accession = REGEXP_EXTRACT(T1.file_submitted, r'^([A-Z0-9]+)_[a-z]+.[a-z]+.gz$')
