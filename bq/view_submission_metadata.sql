CREATE VIEW prj-int-dev-covid19-nf-gls.sarscov2_metadata.submission_metadata
AS SELECT
  T1.analysis_accession,
  T1.sample_accession,
  T1.submitted_bytes,
  T1.analysis_date,
  T2.run_accession,
  T2.instrument_platform AS platform,
  T2.instrument_model AS model,
  T2.first_public,
  T2.first_created,
  T2.country,
  T2.collection_date,
  (CASE
    WHEN T1.analysis_date IS NULL THEN NULL
    WHEN T1.analysis_date < '2022-03-22' THEN '2022-03-22'
    WHEN T1.analysis_date >= '2022-03-22' AND T1.analysis_date < '2022-04-12' THEN '2022-04-12'
    WHEN T1.analysis_date >= '2022-04-12' AND T1.analysis_date < '2022-05-23' THEN '2022-05-23'
    WHEN T1.analysis_date >= '2022-05-23' AND T1.analysis_date < '2022-06-27' THEN '2022-06-27'
    WHEN T1.analysis_date >= '2022-06-27' AND T1.analysis_date < '2022-07-25' THEN '2022-07-25'
    WHEN T1.analysis_date >= '2022-07-25' AND T1.analysis_date < '2022-08-22' THEN '2022-08-22'
    WHEN T1.analysis_date >= '2022-08-22' AND T1.analysis_date < '2022-09-26' THEN '2022-09-26'
    WHEN T1.analysis_date >= '2022-09-26' AND T1.analysis_date < '2022-10-24' THEN '2022-10-24'
    WHEN T1.analysis_date >= '2022-10-24' AND T1.analysis_date < '2022-11-21' THEN '2022-11-21'
    WHEN T1.analysis_date >= '2022-11-21' AND T1.analysis_date < '2022-12-19' THEN '2022-12-19'
    ELSE 'TBD'
  END
    ) AS snapshot_date,
  T3.time_submitted
FROM
  `prj-int-dev-covid19-nf-gls.sarscov2_metadata.analysis_archived` T1
LEFT JOIN
  `prj-int-dev-covid19-nf-gls.sarscov2_metadata.sra_index` T2
ON
  T2.run_accession = T1.run_ref
LEFT JOIN
  `prj-int-dev-covid19-nf-gls.sarscov2_metadata.submission_receipts` T3
ON
  T3.analysis_accession = T1.analysis_accession;
