CREATE VIEW prj-int-dev-covid19-nf-gls.sarscov2_metadata.submission_metadata
AS SELECT
    T1.analysis_accession,
    T1.sample_accession,
    T1.submitted_bytes,
    T1.analysis_date,
    T2.run_accession,
    T2.instrument_platform AS platform,
    T2.instrument_model AS MODEL,
    T2.first_public,
    T2.first_created,
    T2.country,
    T2.collection_date,
    (CASE
         WHEN T1.analysis_date IS NULL THEN NULL
         WHEN T1.analysis_date < '2021-03-22' THEN '2021-03-22'
         WHEN T1.analysis_date >= '2021-03-22' AND T1.analysis_date < '2021-04-19' THEN '2021-04-19'
         WHEN T1.analysis_date >= '2021-04-19'
             AND T1.analysis_date < '2021-05-17' THEN '2021-05-17'
         WHEN T1.analysis_date >= '2021-05-17' AND T1.analysis_date < '2021-06-21' THEN '2021-06-21'
         WHEN T1.analysis_date >= '2021-06-21'
             AND T1.analysis_date < '2021-07-19' THEN '2021-07-19'
         WHEN T1.analysis_date >= '2021-07-19' AND T1.analysis_date < '2021-08-23' THEN '2021-08-23'
         WHEN T1.analysis_date >= '2021-08-23'
             AND T1.analysis_date < '2021-09-20' THEN '2021-09-20'
         WHEN T1.analysis_date >= '2021-09-20' AND T1.analysis_date < '2021-10-18' THEN '2021-10-18'
         WHEN T1.analysis_date >= '2021-10-18'
             AND T1.analysis_date < '2021-11-22' THEN '2021-11-22'
         WHEN T1.analysis_date >= '2021-11-22' AND T1.analysis_date < '2021-12-20' THEN '2021-12-20'
         WHEN T1.analysis_date >= '2021-12-20'
             AND T1.analysis_date < '2022-01-17' THEN '2022-01-17'
         WHEN T1.analysis_date >= '2022-01-17' AND T1.analysis_date < '2022-02-21' THEN '2022-02-21'
         WHEN T1.analysis_date >= '2022-02-21'
             AND T1.analysis_date < '2022-03-22' THEN '2022-03-22'
         WHEN T1.analysis_date >= '2022-03-22' AND T1.analysis_date < '2022-04-12' THEN '2022-04-12'
         WHEN T1.analysis_date >= '2022-04-12'
             AND T1.analysis_date < '2022-05-23' THEN '2022-05-23'
         WHEN T1.analysis_date >= '2022-05-23' AND T1.analysis_date < '2022-06-27' THEN '2022-06-27'
         WHEN T1.analysis_date >= '2022-06-27'
             AND T1.analysis_date < '2022-07-25' THEN '2022-07-25'
         WHEN T1.analysis_date >= '2022-07-25' AND T1.analysis_date < '2022-08-22' THEN '2022-08-22'
         WHEN T1.analysis_date >= '2022-08-22'
             AND T1.analysis_date < '2022-09-26' THEN '2022-09-26'
         WHEN T1.analysis_date >= '2022-09-26' AND T1.analysis_date < '2022-10-24' THEN '2022-10-24'
         WHEN T1.analysis_date >= '2022-10-24'
             AND T1.analysis_date < '2022-11-21' THEN '2022-11-21'
         WHEN T1.analysis_date >= '2022-11-21' AND T1.analysis_date < '2022-12-19' THEN '2022-12-19'
         WHEN T1.analysis_date >= '2022-12-19'
             AND T1.analysis_date < '2023-01-23' THEN '2023-01-23'
         WHEN T1.analysis_date >= '2023-01-23' AND T1.analysis_date < '2023-02-20' THEN '2023-02-20'
         WHEN T1.analysis_date >= '2023-02-20'
             AND T1.analysis_date < '2023-03-20' THEN '2023-03-20'
         WHEN T1.analysis_date >= '2023-03-20' AND T1.analysis_date < '2023-04-17' THEN '2023-04-17'
         WHEN T1.analysis_date >= '2023-04-17'
             AND T1.analysis_date < '2023-05-22' THEN '2023-05-22'
         WHEN T1.analysis_date >= '2023-05-22' AND T1.analysis_date < '2023-06-19' THEN '2023-06-19'
         WHEN T1.analysis_date >= '2023-06-19'
             AND T1.analysis_date < '2023-07-24' THEN '2023-07-24'
         WHEN T1.analysis_date >= '2023-07-24' AND T1.analysis_date < '2023-08-21' THEN '2023-08-21'
         WHEN T1.analysis_date >= '2023-08-21'
             AND T1.analysis_date < '2023-09-18' THEN '2023-09-18'
         WHEN T1.analysis_date >= '2023-09-18' AND T1.analysis_date < '2023-10-23' THEN '2023-10-23'
         WHEN T1.analysis_date >= '2023-10-23'
             AND T1.analysis_date < '2023-11-20' THEN '2023-11-20'
         WHEN T1.analysis_date >= '2023-11-20' AND T1.analysis_date < '2023-12-18' THEN '2023-12-18'
         ELSE
             'TBD'
        END
        ) AS snapshot_date,
    T3.time_submitted
FROM
    `prj-int-dev-covid19-nf-gls.sarscov2_metadata.analysis_archived` T1
        LEFT JOIN
    `prj-int-dev-covid19-nf-gls.sarscov2_metadata.sra_index` T2
    ON
            T2.run_accession = T1.run_accession
        LEFT JOIN
    `prj-int-dev-covid19-nf-gls.sarscov2_metadata.submission_receipts` T3
    ON
            T3.analysis_accession = T1.analysis_accession;