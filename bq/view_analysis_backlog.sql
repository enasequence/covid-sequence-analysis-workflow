CREATE VIEW prj-int-dev-covid19-nf-gls.sarscov2_metadata.analysis_backlog
AS WITH
    tmp_analysed AS (
        SELECT
            analysis_date,
            COUNT(DISTINCT run_accession) AS analysed_runs
        FROM
            `prj-int-dev-covid19-nf-gls.sarscov2_metadata.analysis_archived`
        GROUP BY
            analysis_date),
    tmp_public AS (
        SELECT
            first_public,
            COUNT(DISTINCT run_accession) AS public_runs
        FROM
            `prj-int-dev-covid19-nf-gls.sarscov2_metadata.sra_index`
        GROUP BY
            first_public)
SELECT
    T1.first_public AS dates,
    IFNULL(T1.public_runs, 0) AS public_runs,
    IFNULL(T2.analysed_runs, 0) AS analysed_runs,
    (IFNULL(T1.public_runs, 0) - IFNULL(T2.analysed_runs, 0)) AS analysis_backlog
FROM
    tmp_public T1
        FULL JOIN
    tmp_analysed T2
    ON
            T1.first_public = T2.analysis_date
ORDER BY
    dates asc