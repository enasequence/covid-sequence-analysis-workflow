CREATE VIEW prj-int-dev-covid19-nf-gls.datahub_metadata.all_analyses_storage
AS SELECT
  * EXCEPT(submitted_bytes),
  (
  SELECT
    SUM(PARSE_BIGNUMERIC(val)/1073741824)
  FROM
    UNNEST(SPLIT(submitted_bytes, ';')) AS val) AS storage_size
FROM
  `prj-int-dev-covid19-nf-gls.datahub_metadata.all_analyses`