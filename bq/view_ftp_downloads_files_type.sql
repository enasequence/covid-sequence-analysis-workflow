CREATE VIEW prj-int-dev-covid19-nf-gls.datahub_metadata.ftp_downloads_files_type
AS SELECT
  * EXCEPT(file_size__Descending),
  (file_size__Descending / 1073741824) AS file_size,
  REGEXP_EXTRACT(Resource, '^/vol1/([a-zA-Z0-9]+)/') AS type
FROM
  `prj-int-dev-covid19-nf-gls.datahub_metadata.ftp_downloads_files`