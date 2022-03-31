CREATE VIEW prj-int-dev-covid19-nf-gls.sarscov2_metadata.submission_metadata
AS SELECT DISTINCT
     T2.run_accession AS run_id,
     T2.instrument_platform AS platform,
     T2.instrument_model AS model,
     T2.first_public,
     T2.first_created,
     T2.country,
     T2.collection_date,
     T1.snapshot_date
   FROM
     `prj-int-dev-covid19-nf-gls.sarscov2_metadata.submission_receipts` T1,
     `prj-int-dev-covid19-nf-gls.sarscov2_metadata.sra_index` T2
   WHERE
     T2.run_accession = REGEXP_EXTRACT(T1.file_submitted, r'^([A-Z0-9]+)_[a-z]+.[a-z]+.gz$')
