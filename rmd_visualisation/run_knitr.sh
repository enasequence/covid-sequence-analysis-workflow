R -e "rmarkdown::render('visualise.Rmd', params = list(
stat_folder = './output/data/',
gff_folder = './data/ref/sars2/NC_045512.2.gff3'
))"