# GDA exam

This repo contains a script to automatically generate a report using as input a tsv file with annotated CNVS.
This input file needs as mandatory columns 
> "Sample ID" "Chr" "Start" "End" "CNV Value" "CNV Conf" "Comment"  

The script integrates information from 4 databases:
* DDG2P = CNV associated with developmental pathologies
* OMIM = genes implicated in deseases 
* CLINGEN = dosage sensitivity from ClinGen
* SYNDROME = CNVs associated with syndromes (DECIPHER)

### Basic Usage

First of all check the dependencies. 
To list them
`Rscript install_dependencies.R show `
To install all the dependecies
`Rscript install_dependencies.R`

Next create a DB on postgresql (ex on Linux) with user postgres and no password on our localhost and port 5432

`createdb CNV`

Then
`Rscript generate_report.R -db "CNV2" -dbh "localhost" -dbp 5432 -dbu "postgres" -dbpw ""  -f 20200408_standard_cnv_report.txt`

Will produce two files names output_summary_plots.pdf with some summary plots and output_CNV_report.html with the actual report.

The script will automatically recreate the 4 tables every time, if you want just to add the necessary tables set the argument `--do-not-force-annot`, on the other hand if you have everything (also the input file) already on the Database set `--do-not-create-df` and write the table name `-stn input_table_name_in_the_db`. (**remember tha column names need to be lowercase**)

The script filters out automatically any CNVs equal to 2. To deactivate this beahviour for the X chromosome set `--no-filter-diploid-X`.

For more information
`Rscript generate_report.R --help`

In the repo are already the output files from the exam input file.
