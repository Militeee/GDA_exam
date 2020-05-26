# GDA exam

This repo contains a script to automatically generate a report using as input a tsv file with annotated CNVs.
This input file needs as mandatory columns 
> "Sample ID" "Chr" "Start" "End" "CNV Value" "CNV Conf" "Comment"  

The script integrates information from 4 databases:
* DDG2P = CNV associated with developmental pathologies
* OMIM = genes implicated in deseases 
* CLINGEN = dosage sensitivity from ClinGen
* SYNDROME = CNVs associated with syndromes (DECIPHER)

Annotation table for genes comes from BioMart GRCh37 for hg19 and GRCh38 for hg38

Brief note on overlapping:

With syndrome the overlaps are calculated as the % of the input CNV covered by the syndrome region

With genes the overlaps are calculated as the % of the gene covered by the input CNV

### Basic Usage

Clone the repo

First of all check the dependencies. 
To list them

`Rscript install_dependencies.R show `

To install all the dependecies

`Rscript install_dependencies.R`

Next create a DB on postgresql (ex on Linux) with user postgres and no password on our localhost and port 5432

`createdb CNV`

Then

`Rscript generate_report.R -db "CNV" -dbh "localhost" -dbp 5432 -dbu "postgres" -dbpw ""  -f 20200408_standard_cnv_report.txt`

Will produce two files names output_summary_plots.pdf with some summary plots and output_CNV_report.html with the actual report.

The script will automatically recreate the 4 tables every time, if you want just to add the necessary tables set the argument `--do-not-force-annot`, on the other hand if you have everything (also the input file) already on the database set `--do-not-create-df` and write the table name `-stn input_table_name_in_the_db`. (**remember that column names of the tables need to be lowercase**)

The script filters out automatically any CNVs equal to 2. To deactivate this beahviour for the X chromosome set `--no-filter-diploid-X`. If you want more complex filters (ex on non mandatory columns)  you can write a custom filter after `--custom-filter`. **remember to escape ' with \' **

For more information
`Rscript generate_report.R --help`


In the repo are already the output files from the exam input file.
