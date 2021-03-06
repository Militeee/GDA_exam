require(argparse, quietly = T)

## SQL queries

alter_range = "alter table {table}
               add column range int8range"

range_update = "update {table} set \"range\"= subquery.range
                from (
                  select (\"chr\", \"start\", \"end\") as tp,
                  int8range(cast(\"start\" as INT),
                  cast(\"end\" as INT), '[]') as range
                  from {table}
                  ) as subquery
                where (\"chr\", \"start\", \"end\") = subquery.tp"

overlap_syndrome = "select \"sample id\", \"sindrome\", st.\"chr\" , st.\"start\", st.\"end\", s.\"start\" as syn_start, s.\"end\" as syn_end,
                    ((upper(st.range * s.range) - lower(st.range * s.range))::numeric / (upper(st.range) - lower(st.range)) * 100) as \"perc overlap\",
                    \"classe\", \"cnv value\", \"comment\"
                    from syndrome as s inner join {sample_table} as st on s.chr = st.chr
                    where st.range && s.range and st.\"cnv conf\" >= {qcoff}  and ({filter})
                  "

gene_annotation = "
                   create table gene_snp_annot as (
                	 select distinct \"sample id\",st.\"chr\", st.\"start\", st.\"end\",gn.\"gene name\" ,gn.\"start\" as gene_start, gn.\"end\" as gene_end, st.\"cnv value\",
                	 ((upper(st.range * gn.range) - lower(st.range * gn.range))::numeric / (upper(gn.range) - lower(gn.range)) * 100) as \"perc overlap\",
                	 st.\"comment\"
                	 from {sample_table} as st inner join gene_annotations as gn on st.\"chr\"  = gn.\"chr\"
                	 where st.range && gn.range and st.\"cnv conf\" >= {qcoff}  and ({filter})
                	)"

overlap_ddg2p = " select dd.\"gene_symbol\", st.\"chr\",st.\"sample id\", st.\"perc overlap\", dd.\"allelic_requirement\",dd.\"ddd.category\", st.\"cnv value\", dd.\"organ.specificity.list\", st.\"comment\"
                  from  gene_snp_annot as st inner join ddg2p as dd on st.\"gene name\" = dd.\"gene_symbol\"
                "

overlap_clingen = "
                  select cg.\"symbol\", st.\"chr\",st.\"sample id\", st.\"perc overlap\", cg.\"haploinsufficiency\",cg.\"triplosensitivity\", st.\"cnv value\", cg.\"online_report\", st.\"comment\"
                  from  gene_snp_annot as st inner join clingen_dosage as cg on st.\"gene name\" = cg.\"symbol\"
                  "
overlap_OMIM = "
                  select ph.\"symbol\", st.\"chr\",st.\"start\", st.\"end\", st.\"gene_start\", st.\"gene_end\",st.\"sample id\", st.\"perc overlap\", ph.\"phenotype\", st.\"cnv value\", st.\"comment\"
                  from  gene_snp_annot as st inner join morbidmap_ok as ph on st.\"gene name\" = ph.\"symbol\"
               "

get_sample = "
              select *
              from {sample_table}
            "

get_sample_filt = "
                    select *
                    from {sample_table} as st
                    where st.\"cnv conf\" >= {qcoff}  and ({filter})
                  "

## Parsing arguments

parser <- ArgumentParser()

parser$add_argument("-db", "--database", type= "character", help = "Name of the Postgres database where you want to insert/read annot tables")
parser$add_argument("-dbh", "--database-host",type= "character")
parser$add_argument("-dbp", "--database-port",type= "character")
parser$add_argument("-dbu", "--database-user",type= "character")
parser$add_argument("-dbpw", "--database-password", type= "character")
parser$add_argument("--quality-cutoff", type="integer", default=30)
parser$add_argument("--annotation-dir", type= "character", default ="annot", help = "Directory of annotation table")
parser$add_argument("--delim", type= "character", default =",", help = "Separator for annotation tables")
parser$add_argument("-f","--filename", type= "character", help= "Input file, needs to be tab separated")
parser$add_argument("-o","--out-prefix", type= "character", default = "output", help = "Output prefix to be appendend to the two output files")
parser$add_argument("--do-not-force-annot", action="store_true", default = FALSE, help = "If selected the script will not insert all the tables but only the one different from those present in the Postgres DB")
parser$add_argument("-stn","--sample-table-name", type = "character", default = "sample_table", help = "Name of sample table in the df")
parser$add_argument("--do-not-create-df", action="store_true", default = FALSE , help = "Don't write or create any table on the database, just read")
parser$add_argument("--no-filter-diploid-X", action="store_true", default = FALSE , help = "Do not filter samples with CNV value = 2 on the X chromosome")
parser$add_argument("-gv","--genome-version", type= "character", default ="hg19", help = "Version of the ref [hg38 or hg19] genome to be used for gene position")
parser$add_argument("--custom-filter", type= "character", default =NULL, help = "A custom filter for the sample file")
parser$add_argument("--save-tables", action= "store_true", default =FALSE, help = "Save an RData object with the results of the queries")
parser$add_argument("--verbose", action= "store_true", default =FALSE, help ="Show output and warnings of functions")

args <- parser$parse_args()


suppressPackageStartupMessages({
require(RPostgreSQL, quietly = T)
require(tidyverse, quietly = T)
require(glue, quietly = T)
require(kableExtra, quietly = T)
require(rmarkdown, quietly = T)
require(cowplot, quietly = T)
require(gtools, quietly = T)
require(scales, quietly = T)
require(formattable, quietly = T)})

## Connecting to the db
if(!args$verbose){
  options(warn=-1)
  col <- cols()
} else{
  col <- NULL
}




cat("Connecting to DB")
cat("\n")

con <- dbConnect(RPostgres::Postgres(), dbname = args$database, host=args$database_host, port=args$database_port, user=args$database_user, password= args$database_password)


## Make sure to not reload all the tables if the database when not neaded

files <- dir(args$annotation_dir, full.names = T)

if(!(args$genome_version %in% c("hg19", "hg38"))) stop("Genome version should be hg19 or hg38")

files_no_hg <- grep(files, pattern = "hg", value = T, invert = T)
files_hg <- grep(files, pattern = args$genome_version, value = T)
files <- c(files_hg, files_no_hg)
nms <- sub("(.*/)","",files)
nms <- sub("\\..*","",nms)
nms <- sub(paste0("_",args$genome_version),"",nms)

tables <- dbListTables(con)

tables1 <- setdiff(nms, tables)


if(!args$do_not_create_df & (length(tables1) != 0 | !args$do_not_force_annot)){

  cat("Reading and inserting annotation tables")
  cat("\n")

  # create tables for the annotations
  for(i in seq_along(files)){

    annot_df <- read_delim(files[i], delim = args$delim,col_types = col)
    colnames(annot_df) <- tolower(colnames(annot_df))
    if(args$verbose)
      cat(paste0("Reading ", nms[i], " \n"))
    dbWriteTable(con, nms[i], annot_df, overwrite = T)


    # add range if the table has the right columns
    if(all(c("start", "end", "chr")  %in% colnames(annot_df)) ){
      dbExecute(con, glue(alter_range, table = nms[i]))
      dbExecute(con, glue(range_update, table = nms[i]))
    }
  }


}

if(is.null(args$custom_filter)){
  if(args$no_filter_diploid_X){

      filter = "st.\"cnv value\" <> 2  or st.chr = 'X'"

    } else {

      filter = " st.\"cnv value\" <> 2 "
    }
} else {
  filter = args$custom_filter
}
quality_cutoff <- paste(args$quality_cutoff)

## Processing sample file


cat("Processing sample table")
cat("\n")

if(!args$do_not_create_df){

  sample_file <- read_delim(args$filename, delim = "\t", na = c("","NA"), col_types = col)
  colnames(sample_file) <- tolower(colnames(sample_file))

  dbWriteTable(con, args$sample_table_name, sample_file, overwrite = T)

  k_ <- dbExecute(con, glue(alter_range, table = args$sample_table_name))
  k_ <- dbExecute(con, glue(range_update, table = args$sample_table_name))
} else {
  sample_file <- dbGetQuery(con, glue(get_sample, sample_table = args$sample_table_name))
}

## Calculating overlappings

cat("Calculating overlaps")
cat("\n")

syndrome_overlaps <- dbGetQuery(con, glue(overlap_syndrome, sample_table = args$sample_table_name, filter = filter, qcoff = quality_cutoff))

if(!args$do_not_create_df){
  k_ <- dbExecute(con, "drop table if exists gene_snp_annot")
  k_ <- dbExecute(con, glue(gene_annotation,sample_table = args$sample_table_name, qcoff = quality_cutoff, filter = filter))
}


annotated_ddg2p <- dbGetQuery(con, glue(overlap_ddg2p))

annotated_clingen <- dbGetQuery(con, glue(overlap_clingen))

annotated_morbidmap <- dbGetQuery(con, glue(overlap_OMIM))

sample_file_filtered <- dbGetQuery(con, glue(get_sample_filt,
                                             sample_table = args$sample_table_name,
                                             qcoff = quality_cutoff, filter = filter))

## Save datasets temporarly for rendering the report

cat("Saving temporary file")
cat("\n")

save(sample_file, sample_file_filtered, syndrome_overlaps, annotated_ddg2p,
      annotated_clingen, annotated_morbidmap,file = "tables.RData")

## Prepare datasets for the summary plot

cat("Plotting")
cat("\n")

# DF for plotting the percentage of the different CNVs

plt_data1 <- sample_file_filtered %>% group_by(`cnv value`) %>%
              summarize(counts = n()) %>% arrange(counts) %>% mutate(freq = counts/sum(counts), y_pos = cumsum(freq) - 0.4*freq)

# DF for plotting the percentage of amplifications and deletions

plt_data2 <- sample_file_filtered %>% mutate(del_dup = if_else(`cnv value` > 2, "ampl", "del")) %>%
                  mutate(del_dup = if_else(`cnv value` == 2, "norm", del_dup)) %>% group_by(del_dup) %>%
                    summarize(counts = n()) %>% arrange(counts) %>%
                     mutate(freq = counts/sum(counts), y_pos = cumsum(freq) - 0.5*cumsum(freq))

# DF for plotting the chrs prevalence of CNVs

plt_data3 <- sample_file_filtered %>% group_by(chr) %>% summarize(counts = n()) %>%  mutate(chr = factor(.$chr, levels = mixedsort(.$chr)))




# DF for plotting the overlap sample ids and datasets

plt_data4 <- data.frame(syndrome = length(intersect(sample_file_filtered$`sample id`, syndrome_overlaps$`sample id`)),
                        ddg2p = length(intersect(sample_file_filtered$`sample id`, annotated_ddg2p$`sample id`)),
                        omim = length(intersect(sample_file_filtered$`sample id`, annotated_morbidmap$`sample id`)),
                        clingen = length(intersect(sample_file_filtered$`sample id`, annotated_clingen$`sample id`)))

plt_data4 <- suppressMessages(plt_data4 %>% reshape2::melt())




## Generate summary plots


pie_1 <- ggplot(data = plt_data1, aes(x = "", y = freq, fill = paste(`cnv value`))) + geom_bar(width = 1, stat = "identity") +
          coord_polar("y", start=0)+ ggtitle("Perc of CNVs in the table") + scale_fill_discrete( "CNV value") + ylab("") + xlab("") +
            theme_minimal() + theme(axis.text.x = element_blank())

pie_2 <- ggplot(data = plt_data2, aes(x = "", y = freq, fill = paste(del_dup))) + geom_bar(width = 1, stat = "identity") +
          coord_polar("y", start=0)+ ggtitle("Perc of amp/del in the table") + scale_fill_discrete( "CNV value") + ylab("") + xlab("") +
            theme_minimal() + theme(axis.text.x = element_blank())


hist_1 <- ggplot(data = plt_data3, aes(x = chr, y = counts/sum(counts) * 100, fill = chr)) +
            geom_bar(stat = "identity", show.legend = FALSE)+ ggtitle("Distribution of CNVs in the chromosome") + ylab("") +
              theme_minimal()

hist_2 <- ggplot(data = plt_data4, aes(x = variable, y = value, fill = variable)) +
            geom_bar(stat = "identity", show.legend = FALSE) + ggtitle("Number of annotated samples for db") + ylab("") +
              theme_minimal() + theme(axis.text.x = element_text(angle = 90))


## Save summary plots

cowplot::plot_grid(
  pie_1,
  pie_2,
  hist_1,
  hist_2,
  nrow = 2, ncol = 2, align = "h"
) %>% ggsave(filename = paste0(args$out_prefix,"_summary_plots.pdf"), device = "pdf", width = 10, height = 10)

## Render the report

cat("Generating report")
cat("\n")

rmarkdown::render("print_report.Rmd", output_file = paste0(args$out_prefix,"_CNV_report.html"), quiet = T)


cat("Remove temporary files")
cat("\n")

## Remove the temporary datasets
if(!args$save_tables)
  system("rm tables.RData")

cat("BYE!")
cat("\n")

quit()
