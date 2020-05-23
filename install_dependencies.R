## Dependencies

dependencies <- c("argparse", "RPostgreSQL", "tidyverse", "glue", "kableExtra",
                  "rmarkdown", "cowplot", "gtools", "scales", "formattable", "knitr")

args = commandArgs(trailingOnly=TRUE)

if(any(grepl(pattern = "show", args))){
  print(dependencies)
  quit()
}

for(dep in dependencies)
    install.packages(dep, repos='http://cran.us.r-project.org')
