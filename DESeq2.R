#
# +-----------------------------------------------------------------------------------------+
# +-----------------------------------------------------------------------------------------+
# RNA-Seq Workflow by @furkanmtorun 
    # E-mail: furkanmtorun@gmail.com  
    # GitHub: https://github.com/furkanmtorun 
    # Google Scholar: https://scholar.google.com/citations?user=d5ZyOZ4AAAAJ 
    # Personal Website: https://furkanmtorun.github.io/    
# +-----------------------------------------------------------------------------------------+
# +-----------------------------------------------------------------------------------------+
#


# +--------------------------------------------------+
# Install and import required libraries & packages
# +--------------------------------------------------+
install.packages("dplyr")
BiocManager::install("DESeq2")
BiocManager::install("gage")
BiocManager::install("pathview")
BiocManager::install("gageData")
library(dplyr)
library(DESeq2)
library(gage)
library(pathview)
library(gageData)


# +--------------------------------------------------+
# Argument parsing
# +--------------------------------------------------+

options(echo=FALSE)
args<-commandArgs(TRUE)

if (length(args)==0) {
    print("At least one argument must be supplied!")
} else if (length(args)==1) {
    print("You have only one argument!")
} else {
    print("You have multiple arguments: ")
    print(args)
}


# +--------------------------------------------------+
# Set the environments and files
# +--------------------------------------------------+
working_dir <- "argumentWorkingDirectory"
setwd(working_dir)






