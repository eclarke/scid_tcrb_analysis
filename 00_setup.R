
## To be sourced before main report generation

# Libraries ---------------------------------------------------------------
library(plyr)
library(dplyr)
library(ggplot2)
library(reshape2)
library(xtable)
library(RMySQL)
library(stringr)
library(data.table)
library(pander)
library(ggthemes)
library(gridExtra)
library(ggbeeswarm)
library(magrittr)
library(ineq)
library(eclthemes)
library(lymphclon)


# Global options ---------------------------------------------------------

## Set options for Knitr
knitr::opts_chunk$set(
  echo=FALSE,
  warning=FALSE,
  message=FALSE,
  cache=TRUE
)

## Set theme defaults
theme_set(theme_clarity(base_size = 14))

theme_clarity_bottomlegend <- theme_clarity() +
  theme(legend.position="bottom", legend.direction="horizontal")

## Report options
# Set to TRUE to use saved .RData version of data 
USE_CACHED_DATA <- FALSE
# Set to TRUE to only use US patients
US_PATIENTS_ONLY <- FALSE



# File paths --------------------------------------------------------------

base_fp <- "."
data_fp <- file.path("/Volumes/THING2/erik/101_TCRB")
saved_seq_data_fp <- file.path(base_fp, "SequenceData.RData")
source(file.path(base_fp, "lib/general_functions.R"))
source(file.path(base_fp, "lib/data_loading_functions.R"))
source(file.path(base_fp, "lib/population_diversity.R"))

# Username, password, and host are saved in ~/.my.cnf under group [bushman_mysql]
db_settings <- list(
  defaults = "~/.my.cnf",
  group = "bushman_mysql",
  db = "specimen_management",
  table = "GTSP"
)

# This directory should point to the folder with all the GTSP*.tsv files
seq_files <- file.path(data_fp, "tsv")
# This file should have the same number of rows as the number of files in seq_files
seq_stats_file <- file.path(data_fp, "sampleStats.all.tsv")
# This should contain the ng of DNA used for each sample
sample_dna_file <- file.path(data_fp, "sampleDNA.txt")



# Loading sequence data ------------------------------------------------------------

## Reading in sequence data
if (!USE_CACHED_DATA) {
  seqs  <- read_adaptive_tsv(list.files(seq_files, full.names = TRUE))
  save(seqs, file = saved_seq_data_fp)
} else {
  message("Loading saved sequence data...")
  load(file = saved_seq_data_fp)
}


# Load metadata -----------------------------------------------------------

# DB connection for sample metadata
con <- with(db_settings, {
  dbConnect(MySQL(), default.file=defaults, group=group, dbname=db)
})

mdata <- read_sample_metadata(list.files(seq_files), con, "GTSP")
stats <- read_sequence_stats_file(seq_stats_file)
dna.amts <- read_sample_dna_file(sample_dna_file)
metadata <- merge(merge(mdata, dna.amts, by=c("accn", "replicate")),
                  stats, by=c("accn", "replicate")) 

# Calculate the Unique/Total ratio and clean up factor ordering
metadata %<>%
  arrange(desc(Group), patient, timepoint, cell.type) %>%
  mutate(
    Unique.v.Total = Unique/Total,
    patient_at_timept=reorder_levels_by_uniq(patient_at_timept),
    patient_at_timept_ctype=reorder_levels_by_uniq(patient_at_timept_ctype))

## Write metadata to CSV for any outside manipulation
write.table(
  metadata, 
  file = file.path(data_fp, "metadata.tsv"),
  sep = "\t",
  row.names = FALSE)

# Cleanup
if (!dbDisconnect(con)) warning("Could not close db connection...") else rm(con)
rm(stats, dna.amts)


# (Optional) US patients only ----------------------------------------------

if (US_PATIENTS_ONLY) {
  metadata %<>% 
    filter(trial %in% c("Control") | (trial == "SCIDn2" & grepl("SCID", patient)))
}
