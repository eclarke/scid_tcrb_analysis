# 0 - Data Loading Functions --------------------------------------------------

#' Reads a list of Adaptive .tsv files quickly and combines them, filling in
#' columns that are not shared with NAs. Also adds accn and replicate fields.
read_adaptive_tsv <- function(files) {
  message("Reading ", length(files), " files...")
  
  contents <- lapply(files, function(fp) {
    fname <- basename(fp)
    # data.table::fread is orders of magnitude faster than read.table
    d <- suppressWarnings(data.table::fread(fp, na.strings=c("(undefined)")))
    d$accn <- .split_fname(fname)$accn
    d$replicate <- .split_fname(fname)$replicate
    data.table::setnames(d, "frequencyCount (%)", "frequency")
    message("Read ", fname)
    return(d)
  })
  contents <- dplyr::tbl_dt(data.table::rbindlist(contents, use.names = TRUE, fill = TRUE))
}


#' Retrieves metadata from the sample database 
read_sample_metadata <- function(files, con, sql_table) {
  message("Reading sample information from database...")
  metadata <- lapply(files, function(fp) {
    sample <- .split_fname(basename(fp))
    query <- sprintf("select SamplePatientCode, CellType, Timepoint, Trial from %s where SpecimenAccNum like '%s' or SamplePatientCode like '%s';",
                     sql_table, sample$accn, sample$accn)
    m <- DBI::dbGetQuery(con, query)
    with(m, {
      return(data.frame(accn=sample$accn, 
                        patient=SamplePatientCode, 
                        cell.type=CellType,
                        timepoint=Timepoint, 
                        trial=Trial,
                        replicate=sample$replicate))
    })
  })
  metadata <- suppressWarnings(dplyr::tbl_df(dplyr::rbind_all(metadata)))
  metadata$group <- "patient"
  within(metadata, {
    # Convert NA strings to NA values
    timepoint[timepoint == "NA"] <- NA
    # A row missing a timepoint is a control
    trial[is.na(timepoint)] <- "Control"
    # Specify what kind of control (adult/child)
    group[trial == "Control" & grepl("ND", patient)] <- "ctrl_adult"
    group[trial == "Control" & group != "ctrl_adult"] <- "ctrl_child"
    
    # Re-add healthy baby ages
    timepoint[patient == "060"] <- "m21"
    timepoint[patient == "062"] <- "m43"
    timepoint[patient == "063"] <- "m37"
    timepoint[patient == "064"] <- "m30"
    timepoint[patient == "065"] <- "m35"
    
    # More informative group notation, for use in display
    Group <- ifelse(group=="patient", trial, group)
    # Standardize timepoints to help convert into continuous nmonths variable
    timepoint <- plyr::revalue(timepoint, c(d180="m6", d365="m12"))
    timepoint <- stringr::str_replace(timepoint, "_", ".")
    nmonths <- as.numeric(stringr::str_replace(timepoint, "m", ""))
    # More columns for display/grouping
    #  Patient @ timepoint
    patient_at_timept = sprintf("%s @ %im", patient, nmonths)
    patient_at_timept[is.na(nmonths)] <- patient[is.na(nmonths)]
    #  Patient @ timepoint (cell type )
    patient_at_timept_ctype = sprintf("%s @ %im (%s)", patient, nmonths, cell.type)
    patient_at_timept_ctype[is.na(nmonths)] <- sprintf(
      "%s (%s)", patient[is.na(nmonths)], cell.type[is.na(nmonths)])
  })
}

find_additional_samples <- function(files, con, sql_table) {
  message("Reading sample information from database...")
  metadata <- lapply(files, function(fp) {
    sample <- .split_fname(basename(fp))
    query <- sprintf("select * from %s where SpecimenAccNum like '%s' or SamplePatientCode like '%s';",
                     sql_table, sample$accn, sample$accn)
    m <- DBI::dbGetQuery(con, query)
    return(data.frame(m))
  })
}


#' Reads the sample summary file from Adaptive
read_sequence_stats_file <- function(stats.file) {
  message("Reading sequence statistics file...")
  s <- dplyr::tbl_df(read.table(stats.file, header = TRUE, sep="\t"))
  accn.repl <- plyr::adply(as.character(s$Sample.Name), 1, 
                     function(i) data.frame(.split_fname(i)))[2:3]
  cbind(accn.repl, s[, -1])
}

read_sample_dna_file <- function(dna.file) {
  message("Reading sample DNA amounts file...")
  s <- dplyr::tbl_df(read.table(dna.file, header=TRUE, sep="\t"))
  accn.repl <- plyr::adply(as.character(s$Sample.Name), 1, 
                           function(i) data.frame(.split_fname(i)))[2:3]
  cbind(accn.repl, s[, -1])
}


#' Splits a [GTSP####[a-z].tsv] filename into its accession and replicate,
#' and converts character replicate to numerical (a=1, b=2, etc)
.split_fname <- function(fname) {
  # pretend it's a filename if it's missing an extension
  if (!stringr::str_detect(fname, "\\.tsv$")) fname <- paste(c(fname, ".tsv"), collapse="")
  accn = stringr::str_extract(fname, stringr::perl("[A-Z]*\\d*(?=[a-z]?.tsv)"))
  repl = stringr::str_extract(fname, stringr::perl("(?<=\\d)[a-z](?=\\.tsv)"))
  # if we can't find a replicate letter, assume it's the first replicate
  if (is.na(repl)) repl <- 1 else repl <- which(repl == letters)
  return(list(accn=accn, replicate=repl))
}