#' simple_dada
#'
#' Runs through the DADA2 routine with dynamic paramter
#' assignment. In situations where the data is clean this
#' will work well. If the data is not clean, it will be 
#' bet to either run the DADA2 routine. or each
#' sub-function.
#' 
#' @usage simple_dada(path, paired = TRUE, cores = 0)
#' @param path Directory containing raw fastq files.
#' @param paired Whether or not the program should look
#' for forward (R1) and reverse (R2) fastq files.
#' @param cores The number of CPU cores/threads to use.
#' @import dada2
#' @import data.table
#' @import doParallel
#' @import parallel
#' @export
#' @return data.table

simple_dada <- function(path, paired = TRUE, cores = 0){
  if(cores != 1){requireNamespace('doParallel')}
  if(cores == 0){cores <- detectCores()-1}
  if(paired){
    forward_files <- sort(list.files(file.path(path), pattern = "_R1_001.fastq", full.names = TRUE))
    reverse_files <- sort(list.files(file.path(path), pattern = "_R2_001.fastq", full.names = TRUE))
    fwd <- read_quality_report(forward_files, cores = cores)
    rev <- read_quality_report(reverse_files, cores = cores)
    fwd <- fwd[count > 0.1*median(count)]
    rev <- rev[count > 0.1*median(count)]
    if(cores == 1){
      out <- for(i in seq_along(fwd$sample)){
        filtered_f <- file.path(path, "filtered", paste0(fwd[i, 'sample'], "_F.filtered.fastq.gz"))
        filtered_r <- file.path(path, "filtered", paste0(rev[i, 'sample'], "_R.filtered.fastq.gz"))
        filterAndTrim(fwd[['file']][i], filtered_f, rev[['file']][i], filtered_r,
                      truncLen=c(fwd[['quality_length']][i], rev[['quality_length']][i]),
                      maxN = 0, maxEE = c(2,2), truncQ = 2, rm.phix = TRUE,
                      compress = TRUE, multithread = FALSE)
      }
    } else {
      cl <- makeCluster(cores, type="FORK")  
      registerDoParallel(cl)
      out <- foreach(i = seq_along(fwd$sample)) %dopar% {
        filtered_f <- file.path(path, "filtered", paste0(fwd[i, 'sample'], "_F.filtered.fastq.gz"))
        filtered_r <- file.path(path, "filtered", paste0(rev[i, 'sample'], "_R.filtered.fastq.gz"))
        filterAndTrim(fwd[['file']][i], filtered_f, rev[['file']][i], filtered_r,
                      truncLen = c(fwd[['quality_length']][i], rev[['quality_length']][i]),
                      maxN = 0, maxEE = c(2,2), truncQ = 2, rm.phix = TRUE,
                      compress = TRUE, multithread = FALSE)
      }
      stopCluster(cl)
    }
    filtered_f <- sort(list.files(file.path(path, 'filtered'), pattern = "_F.filtered", full.names = TRUE))
    filtered_r <- sort(list.files(file.path(path, 'filtered'), pattern = "_R.filtered", full.names = TRUE))
    errF <- learnErrors(filtered_f, randomize = TRUE, multithread = cores)
    errR <- learnErrors(filtered_r, randomize = TRUE, multithread = cores)
    
    sample_ids <- sapply(strsplit(basename(filtered_f), "_F.filtered"), `[`, 1)
    dadas <- vector("list", length(sample_ids))
    names(dadas) <- sample_ids
    if(cores == 1){
      for(i in seq_along(dadas)){
        derep_forward <- derepFastq(filtered_f[i], verbose = FALSE)
        derep_reverse <- derepFastq(filtered_r[i], verbose = FALSE)
        dada_forward <- dada(derep_forward, err = errF, multithread = FALSE)
        dada_reverse <- dada(derep_reverse, err = errR, multithread = FALSE) 
        dadas[[i]] <-  mergePairs(dada_forward, derep_forward, 
                                  dada_reverse, derep_reverse)
      }
    } else {
      cl <- makeCluster(cores, type = "FORK")  
      registerDoParallel(cl)
      dadas <- foreach(i = seq_along(dadas)) %dopar% {
        derep_forward <- derepFastq(filtered_f[i], verbose = FALSE)
        derep_reverse <- derepFastq(filtered_r[i], verbose = FALSE)
        dada_forward <- dada(derep_forward, err = errF, multithread = FALSE)
        dada_reverse <- dada(derep_reverse, err = errR, multithread = FALSE) 
        merged_reads <- mergePairs(dada_forward, derep_forward, 
                                   dada_reverse, derep_reverse,
                                   verbose = FALSE)
      }
      stopCluster(cl)
    }
  }
  
  if(!paired){
    files <- sort(list.files(file.path(path), pattern = ".fastq", full.names = TRUE))
    reads <- read_quality_report(files, cores = cores)
    reads <- reads[count > 0.1*median(count)]
    if(cores == 1){
      out <- for(i in seq_along(reads$sample)){
        filtered <- file.path(path, "filtered", paste0(reads[i, 'sample'], ".filtered.fastq.gz"))
        filterAndTrim(reads[['file']][i], filtered,
                      truncLen = c(reads[['quality_length']][i]),
                      maxN = 0, maxEE = c(2), truncQ = 2, rm.phix = TRUE,
                      compress = TRUE, multithread = FALSE)
      }
    } else {
      cl <- makeCluster(cores, type = "FORK")  
      registerDoParallel(cl)
      out <- foreach(i = seq_along(reads$sample)) %dopar% {
        filtered <- file.path(path, "filtered", paste0(reads[i, 'sample'], ".filtered.fastq.gz"))
        filterAndTrim(reads[['file']][i], filtered,
                      truncLen=c(reads[['quality_length']][i]),
                      maxN = 0, maxEE = c(2), truncQ = 2, rm.phix = TRUE,
                      compress = TRUE, multithread = FALSE)
      }
      stopCluster(cl)
    }
    
    filtered <- sort(list.files(file.path(path, 'filtered'), pattern = "filtered.fastq", full.names = TRUE))
    err <- learnErrors(filtered, randomize = TRUE, multithread = cores)
    
    sample_ids <- sapply(strsplit(basename(filtered), ".filtered"), `[`, 1)
    dadas <- vector("list", length(sample_ids))
    names(dadas) <- sample_ids
    if(cores == 1){
      for(i in seq_along(dadas)){
        derep <- derepFastq(filtered[i], verbose = FALSE)
        dadas[[i]] <- dada(derep, err = err, multithread = FALSE)
      }
    } else {
      cl <- makeCluster(cores, type="FORK")  
      registerDoParallel(cl)
      dadas <- foreach(i = seq_along(dadas)) %dopar% {
        derep <- derepFastq(filtered[i], verbose = FALSE)
        dada(derep, err = err, multithread = FALSE)
      }
      stopCluster(cl)
    }
  }
  
  seq_table <- makeSequenceTable(dadas)
  rownames(seq_table) <- sample_ids
  saveRDS(t(seq_table), file.path(path, 'seq_table.RDS'))
  seq_table <- removeBimeraDenovo(seq_table, method = "consensus", multithread = TRUE)
  rownames(seq_table) <- sample_ids
  saveRDS(t(seq_table), file.path(path, 'seq_table.RDS'))
  tax_table <- assignTaxonomy(seq_table, system.file("extdata", "silva_nr_v132_train_set.fa.gz", package = "simple.dada"), multithread = TRUE)
  saveRDS(tax_table, file.path(path, 'tax_table.RDS'))
}


