#' read_quality_report
#'
#' Processes fastq files to look at read counts, read lengths,
#' and at what read cycle the quality drops below the \code{q}
#' quality threshold.
#' @usage read_quality_report(path, q = 20, k = 3, n = 5e+06, cores = 1)
#' @param path File path(s) to fastq or fastq.gz file(s).
#' @param q Quality score cutoff for the read. Will look at the mean average score for 
#' \code{k} bases beyond where the potential read length cutoff would be recommended.
#' @param k number of bases beyond the current to consider for quality cutoff.
#' @param n The number of records to sample from the fastq file.
#' @param cores The number of CPU cores/threads to use.
#' @import dada2
#' @import data.table
#' @import doParallel
#' @import parallel
#' @importFrom ShortRead qa
#' @importFrom stats median
#' @export
#' @return data.table

read_quality_report <- function(path, q = 20, k = 3, n = 5e+06, cores = 1){
  if(cores != 1){requireNamespace('doParallel')}
  if(cores == 0){cores <- detectCores()-1}
  
  if(length(path) == 1 && dir.exists(path)){ 
    path <- gsub('/$','',path)
    files <- vector()
    for(extension in c(".fastq.gz$", ".fastq.bz2$", ".fastq$")){
      files <- append(files, dir(path, extension, full.names = TRUE))
    }
  } else {files <- path}
  files <- normalizePath(files)
  
  if(cores == 1){
    read_report <- data.table(file = character(), sample = character(), count = numeric(), 
                              length = numeric(), quality_length = numeric())
    for(file in files){
      srqa <- tryCatch(
        {qa(file, n = n)
          df <- srqa[["perCycle"]]$quality
          read_counts <- sum(srqa[["readCounts"]]$read)
          lengths <- as.vector(by(df, df$Cycle, function(cycle){
            sum(cycle$Count)
          }, simplify = TRUE))
          averages <- as.vector(by(df, df$Cycle, function(cycle){
            cycle$Score[min(which(cumsum(cycle$Score) >= sum(cycle$Score)/2))]
          }, simplify = TRUE))[which(lengths/read_counts > .10)]
          # averages <- rowsum(df$Score * df$Count, df$Cycle)/
          # rowsum(df$Count, df$Cycle)
          q_length <- length(averages)
          for(cycle in seq_along(averages)-1){
            if(mean(averages[(cycle+1):(cycle+(k-1))], na.rm = T) < q){
              q_length <- cycle
              break
            }
          }
          read_report <- rbind(read_report, list(file, gsub('\\..*','',basename(file)), read_counts, 
                                                 length(averages), q_length))
        },
        error = function(e){
          read_report <- rbind(read_report, 
                               list(file, gsub('\\..*','',basename(file)), 0, 0, 0))
        })
    }
  } else {
    cl <- makeCluster(cores, type="FORK")  
    registerDoParallel(cl)
    on.exit(stopCluster(cl))
    read_report <- foreach(i = seq_along(files), .combine = 'rbind') %dopar% {
      file = files[i]
      srqa <- tryCatch(
        {ShortRead::qa(file, n = n)
          df <- srqa[["perCycle"]]$quality
          read_counts <- sum(srqa[["readCounts"]]$read)
          lengths <- as.vector(by(df, df$Cycle, function(cycle){
            sum(cycle$Count)
          }, simplify = TRUE))
          lengths <- lengths[which(lengths/read_counts > .10)]
          averages <- as.vector(by(df, df$Cycle, function(cycle){
            cycle$Score[min(which(cumsum(cycle$Score) >= sum(cycle$Score)/2))]
          }, simplify = TRUE))
          # averages <- rowsum(df$Score * df$Count, df$Cycle)/
          # rowsum(df$Count, df$Cycle)
          q_length <- length(averages)
          for(cycle in seq_along(averages)-1){
            if(mean(averages[(cycle+1):(cycle+(k-1))], na.rm = T) < q){
              q_length <- cycle
              break
            }
          }
          return(data.table(file = file, sample = gsub('\\..*','',basename(file)), count = read_counts, length = length(averages), quality_length = q_length))
        },
        error = function(e){
          return(data.table(file = file, sample = gsub('\\..*','',basename(file)), count = 0, length = 0, quality_length = 0))
        })
}  
  }
  return(read_report)
}
