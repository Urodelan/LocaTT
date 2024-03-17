#' Truncate and Merge Forward and Reverse DNA Sequence Reads
#'
#' @description Removes DNA read pairs containing ambiguous nucleotides, truncates reads by length and quality score, and merges forward and reverse reads.
#' @details For each pair of input FASTQ files, removes DNA read pairs containing ambiguous nucleotides, truncates reads by length, quality score threshold, and probability threshold (in that order), and then merges forward and reverse reads. Merged reads are summarized by frequency of occurrence and written to a FASTA file. See [`contains_wildcards`][contains_wildcards()], [`truncate_sequences.length`][truncate_sequences.length()], [`truncate_sequences.quality_score`][truncate_sequences.quality_score()], [`truncate_sequences.probability`][truncate_sequences.probability()], and [`merge_pairs`][merge_pairs()] for methods. Quality scores are assumed to be encoded in Sanger format. Forward and reverse reads can be truncated by different thresholds (see `truncation_length`, `threshold.quality_score`, and `threshold.probability` arguments). \cr \cr
#' Multicore parallel processing is supported on Mac and Linux operating systems (not available on Windows). When `cores > 1` (parallel processing enabled), warnings and errors are printed to the console in addition to being invisibly returned as a list (see the return value section), and errors produced while processing a pair of FASTQ files will not interrupt the processing of other FASTQ file pairs. When `cores = 1`, FASTQ file pairs are processed sequentially on a single core, and errors will prevent the processing of subsequent FASTQ file pairs (but warnings will not).
#' @param forward_files A character vector of file paths to FASTQ files containing forward DNA sequence reads.
#' @param reverse_files A character vector of file paths to FASTQ files containing reverse DNA sequence reads.
#' @param output_files A character vector of file paths to output FASTA files.
#' @param truncation_length Numeric. The length to truncate DNA sequences to (passed to the `length` argument of [`truncate_sequences.length`][truncate_sequences.length()]). If `NA` (the default), then DNA sequences are not truncated by length. If a single value is supplied, then both forward and reverse reads are truncated to the same length. If two values are supplied in a numeric vector, then the first value is used to truncate the forward reads, and the second value is used to truncate the reverse reads. `NA` can also be supplied as either the first or second element of the numeric vector to prevent length truncation of the respective read direction while allowing the other read direction to be length truncated.
#' @param threshold.quality_score Numeric. The Phred quality score threshold used for truncation (passed to the `threshold` argument of [`truncate_sequences.quality_score`][truncate_sequences.quality_score()]). The default is `3` (*i.e.*, each base in a truncated sequence has a greater than 50% probability of having been called correctly). If `NA`, then DNA sequences are not truncated by quality score threshold. If a single value is supplied, then both forward and reverse reads are truncated by the same quality score threshold. If two values are supplied in a numeric vector, then the first value is used to truncate the forward reads, and the second value is used to truncate the reverse reads. `NA` can also be supplied as either the first or second element of the numeric vector to prevent quality-score-threshold truncation of the respective read direction while allowing the other read direction to be quality-score-threshold truncated.
#' @param threshold.probability Numeric. The probability threshold used for truncation (passed to the `threshold` argument of [`truncate_sequences.probability`][truncate_sequences.probability()]). The default is `0.5` (*i.e.*, each truncated sequence has a greater than 50% probability that all bases were called correctly). If `NA`, then DNA sequences are not truncated by probability threshold. If a single value is supplied, then both forward and reverse reads are truncated by the same probability threshold. If two values are supplied in a numeric vector, then the first value is used to truncate the forward reads, and the second value is used to truncate the reverse reads. `NA` can also be supplied as either the first or second element of the numeric vector to prevent probability-threshold truncation of the respective read direction while allowing the other read direction to be probability-threshold truncated.
#' @param minimum_overlap Numeric. The minimum length of an overlap that must be found between the end of the forward read and the start of the reverse complement of the reverse read in order for a read pair to be merged (passed to [`merge_pairs`][merge_pairs()]). The default is `10`.
#' @param cores Numeric. If `1` (the default), then FASTQ file pairs are processed sequentially on a single core. If greater than `1`, then FASTQ file pairs are processed in parallel across the specified number of cores. Parallel processing is not supported on Windows.
#' @param progress Logical. If `TRUE`, then a progress indicator is printed to the console. Ignored if `cores > 1`. If `FALSE` (the default), then no progress indicator is displayed.
#' @returns If `cores = 1`, then no return value. Writes a FASTA file for each pair of input FASTQ files with DNA sequence counts stored in the header lines. If `cores > 1`, then also invisibly returns a list where each element contains warning or error messages associated with processing each pair of input FASTQ files. A `NULL` value in the returned list means that no warnings or errors were generated from processing the respective pair of FASTQ files.
#' @seealso
#' [`contains_wildcards`][contains_wildcards()] for detecting ambiguous nucleotides in DNA sequences. \cr \cr
#' [`truncate_sequences.length`][truncate_sequences.length()] for truncating DNA sequences to a specified length. \cr \cr
#' [`truncate_sequences.quality_score`][truncate_sequences.quality_score()] for truncating DNA sequences by Phred quality score. \cr \cr
#' [`truncate_sequences.probability`][truncate_sequences.probability()] for truncating DNA sequences by cumulative probability that all bases were called correctly. \cr \cr
#' [`merge_pairs`][merge_pairs()] for merging forward and reverse DNA sequence reads. \cr \cr
#' [`filter_sequences`][filter_sequences()] for filtering merged read pairs by PCR replicate.
#' @references A manuscript describing these methods is in preparation.
#' @examples
#' print("Insert example here.")
#' @export
truncate_and_merge_pairs<-function(forward_files,reverse_files,output_files,truncation_length=NA,threshold.quality_score=3,threshold.probability=0.5,minimum_overlap=10,cores=1,progress=FALSE){
  
  # Check forward files input.
  ## Throw an error if the forward files are not a character vector.
  if(!is.character(forward_files)) stop("Forward files must be a character vector.")
  ## Throw an error if the forward files contain any NAs.
  if(any(is.na(forward_files))) stop("Forward files cannot contain NAs.")
  
  # Check reverse files input.
  ## Throw an error if the reverse files are not a character vector.
  if(!is.character(reverse_files)) stop("Reverse files must be a character vector.")
  ## Throw an error if the reverse files contain any NAs.
  if(any(is.na(reverse_files))) stop("Reverse files cannot contain NAs.")
  
  # Check output files input.
  ## Throw an error if the output files are not a character vector.
  if(!is.character(output_files)) stop("Output files must be a character vector.")
  ## Throw an error if the output files contain any NAs.
  if(any(is.na(output_files))) stop("Output files cannot contain NAs.")
  
  # Throw an error if the number of forward and reverse files are not the same.
  if(length(forward_files)!=length(reverse_files)) stop("There must be the same number of forward and reverse files.")
  
  # Throw an error if the number of output files is not
  # the same as the number of paired input files.
  if(length(forward_files)!=length(output_files)) stop("There must be the same number of output files as paired input files.")
  
  # Check truncation length input.
  ## Throw an error if truncation length is not a vector.
  if(!is.vector(truncation_length)) stop("Truncation length must be a vector.")
  ## Throw an error if truncation length has multiple elements.
  if(length(truncation_length) > 2) stop("Truncation length cannot have more than two elements.")
  ## If truncation length has a single element, duplicate the value
  ## for forward and reverse reads.
  if(length(truncation_length)==1){
    ## If truncation length is NA.
    if(is.na(truncation_length)){
      ## Format NAs as numeric.
      truncation_length<-as.numeric(rep(x=truncation_length,times=2))
    } else { ## If truncation length is not NA.
      ## Leave original class.
      truncation_length<-rep(x=truncation_length,times=2)
    }
  } else { ## If truncation length has two elements.
    ## If both elements of truncation length are NA.
    if(all(is.na(truncation_length))){
      ## Format NAs as numeric.
      truncation_length<-as.numeric(truncation_length)
    }
  }
  ## Throw an error if truncation length is not a numeric vector.
  if(!is.numeric(truncation_length)) stop("Truncation length must be class numeric.")
  ## Throw an error if either truncation length is not an integer value.
  if(!all(round(x=truncation_length,digits=0)==truncation_length,na.rm=TRUE)) stop("Truncation length must be an integer value.")
  ## Throw an error if either truncation length is less than one.
  if(any(truncation_length < 1,na.rm=TRUE)) stop("Truncation length must be greater than zero.")
  
  # Check quality score threshold input.
  ## Throw an error if quality score threshold is not a vector.
  if(!is.vector(threshold.quality_score)) stop("Quality score threshold must be a vector.")
  ## Throw an error if quality score threshold has multiple elements.
  if(length(threshold.quality_score) > 2) stop("Quality score threshold cannot have more than two elements.")
  ## If quality score threshold has a single element, duplicate the value
  ## for forward and reverse reads.
  if(length(threshold.quality_score)==1){
    ## If quality score threshold is NA.
    if(is.na(threshold.quality_score)){
      ## Format NAs as numeric.
      threshold.quality_score<-as.numeric(rep(x=threshold.quality_score,times=2))
    } else { ## If quality score threshold is not NA.
      ## Leave original class.
      threshold.quality_score<-rep(x=threshold.quality_score,times=2)
    }
  } else { ## If quality score threshold has two elements.
    ## If both elements of quality score threshold are NA.
    if(all(is.na(threshold.quality_score))){
      ## Format NAs as numeric.
      threshold.quality_score<-as.numeric(threshold.quality_score)
    }
  }
  ## Throw an error if quality score threshold is not a numeric vector.
  if(!is.numeric(threshold.quality_score)) stop("Quality score threshold length must be class numeric.")
  ## Throw an error if either quality score threshold is not an integer value.
  if(!all(round(x=threshold.quality_score,digits=0)==threshold.quality_score,na.rm=TRUE)) stop("Quality score threshold must be an integer value.")
  ## Throw an error if either quality score threshold is less than zero.
  if(any(threshold.quality_score < 0,na.rm=TRUE)) stop("Quality score threshold cannot be less than zero.")
  
  # Check probability threshold input.
  ## Throw an error if probability threshold is not a vector.
  if(!is.vector(threshold.probability)) stop("Probability threshold must be a vector.")
  ## Throw an error if probability threshold has multiple elements.
  if(length(threshold.probability) > 2) stop("Probability threshold cannot have more than two elements.")
  ## If probability threshold has a single element, duplicate the value
  ## for forward and reverse reads.
  if(length(threshold.probability)==1){
    ## If probability threshold is NA.
    if(is.na(threshold.probability)){
      ## Format NAs as numeric.
      threshold.probability<-as.numeric(rep(x=threshold.probability,times=2))
    } else { ## If probability threshold is not NA.
      ## Leave original class.
      threshold.probability<-rep(x=threshold.probability,times=2)
    }
  } else { ## If probability threshold has two elements.
    ## If both elements of probability threshold are NA.
    if(all(is.na(threshold.probability))){
      ## Format NAs as numeric.
      threshold.probability<-as.numeric(threshold.probability)
    }
  }
  ## Throw an error if probability threshold is not a numeric vector.
  if(!is.numeric(threshold.probability)) stop("Probability threshold length must be class numeric.")
  ## Throw an error if either probability threshold is less than zero.
  if(any(threshold.probability < 0,na.rm=TRUE)) stop("Probability threshold cannot be less than zero.")
  ## Throw an error if either probability threshold is greater than one.
  if(any(threshold.probability > 1,na.rm=TRUE)) stop("Probability threshold cannot be greater than one.")
  
  # Check minimum overlap input.
  ## Throw an error if minimum overlap is not a numeric vector.
  if(!is.numeric(minimum_overlap)) stop("Minimum overlap must be class numeric.")
  ## Throw an error if minimum overlap has multiple elements.
  if(length(minimum_overlap) > 1) stop("Minimum overlap cannot have multiple elements.")
  ## Throw an error if minimum overlap is not an integer value.
  if(round(x=minimum_overlap,digits=0)!=minimum_overlap) stop("Minimum overlap must be an integer value.")
  ## Throw an error if minimum overlap is less than one.
  if(minimum_overlap < 1) stop("Minimum overlap must be greater than zero.")
  
  # Check cores input.
  ## Throw an error if cores is not a numeric vector.
  if(!is.numeric(cores)) stop("Cores must be class numeric.")
  ## Throw an error if cores has multiple elements.
  if(length(cores) > 1) stop("Cores cannot have multiple elements.")
  ## Throw an error if cores is not an integer value.
  if(round(x=cores,digits=0)!=cores) stop("Cores must be an integer value.")
  ## Throw an error if cores is less than one.
  if(cores < 1) stop("Cores must be greater than zero.")
  
  # Check progress input.
  ## Throw an error if the progress argument is not TRUE or FALSE.
  if(!is.logical(progress) | is.na(progress)) stop("Progress must be TRUE or FALSE.")
  
  # Define internal function for truncating and
  # merging read pairs for a single file pair.
  truncate_and_merge<-function(forward_file,reverse_file,output_file,truncation_length,threshold.quality_score,threshold.probability,minimum_overlap){
    
    # Read in forward fastq file.
    fwd<-read.fastq(file=forward_file)
    
    # Read in reverse fastq file.
    rev<-read.fastq(file=reverse_file)
    
    # Throw an error if the lengths of reads and quality scores are not the same.
    ## Forward reads and quality scores.
    if(!identical(nchar(fwd$Sequence),nchar(fwd$Quality_scores))) stop(paste0("The lengths of forward reads and forward quality scores are not the same for foward input file: ",forward_file," (unwritten output file: ",output_file,")"))
    ## Reverse reads and quality scores.
    if(!identical(nchar(rev$Sequence),nchar(rev$Quality_scores))) stop(paste0("The lengths of reverse reads and reverse quality scores are not the same for reverse input file: ",reverse_file," (unwritten output file: ",output_file,")"))
    
    # Remove the trailing portion of read names denoting forward or reverse,
    # which is assumed to be delimited by a space.
    ## Forward.
    fwd$Name.simple<-sub(pattern=" .*$",replacement="",x=fwd$Name)
    ## Reverse.
    rev$Name.simple<-sub(pattern=" .*$",replacement="",x=rev$Name)
    
    # Throw an error if forward and reverse read names do not match.
    if(!identical(fwd$Name.simple,rev$Name.simple)) stop(paste0("Forward and reverse read names do not match for input file pair: ",forward_file," & ",reverse_file," (unwritten output file: ",output_file,")"))
    
    # Throw an error if not all read pair names are unique.
    if(any(duplicated(fwd$Name.simple))) stop(paste0("Duplicate read pair names exist in input file pair: ",forward_file," & ",reverse_file," (unwritten output file: ",output_file,")"))
    
    # Rename fields to denote which are associated with forward or reverse reads.
    ## Forward.
    colnames(fwd)[colnames(fwd) %in% c("Name","Sequence","Comment","Quality_scores")]<-paste0(colnames(fwd)[colnames(fwd) %in% c("Name","Sequence","Comment","Quality_scores")],".fwd")
    ## Reverse.
    colnames(rev)[colnames(rev) %in% c("Name","Sequence","Comment","Quality_scores")]<-paste0(colnames(rev)[colnames(rev) %in% c("Name","Sequence","Comment","Quality_scores")],".rev")
    
    # Combine forward and reverse reads into a single data frame.
    df<-merge(x=fwd,y=rev)
    
    # Remove read pairs which contain ambiguous nucleotides.
    df<-df[!(contains_wildcards(sequences=df$Sequence.fwd) | contains_wildcards(sequences=df$Sequence.rev)),]
    
    # If at least some sequences remain after removing
    # those containing ambiguous nucleotides.
    if(nrow(df)!=0){
      
      # Create initial lists of sequences and quality scores.
      ## Forward reads.
      fwd.trunc<-list(sequences=df$Sequence.fwd,quality_scores=df$Quality_scores.fwd)
      ## Reverse reads.
      rev.trunc<-list(sequences=df$Sequence.rev,quality_scores=df$Quality_scores.rev)
      
      # Truncate reads by length.
      ## Forward reads.
      ### If the forward truncation length is not set to NA.
      if(!is.na(truncation_length[1])){
        ### Truncate forward reads to the user-specified length.
        fwd.trunc<-truncate_sequences.length(sequences=fwd.trunc$sequences,
                                             quality_scores=fwd.trunc$quality_scores,
                                             length=truncation_length[1])
      }
      ## Reverse reads.
      ### If the reverse truncation length is not set to NA.
      if(!is.na(truncation_length[2])){
        ### Truncate reverse reads to the user-specified length.
        rev.trunc<-truncate_sequences.length(sequences=rev.trunc$sequences,
                                             quality_scores=rev.trunc$quality_scores,
                                             length=truncation_length[2])
      }
      
      # Truncate reads by quality score.
      ## Forward reads.
      ### If the forward quality score threshold is not set to NA.
      if(!is.na(threshold.quality_score[1])){
        # Truncate forward reads such that each base in a truncated read has a
        # quality score higher than the user-defined quality score.
        fwd.trunc<-truncate_sequences.quality_score(sequences=fwd.trunc$sequences,
                                                    quality_scores=fwd.trunc$quality_scores,
                                                    threshold=threshold.quality_score[1])
      }
      ## Reverse reads.
      ### If the reverse quality score threshold is not set to NA.
      if(!is.na(threshold.quality_score[2])){
        # Truncate reverse reads such that each base in a truncated read has a
        # quality score higher than the user-defined quality score.
        rev.trunc<-truncate_sequences.quality_score(sequences=rev.trunc$sequences,
                                                    quality_scores=rev.trunc$quality_scores,
                                                    threshold=threshold.quality_score[2])
      }
      
      # Truncate reads by probability.
      ## Forward reads.
      ### If the forward probability threshold is not set to NA.
      if(!is.na(threshold.probability[1])){
        # Truncate forward reads such that each truncated read has a greater
        # than the user-defined probability that all bases were called correctly.
        fwd.trunc<-truncate_sequences.probability(sequences=fwd.trunc$sequences,
                                                  quality_scores=fwd.trunc$quality_scores,
                                                  threshold=threshold.probability[1])
      }
      ## Reverse reads.
      ### If the reverse probability threshold is not set to NA.
      if(!is.na(threshold.probability[2])){
        # Truncate reverse reads such that each truncated read has a greater
        # than the user-defined probability that all bases were called correctly.
        rev.trunc<-truncate_sequences.probability(sequences=rev.trunc$sequences,
                                                  quality_scores=rev.trunc$quality_scores,
                                                  threshold=threshold.probability[2])
      }
      
      # Add truncated reads and quality scores back to the primary data frame.
      ## Forward.
      ### Reads.
      df$Sequence.fwd<-fwd.trunc$sequences
      ### Quality scores.
      df$Quality_scores.fwd<-fwd.trunc$quality_scores
      ## Reverse.
      ### Reads.
      df$Sequence.rev<-rev.trunc$sequences
      ### Quality scores.
      df$Quality_scores.rev<-rev.trunc$quality_scores
      
      # Reduce primary data frame to remaining fields of interest.
      df<-df[,c("Name.simple","Sequence.fwd","Sequence.rev")]
      
      # Get unique read pairs.
      unique_read_pairs<-unique(df[,c("Sequence.fwd","Sequence.rev")])
      
      # Merge unique read pairs. NAs will be returned for read pairs which cannot be merged,
      # or for read pairs whose forward or reverse reads are blank (i.e., ""). Merging only
      # unique read pairs reduces computation load.
      unique_read_pairs$Sequence.merged<-merge_pairs(
        forward_reads=unique_read_pairs$Sequence.fwd,
        reverse_reads=unique_read_pairs$Sequence.rev,
        minimum_overlap=minimum_overlap)
      
      # Add merged read pairs to the primary data frame.
      df<-merge(x=df,y=unique_read_pairs)
      
      # Remove read pairs from the primary data frame which could not be merged
      # (i.e., read pairs which contain NAs in the merged field).
      df<-df[!is.na(df$Sequence.merged),]
      
      # If at least some read pairs could be merged.
      if(nrow(df)!=0){
        
        # Create a contingency table of sequence counts.
        counts<-as.data.frame(table(df$Sequence.merged),stringsAsFactors=FALSE)
        
        # Provide field names for the contingency table of sequence counts.
        colnames(counts)<-c("Sequence","Count")
        
        # Write out sequence information as a fasta file
        # with sequence counts stored in the header lines.
        write.fasta(names=paste0("Frequency: ",counts$Count),
                    sequences=counts$Sequence,
                    file=output_file)
        
      } else { # If no read pairs could be merged.
        
        # Throw a warning stating that no read pairs could be merged for the file pair.
        warning(paste0("No read pairs could be merged for input file pair: ",forward_file," & ",reverse_file," (unwritten output file: ",output_file,")"))
        
      }
      
    } else { # If no sequences remain after removing those containing ambiguous nucleotides.
      
      # Throw a warning stating that no read pairs could be merged for the file pair.
      warning(paste0("No read pairs could be merged for input file pair: ",forward_file," & ",reverse_file," (unwritten output file: ",output_file,")"))
      
    }
    
  }
  
  # If the use of a single core is requested.
  if(cores==1){
    
    # Report initial progress.
    if(progress){
      print("Reporting progress...")
      print(paste0("- Progress: 0 of ",length(forward_files)," (0%)"))
    }
    
    # Loop through each pair of fastq files.
    for(i in 1:length(forward_files)){
      
      # Get the file paths.
      ## Forward file.
      forward_file<-forward_files[i]
      ## Reverse file.
      reverse_file<-reverse_files[i]
      ## Output file.
      output_file<-output_files[i]
      
      # Truncate and merge read pairs for the pair of fastq files.
      truncate_and_merge(forward_file=forward_file,
                         reverse_file=reverse_file,
                         output_file=output_file,
                         truncation_length=truncation_length,
                         threshold.quality_score=threshold.quality_score,
                         threshold.probability=threshold.probability,
                         minimum_overlap=minimum_overlap)
      
      # Report progress.
      if(progress){
        # If this is the last iteration.
        if(i==length(forward_files)){
          # Report final progress percentage without a decimal place.
          print(paste0("- Progress: ",i," of ",length(forward_files)," (100%)"))
        } else { # If this is not the last iteration.
          # Report progress percentage with a single decimal place.
          print(paste0("- Progress: ",i," of ",length(forward_files),
                       " (",formatC(x=i/length(forward_files)*100,
                                    digits=1,format="f"),"%)"))
        }
      }
      
    }
    
  } else { # If the use of multiple cores is requested.
    
    # Truncate and merge read pairs in parallel using the user-specified number of cores.
    messages<-parallel::mclapply(X=1:length(forward_files),
                                 FUN=function(i) truncate_and_merge(
                                   forward_file=forward_files[i],
                                   reverse_file=reverse_files[i],
                                   output_file=output_files[i],
                                   truncation_length=truncation_length,
                                   threshold.quality_score=threshold.quality_score,
                                   threshold.probability=threshold.probability,
                                   minimum_overlap=minimum_overlap),
                                 mc.cores=cores)
    
    # Multicore error handling.
    ## Locate error messages in the returned message list.
    contains.errors<-sapply(X=messages,FUN=inherits,what="try-error")
    ## Locate warning messages in the returned message list.
    contains.warnings<-!(contains.errors | sapply(X=messages,FUN=is.null))
    ## Create storage vector for formatted error and warning messages.
    messages.print<-c()
    ## If the returned message list contains any error messages.
    if(any(contains.errors)){
      ## Get the list elements which are error messages.
      messages.errors<-messages[contains.errors]
      ## Extract the message attribute of the error messages.
      messages.errors<-sapply(X=messages.errors,
                              FUN=function(x) attributes(x)$condition$message)
      ## Append "Error: " to the front of the error messages.
      messages.errors<-paste0("Error: ",messages.errors)
      ## Add the error messages to the messages storage vector.
      messages.print<-c(messages.print,messages.errors)
    }
    ## If the returned message list contains any warning messages.
    if(any(contains.warnings)){
      ## Get the list elements which are warning messages.
      messages.warnings<-unlist(messages[contains.warnings])
      ## Append "Warning: " to the front of the warning messages.
      messages.warnings<-paste0("Warning: ",messages.warnings)
      ## Add the warning messages to the messages storage vector.
      messages.print<-c(messages.print,messages.warnings)
    }
    ## If any error or warning messages exist.
    if(length(messages.print) > 0){
      ## Loop through each error or warning message.
      for(i in 1:length(messages.print)){
        ## Produce a warning containing the error or warning message.
        warning(messages.print[i])
      }
    }
    
    # Invisibly return list of warning and error messages generated by parallel::mclapply.
    return(invisible(messages))
    
  }
  
}