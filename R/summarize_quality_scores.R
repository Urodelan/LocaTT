#' Summarize Quality Scores
#'
#' @description For each base pair position, summarizes read length, Phred quality score, and the cumulative probability that all bases were called correctly.
#' @details For each combination of base pair position and read direction, calculates summary statistics of read length, Phred quality score, and the cumulative probability that all bases were called correctly. The cumulative probability is calculated from the first base pair up to the current position. Quality scores are assumed to be encoded in Sanger format. Read pairs are selected by randomly sampling up to `n.each` read pairs from each pair of input FASTQ files. By default, `n.each` is derived from `n.total`, and `n.total` will be ignored if `n.each` is provided. By default, [`mean`][mean()] is used to compute the summary statistics, but the user may provide another summary function instead (*e.g.*, [`median`][stats::median()]). Functions which return multiple summary statistics are also supported (*e.g.*, [`summary`][summary()] and [`quantile`][stats::quantile()]). Arguments in `...` are passed to the summary function.
#' @param forward_files A character vector of file paths to FASTQ files containing forward DNA sequence reads.
#' @param reverse_files A character vector of file paths to FASTQ files containing reverse DNA sequence reads.
#' @param n.total Numeric. The number of read pairs to randomly sample from the input FASTQ files. Ignored if `n.each` is specified. The default is `10000`.
#' @param n.each Numeric. The number of read pairs to randomly sample from each pair of input FASTQ files. The default is `ceiling(n.total/length(forward_files))`.
#' @param seed Numeric. The seed for randomly sampling read pairs. If `NULL` (the default), then a random seed is used.
#' @param FUN A function to compute summary statistics of the quality scores. The default is [`mean`][mean()].
#' @param ... Additional arguments passed to `FUN`.
#' @returns Returns a data frame containing summary statistics of read length and quality score at each base pair position. The returned data frame contains the following fields:
#' * Direction: The read direction (_i.e._, "Forward" or "Reverse").
#' * Position: The base pair position.
#' * Length: The summary statistic(s) of read lengths. If `FUN` returns multiple summary statistics, then a matrix of the summary statistics will be stored in this field, which can be accessed with [`$Length`][Extract()].
#' * Score: The summary statistic(s) of Phred quality scores. If `FUN` returns multiple summary statistics, then a matrix of the summary statistics will be stored in this field, which can be accessed with [`$Score`][Extract()].
#' * Probability: The summary statistic(s) of the cumulative probability that all bases were called correctly. If `FUN` returns multiple summary statistics, then a matrix of the summary statistics will be stored in this field, which can be accessed with [`$Probability`][Extract()].
#' @seealso
#' [`decode_quality_scores`][decode_quality_scores()] for decoding quality scores.
#' @examples
#' print("Insert example here.")
#' @export
summarize_quality_scores<-function(forward_files,reverse_files,n.total=10000,n.each=ceiling(n.total/length(forward_files)),seed=NULL,FUN=mean,...){
  
  ### Check arguments.
  
  # If both n.total and n.each are provided, provide a warning that n.total will be ignored.
  if(!missing(n.total) & !missing(n.each)) warning("The n.total argument is ignored if n.each is provided.")
  
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
  
  # Throw an error if the number of forward and reverse files are not the same.
  if(length(forward_files)!=length(reverse_files)) stop("There must be the same number of forward and reverse files.")
  
  # Check n.total.
  ## Throw an error if n.total is not a numeric vector.
  if(!is.numeric(n.total)) stop("n.total must be class numeric.")
  ## Throw an error if n.total has multiple elements.
  if(length(n.total) > 1) stop("n.total cannot have multiple elements.")
  ## Throw an error if n.total is not an integer value.
  if(round(x=n.total,digits=0)!=n.total) stop("n.total must be an integer value.")
  ## Throw an error if n.total is less than one.
  if(n.total < 1) stop("n.total must be greater than zero.")
  
  # Check n.each.
  ## Throw an error if n.each is not a numeric vector.
  if(!is.numeric(n.each)) stop("n.each must be class numeric.")
  ## Throw an error if n.each has multiple elements.
  if(length(n.each) > 1) stop("n.each cannot have multiple elements.")
  ## Throw an error if n.each is not an integer value.
  if(round(x=n.each,digits=0)!=n.each) stop("n.each must be an integer value.")
  ## Throw an error if n.each is less than one.
  if(n.each < 1) stop("n.each must be greater than zero.")
  
  # Check seed.
  ## If the seed is not NULL.
  if(!is.null(seed)){
    ## Throw an error if seed is not a numeric vector.
    if(!is.numeric(seed)) stop("seed must be class numeric.")
    ## Throw an error if seed has multiple elements.
    if(length(seed) > 1) stop("seed cannot have multiple elements.")
    ## Throw an error if seed is not an integer value.
    if(round(x=seed,digits=0)!=seed) stop("seed must be an integer value.")
  }
  
  # Check FUN.
  ## Throw an error if FUN is not a function.
  if(!is.function(FUN)) stop("FUN must be a function.")
  
  ### Begin operations.
  
  # Set seed.
  set.seed(seed=seed)
  
  # Create empty storage data frame.
  scores<-data.frame(NULL)
  
  # Loop through each pair of FASTQ files.
  for(i in 1:length(forward_files)){
    
    # Get the FASTQ file pair.
    ## Forward.
    forward_file<-forward_files[i]
    ## Reverse.
    reverse_file<-reverse_files[i]
    
    # Read in FASTQ files.
    ## Forward.
    fwd<-read.fastq(file=forward_file)
    ## Reverse.
    rev<-read.fastq(file=reverse_file)
    
    # Throw an error if the lengths of reads and quality scores are not the same.
    ## Forward reads and quality scores.
    if(!identical(nchar(fwd$Sequence),nchar(fwd$Quality_scores))) stop(paste0("The lengths of forward reads and forward quality scores are not the same for foward input file: ",forward_file))
    ## Reverse reads and quality scores.
    if(!identical(nchar(rev$Sequence),nchar(rev$Quality_scores))) stop(paste0("The lengths of reverse reads and reverse quality scores are not the same for reverse input file: ",reverse_file))
    
    # Remove the trailing portion of read names denoting forward or reverse,
    # which is assumed to be delimited by a space.
    ## Forward.
    fwd$Name.simple<-sub(pattern=" .*$",replacement="",x=fwd$Name)
    ## Reverse.
    rev$Name.simple<-sub(pattern=" .*$",replacement="",x=rev$Name)
    
    # Throw an error if forward and reverse read names do not match.
    if(!identical(fwd$Name.simple,rev$Name.simple)) stop(paste0("Forward and reverse read names do not match for input file pair: ",forward_file," & ",reverse_file))
    
    # Throw an error if not all read pair names are unique.
    if(any(duplicated(fwd$Name.simple))) stop(paste0("Duplicate read pair names exist in input file pair: ",forward_file," & ",reverse_file))
    
    # Rename fields to denote which are associated with forward or reverse reads.
    ## Forward.
    colnames(fwd)[colnames(fwd) %in% c("Name","Sequence","Comment","Quality_scores")]<-paste0(colnames(fwd)[colnames(fwd) %in% c("Name","Sequence","Comment","Quality_scores")],".fwd")
    ## Reverse.
    colnames(rev)[colnames(rev) %in% c("Name","Sequence","Comment","Quality_scores")]<-paste0(colnames(rev)[colnames(rev) %in% c("Name","Sequence","Comment","Quality_scores")],".rev")
    
    # Combine forward and reverse reads into a single data frame.
    df<-merge(x=fwd,y=rev)
    
    # If there are more sequences than requested.
    if(nrow(df) > n.each){
      # Randomly subset sequences to the requested number.
      df<-df[sample.int(n=nrow(df),size=n.each),]
    }
    
    # Get just the quality score fields.
    df<-df[,c("Quality_scores.fwd","Quality_scores.rev")]
    
    # Reset row names.
    row.names(df)<-1:nrow(df)
    
    # Add quality scores to the storage data frame.
    scores<-rbind(scores,df)
    
  }
  
  # Decode quality scores.
  ## Forward.
  decoded.fwd<-sapply(X=scores$Quality_scores.fwd,FUN=decode_quality_scores,simplify=FALSE)
  ## Reverse.
  decoded.rev<-sapply(X=scores$Quality_scores.rev,FUN=decode_quality_scores,simplify=FALSE)
  
  # Get length of quality scores.
  ## Forward.
  len.fwd<-sapply(X=decoded.fwd,FUN=length)
  ## Reverse.
  len.rev<-sapply(X=decoded.rev,FUN=length)
  
  # Define cumulative probability function using decoded quality scores.
  cumprob<-function(decoded){
    # Convert quality scores into probabilities that each base call was incorrect.
    prob_error<-10^(-decoded/10)
    # Get the probabilities that each base call was correct.
    prob_correct<-1-prob_error
    # Get the probability that all base calls were correct up to each base in the sequence.
    prob_all_correct<-cumprod(x=prob_correct)
    # Return the probability that all base calls were correct up to each base in the sequence.
    return(prob_all_correct)
  }
  
  # Calculate the probability that all base calls were
  # correct up to each base in the sequence.
  ## Forward.
  prob.fwd<-sapply(X=decoded.fwd,FUN=cumprob,simplify=FALSE)
  ## Reverse.
  prob.rev<-sapply(X=decoded.rev,FUN=cumprob,simplify=FALSE)
  
  # Create an empty storage data frame.
  df.full<-data.frame(NULL)
  
  # Loop through each position.
  ## Forward.
  for(i in 1:max(len.fwd)){
    ## Set the position.
    position<-i
    ## Check if each read includes position.
    len.pos<-as.numeric(position <= unname(len.fwd))
    ## Get the decoded quality scores at the position.
    decoded<-unname(sapply(X=decoded.fwd,FUN="[",position,simplify=TRUE))
    ## Get the cumulative probability at the position.
    probability<-unname(sapply(X=prob.fwd,FUN="[",position,simplify=TRUE))
    ## Collect the information for the position into a data frame.
    df.position<-data.frame(Direction="Forward",
                            Position=position,
                            Length=len.pos,
                            Score=decoded,
                            Probability=probability,
                            stringsAsFactors=FALSE)
    ## Add the information to the storage data frame.
    df.full<-rbind(df.full,df.position)
  }
  ## Reverse.
  for(i in 1:max(len.rev)){
    ## Set the position.
    position<-i
    ## Check if each read includes position.
    len.pos<-as.numeric(position <= unname(len.rev))
    ## Get the decoded quality scores at the position.
    decoded<-unname(sapply(X=decoded.rev,FUN="[",position,simplify=TRUE))
    ## Get the cumulative probability at the position.
    probability<-unname(sapply(X=prob.rev,FUN="[",position,simplify=TRUE))
    ## Collect the information for the position into a data frame.
    df.position<-data.frame(Direction="Reverse",
                            Position=position,
                            Length=len.pos,
                            Score=decoded,
                            Probability=probability,
                            stringsAsFactors=FALSE)
    ## Add the information to the storage data frame.
    df.full<-rbind(df.full,df.position)
  }
  
  # Summarize quality scores.
  ## Length.
  sum.length<-stats::aggregate(Length~Direction+Position,data=df.full,FUN=FUN,...)
  ## Score.
  sum.score<-stats::aggregate(Score~Direction+Position,data=df.full,FUN=FUN,...)
  ## Probability.
  sum.prob<-stats::aggregate(Probability~Direction+Position,data=df.full,FUN=FUN,...)
  ## Combine length and score information into a single data frame.
  summarized<-merge(x=sum.length,y=sum.score)
  ## Combine length, score, and probability information into a single data frame.
  summarized<-merge(x=summarized,y=sum.prob)
  ## Re-order records by direction and position.
  summarized<-summarized[order(summarized$Direction,summarized$Position),]
  ## Reset row names.
  row.names(summarized)<-1:nrow(summarized)
  
  # Return summarized quality scores.
  return(summarized)
  
}