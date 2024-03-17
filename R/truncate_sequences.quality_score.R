#' Truncate DNA Sequences at Specified Quality Score
#'
#' @description Truncates DNA sequences immediately prior to the first occurrence of a Phred quality score being equal to or less than a specified value.
#' @param sequences A character vector of DNA sequences to truncate.
#' @param quality_scores A character vector of DNA sequence quality scores encoded in Sanger format.
#' @param threshold Numeric. The Phred quality score threshold used for truncation. The default is `3` (*i.e.*, each base in a truncated sequence has a greater than 50% probability of having been called correctly).
#' @returns A list containing two elements. The first element is a character vector of truncated DNA sequences, and the second element is a character vector of quality scores which have been truncated to their corresponding truncated DNA sequences.
#' @seealso
#' [`truncate_sequences.length`][truncate_sequences.length()] for truncating DNA sequences to a specified length. \cr \cr
#' [`truncate_sequences.probability`][truncate_sequences.probability()] for truncating DNA sequences by cumulative probability that all bases were called correctly. \cr \cr
#' [`truncate_and_merge_pairs`][truncate_and_merge_pairs()] for truncating and merging forward and reverse DNA sequence reads.
#' @examples
#' truncate_sequences.quality_score(sequences=c("ATATAGCGCG","TGCCGATATA","ATCTATCACCGC"),
#'                                  quality_scores=c("989!.C;F@\"","A((#-#;,2F","HD8I/+67=1>?"),
#'                                  threshold=3)
#' @export
truncate_sequences.quality_score<-function(sequences,quality_scores,threshold=3){
  
  # Throw an error if the sequences are not a character vector.
  if(!is.character(sequences)) stop("The sequences must be a character vector.")
  
  # Throw an error if there are NAs in the sequences vector.
  if(any(is.na(sequences))) stop("There are NAs in the sequences vector.")
  
  # Remove sequence names, if present.
  sequences<-unname(sequences)
  
  # Throw an error if the quality scores are not a character vector.
  if(!is.character(quality_scores)) stop("The quality scores must be a character vector.")
  
  # Throw an error if there are NAs in the quality scores vector.
  if(any(is.na(quality_scores))) stop("There are NAs in the quality scores vector.")
  
  # Remove quality score names, if present.
  quality_scores<-unname(quality_scores)
  
  # Throw an error if the sequence and quality scores vectors are not the same length.
  if(length(sequences)!=length(quality_scores)) stop("The sequence and quality scores vectors must have the same length.")
  
  # Throw an error if the number of characters differs between sequences and quality scores.
  if(!identical(nchar(quality_scores),nchar(sequences))) stop("Sequences and quality scores must have the same number of characters for each vector element.")
  
  # Throw an error if threshold is not of class numeric.
  if(!is.numeric(threshold)) stop("Threshold must be class numeric.")
  
  # Throw an error if threshold contains multiple elements.
  if(length(threshold) > 1) stop("Threshold cannot have multiple elements.")
  
  # Throw an error if threshold is not an integer value.
  if(round(x=threshold,digits=0)!=threshold) stop("Threshold must be an integer value.")
  
  # Throw an error if threshold is less than zero.
  if(threshold < 0) stop("Threshold cannot be less than zero.")
  
  # Define function for getting the length to trim sequences to.
  get_trim_length<-function(Qscores,limit){
    # Decode the quality scores.
    decoded<-decode_quality_scores(symbols=Qscores)
    # If any quality scores are equal to or less than the limit.
    if(any(decoded <= limit)){
      # Set trimmed length to immediately before the first occurrence
      # of the quality score equal to or less than the limit.
      trim_length<-min(which(decoded <= limit))-1
    } else { # If all quality scores are greater than the limit.
      # Set the trimmed length to the length of the sequence.
      trim_length<-length(decoded)
    }
    # Return the trimmed length.
    return(trim_length)
  }
  
  # Get the trim lengths for each sequence.
  trim_lengths<-sapply(X=quality_scores,
                       FUN=get_trim_length,
                       limit=threshold,
                       USE.NAMES=FALSE)
  
  # Trim sequences.
  trimmed_sequences<-substr(x=sequences,start=1,stop=trim_lengths)
  
  # Trim quality scores.
  trimmed_scores<-substr(x=quality_scores,start=1,stop=trim_lengths)
  
  # Create a list containing the trimmed sequences and quality scores.
  trimmed<-list(sequences=trimmed_sequences,quality_scores=trimmed_scores)
  
  # Return trimmed sequences.
  return(trimmed)
  
}