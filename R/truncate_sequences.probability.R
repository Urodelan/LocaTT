#' Truncate DNA Sequences at Specified Probability that All Bases were Called Correctly
#'
#' @description Calculates the cumulative probability that all bases were called correctly along each DNA sequence and truncates the DNA sequence immediately prior to the first occurrence of a probability being equal to or less than a specified value.
#' @param sequences A character vector of DNA sequences to truncate.
#' @param quality_scores A character vector of DNA sequence quality scores encoded in Sanger format.
#' @param threshold Numeric. The probability threshold used for truncation. The default is `0.5` (*i.e.*, each trimmed sequence has a greater than 50% probability that all bases were called correctly).
#' @returns A list containing two elements. The first element is a character vector of truncated DNA sequences, and the second element is a character vector of quality scores which have been truncated to their corresponding truncated DNA sequences.
#' @examples
#' truncate_sequences.probability(sequences=c("ATATAGCGCG","TGCCGATATA","ATCTATCACCGC"),
#'                                quality_scores=c("989!.C;F@\"","A((#-#;,2F","HD8I/+67=1>?"),
#'                                threshold=0.5)
#' @export
truncate_sequences.probability<-function(sequences,quality_scores,threshold=0.5){
  
  # Throw an error if the sequences are not a character vector.
  if(!is.character(sequences)) stop("The sequences must be a character vector.")
  
  # Remove sequence names, if present.
  sequences<-unname(sequences)
  
  # Throw an error if the quality scores are not a character vector.
  if(!is.character(quality_scores)) stop("The quality scores must be a character vector.")
  
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
  
  # Throw an error if threshold is less than zero.
  if(threshold < 0) stop("Threshold cannot be less than zero.")
  
  # Throw an error if threshold is greater than one.
  if(threshold > 1) stop("Threshold cannot be greater than one.")
  
  # Define function for getting the length to trim sequences to.
  get_trim_length<-function(Qscores,limit){
    # Decode the quality scores.
    decoded<-decode_quality_scores(symbols=Qscores)
    # Convert quality scores into probabilities that each base call was incorrect.
    prob_error<-10^(-decoded/10)
    # Get the probabilities that each base call was correct.
    prob_correct<-1-prob_error
    # Get the probability that all base calls were correct up to each base in the sequence.
    prob_all_correct<-cumprod(x=prob_correct)
    # If the probability that all base calls were correct up to each base
    # in the sequence is equal to or less than the limit for any base.
    if(any(prob_all_correct <= limit)){
      # Set trimmed length to immediately before the first occurrence
      # of the probability which is equal to or less than the limit.
      trim_length<-min(which(prob_all_correct <= limit))-1
    } else {
      # If the probability that any base calls were incorrect up to each base
      # in the sequence does not equal or exceed the limit for any base.
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