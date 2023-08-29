#' Truncate DNA Sequences to Specified Length
#'
#' @description Truncates DNA sequences to a specified length.
#' @param sequences A character vector of DNA sequences to truncate.
#' @param length Numeric. The length to truncate DNA sequences to.
#' @param quality_scores An optional character vector of DNA sequence quality scores. If supplied, these will be truncated to their corresponding truncated DNA sequences.
#' @returns If quality scores are not provided, then a character vector of truncated DNA sequences is returned. If quality scores are provided, then a list containing two elements is returned. The first element is a character vector of truncated DNA sequences, and the second element is a character vector of quality scores which have been truncated to their corresponding truncated DNA sequences.
#' @examples
#' truncate_sequences.length(sequences=c("ATATAGCGCG","TGCCGATATA","ATCTATCACCGC"),
#'                           length=5,
#'                           quality_scores=c("989!.C;F@\"","A((#-#;,2F","HD8I/+67=1>?"))
#' @export
truncate_sequences.length<-function(sequences,length,quality_scores){
  
  # Throw an error if the sequences are not a character vector.
  if(!is.character(sequences)) stop("The sequences must be a character vector.")
  
  # Throw an error if there are NAs in the sequences vector.
  if(any(is.na(sequences))) stop("There are NAs in the sequences vector.")
  
  # Remove sequence names, if present.
  sequences<-unname(sequences)
  
  # Throw an error if length is not of class numeric.
  if(!is.numeric(length)) stop("Length must be class numeric.")
  
  # Throw an error if length contains multiple elements.
  if(length(length) > 1) stop("Length cannot have multiple elements.")
  
  # Throw an error if length is not an integer value.
  if(round(x=length,digits=0)!=length) stop("Length must be an integer value.")
  
  # Throw an error if length is less than one.
  if(length < 1) stop("Length must be greater than zero.")
  
  # If a vector of quality scores is provided.
  if(!missing(quality_scores)){
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
  }
  
  # Trim the sequences by length.
  trimmed<-substr(x=sequences,start=1,stop=length)
  
  # If a vector of quality scores is provided.
  if(!missing(quality_scores)){
    # Trim the quality scores by length.
    trimmed_scores<-substr(x=quality_scores,start=1,stop=length)
    # Create a list containing the trimmed sequences and quality scores.
    trimmed<-list(sequences=trimmed,quality_scores=trimmed_scores)
  }
  
  # Return trimmed sequences.
  return(trimmed)
  
}