#' Trim Target Nucleotide Sequence from DNA Sequences
#'
#' @description Trims a target nucleotide sequence from the front or back of DNA sequences. Ambiguous nucleotides in the target nucleotide sequence are supported.
#' @details For each DNA sequence, the target nucleotide sequence is searched for at either the front or back of the DNA sequence, depending on the value of the anchor argument. If the target nucleotide sequence is found, then it is removed from the DNA sequence. If the required argument is set to `TRUE`, then DNA sequences in which the target nucleotide sequence was not found will be returned as `NA`s. If the required argument is set to `FALSE`, then untrimmed DNA sequences will be returned along with DNA sequences for which trimming was successful. Ambiguous nucleotides in the target nucleotide sequence are supported through the internal use of the `substitute_wildcards` function on the target nucleotide sequence, and a regular expression with a leading or ending anchor is used to search for the target nucleotide sequence in the DNA sequences. Trimming will fail for DNA sequences which contain ambiguous nucleotides (*e.g.*, Ns) in their target nucleotide sequence region, resulting in `NA`s for those sequences if the required argument is set to `TRUE`.
#' @param sequences A character vector of DNA sequences to trim.
#' @param target A string specifying the target nucleotide sequence.
#' @param anchor A string specifying whether the target nucleotide sequence should be trimmed from the start or end of the DNA sequences. Allowable values are "start" (the default) and "end".
#' @param required A logical value specifying whether trimming is required. If `TRUE` (the default), then sequences which could not be trimmed are returned as `NA`s. If `FALSE`, then untrimmed sequences are returned along with DNA sequences for which trimming was successful.
#' @param quality_scores An optional character vector of DNA sequence quality scores. If supplied, these will be trimmed to their corresponding trimmed DNA sequences.
#' @returns If quality scores are not provided, then a character vector of trimmed DNA sequences is returned. If quality scores are provided, then a list containing two elements is returned. The first element is a character vector of trimmed DNA sequences, and the second element is a character vector of quality scores which have been trimmed to their corresponding trimmed DNA sequences.
#' @examples
#' trim_sequences(sequences=c("ATATAGCGCG","TGCCGATATA","ATCTATCACCGC"),
#'                target="ATMTA",
#'                anchor="start",
#'                required=TRUE,
#'                quality_scores=c("989!.C;F@\"","A((#-#;,2F","HD8I/+67=1>?"))
#' @export
trim_sequences<-function(sequences,target,anchor="start",required=TRUE,quality_scores){
  
  # Create a vector of supported nucleotides.
  supported_nucleotides<-c("A","T","G","C","R","Y","S","W","K","M","B","D","H","V","N")
  
  # Throw an error if the sequences are not a character vector.
  if(!is.character(sequences)) stop("The sequences must be a character vector.")
  
  # Remove sequence names, if present.
  sequences<-unname(sequences)
  
  # Throw an error if the target is not a character string.
  if(!is.character(target)) stop("The target must be a character string.")
  # Throw an error if the length of the target character string is not one.
  if(length(target)!=1) stop("The target must be a character string of length 1.")
  # Split the target into individual nucleotides.
  target_split<-strsplit(x=target,split="")[[1]]
  # If any nucleotides in the target are not supported.
  if(!all(target_split %in% supported_nucleotides)){
    # Get the unsupported nucleotides from the target.
    target_unsupported<-sort(unique(target_split[!(target_split %in% supported_nucleotides)]))
    # Throw an error listing the unsupported nucleotides from the target.
    stop(paste0("The target contains the following unsupported nucleotides: ",paste(target_unsupported,collapse=", ")),". Supported nucleotides are: ",paste(supported_nucleotides,collapse=", "),".")
  }
  
  # Throw an error if anchor is not a character string.
  if(!is.character(anchor)) stop("Anchor must be a character string.")
  # Throw an error if the length of the anchor character string is not one.
  if(length(anchor)!=1) stop("Anchor must be a character string of length 1.")
  # Throw an error if anchor is not "start" or "end".
  if(!(anchor %in% c("start","end"))) stop("Anchor must be 'start' or 'end'.")
  
  # Throw an error if required is not a logical string.
  if(!is.logical(required)) stop("Required must be a logical value.")
  # Throw an error if the length of required is not one.
  if(length(required)!=1) stop("Required must be a logical value of length 1.")
  # Throw an error if required is not TRUE or FALSE (e.g., NA).
  if(!(required %in% c(TRUE,FALSE))) stop("Required must be TRUE or FALSE.")
  
  # If a vector of quality scores is provided.
  if(!missing(quality_scores)){
    # Throw an error if the quality scores are not a character vector.
    if(!is.character(quality_scores)) stop("The quality scores must be a character vector.")
    # Remove quality score names, if present.
    quality_scores<-unname(quality_scores)
    # Throw an error if the sequence and quality scores vectors are not the same length.
    if(length(sequences)!=length(quality_scores)) stop("The sequence and quality scores vectors must have the same length.")
    # Throw an error if the number of characters differs between sequences and quality scores.
    if(!identical(nchar(quality_scores),nchar(sequences))) stop("Sequences and quality scores must have the same number of characters for each vector element.")
  }
  
  # Substitute wildcard characters in the target with their respective nucleotides.
  target<-substitute_wildcards(sequence=target)
  
  # If the anchor is at the start.
  if(anchor=="start"){
    # Anchor the target in front.
    target<-paste0("^",target)
  } else {
    # Anchor the target in back.
    target<-paste0(target,"$")
  }
  
  # If trimming is required.
  if(required){
    # Get which sequences the target is found in.
    matches<-grepl(pattern=target,x=sequences)
    # Create an empty character vector for trimmed sequences.
    trimmed<-vector(mode="character",length=length(sequences))
    # Trim the sequences.
    trimmed[matches]<-sub(pattern=target,replacement="",x=sequences[matches])
    # Replace sequences which could not be trimmed with NAs.
    trimmed[!matches]<-NA
  } else { # If trimming is not required.
    # Trim the sequences, leaving in untrimmed sequences.
    trimmed<-sub(pattern=target,replacement="",x=sequences)
  }
  
  # If a vector of quality scores is provided.
  if(!missing(quality_scores)){
    # Get the number of characters remaining in each trimmed sequence.
    trimmed_lengths<-nchar(trimmed)
    # Get the number of characters in the original sequences.
    original_lengths<-nchar(sequences)
    # If the anchor is at the start of the sequence.
    if(anchor=="start"){
      # Get start positions for the trimmed quality scores.
      start_positions<-original_lengths-trimmed_lengths+1  
      # Trim the quality scores to the trimmed sequences.
      trimmed_scores<-substr(x=quality_scores,start=start_positions,stop=original_lengths)
      
    } else { # If the anchor is at the end of the sequence.
      # Trim the quality scores to the trimmed sequences.
      trimmed_scores<-substr(x=quality_scores,start=1,stop=trimmed_lengths)
    }
    # Create a list containing the trimmed sequences and quality scores.
    trimmed<-list(sequences=trimmed,quality_scores=trimmed_scores)
  }
  
  # Return trimmed sequences.
  return(trimmed)
  
}