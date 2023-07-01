#' Check Whether DNA Sequences Contain Wildcard Characters
#'
#' Checks whether DNA sequences contain wildcard characters.
#' @param sequences A character vector of DNA sequences.
#' @returns A logical vector indicating whether each DNA sequence contains wildcard characters.
#' @examples
#' contains_wildcards(sequences=c("TKCTAGGTGW","CATAATTAGG","ATYGGCTATG"))
#' @export
contains_wildcards<-function(sequences){
  
  # Throw an error if the sequences is not a character vector.
  if(!is.character(sequences)) stop("The sequences must be a character vector.")
  
  # Throw an error if there are NAs in the sequences vector.
  if(any(is.na(sequences))) stop("There are NAs in the sequences vector.")
  
  # Remove sequence names, if present.
  sequences<-unname(sequences)
  
  # Define vector of wildcard characters.
  wildcards<-c("Y","R","W","S","K","M","D","V","H","B","N")
  
  # Define vector of supported nucleotides.
  supported<-c("A","G","C","T",wildcards)
  
  # Split sequences into individual bases.
  seqs_split<-strsplit(x=sequences,split="")
  
  # Get all unique nucleotides included in the sequences.
  all_seq_nucleotides<-unique(unlist(seqs_split))
  
  # If any nucleotides in the sequences are not supported.
  if(!all(all_seq_nucleotides %in% supported)){
    # Get the unsupported nucleotides.
    unsupported<-sort(all_seq_nucleotides[!(all_seq_nucleotides %in% supported)])
    # Throw an error listing the unsupported nucleotides.
    stop(paste0("The sequences contain the following unsupported nucleotides: ",paste(unsupported,collapse=", ")),". Supported nucleotides are: ",paste(supported,collapse=", "),".")
      
  }
  
  # Check whether each sequence contains wildcard characters.
  result<-sapply(X=seqs_split,FUN=function(x) any(x %in% wildcards))
  
  # Return the results.
  return(result)
  
}