#' Trim DNA Sequences to an Amplicon Region Using Forward and Reverse Primer Sequences
#'
#' @description Trims DNA sequences to an amplicon region using forward and reverse primer sequences. Ambiguous nucleotides in forward and reverse primers are supported.
#' @details For each DNA sequence, nucleotides matching and preceding the forward primer are removed, and nucleotides matching and following the reverse complement of the reverse primer are removed. The reverse complement of the reverse primer is internally derived from the reverse primer using the `reverse_complement` function. Ambiguous nucleotides in primers (*i.e.*, the forward and reverse primer arguments) are supported through the internal use of the `substitute_wildcards` function on the forward primer and the reverse complement of the reverse primer, and primer regions in DNA sequences are located using regular expressions. Trimming will fail for DNA sequences which contain ambiguous nucleotides in their primer regions (*e.g.*, Ns), resulting in `NA`s for those sequences.
#' @param sequences A character vector of DNA sequences to trim to the amplicon region.
#' @param forward_primer A string specifying the forward primer sequence. Can contain ambiguous nucleotides.
#' @param reverse_primer A string specifying the reverse primer sequence. Can contain ambiguous nucletodies.
#' @returns A character vector of DNA sequences trimmed to the amplicon region. `NA`s are returned for DNA sequences which could not be trimmed, which occurs when either primer region is missing from the DNA sequence or when the forward primer region occurs after a region matching the reverse complement of the reverse primer.
#' @examples
#' isolate_amplicon(sequences=c("ACACAATCGTGTTTATATTAACTTCAAGAGTGGGCATAGG",
#'                              "CGTGACAATCATGTTTGTGATTCGTACAAAAGTGCGTCCT"),
#'                  forward_primer="AATCRTGTTT",
#'                  reverse_primer="CSCACTHTTG")
#' @export
isolate_amplicon<-function(sequences,forward_primer,reverse_primer){
  
  # Create a vector of supported nucleotides.
  supported_nucleotides<-c("A","T","G","C","R","Y","S","W","K","M","B","D","H","V","N")
  
  # Throw an error if the sequences are not a character vector.
  if(!is.character(sequences)) stop("The sequences must be a character vector.")
  
  # Throw an error if there are NAs in the sequences vector.
  if(any(is.na(sequences))) stop("There are NAs in the sequences vector.")
  
  # Remove sequence names, if present.
  sequences<-unname(sequences)
  
  # Throw an error if the forward primer is not a character string.
  if(!is.character(forward_primer)) stop("The forward primer must be a character string.")
  # Throw an error if the length of the forward primer character string is not one.
  if(length(forward_primer)!=1) stop("The forward primer must be a character string of length 1.")
  # Split the forward primer into individual nucleotides.
  fwd_split<-strsplit(x=forward_primer,split="")[[1]]
  # If any nucleotides in the forward primer are not supported.
  if(!all(fwd_split %in% supported_nucleotides)){
    # Get the unsupported nucleotides from the forward primer.
    fwd_unsupported<-sort(unique(fwd_split[!(fwd_split %in% supported_nucleotides)]))
    # Throw an error listing the unsupported nucleotides from the forward primer.
    stop(paste0("The forward primer contains the following unsupported nucleotides: ",paste(fwd_unsupported,collapse=", ")),". Supported nucleotides are: ",paste(supported_nucleotides,collapse=", "),".")
  }
  
  # Throw an error if the reverse primer is not a character string.
  if(!is.character(reverse_primer)) stop("The reverse primer must be a character string.")
  # Throw an error if the length of the reverse primer character string is not one.
  if(length(reverse_primer)!=1) stop("The reverse primer must be a character string of length 1.")
  # Split the reverse primer into individual nucleotides.
  rev_split<-strsplit(x=reverse_primer,split="")[[1]]
  # If any nucleotides in the reverse primer are not supported.
  if(!all(rev_split %in% supported_nucleotides)){
    # Get the unsupported nucleotides from the reverse primer.
    rev_unsupported<-sort(unique(rev_split[!(rev_split %in% supported_nucleotides)]))
    # Throw an error listing the unsupported nucleotides from the reverse primer.
    stop(paste0("The reverse primer contains the following unsupported nucleotides: ",paste(rev_unsupported,collapse=", ")),". Supported nucleotides are: ",paste(supported_nucleotides,collapse=", "),".")
  }
  
  # Get the reverse complement of the reverse primer.
  reverse_primer_reverse_complement<-reverse_complement(sequence=reverse_primer)
  
  # Substitute wildcard characters in the forward primer
  # with their respective nucleotides.
  forward_primer<-substitute_wildcards(sequence=forward_primer)
  
  # Substitute wildcard characters in the reverse complement of the reverse primer
  # with their respective nucleotides.
  reverse_primer_reverse_complement<-substitute_wildcards(sequence=reverse_primer_reverse_complement)
  
  # Define function for trimming a single sequence.
  isolate_amp<-function(sequence,forward,reverse_complement_of_reverse){
    # If the forward primer is found in the sequence.
    if(grepl(pattern=forward,x=sequence)){
      # Trim up to and including the forward primer.
      trimmed<-gsub(pattern=paste0("^.*",forward),replacement="",x=sequence)
      # If the reverse complement of the reverse primer is found in the sequence.
      if(grepl(pattern=reverse_complement_of_reverse,x=trimmed)){
        # Trim behind and including the reverse complement of the reverse primer.
        trimmed<-gsub(pattern=paste0(reverse_complement_of_reverse,".*$"),replacement="",x=trimmed)
      } else { # If the reverse complement of the reverse primer is not found in the sequence.
        # Set the trimmed sequence to NA.
        trimmed<-NA
      }
    } else { # If the forward primer is not found in the sequence.
      # Set the trimmed sequence to NA.
      trimmed<-NA
    }
    # Return the trimmed sequence.
    return(trimmed)
  }
  
  # Trim sequences.
  isolated<-sapply(X=sequences,
                   FUN=isolate_amp,
                   forward=forward_primer,
                   reverse_complement_of_reverse=reverse_primer_reverse_complement,
                   USE.NAMES=FALSE)
  
  # Return trimmed sequences.
  return(isolated)
  
}