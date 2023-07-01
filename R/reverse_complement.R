#' Get the Reverse Complement of a DNA Sequence
#'
#' Gets the reverse complement of a DNA Sequence. Ambiguous nucleotides are supported.
#' @param sequence A string specifying the DNA sequence. Can contain ambiguous nucleotides.
#' @returns A string of the reverse complement of the DNA sequence.
#' @examples
#' reverse_complement(sequence="TTCTCCASCCGCGGATHTTG")
#' @export
reverse_complement<-function(sequence){
  # Create a data frame defining nucleotide complements.
  complements<-data.frame(Nucleotide=c("A","G","C","T","Y","R","W","S","K","M","D","V","H","B","N"),
                          Complement=c("T","C","G","A","R","Y","W","S","M","K","H","B","D","V","N"),
                          stringsAsFactors=F)
  # Throw an error if the sequence is not a character string.
  if(!is.character(sequence)) stop("The sequence must be a character string.")
  # Throw an error if the length of the sequence character string is not one.
  if(length(sequence)!=1) stop("The sequence must be a character string of length 1.")
  # Split the sequence into individual nucleotides.
  tmp<-strsplit(x=sequence,split="")[[1]]
  # If any nucleotides in the sequence are not supported.
  if(!all(tmp %in% complements$Nucleotide)){
    # Get the unsupported nucleotides.
    unsupported<-sort(unique(tmp[!(tmp %in% complements$Nucleotide)]))
    # Throw an error listing the unsupported nucleotides.
    stop(paste0("The sequence contains the following unsupported nucleotides: ",paste(unsupported,collapse=", ")),". Supported nucleotides are: ",paste(complements$Nucleotide,collapse=", "),".")
  }
  # Reverse the orientation of the sequence.
  tmp<-tmp[length(tmp):1]
  # Translate each nucleotide to its complement.
  r_cmp<-complements$Complement[match(tmp,complements$Nucleotide)]
  # Collapse the individual nucleotides into a sequence.
  r_cmp<-paste(r_cmp,collapse="")
  # Return the reverse complement.
  return(r_cmp)
}