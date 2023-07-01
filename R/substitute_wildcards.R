#' Substitute Wildcard Characters in a DNA Sequence
#'
#' Substitutes wildcard characters in a DNA sequence with their associated nucleotides surrounded by square brackets. The output is useful for matching in regular expressions.
#' @param sequence A string specifying the DNA sequence containing wildcard characters.
#' @returns A string of the DNA sequence in which wildcard characters are replaced with their associated nucleotides surrounded by square brackets.
#' @examples
#' substitute_wildcards(sequence="CAADATCCGCGGSTGGAGAA")
#' @export
substitute_wildcards<-function(sequence){
  
  # Throw an error if the sequence is not a character string.
  if(!is.character(sequence)) stop("The sequence must be a character string.")
  
  # Throw an error if the length of the sequence character string is not one.
  if(length(sequence)!=1) stop("The sequence must be a character string of length 1.")
  
  # Define IUPAC wildcard characters.
  IUPAC_wildcards<-list(R=c("A","G"),
                        Y=c("C","T"),
                        S=c("G","C"),
                        W=c("A","T"),
                        K=c("G","T"),
                        M=c("A","C"),
                        B=c("C","G","T"),
                        D=c("A","G","T"),
                        H=c("A","C","T"),
                        V=c("A","C","G"),
                        N=c("A","T","G","C"))
  
  # Create a lookup table for translating wildcard characters to their associated nucleotides.
  translation<-data.frame(Character=names(IUPAC_wildcards),
                          Substitute=paste0("[",
                                            sapply(X=IUPAC_wildcards,
                                                   FUN=paste,collapse=""),
                                            "]"),
                          stringsAsFactors=FALSE)
  
  # Create a translation table for non-ambiguous nucleotides which keeps them the same.
  nucleotides<-data.frame(Character=c("A","T","G","C"),
                          Substitute=c("A","T","G","C"),
                          stringsAsFactors=FALSE)
  
  # Add the non-ambiguous nucleotide translations to the lookup table.
  translation<-rbind(nucleotides,translation)
  
  # Split the sequence into individual nucleotides.
  seq_split<-strsplit(x=sequence,split="")[[1]]
  
  # If any nucleotides in the sequence are not supported.
  if(!all(seq_split %in% translation$Character)){
    # Get the unsupported nucleotides.
    unsupported<-sort(unique(seq_split[!(seq_split %in% translation$Character)]))
    # Throw an error listing the unsupported nucleotides.
    stop(paste0("The sequence contains the following unsupported nucleotides: ",paste(unsupported,collapse=", ")),". Supported nucleotides are: ",paste(translation$Character,collapse=", "),".")
  }
  
  # Substitute individual nucleotides in the sequence according to the lookup table.
  seq_sub<-translation$Substitute[match(seq_split,translation$Character)]
  
  # Collapse the individual nucleotides back into a sequence.
  seq_sub<-paste(seq_sub,collapse="")
  
  # Return the sequence with substituted ambiguous nucleotides.
  return(seq_sub)
  
}