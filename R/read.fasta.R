#' Read FASTA Files
#'
#' Reads FASTA files.
#' @param file A string specifying the path to a FASTA file (with '.fasta' extension) to read.
#' @returns A data frame with fields for sequence names and sequences.
#' @examples
#' # Get path to example FASTA file.
#' path_to_fasta_file<-system.file("extdata",
#'                                 "example_query_sequences.fasta",
#'                                 package="LocaTT",
#'                                 mustWork=TRUE)
#' 
#' # Read the example FASTA file.
#' read.fasta(file=path_to_fasta_file)
#' @export
read.fasta<-function(file){
  # Read in fasta file.
  fasta<-utils::read.delim(file=file,header=FALSE,stringsAsFactors=FALSE)
  # Get the indices of rows which contain sequence header information.
  indices_of_sequence_headers<-which(grepl(pattern="^>",x=fasta$V1))
  # Get the number of rows associated with each sequence header within the file.
  number_of_row_repeats_per_sequence<-diff(x=c(indices_of_sequence_headers,
                                               nrow(fasta)+1),lag=1)
  # Create a vector of integer values representing which records belong to each
  # sequence within a sample.
  sequence_split_groups<-rep(x=1:length(indices_of_sequence_headers),
                             times=number_of_row_repeats_per_sequence)
  # Get the sequence segments associated with each sequence.
  sequences_split<-split(x=fasta$V1,f=sequence_split_groups)
  # Add the first element information (sequence header) as a name to the lists.
  names(sequences_split)<-sub(pattern="^>",replacement="",
                              x=sapply(X=sequences_split,FUN="[[",1))
  # Remove the first element information (sequence header) from the lists.
  sequences_split<-sapply(X=sequences_split,FUN="[",-1,simplify=F)
  # Collapse the sequence segments together.
  sequences<-sapply(X=sequences_split,FUN=paste,collapse="")
  # Create data frame with names and sequences.
  df<-data.frame(Name=names(sequences),Sequence=unname(sequences),stringsAsFactors=FALSE)
  # Return the data frame.
  return(df)
}