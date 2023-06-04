#' Write FASTA Files
#'
#' Writes FASTA files.
#' @param names A character vector of sequence names.
#' @param sequences A character vector of sequences.
#' @param file A string specifying the path to a FASTA file (with '.fasta' extension) to write.
#' @examples
#' # Get path to example sequences CSV file.
#' path_to_CSV_file<-system.file("extdata",
#'                               "example_query_sequences.csv",
#'                               package="LocaTT",
#'                               mustWork=TRUE)
#' 
#' # Read the example sequences CSV file.
#' df<-read.csv(file=path_to_CSV_file)
#' 
#' # Create a temporary file path for the FASTA file to write.
#' path_to_FASTA_file<-tempfile(fileext=".fasta")
#' 
#' # Write the example sequences as a FASTA file.
#' write.fasta(names=df$Name,
#'             sequences=df$Sequence,
#'             file=path_to_FASTA_file)
#' @export
write.fasta<-function(names,sequences,file){
  # Check that names and sequences are of the same length.
  if(length(names)!=length(sequences)) stop("Names and sequences are not of the same length.")
  # Append '>' to the start of sequence names.
  names<-paste0(">",names)
  # Create an empty data frame to store what will be written to a text file.
  df<-as.data.frame(matrix(NA,nrow=2*length(names),ncol=1))
  # Add names to the text file data frame.
  df$V1[seq(from=1,to=nrow(df),by=2)]<-names
  # Add sequences to the text file data frame.
  df$V1[seq(from=2,to=nrow(df),by=2)]<-as.character(sequences)
  # Write out fasta file.
  utils::write.table(x=df,file=file,row.names=FALSE,col.names=FALSE,quote=FALSE)
}