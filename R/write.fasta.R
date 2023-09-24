#' Write FASTA Files
#'
#' Writes FASTA files.
#' @param names A character vector of sequence names.
#' @param sequences A character vector of sequences.
#' @param file A string specifying the path to a FASTA file to write.
#' @seealso
#' [`read.fasta`][read.fasta()] for reading FASTA files. \cr
#' [`write.fastq`][write.fastq()] for writing FASTQ files. \cr
#' [`read.fastq`][read.fastq()] for reading FASTQ files.
#' @examples
#' # Get path to example sequences CSV file.
#' path_to_CSV_file<-system.file("extdata",
#'                               "example_query_sequences.csv",
#'                               package="LocaTT",
#'                               mustWork=TRUE)
#' 
#' # Read the example sequences CSV file.
#' df<-read.csv(file=path_to_CSV_file,stringsAsFactors=FALSE)
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
  # Check that sequence names is a character vector.
  if(!is.character(names)) stop("Names must be a character vector.")
  # Check that sequences is a character vector.
  if(!is.character(sequences)) stop("Sequences must be a character vector.")
  # Check that names and sequences are of the same length.
  if(length(names)!=length(sequences)) stop("Names and sequences are not of the same length.")
  # Append '>' to the start of sequence names.
  names<-paste0(">",names)
  # Create an empty vector to store what will be written to a text file.
  fasta<-vector(mode="character",length=2*length(names))
  # Add names to the text file vector.
  fasta[seq(from=1,to=length(fasta),by=2)]<-names
  # Add sequences to the text file vector.
  fasta[seq(from=2,to=length(fasta),by=2)]<-sequences
  # Write out fasta file.
  writeLines(text=fasta,con=file)
}