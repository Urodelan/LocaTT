#' Write FASTA Files
#'
#' Writes FASTA files.
#' @param names A character vector of sequence names.
#' @param sequences A character vector of sequences.
#' @param file A string specifying the path to a FASTA file (with '.fasta' extension) to write.
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
  write.table(x=df,file=file,row.names=F,col.names=F,quote=F)
}