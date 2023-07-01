#' Read FASTQ Files
#'
#' Reads FASTQ files. Does not support the reading of FASTQ files with sequences or quality scores wrapping multiple lines.
#' @param file A string specifying the path to a FASTQ file to read.
#' @returns A data frame with fields for sequence names, sequences, comments, and quality scores.
#' @examples
#' # Get path to example FASTQ file.
#' path_to_fastq_file<-system.file("extdata",
#'                                 "example_query_sequences.fastq",
#'                                 package="LocaTT",
#'                                 mustWork=TRUE)
#' 
#' # Read the example FASTQ file.
#' read.fastq(file=path_to_fastq_file)
#' @export
read.fastq<-function(file){
  # Read in the lines of the fastq file.
  fastq<-readLines(con=file)
  # Get sequence names.
  names<-sub(pattern="^@",replacement="",x=fastq[seq(from=1,to=length(fastq),by=4)])
  # Get sequences.
  sequences<-fastq[seq(from=2,to=length(fastq),by=4)]
  # Get comments.
  comments<-sub(pattern="^\\+",replacement="",x=fastq[seq(from=3,to=length(fastq),by=4)])
  # Get quality scores.
  quality_scores<-fastq[seq(from=4,to=length(fastq),by=4)]
  # Create data frame with names, sequences, and quality scores.
  df<-data.frame(Name=names,
                 Sequence=sequences,
                 Comment=comments,
                 Quality_scores=quality_scores,
                 stringsAsFactors=FALSE)
  # Return the data frame.
  return(df)
}