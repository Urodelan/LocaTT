#' Write FASTQ Files
#'
#' Writes FASTQ files.
#' @param names A character vector of sequence names.
#' @param sequences A character vector of sequences.
#' @param quality_scores A character vector of quality scores.
#' @param file A string specifying the path to a FASTQ file to write.
#' @param comments An optional character vector of sequence comments.
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
#' # Create a temporary file path for the FASTQ file to write.
#' path_to_FASTQ_file<-tempfile(fileext=".fastq")
#' 
#' # Write the example sequences as a FASTQ file.
#' write.fastq(names=df$Name,
#'             sequences=df$Sequence,
#'             quality_scores=df$Quality_score,
#'             file=path_to_FASTQ_file,
#'             comments=df$Comment)
#' @export
write.fastq<-function(names,sequences,quality_scores,file,comments){
  
  # Check that arguments are character vectors.
  ## Sequence names.
  if(!is.character(names)) stop("Names must be a character vector.")
  ## Sequences.
  if(!is.character(sequences)) stop("Sequences must be a character vector.")
  ## Quality scores.
  if(!is.character(quality_scores)) stop("Quality scores must be a character vector.")
  ## Comments.
  if(!missing(comments)){
    if(!is.character(comments)) stop("Comments must be a character vector.")
  }
  
  # If comments were not provided.
  if(missing(comments)){
    # Check that names, sequences, and quality scores are of the same length.
    if(!all(sapply(list(length(names),length(sequences),length(quality_scores)),
                   function(x) x==length(names)))){
      # Throw an error if names, sequences, and quality scores
      # are not of the same length.
      stop("Names, sequences, and quality scores are not of the same length.")
    }
  } else { # If comments were provided.
    # Check that names, sequences, comments, and quality scores are of the same length.
    if(!all(sapply(list(length(names),length(sequences),length(comments),length(quality_scores)),
                   function(x) x==length(names)))){
      # Throw an error if names, sequences, comments, and quality scores
      # are not of the same length.
      stop("Names, sequences, comments, and quality scores are not of the same length.")
    }
  }
  
  # Append '@' to the start of sequence names.
  names<-paste0("@",names)
  # If comments were not provided.
  if(missing(comments)){
    # Set the sequence comments to '+'.
    comments<-"+"
  } else {
    # Append '+' to the start of sequence comments.
    comments<-paste0("+",comments)
  }
  
  # Create an empty vector to store what will be written to a text file.
  fastq<-vector(mode="character",length=4*length(names))
  # Add names to the text file vector.
  fastq[seq(from=1,to=length(fastq),by=4)]<-names
  # Add sequences to the text file vector.
  fastq[seq(from=2,to=length(fastq),by=4)]<-sequences
  # Add comments to the text file vector.
  fastq[seq(from=3,to=length(fastq),by=4)]<-comments
  # Add quality scores to the text file data frame.
  fastq[seq(from=4,to=length(fastq),by=4)]<-quality_scores
  
  # Write out fastq file.
  writeLines(text=fastq,con=file)
  
}