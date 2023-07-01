#' Decode DNA Sequence Quality Scores
#'
#' Decodes Phred quality scores in Sanger format from symbols to numeric values.
#' @param symbols A string containing quality scores encoded as symbols in Sanger format.
#' @returns A numeric vector of Phred quality scores.
#' @examples
#' decode_quality_scores(symbols="989!.C;F@\"")
#' @export
decode_quality_scores<-function(symbols){
  # Throw an error if symbols is not a character string.
  if(!is.character(symbols)) stop("Symbols must be a character string.")
  # Throw an error if the length of the symbols character string is not one.
  if(length(symbols)!=1) stop("Symbols must be a character string of length 1.")
  # Decode quality scores from symbols to numeric values.
  scores_numeric<-as.numeric(charToRaw(x=symbols))-33
  # Return the numeric quality scores.
  return(scores_numeric)
}