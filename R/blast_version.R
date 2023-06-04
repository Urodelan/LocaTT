#' Get BLAST Version
#' 
#' Gets the version of a BLAST program.
#' @param blast_command String specifying the path to a BLAST program. The default (`'blastn'`) should return the version of the blastn program for standard BLAST installations. The user can provide a path to a BLAST program for non-standard BLAST installations.
#' @returns Returns a string of the version of the BLAST program.
#' @export
#' @examples
#' blast_version()
blast_version<-function(blast_command="blastn"){
  
  # Run code even if it throws errors or warnings.
  tryCatch(
    # Multiple lines go in brackets.
    {
      # Check the BLAST command version.
      blast_version<-system2(command=blast_command,args="-version",stdout=TRUE)
      # Collapse the BLAST command version vector.
      blast_version<-paste(blast_version,collapse="")
      # Return the BLAST command version.
      return(blast_version)
    },
    # Replace the default error message with a custom one.
    ## If an error was encountered.
    error=function(e){
      message(paste0("The ",blast_command," program could not be found. If using a non-standard installation of BLAST, set the path to the BLAST program using the blast_command argument."))
    },
    ## If a warning was encountered.
    warning=function(w){
      message(paste0("The ",blast_command," program could not be found. If using a non-standard installation of BLAST, set the path to the BLAST program using the blast_command argument."))
    }
  )
  
}