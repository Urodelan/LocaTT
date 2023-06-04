#' Check BLAST Installation
#'
#' Checks whether a BLAST program can be found.
#' @param blast_command String specifying the path to a BLAST program.
#' @returns Logical. Returns `TRUE` if the BLAST program could be found.
#' @export
#' @examples
#' blast_command_found(blast_command="blastn")
blast_command_found<-function(blast_command){
  
  # Throw an error if the blast_command argument is missing.
  if(missing(blast_command)) stop("The blast_command argument is missing. Please provide this argument.")
  
  # Set the initial command_was_found value to FALSE.
  command_was_found<-FALSE
  
  # Run code even if it throws errors or warnings.
  tryCatch(
    # Multiple lines go in brackets.
    {
      # Check that the BLAST command can be found.
      blast_test<-system2(command=blast_command,args="-version",stdout=TRUE)
      # If the command was found, update command_was_found to TRUE.
      command_was_found<-TRUE
      # Return whether the command was found.
      return(command_was_found)
    },
    # Return whether the command was found.
    ## If an error was encountered.
    error=function(e){
      return(command_was_found)
    },
    ## If a warning was encountered.
    warning=function(w){
      return(command_was_found)
    }
  )
  
}