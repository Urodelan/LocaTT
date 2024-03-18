#' Get Specified Taxonomic Level from Taxonomic Strings
#'
#' Gets the specified taxonomic level from a vector of taxonomic strings.
#' @param taxonomies A character vector of taxonomic strings.
#' @param level A numeric value representing the taxonomic level to be extracted. A value of `1` retrieves the highest taxonomic level (*e.g.*, domain) from the input taxonomies, with each sequentially higher value retrieving sequentially lower taxonomic levels. `0` is a special value which retrieves the lowest taxonomic level available in the input taxonomies.
#' @param full_names Logical. If `TRUE` (the default), then full taxonomies are returned down to the requested taxonomic level. If `FALSE`, then only the requested taxonomic level is returned.
#' @param delimiter A character string of the delimiter between taxonomic levels in the input taxonomies. The default is `";"`.
#' @returns A character vector containing the requested taxonomic level for each element of the input taxonomies.
#' @seealso
#' [`expand_taxonomies`][expand_taxonomies()] for extracting each taxonomic level from a vector of taxonomic strings. \cr \cr
#' [`get_consensus_taxonomy`][get_consensus_taxonomy()] for generating a consensus taxonomy from taxonomic strings.
#' @examples
#' get_taxonomic_level(taxonomies=
#'    c("Eukaryota;Chordata;Amphibia;Caudata;Ambystomatidae;Ambystoma;Ambystoma_mavortium",
#'      "Eukaryota;Chordata;Amphibia;Anura;Bufonidae;Anaxyrus;Anaxyrus_boreas",
#'      "Eukaryota;Chordata;Amphibia;Anura;Ranidae;Rana;Rana_luteiventris"),
#'    level=5,
#'    full_names=TRUE,
#'    delimiter=";")
#' @export
get_taxonomic_level<-function(taxonomies,level,full_names=TRUE,delimiter=";"){
  
  # Throw an error if taxonomies is not a character vector.
  if(!is.character(taxonomies)) stop("Taxonomies must be a character vector.")
  
  # Remove names from taxonomies vector, if present.
  taxonomies<-unname(taxonomies)
  
  # Throw an error if level is not numeric.
  if(!is.numeric(level)) stop("Level must be numeric.")
  
  # Throw an error if level has multiple elements.
  if(length(level) > 1) stop("Level cannot have multiple elements.")
  
  # Throw an error if level is not an integer.
  if(round(x=level,digits=0)!=level) stop("Level must be an integer.")
  
  # Throw an error if level is less than zero.
  if(level < 0) stop("Level must be at least zero.")
  
  # Throw an error if the full names argument is not logical.
  if(!is.logical(full_names)) stop("The full names argument must be logical.")
  
  # Throw an error if the full names argument has multiple elements.
  if(length(full_names) > 1) stop("The full names argument cannot have multiple elements.")
  
  # Throw an error if the full names argument is not TRUE or FALSE.
  if(!(full_names %in% c(TRUE,FALSE))) stop("The full names argument must be TRUE or FALSE.")
  
  # Throw an error if the delimiter is not a character string.
  if(!is.character(delimiter)) stop("The delimiter must be a character string.")
  
  # Throw an error if the full names argument has multiple elements.
  if(length(delimiter) > 1) stop("The delimiter cannot have multiple elements.")
  
  # Split taxonomies using the specified delimiter.
  split_taxonomies<-strsplit(x=taxonomies,split=delimiter)
  
  # Get the number of levels included in the taxonomies.
  length_taxonomies<-sapply(X=split_taxonomies,FUN=length)
  
  # Throw an error if not all taxonomies have the same number of levels.
  if(!all(length_taxonomies==length_taxonomies[1])) stop("Not all taxonomies have the same number of levels.")
  
  # Throw an error if level exceeds the number of taxonomic
  # levels in the elements of the taxonomies vector.
  if(level > length_taxonomies[1]) stop("Level exceeds the number of taxonomic levels in the elements of the taxonomies vector.")
  
  # If the level argument is set to zero.
  if(level==0){
    
    # If full names are requested.
    if(full_names){
      
      # Return the full taxonomies.
      taxonomic_level<-taxonomies
      
    } else { # If full names are not requested.
      
      # Get the finest taxonomic level without the higher levels.
      taxonomic_level<-sapply(X=split_taxonomies,FUN="[[",length_taxonomies[1])
      
    }
    
  } else { # If the level argument is not set to zero.
    
    # If full names are requested.
    if(full_names){
      
      # Define internal function for collapsing split
      # taxonomies at the desired taxonomic level.
      collapse_split_taxonomies<-function(split_name,lev,delim){
        paste(split_name[1:lev],collapse=delim)
      }
      
      # Apply the collapse taxonomies function
      # across the list of split taxonomies.
      taxonomic_level<-sapply(X=split_taxonomies,
                              FUN=collapse_split_taxonomies,
                              lev=level,
                              delim=delimiter)
      
    } else { # If full names are not requested.
      
      # Get the desired taxonomic level without the higher levels.
      taxonomic_level<-sapply(X=split_taxonomies,FUN="[[",level)
      
    }
    
  }
  
  # Return the requested taxonomic level.
  return(taxonomic_level)
  
}