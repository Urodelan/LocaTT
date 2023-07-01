#' Get Consensus Taxonomy from Taxonomic Strings
#'
#' Gets the consensus taxonomy from a vector of taxonomic strings.
#' @param taxonomies A character vector of taxonomic strings.
#' @param full_names Logical. If `TRUE` (the default), then the full consensus taxonomy is returned. If `FALSE`, then only the lowest taxonomic level of the consensus taxonomy is returned.
#' @param delimiter A character string of the delimiter between taxonomic levels in the input taxonomies. The default is `";"`.
#' @returns A character string containing the taxonomy agreed upon by all input taxonomies. If the input taxonomies are not the same at any taxonomic level, then `NA` is returned.
#' @examples
#' get_consensus_taxonomy(taxonomies=
#'    c("Eukaryota;Chordata;Amphibia;Caudata;Ambystomatidae;Ambystoma;Ambystoma_mavortium",
#'      "Eukaryota;Chordata;Amphibia;Anura;Bufonidae;Anaxyrus;Anaxyrus_boreas",
#'      "Eukaryota;Chordata;Amphibia;Anura;Ranidae;Rana;Rana_luteiventris"),
#'                        full_names=TRUE,
#'                        delimiter=";")
#' @export
get_consensus_taxonomy<-function(taxonomies,full_names=TRUE,delimiter=";"){
  
  # Throw an error if taxonomies is not a character vector.
  if(!is.character(taxonomies)) stop("Taxonomies must be a character vector.")
  
  # Remove names from taxonomies vector, if present.
  taxonomies<-unname(taxonomies)
  
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
  
  # Define internal function for checking whether
  # taxonomies match at a given taxonomic level.
  taxonomies_match<-function(tax,lev,delim){
    # Get the taxonomic level from the taxonomies.
    taxonomic_level<-get_taxonomic_level(taxonomies=tax,
                                         level=lev,
                                         full_names=TRUE,
                                         delimiter=delim)
    # Check if all taxonomies are the same at the taxonomic level.
    all_same<-all(taxonomic_level==taxonomic_level[1])
    # Return whether all taxonomies are the same at the taxonomic level.
    return(all_same)
  }
  
  # Get whether taxonomies match at each taxonomic
  # level available in the input taxonomies.
  matches<-sapply(X=1:length_taxonomies[1],
                  FUN=taxonomies_match,
                  tax=taxonomies,
                  delim=delimiter)
  
  # If the input taxonomies match at any taxonomic level.
  if(any(matches)){
    # Get the maximum matching taxonomic levels.
    max_level<-max(which(matches))
    # Get the consensus taxonomy.
    consensus<-get_taxonomic_level(taxonomies=taxonomies[1],
                                   level=max_level,
                                   full_names=full_names,
                                   delimiter=delimiter)
  } else { # If the input taxonomies do not match at any taxonomic level.
    # Set the consensus taxonomy to NA.
    consensus<-NA
  }
  
  # Return the consensus taxonomy.
  return(consensus)
  
}