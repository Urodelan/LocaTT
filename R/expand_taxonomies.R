#' Expand Taxonomies
#'
#' Extracts each taxonomic level from a vector of taxonomic strings.
#' @param taxonomies A character vector of taxonomic strings.
#' @param levels A character vector of taxonomic level names. The length of `levels` determines the number of taxonomic levels to extract from the taxonomies, and the order of elements in `levels` is assumed to match the order of taxonomic levels in the taxonomies. `levels` is also used as the field names of the returned data frame (see return value section). The default vector includes: `"Domain"`, `"Phylum"`, `"Class"`, `"Order"`, `"Family"`, `"Genus"`, and `"Species"`.
#' @param full_names Logical. If `TRUE` (the default), then full taxonomies are returned down to the extracted taxonomic level. If `FALSE`, then only the extracted taxonomic level is returned.
#' @param delimiter A character string of the delimiter between taxonomic levels in the input taxonomies. The default is `";"`.
#' @param ignore An optional character vector of taxonomic strings for which taxonomic expansion will be skipped. In the returned data frame (see return value section), the record for each skipped taxonomic string will be filled with `NA`s.
#' @returns Returns a data frame of extracted taxonomic levels. One record for each element of `taxonomies`, and one field for each element of `levels`. Field names are inherited from `levels`. If a taxonomic level is not present in a taxonomic string, then the respective cell in the returned data frame will contain `NA`.
#' @seealso
#' [`get_taxonomic_level`][get_taxonomic_level()] for extracting a taxonomic level from taxonomic strings. \cr \cr
#' [`get_consensus_taxonomy`][get_consensus_taxonomy()] for generating a consensus taxonomy from taxonomic strings.
#' @examples
#' expand_taxonomies(taxonomies=
#'    c("Eukaryota;Chordata;Amphibia;Caudata;Ambystomatidae;Ambystoma;Ambystoma mavortium",
#'      "Eukaryota;Chordata;Amphibia;Anura;Bufonidae;Anaxyrus;Anaxyrus boreas",
#'      "Eukaryota;Chordata;Amphibia;Anura;Ranidae;Rana;Rana luteiventris"),
#'                        full_names=FALSE,
#'                        delimiter=";")
#' @export
expand_taxonomies<-function(taxonomies,levels=c("Domain","Phylum","Class","Order","Family","Genus","Species"),full_names=TRUE,delimiter=";",ignore){
  
  # Check arguments.
  
  # Taxonomies.
  ## Throw an error if the input taxonomies are not a character vector.
  if(!is.character(taxonomies)) stop("The input taxonomies must be a character vector.")
  
  # Taxonomic levels.
  ## Throw an error if the input taxonomic levels are not a character vector.
  if(!is.character(levels)) stop("The levels argument must be a character vector.")
  ## Throw an error if any input taxonomic levels are NA.
  if(any(is.na(levels))) stop("There cannot be NAs in the levels argument.")
  
  # Full names.
  ## Throw an error if the full names argument is not logical.
  if(!is.logical(full_names)) stop("The full names argument must be logical.")
  ## Throw an error if the full names argument has multiple elements.
  if(length(full_names) > 1) stop("The full names argument cannot have multiple elements.")
  ## Throw an error if the full names argument is NA.
  if(is.na(full_names)) stop("The full names argument cannot be NA.")
  
  # Taxonomic delimiter.
  ## Throw an error if the taxonomic delimiter is not a character string.
  if(!is.character(delimiter)) stop("The taxonomic delimiter must be a character string.")
  ## Throw an error if the taxonomic delimiter has multiple elements.
  if(length(delimiter) > 1) stop("The taxonomic delimiter cannot have multiple elements.")
  ## Throw an error if the taxonomic delimiter is NA.
  if(is.na(delimiter)) stop("The taxonomic delimiter cannot be NA.")
  
  # Ignore.
  ## If the ignore argument is provided.
  if(!missing(ignore)){
    ## Throw an error if the ignore argument is not a character vector.
    if(!is.character(ignore)) stop("The ignore argument must be a character vector.")
  }
  
  # Begin operations.
  
  # Get the number of taxonomic levels.
  n_levels<-length(levels)
  
  # Create a data frame for storing levels of the expanded taxonomies.
  df<-as.data.frame(matrix(data=NA,nrow=length(taxonomies),ncol=n_levels))
  
  # Add field names to the taxonomy data frame.
  colnames(df)<-levels
  
  # Loop through each input taxonomy.
  for(i in 1:length(taxonomies)){
    
    # Get element of the input taxonomy.
    taxonomy<-taxonomies[i]
    
    # Set expand flag to TRUE.
    expand<-TRUE
    
    # If there are any taxonomies to ignore.
    if(!missing(ignore)){
      # If the taxonomy matches any elements of ignore.
      if(taxonomy %in% ignore){
        # Change the expand flag to FALSE.
        expand<-FALSE
      }
    }
    
    # If the expand flag is TRUE.
    if(expand){
      # Loop through each taxonomic level.
      for(taxonomic_level in 1:n_levels){
        # If the taxonomy contains the desired taxonomic level.
        if(length(strsplit(x=taxonomy,split=delimiter)[[1]]) >= taxonomic_level){
          # Get the desired taxonomic level from the consensus taxonomy.
          df[i,taxonomic_level]<-get_taxonomic_level(taxonomies=taxonomy,
                                                     level=taxonomic_level,
                                                     full_names=full_names,
                                                     delimiter=delimiter)
        }
      }
    }
    
  }
  
  # Return the data frame of expanded taxonomies.
  return(df)
  
}