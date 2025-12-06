#' Get NCBI Taxonomies from Species Binomials
#'
#' @description Remotely fetches taxonomies from the NCBI taxonomy database for a list of species binomials. Installation of the `taxize` package is required to use this function.
#' @param path_to_species_binomials String specifying path to input species list with common and scientific names. The file should be in CSV format and contain the following fields: 'Common_Name', 'Scientific_Name'. Values in the 'Common_Name' field are optional. Values in the 'Scientific_Name' field are required.
#' @param path_to_output_file String specifying path to output species list with added NCBI taxonomies. The output file will be in CSV format.
#' @param path_to_taxonomy_edits String specifying path to taxonomy edits file in CSV format. The file must contain the following fields: 'Old_Taxonomy', 'New_Taxonomy', 'Notes'. Old taxonomies are replaced with new taxonomies in the order the records appear in the file. The taxonomic levels in the 'Old_Taxonomy' and 'New_Taxonomy' fields should be delimited by a semi-colon. If no taxonomy edits are desired, then set this variable to `NA` (the default).
#' @param print_queries Logical. Whether taxa queries should be printed. The default is `TRUE`.
#' @param ... Accepts former argument names for backwards compatibility.
#' @returns No return value. Writes an output CSV file with added taxonomies. Species which could not be found in the NCBI taxonomy database appear in the top records of the output file.
#' @seealso
#' [`get_taxonomies.IUCN`][get_taxonomies.IUCN()] for formatting taxonomies from the IUCN Red List. \cr \cr
#' [`adjust_taxonomies`][adjust_taxonomies()] for adjusting a taxonomy system.
#' @examplesIf interactive()
#' # Get path to example input species binomials CSV file.
#' path_to_species_binomials<-system.file("extdata",
#'                                        "example_species_binomials.csv",
#'                                        package="LocaTT",
#'                                        mustWork=TRUE)
#' 
#' # Create a temporary file path for the output CSV file.
#' path_to_output_file<-tempfile(fileext=".csv")
#' 
#' # Fetch taxonomies from species binomials.
#' get_taxonomies.species_binomials(path_to_species_binomials=path_to_species_binomials,
#'                                  path_to_output_file=path_to_output_file,
#'                                  print_queries=FALSE)
#' @export
get_taxonomies.species_binomials<-function(path_to_species_binomials,path_to_output_file,path_to_taxonomy_edits=NA,print_queries=TRUE,...){
  
  ### Check dependencies.
  
  # Ensure that the taxize package is installed.
  if(!requireNamespace(package="taxize",quietly=TRUE)){
    stop("Please install package 'taxize' to use this function.")
  }
  
  ### Ensure backwards compatibility.
  
  # Handle changes in argument names.
  ## Define data frame relating former to current argument names.
  back.compat<-data.frame(arg.former=c("path_to_input_species_binomials",
                                       "path_to_output_local_taxa_list",
                                       "print_taxize_queries"),
                          arg.current=c("path_to_species_binomials",
                                        "path_to_output_file",
                                        "print_queries"),
                          stringsAsFactors=FALSE)
  ## Collect extra arguments.
  dots<-list(...)
  ## If there are extra arguments.
  if(!missing(...)){
    ## Loop through each extra argument.
    for(i in 1:length(dots)){
      ## Get the extra argument.
      dot<-dots[i]
      ## If the extra argument is a former argument.
      if(names(dot) %in% back.compat$arg.former){
        ## Get the argument translation record.
        arg<-back.compat[back.compat$arg.former==names(dot),]
        ## Throw an error if the former argument name matches multiple current argument names.
        if(nrow(arg) > 1) stop("Former argument name matches multiple current argument names.")
        ## If the current argument is missing.
        if(do.call(missing,list(arg$arg.current))){
          ## Pass the former argument to the current argument.
          assign(x=arg$arg.current,value=dot[[1]])
          ## Provide a warning that the former argument has been
          ## renamed to the current argument.
          warning(paste0("Former argument '",names(dot),
                         "' has been renamed to current argument '",
                         arg$arg.current,"'."))
        } else { ## If the current argument is not missing.
          ## If both old and current arguments are provided, throw an error.
          stop(paste0("Multiple equivalent arguments found for current argument '",
                      arg$arg.current,"' (former argument '",arg$arg.former,"')."))
        }
      } else { # If the extra argument is not a former argument.
        # Throw an error for the unused argument.
        stop(paste0("unused argument (",names(dot)," = ",unname(unlist(dot)),")"))
      }
    }
  }
  
  ### Check arguments.
  
  # Path to species binomials.
  ## Throw an error if the path to species binomials is not a character string.
  if(!is.character(path_to_species_binomials)) stop("The path to species binomials must be a character string.")
  ## Throw an error if the path to species binomials has multiple elements.
  if(length(path_to_species_binomials) > 1) stop("The path to species binomials cannot have multiple elements.")
  ## Throw an error if the path to species binomials is NA.
  if(is.na(path_to_species_binomials)) stop("The path to species binomials cannot be NA.")
  
  # Path to output file.
  ## Throw an error if the path to output file is not a character string.
  if(!is.character(path_to_output_file)) stop("The path to output file must be a character string.")
  ## Throw an error if the path to output file has multiple elements.
  if(length(path_to_output_file) > 1) stop("The path to output file cannot have multiple elements.")
  ## Throw an error if the path to output file is NA.
  if(is.na(path_to_output_file)) stop("The path to output file cannot be NA.")
  
  # Path to taxonomy edits.
  ## Throw an error if the path to taxonomy edits is not a vector.
  if(!is.vector(path_to_taxonomy_edits)) stop("The path to taxonomy edits must be a vector.")
  ## Throw an error if the path to taxonomy edits has multiple elements.
  if(length(path_to_taxonomy_edits) > 1) stop("The path to taxonomy edits cannot have multiple elements.")
  ## Throw an error if the path to taxonomy edits is not a character string, and not NA.
  if(!(is.na(path_to_taxonomy_edits) | is.character(path_to_taxonomy_edits))) stop("The path to taxonomy edits must be a character string or NA.")
  
  # Print queries argument.
  ## Throw an error if the print queries argument is not logical.
  if(!is.logical(print_queries)) stop("The print queries argument must be logical.")
  ## Throw an error if the print queries argument has multiple elements.
  if(length(print_queries) > 1) stop("The print queries argument cannot have multiple elements.")
  ## Throw an error if the print queries argument is NA.
  if(is.na(print_queries)) stop("The print queries argument cannot be NA.")
  
  ### Begin operations.
  
  # Read in input csv file.
  taxa<-utils::read.csv(file=path_to_species_binomials,stringsAsFactors=FALSE)
  
  # Check that field names are right.
  if(!identical(colnames(taxa),c("Common_Name","Scientific_Name"))) stop('Fields in the input csv file should be "Common_Name" and "Scientific_Name".')
  
  # Throw an error if any blanks or NAs exist in the species binomial field
  # of the input csv file.
  if(any(taxa$Scientific_Name=="" | is.na(taxa$Scientific_Name))) stop("There are blanks or NAs in the species binomials field of the input csv file.")
  
  # Check that there are no duplicates in the input csv file.
  if(any(duplicated(taxa$Scientific_Name))) stop("There are duplicated species binomials in the input csv file.")
  
  # Get NCBI taxonomies from the scientific names.
  # Species synonyms are accounted for (if the synonyms are present in NCBI),
  # but mispellings are not.
  taxonomies<-taxize::tax_name(sci=taxa$Scientific_Name,get=c("domain","phylum","class","order","family","genus","species"),db="ncbi",messages=print_queries)
  
  # Add common name to the taxonomies.
  taxonomies$Common_Name<-taxa$Common_Name
  
  # Subset to just desired fields.
  taxa<-taxonomies[,c("Common_Name","query","domain","phylum","class","order","family","genus","species")]
  
  # Rename fields.
  colnames(taxa)<-c("Common_Name","Query","Domain","Phylum","Class","Order","Family","Genus","Species")
  
  # Throw an error if NCBI taxonomies were not found for any local taxa.
  if(all(is.na(taxa$Species))) stop("NCBI taxonomies were not found for any species.")
  
  # Check whether any taxa did not receive NCBI taxonomies.
  trouble_taxa_present<-any(is.na(taxa$Species))
  
  # If there are taxa which did not receive NCBI taxonomies.
  if(trouble_taxa_present){
    
    # Get taxa which did not receive NCBI taxonomies.
    trouble_taxa<-taxa[is.na(taxa$Species),]
    # Provide the query taxa names of the taxa which did not recieve NCBI
    # taxonomies in the species field.
    trouble_taxa$Species<-trouble_taxa$Query
    # Get taxa which did receive NCBI taxonomies.
    taxa<-taxa[!is.na(taxa$Species),]
    
  }
  
  # If there are taxa which received incomplete NCBI taxonomies.
  if(any(is.na(taxa[,3:8]))){
    
    # Get row indices of taxa which received incomplete NCBI taxonomies.
    partial_indices<-which(apply(X=taxa[,3:8],MARGIN=1,FUN=function(x) any(is.na(x))))
    # Get rows of taxa which received incomplete NCBI taxonomies.
    partial<-taxa[partial_indices,]
    
    # Loop through each taxa which received incomplete NCBI taxonomies.
    for(i in 1:nrow(partial)){
      
      # Get the taxon.
      partial_row<-partial[i,]
      
      # Loop through each taxonomic level above species.
      for(j in 8:3){
        
        # If the taxon is NA.
        if(is.na(partial_row[,j])){
          
          # Get the lower taxon.
          partial_lower<-partial_row[,j+1]
          # Add prefix for the current taxonomic level.
          partial_current<-paste(tolower(colnames(partial)[j]),partial_lower)
          # Set the taxon to the lower taxon with a prefix added for the current taxonomic level.
          partial_row[,j]<-partial_current
          
        }
        
      }
      
      # Add the updated row back to the partial data frame.
      partial[i,]<-partial_row
      
    }
    
    # Replace incomplete with complete NCBI taxonomies.
    taxa[partial_indices,]<-partial
    
  }
  
  # If a path to a taxonomy edit file is provided.
  if(!is.na(path_to_taxonomy_edits)){
    
    # Read in edits to reference taxonomies.
    taxonomy_edits<-utils::read.csv(file=path_to_taxonomy_edits,stringsAsFactors=FALSE)
    
    # Throw an error if the fields of the taxonomy edits file are not
    # Old_Taxonomy, New_Taxonomy, Notes.
    if(!identical(colnames(taxonomy_edits),c("Old_Taxonomy","New_Taxonomy","Notes"))) stop("The fields of the taxonomy edits file must be: 'Old_Taxonomy', 'New_Taxonomy' ,'Notes'.")
    
    # Check that there are no NAs in the taxonomy edits fields.
    if(any(is.na(taxonomy_edits[,c("Old_Taxonomy","New_Taxonomy")]) | taxonomy_edits[,c("Old_Taxonomy","New_Taxonomy")]=="")) stop("There are NAs or blanks in the 'Old_Taxonomy' or 'New_Taxonomy' fields of the taxonomy edits file. Please ensure that these fields have entries for all records.")
    
    # Check that there are no spaces in the taxonomy edits fields.
    if(any(t(apply(X=taxonomy_edits[,c("Old_Taxonomy","New_Taxonomy")],MARGIN=1,FUN=grepl,pattern=" ")))) stop("There cannot be spaces in the 'Old_Taxonomy' or 'New_Taxonomy' fields of the taxonomy edits file.")
    
    # Add a carrot to anchor the start of the old taxonomies field.
    taxonomy_edits$Old_Taxonomy<-paste0("^",taxonomy_edits$Old_Taxonomy)
    
    # Collapse local taxa taxonomies by semi-colons.
    taxa_names<-apply(X=taxa[,c("Domain","Phylum","Class","Order","Family","Genus","Species")],MARGIN=1,FUN=paste,collapse=";")
    
    # Replace spaces with underscores in collapsed local taxa taxonomies.
    taxa_names<-gsub(pattern=" ",replacement="_",x=taxa_names)
    
    # Loop through each taxonomy edit.
    for(i in 1:nrow(taxonomy_edits)){
      
      # Get the taxonomy edit.
      taxonomy_edit<-taxonomy_edits[i,]
      # Apply translation to local taxa taxonomies.
      taxa_names<-sub(pattern=taxonomy_edit$Old_Taxonomy,
                      replacement=taxonomy_edit$New_Taxonomy,
                      x=taxa_names)
      
    }
    
    # Replace underscores with spaces in collapsed local taxa taxonomies.
    taxa_names<-gsub(pattern="_",replacement=" ",x=taxa_names)
    
    # Get the number of taxonomic levels each taxonomy has.
    check_num_taxonomic_levels<-sapply(X=strsplit(x=taxa_names,split=";"),FUN=length)
    
    # If any taxonomies lack 7 levels.
    if(any(check_num_taxonomic_levels!=7)){
      
      # Get the taxonomies which lack 7 levels.
      incorrect_num_taxonomic_levels_taxonomies<-taxa_names[check_num_taxonomic_levels!=7]
      # Get the common names for the taxonomies which lack 7 levels.
      incorrect_num_taxonomic_levels_common_names<-taxa$Common_Name[check_num_taxonomic_levels!=7]
      # Create a field containing information on both the sequence name
      # and nucleotide sequence.
      incorrect_num_taxonomic_levels_message<-paste0(incorrect_num_taxonomic_levels_common_names," (",incorrect_num_taxonomic_levels_taxonomies,")")
      # Throw a error message mentioning the taxonomies which lack 7-levels.
      stop(paste0("Not all taxonomies have 7 taxonomic levels. If applying taxonomy edits, please ensure that these edits preserve the 7-level taxonomy. The following taxonomies lack 7 levels: ",paste(incorrect_num_taxonomic_levels_message,collapse=", ")))
      
    }
    
    # Get updated taxonomies by splitting the character strings by semi-colons.
    updated_taxonomies<-as.data.frame(do.call(rbind,strsplit(x=taxa_names,split=";")))
    
    # Update the taxonomies of local taxa.
    taxa[,c("Domain","Phylum","Class","Order","Family","Genus","Species")]<-updated_taxonomies
    
  }
  
  # If there are taxa which did not receive NCBI taxonomies.
  if(trouble_taxa_present){
    
    # Put taxa which did and did not receive NCBI taxonomies back together,
    # with the taxa which did not receive NCBI taxonomies at the top of the list.
    taxa<-rbind(trouble_taxa,taxa)
    
    # Issue a warning that some species could not be found in the NCBI taxonomy database.
    warning(paste0("NCBI taxonomies could not be found for ",nrow(trouble_taxa)," species."))
    
  }
  
  # Remove the query field.
  taxa<-taxa[,colnames(taxa)!="Query"]
  
  # Replace NA cells with blanks.
  taxa[is.na(taxa)]<-""
  
  # Write out taxonomies.
  utils::write.csv(x=taxa,file=path_to_output_file,row.names=FALSE)
  
}