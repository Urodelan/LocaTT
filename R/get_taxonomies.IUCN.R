#' Get Taxonomies from IUCN Red List Files
#' 
#' Formats taxonomies from IUCN Red List taxonomy.csv and common_names.csv files for use with the [`local_taxa_tool`][local_taxa_tool()] function.
#' @param path_to_taxonomies String specifying path to input IUCN Red List taxonomy.csv file.
#' @param path_to_common_names String specifying path to input IUCN Red List common_names.csv file.
#' @param path_to_output_file String specifying path to output species list (in CSV format) with formatted taxonomies.
#' @param domain String specifying the domain name to use for all species. The IUCN Red List files do not include domain information, so a domain name must be provided. If using a reference database from UNITE, provide a kingdom name here (*e.g.*, `'Fungi'`). The default is `'Eukaryota'`.
#' @param path_to_taxonomy_edits String specifying path to taxonomy edits file in CSV format. The file must contain the following fields: 'Old_Taxonomy', 'New_Taxonomy', 'Notes'. Old taxonomies are replaced with new taxonomies in the order the records appear in the file. The taxonomic levels in the 'Old_Taxonomy' and 'New_Taxonomy' fields should be delimited by a semi-colon. If no taxonomy edits are desired, then set this variable to `NA` (the default).
#' @param ... Accepts former argument names for backwards compatibility.
#' @returns No return value. Writes an output CSV file with formatted taxonomies.
#' @seealso
#' [`get_taxonomies.species_binomials`][get_taxonomies.species_binomials()] for remotely fetching NCBI taxonomies from species binomials. \cr \cr
#' [`adjust_taxonomies`][adjust_taxonomies()] for adjusting a taxonomy system.
#' @examples
#' # Get path to example taxonomy CSV file.
#' path_to_taxonomies<-system.file("extdata",
#'                                 "example_taxonomy.csv",
#'                                 package="LocaTT",
#'                                 mustWork=TRUE)
#' 
#' # Get path to example common names CSV file.
#' path_to_common_names<-system.file("extdata",
#'                                   "example_common_names.csv",
#'                                   package="LocaTT",
#'                                   mustWork=TRUE)
#' 
#' # Create a temporary file path for the output CSV file.
#' path_to_output_file<-tempfile(fileext=".csv")
#' 
#' # Format common names and taxonomies.
#' get_taxonomies.IUCN(path_to_taxonomies=path_to_taxonomies,
#'                     path_to_common_names=path_to_common_names,
#'                     path_to_output_file=path_to_output_file)
#' @export
get_taxonomies.IUCN<-function(path_to_taxonomies,path_to_common_names,path_to_output_file,domain="Eukaryota",path_to_taxonomy_edits=NA,...){
  
  ### Ensure backwards compatibility.
  
  # Handle changes in argument names.
  ## Define data frame relating former to current argument names.
  back.compat<-data.frame(arg.former=c("path_to_IUCN_taxonomies",
                                       "path_to_IUCN_common_names",
                                       "path_to_output_local_taxa_list",
                                       "domain_name"),
                          arg.current=c("path_to_taxonomies",
                                        "path_to_common_names",
                                        "path_to_output_file",
                                        "domain"),
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
  
  # Path to taxonomies.
  ## Throw an error if the path to taxonomies is not a character string.
  if(!is.character(path_to_taxonomies)) stop("The path to taxonomies must be a character string.")
  ## Throw an error if the path to taxonomies has multiple elements.
  if(length(path_to_taxonomies) > 1) stop("The path to taxonomies cannot have multiple elements.")
  ## Throw an error if the path to taxonomies is NA.
  if(is.na(path_to_taxonomies)) stop("The path to taxonomies cannot be NA.")
  
  # Path to common names.
  ## Throw an error if the path to common names is not a character string.
  if(!is.character(path_to_common_names)) stop("The path to common names must be a character string.")
  ## Throw an error if the path to common names has multiple elements.
  if(length(path_to_common_names) > 1) stop("The path to common names cannot have multiple elements.")
  ## Throw an error if the path to common names is NA.
  if(is.na(path_to_common_names)) stop("The path to common names cannot be NA.")
  
  # Path to output file.
  ## Throw an error if the path to output file is not a character string.
  if(!is.character(path_to_output_file)) stop("The path to output file must be a character string.")
  ## Throw an error if the path to output file has multiple elements.
  if(length(path_to_output_file) > 1) stop("The path to output file cannot have multiple elements.")
  ## Throw an error if the path to output file is NA.
  if(is.na(path_to_output_file)) stop("The path to output file cannot be NA.")
  
  # Domain argument.
  ## Throw an error if the domain argument is not a character string.
  if(!is.character(domain)) stop("The domain argument must be a character string.")
  ## Throw an error if the domain argument has multiple elements.
  if(length(domain) > 1) stop("The domain argument cannot have multiple elements.")
  ## Throw an error if the domain argument is NA.
  if(is.na(domain)) stop("The domain argument cannot be NA.")
  
  # Path to taxonomy edits.
  ## Throw an error if the path to taxonomy edits is not a vector.
  if(!is.vector(path_to_taxonomy_edits)) stop("The path to taxonomy edits must be a vector.")
  ## Throw an error if the path to taxonomy edits has multiple elements.
  if(length(path_to_taxonomy_edits) > 1) stop("The path to taxonomy edits cannot have multiple elements.")
  ## Throw an error if the path to taxonomy edits is not a character string, and not NA.
  if(!(is.na(path_to_taxonomy_edits) | is.character(path_to_taxonomy_edits))) stop("The path to taxonomy edits must be a character string or NA.")
  
  ### Begin operations.
  
  # Read in IUCN taxonomies.
  taxa<-utils::read.csv(file=path_to_taxonomies,stringsAsFactors=FALSE)
  
  # Read in IUCN common names.
  common<-utils::read.csv(file=path_to_common_names,stringsAsFactors=FALSE)
  
  # Get the main common names.
  common<-common[common$main=="true",]
  
  # Remove duplicate common names for each species.
  common<-common[!duplicated(common$internalTaxonId),]
  
  # Adding common names to taxonomies.
  taxa<-merge(x=taxa,y=common,all.x=TRUE,all.y=FALSE)
  
  # Add a field for domain.
  taxa$domainName<-domain
  
  # Get the common name and taxonomy fields.
  taxa<-taxa[,c("name","domainName","phylumName","className","orderName","familyName","genusName","scientificName")]
  
  # Rename the fields.
  colnames(taxa)<-c("Common_Name","Domain","Phylum","Class","Order","Family","Genus","Species")
  
  # Create a function to capitalize the first letter and the rest to lower case.
  capitalize<-function(x){
    # Conver the whole string to lower case.
    x<-tolower(x)
    # Capitalize the first letter.
    substr(x=x,start=1,stop=1)<-toupper(substr(x=x,start=1,stop=1))
    # Return the capitalized string.
    return(x)
  }
  
  # Convert to lower case and capitalize to first letter of phylum through family.
  taxa[,c("Phylum","Class","Order","Family")]<-sapply(X=taxa[,c("Phylum","Class","Order","Family")],FUN=capitalize)
  
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
  
  # Change NAs in common name field to blanks.
  taxa$Common_Name<-ifelse(is.na(taxa$Common_Name),"",taxa$Common_Name)
  
  # Write formatted IUCN taxonomies.
  utils::write.csv(x=taxa,file=path_to_output_file,row.names=FALSE)
  
}