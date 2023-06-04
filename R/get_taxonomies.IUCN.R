#' Get Taxonomies from IUCN Red List Files
#' 
#' Formats taxonomies from IUCN Red List taxonomy.csv and common_names.csv files for use with the `local_taxa_tool` function.
#' @param path_to_IUCN_taxonomies String specifying path to input IUCN Red List taxonomy.csv file.
#' @param path_to_IUCN_common_names String specifying path to input IUCN Red List common_names.csv file.
#' @param path_to_output_local_taxa_list String specifying path to output species list (in CSV format with '.csv' extension) with formatted taxonomies.
#' @param domain_name String specifying the domain name to use for all species. The IUCN Red List files do not include domain information, so a domain name must be provided. If using a reference database from UNITE, provide a kingdom name here (e.g., `'Fungi'`). The default is `'Eukaryota'`.
#' @param path_to_taxonomy_edits String specifying path to taxonomy edits file in CSV format (with '.csv' extension). The file must contain the following fields: 'Old_Taxonomy', 'New_Taxonomy', 'Notes'. Old taxonomies are replaced with new taxonomies in the order the records appear in the file. The taxonomic levels in the 'Old_Taxonomy' and 'New_Taxonomy' fields should be delimited by a semi-colon. If no taxonomy edits are desired, then set this variable to `NA` (the default).
#' @examples
#' # Get path to example taxonomy CSV file.
#' path_to_taxonomy_file<-system.file("extdata",
#'                                    "example_taxonomy.csv",
#'                                    package="LocaTT",
#'                                    mustWork=TRUE)
#' 
#' # Get path to example common names CSV file.
#' path_to_common_names_file<-system.file("extdata",
#'                                        "example_common_names.csv",
#'                                        package="LocaTT",
#'                                        mustWork=TRUE)
#' 
#' # Create a temporary file path for the output CSV file.
#' path_to_output_file<-tempfile(fileext=".csv")
#' 
#' # Format common names and taxonomies.
#' get_taxonomies.IUCN(path_to_IUCN_taxonomies=path_to_taxonomy_file,
#'                     path_to_IUCN_common_names=path_to_common_names_file,
#'                     path_to_output_local_taxa_list=path_to_output_file)
#' @export
get_taxonomies.IUCN<-function(path_to_IUCN_taxonomies,path_to_IUCN_common_names,path_to_output_local_taxa_list,domain_name="Eukaryota",path_to_taxonomy_edits=NA){
  
  # Read in IUCN taxonomies.
  taxa<-utils::read.csv(file=path_to_IUCN_taxonomies,stringsAsFactors=FALSE)
  
  # Read in IUCN common names.
  common<-utils::read.csv(file=path_to_IUCN_common_names,stringsAsFactors=FALSE)
  
  # Get the main common names.
  common<-common[common$main=="true",]
  
  # Remove duplicate common names for each species.
  common<-common[!duplicated(common$internalTaxonId),]
  
  # Adding common names to taxonomies.
  taxa<-merge(x=taxa,y=common,all.x=TRUE,all.y=FALSE)
  
  # Add a field for domain.
  taxa$domainName<-domain_name
  
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
      stop(paste0("Not all taxonomies for local taxa have 7 taxonomic levels. If applying taxonomy edits, please ensure that these edits preserve the 7-level taxonomy. The following taxonomies lack 7 levels: ",paste(incorrect_num_taxonomic_levels_message,collapse=", ")))
      
    }
    
    # Get updated taxonomies by splitting the character strings by semi-colons.
    updated_taxonomies<-as.data.frame(do.call(rbind,strsplit(x=taxa_names,split=";")))
    
    # Update the taxonomies of local taxa.
    taxa[,c("Domain","Phylum","Class","Order","Family","Genus","Species")]<-updated_taxonomies
    
  }
  
  # Change NAs in common name field to blanks.
  taxa$Common_Name<-ifelse(is.na(taxa$Common_Name),"",taxa$Common_Name)
  
  # Write formatted IUCN taxonomies.
  utils::write.csv(x=taxa,file=path_to_output_local_taxa_list,row.names=FALSE)
  
}