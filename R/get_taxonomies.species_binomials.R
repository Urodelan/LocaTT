#' Get NCBI Taxonomies from Species Binomials
#'
#' Remotely fetches taxonomies from the NCBI taxonomy database for a list of species binomials.
#' @param path_to_input_species_binomials String specifying path to input species list with common and scientific names. The file should be in CSV format with a '.csv' extension. The file must contain the following fields: 'Common_Name', 'Scientific_Name'. Values in the 'Common_Name' field are optional. Values in the 'Scientific_Name' field are required.
#' @param path_to_output_local_taxa_list String specifying path to output species list with added NCBI taxonomies. The output file will be in CSV format, and the file path should have a '.csv' extension.
#' @param path_to_taxonomy_edits String specifying path to taxonomy edits file in CSV format with a '.csv' extension. The file must contain the following fields: 'Old_Taxonomy', 'New_Taxonomy', 'Notes'. Old taxonomies are replaced with new taxonomies in the order the records appear in the file. The taxonomic levels in the 'Old_Taxonomy' and 'New_Taxonomy' fields should be delimited by a semi-colon. If no taxonomy edits are desired, then set this variable to `NA` (the default).
#' @param print_taxize_queries Logical. Whether taxa queries should be printed. The default is `TRUE`.
#' @export
get_taxonomies.species_binomials<-function(path_to_input_species_binomials,path_to_output_local_taxa_list,path_to_taxonomy_edits=NA,print_taxize_queries=TRUE){
  
  # Check that the value supplied to whether taxize messages should be printed is logical.
  if((!is.logical(print_taxize_queries)) | is.na(print_taxize_queries)) stop("The value supplied to whether taxize messages should be printed must be TRUE or FALSE.")
  
  # Read in input csv file.
  taxa<-read.csv(file=path_to_input_species_binomials,stringsAsFactors=F)
  
  # Check that field names are right.
  if(!identical(colnames(taxa),c("Common_Name","Scientific_Name"))) stop('Fields in the input csv file should be "Common_Name" and "Scientific_Name".')
  
  # Throw an error if any blanks or NAs exist in the species binomial field
  # of the input csv file.
  if(any(taxa$Scientific_Name=="" | is.na(taxa$Scientific_Name))) stop("There are blanks or NAs in the species binomials field of the input csv file.")
  
  # Check that there are no duplicates in the input csv file.
  if(sum(duplicated(taxa$Scientific_Name)) > 0) stop("There are duplicated species binomials in the input csv file.")
  
  # Get NCBI taxonomies from the scientific names.
  # Species synonyms are accounted for (if the synonyms are present in NCBI),
  # but mispellings are not.
  taxonomies<-taxize::tax_name(sci=taxa$Scientific_Name,get=c("superkingdom","phylum","class","order","family","genus","species"),db="ncbi",messages=print_taxize_queries)
  
  # Add common name to the taxonomies.
  taxonomies$Common_Name<-taxa$Common_Name
  
  # Subset to just desired fields.
  taxa<-taxonomies[,c("Common_Name","query","superkingdom","phylum","class","order","family","genus","species")]
  
  # Rename fields.
  colnames(taxa)<-c("Common_Name","Query","Domain","Phylum","Class","Order","Family","Genus","Species")
  
  # Throw an error if NCBI taxonomies were not found any local taxa.
  if(all(is.na(taxa$Species))) stop("NCBI taxonomies were not found for any local taxa.")
  
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
  if(sum(is.na(taxa[,3:8])) > 0){
    
    # Get row indices of taxa which received incomplete NCBI taxonomies.
    partial_indices<-which(apply(X=taxa[,3:8],MARGIN=1,FUN=function(x) sum(is.na(x)) > 0))
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
    taxonomy_edits<-read.csv(file=path_to_taxonomy_edits,stringsAsFactors=F)
    
    # Throw an error if the fields of the taxonomy edits file are not
    # Old_Taxonomy, New_Taxonomy, Notes.
    if(!identical(colnames(taxonomy_edits),c("Old_Taxonomy","New_Taxonomy","Notes"))) stop("The fields of the taxonomy edits file must be: 'Old_Taxonomy', 'New_Taxonomy' ,'Notes'.")
    
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
      stop(paste0("Not all taxonomies found for local taxa have 7 taxonomic levels. If applying taxonomy edits, please ensure that these edits preserve the 7-level taxonomy. The following taxonomies lack 7 levels: ",paste(incorrect_num_taxonomic_levels_message,collapse=", ")))
      
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
    
    # Issue a warning to the user to update the unknown taxonomies in the local taxa
    # list before using the Local Taxa Tool.
    warning(paste0("NCBI taxonomies could not be found for ",nrow(trouble_taxa)," local species. Be sure to update unknown taxonomies in the local taxa list before using the Local Taxa Tool. Unknown taxonomies appear at the top of the local taxa list."))
    
  }
  
  # Remove the query field.
  taxa<-taxa[,-which(colnames(taxa)=="Query")]
  
  # Replace NA cells with blanks.
  taxa[is.na(taxa)]<-""
  
  # Write out taxonomies.
  write.csv(x=taxa,file=path_to_output_local_taxa_list,row.names=F)
  
}