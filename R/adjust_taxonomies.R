#' Adjust Taxonomies
#'
#' Performs adjustments to a taxonomy system according to a taxonomy edits file.
#' @param path_to_input_file String specifying path to list of species (in CSV format) whose taxonomies are to be adjusted. The file should contain the following fields: 'Common_Name', 'Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species'. There should be no `NA`s or blanks in the taxonomy fields, and the species field should contain the binomial name. Additional fields may be present in the input file, and fields can be in any order.
#' @param path_to_output_file String specifying path to output species list with adjusted taxonomies. The output file will be in CSV format.
#' @param path_to_taxonomy_edits String specifying path to taxonomy edits file in CSV format. The file must contain the following fields: 'Old_Taxonomy', 'New_Taxonomy', 'Notes'. Old taxonomies are replaced with new taxonomies in the order the records appear in the file. The taxonomic levels in the 'Old_Taxonomy' and 'New_Taxonomy' fields should be delimited by a semi-colon.
#' @returns No return value. Writes an output CSV file with adjusted taxonomies.
#' @seealso
#' [`get_taxonomies.species_binomials`][get_taxonomies.species_binomials()] for remotely fetching NCBI taxonomies from species binomials. \cr \cr
#' [`get_taxonomies.IUCN`][get_taxonomies.IUCN()] for formatting taxonomies from the IUCN Red List.
#' @examples
#' # Get path to input file.
#' path_to_input_file<-system.file("extdata",
#'                                 "example_local_taxa_list.csv",
#'                                 package="LocaTT",
#'                                 mustWork=TRUE)
#' 
#' # Get path to taxonomy edits.
#' path_to_taxonomy_edits<-system.file("extdata",
#'                                     "example_taxonomy_edits.csv",
#'                                     package="LocaTT",
#'                                     mustWork=TRUE)
#' 
#' # Create temporary output file path.
#' path_to_output_file<-tempfile(fileext=".csv")
#' 
#' # Adjust taxonomies.
#' adjust_taxonomies(path_to_input_file=path_to_input_file,
#'                   path_to_output_file=path_to_output_file,
#'                   path_to_taxonomy_edits=path_to_taxonomy_edits)
#' @export
adjust_taxonomies<-function(path_to_input_file,path_to_output_file,path_to_taxonomy_edits){
  
  # Check arguments.
  
  # Input file.
  ## Throw an error if the input file path is not a character string.
  if(!is.character(path_to_input_file)) stop("The input file path must be a character string.")
  ## Throw an error if the input file path has multiple elements.
  if(length(path_to_input_file) > 1) stop("The input file path cannot have multiple elements.")
  ## Throw an error if the input file path is NA.
  if(is.na(path_to_input_file)) stop("The input file path cannot be NA.")
  
  # Output file.
  ## Throw an error if the output file path is not a character string.
  if(!is.character(path_to_output_file)) stop("The output file path must be a character string.")
  ## Throw an error if the output file path has multiple elements.
  if(length(path_to_output_file) > 1) stop("The output file path cannot have multiple elements.")
  ## Throw an error if the output file path is NA.
  if(is.na(path_to_output_file)) stop("The output file path cannot be NA.")
  
  # Taxonomy edits file.
  ## Throw an error if the path to the taxonomy edits file is not a character string.
  if(!is.character(path_to_taxonomy_edits)) stop("The path to the taxonomy edits file must be a character string.")
  ## Throw an error if the path to the taxonomy edits file has multiple elements.
  if(length(path_to_taxonomy_edits) > 1) stop("The path to the taxonomy edits file path cannot have multiple elements.")
  ## Throw an error if the path to the taxonomy edits file is NA.
  if(is.na(path_to_taxonomy_edits)) stop("The path to the taxonomy edits file cannot be NA.")
  
  # Begin operations.
  
  # Read in input csv file.
  df<-utils::read.csv(file=path_to_input_file,stringsAsFactors=FALSE)
  
  # Set vector of required fields.
  required_fields<-c("Common_Name","Domain","Phylum","Class","Order","Family","Genus","Species")
  
  # If the input file is missing any of the required fields.
  if(!all(required_fields %in% colnames(df))){
    # Get the missing field names.
    missing_fields<-required_fields[!(required_fields %in% colnames(df))]
    # Throw an error stating which required fields are missing from the input file.
    stop(paste0('The input file is missing the following required fields: "',paste(missing_fields,collapse='", "'),'".'))
  }
  
  # Remove common name from the required fields vector.
  required_fields<-required_fields[required_fields!="Common_Name"]
  
  # Check that there are no NAs in the taxonomies of the input file.
  if(any(is.na(df[,required_fields]) | df[,required_fields]=="")) stop("There are NAs or blanks in the taxonomy fields of the input file. Please ensure that all taxonomy fields have entries for all records.")
  
  # Check that there are no underscores in the taxonomy fields of the input file.
  if(any(t(apply(X=df[,required_fields],MARGIN=1,FUN=grepl,pattern="_")))) stop("There cannot be underscores in the taxonomy fields of the input file.")
  
  # Collapse taxa names by semi-colons.
  taxa_names<-apply(X=df[,required_fields],MARGIN=1,FUN=paste,collapse=";")
  
  # Replace spaces in the taxa names vector with underscores.
  taxa_names<-gsub(pattern=" ",replacement="_",x=taxa_names)
  
  # Read in taxonomy edits.
  taxonomy_edits<-utils::read.csv(file=path_to_taxonomy_edits,stringsAsFactors=FALSE)
  
  # Throw an error if the fields of the taxonomy edits file are not
  # Old_Taxonomy, New_Taxonomy, Notes.
  if(!identical(colnames(taxonomy_edits),c("Old_Taxonomy","New_Taxonomy","Notes"))) stop('The fields of the taxonomy edits file must be: "Old_Taxonomy", "New_Taxonomy" ,"Notes".')
  
  # Check that there are no NAs in the taxonomy edits fields.
  if(any(is.na(taxonomy_edits[,c("Old_Taxonomy","New_Taxonomy")]) | taxonomy_edits[,c("Old_Taxonomy","New_Taxonomy")]=="")) stop('There are NAs or blanks in the "Old_Taxonomy" or "New_Taxonomy" fields of the taxonomy edits file. Please ensure that these fields have entries for all records.')
  
  # Check that there are no spaces in the taxonomy edits fields.
  if(any(t(apply(X=taxonomy_edits[,c("Old_Taxonomy","New_Taxonomy")],MARGIN=1,FUN=grepl,pattern=" ")))) stop('There cannot be spaces in the "Old_Taxonomy" or "New_Taxonomy" fields of the taxonomy edits file.')
  
  # Add a carrot to anchor the start of the old taxonomies field.
  taxonomy_edits$Old_Taxonomy<-paste0("^",taxonomy_edits$Old_Taxonomy)
  
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
    incorrect_num_taxonomic_levels_common_names<-df$Common_Name[check_num_taxonomic_levels!=7]
    # Create a field containing information on both the sequence name
    # and nucleotide sequence.
    incorrect_num_taxonomic_levels_message<-paste0(incorrect_num_taxonomic_levels_common_names," (",incorrect_num_taxonomic_levels_taxonomies,")")
    # Throw a error message mentioning the taxonomies which lack 7-levels.
    stop(paste0("Not all adjusted taxonomies have 7 taxonomic levels. Please ensure that taxonomy edits preserve the 7-level taxonomy. The following adjusted taxonomies lack 7 levels: ",paste(incorrect_num_taxonomic_levels_message,collapse=", ")))
    
  }
  
  # Get updated taxonomies by splitting the character strings by semi-colons.
  updated_taxonomies<-as.data.frame(do.call(rbind,strsplit(x=taxa_names,split=";")))
  
  # Update the taxonomies of the input file.
  df[,required_fields]<-updated_taxonomies
  
  # Write out adjusted taxonomies.
  utils::write.csv(x=df,file=path_to_output_file,row.names=FALSE)
  
}