#' Format Reference Databases
#'
#' Formats reference databases from MIDORI or UNITE for use with the [`local_taxa_tool`][local_taxa_tool()] function.
#' @param path_to_input_reference_database String specifying path to input reference database in FASTA format.
#' @param path_to_output_BLAST_database String specifying path to output BLAST database in FASTA format. File path cannot contain spaces.
#' @param input_reference_database_source String specifying input reference database source (`'MIDORI'` or `'UNITE'`). The default is `'MIDORI'`.
#' @param path_to_taxonomy_edits String specifying path to taxonomy edits file in CSV format. The file must contain the following fields: 'Old_Taxonomy', 'New_Taxonomy', 'Notes'. Old taxonomies are replaced with new taxonomies in the order the records appear in the file. The taxonomic levels in the 'Old_Taxonomy' and 'New_Taxonomy' fields should be delimited by a semi-colon. If no taxonomy edits are desired, then set this variable to `NA` (the default).
#' @param path_to_sequence_edits String specifying path to sequence edits file in CSV format. The file must contain the following fields: 'Action', 'Common_Name', 'Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species', 'Sequence', 'Notes'. The values in the 'Action' field must be either 'Add' or 'Remove', which will add or remove the respective sequence from the reference database. Values in the 'Common_Name' field are optional. Values should be supplied to all taxonomy fields. If using a reference database from MIDORI, then use NCBI superkingdom names (*e.g.*, 'Eukaryota') in the 'Domain' field. If using a reference database from UNITE, then use kingdom names (*e.g.*, 'Fungi') in the 'Domain' field. The 'Species' field should contain species binomials. Sequence edits are performed after taxonomy edits, if applied. If no sequence edits are desired, then set this variable to `NA` (the default).
#' @param path_to_list_of_local_taxa_to_subset String specifying path to list of species (in CSV format) to subset the reference database to. This option is helpful if the user wants the reference database to include only the sequences of local species. The file should contain the following fields: 'Common_Name', 'Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species'. There should be no 'NA's or blanks in the taxonomy fields. The species field should contain the binomial name without subspecies or other information below the species level. There should be no duplicate species (*i.e.*, multiple records with the same species binomial and taxonomy) in the species list. Subsetting the reference database to the sequences of certain species is performed after taxonomy and sequence edits are applied to the reference database, and species must match at all taxonomic levels in order to be retained in the reference database. If subsetting the reference database to the sequences of certain species is not desired, set this variable to `NA` (the default).
#' @param makeblastdb_command String specifying path to the makeblastdb program, which is a part of BLAST. The default (`'makeblastdb'`) should work for standard BLAST installations. The user can provide a path to the makeblastdb program for non-standard BLAST installations.
#' @returns No return value. Writes formatted BLAST database files.
#' @seealso
#' [`local_taxa_tool`][local_taxa_tool()] for performing geographically-conscious taxonomic assignment. \cr \cr
#' [`adjust_taxonomies`][adjust_taxonomies()] for adjusting a taxonomy system.
#' @examplesIf blast_command_found(blast_command="makeblastdb")
#' # Get path to example reference sequences FASTA file.
#' path_to_input_file<-system.file("extdata",
#'                                 "example_reference_sequences.fasta",
#'                                  package="LocaTT",
#'                                  mustWork=TRUE)
#' 
#' # Create a temporary file path for the output reference database FASTA file.
#' path_to_output_file<-tempfile(fileext=".fasta")
#' 
#' # Format reference database.
#' format_reference_database(path_to_input_reference_database=path_to_input_file,
#'                           path_to_output_BLAST_database=path_to_output_file)
#' @export
format_reference_database<-function(path_to_input_reference_database,path_to_output_BLAST_database,input_reference_database_source="MIDORI",path_to_taxonomy_edits=NA,path_to_sequence_edits=NA,path_to_list_of_local_taxa_to_subset=NA,makeblastdb_command="makeblastdb"){
  
  # Throw an error if the makeblastdb command cannot not be found.
  if(!blast_command_found(blast_command=makeblastdb_command)) stop("The makeblastdb command could not be found. If using a non-standard installation of BLAST, set the path to the makeblastdb command using the makeblastdb_command argument.")
  
  # Throw an error if the user-defined reference database source is not MIDORI or UNITE.
  if(!(input_reference_database_source %in% c("MIDORI","UNITE"))) stop("User-defined input reference database source must be 'MIDORI' or 'UNITE'.")
  
  # Throw an error if the path to the user-defined output BLAST database contains spaces.
  if(grepl(pattern=" ",x=path_to_output_BLAST_database)) stop("There cannot be spaces in path_to_output_BLAST_database.")
  
  # Read in reference database.
  reference<-read.fasta(file=path_to_input_reference_database)
  
  # If the reference database is from MIDORI.
  if(input_reference_database_source=="MIDORI"){
    
    # Remove everything in the reference names before and including ###root;.
    reference$Name<-sub(pattern="^.*###root_1;",replacement="",x=reference$Name)
    # Split up taxonomies by semi-colon.
    reference_taxonomy<-strsplit(x=reference$Name,split=";")
    # Remove everything after and including the last underscore of taxonomic level names.
    reference_taxonomy<-as.data.frame(t(sapply(X=reference_taxonomy,FUN=sub,pattern="_[^_]+$",replacement="")),stringsAsFactors=FALSE)
    # Remove any trailing underscores in the taxonomic level names.
    reference_taxonomy<-as.data.frame(apply(X=reference_taxonomy,MARGIN=2,FUN=gsub,pattern="_*$",replacement=""),stringsAsFactors=FALSE)
    # Replace any sets of repeating underscores with a single underscore.
    reference_taxonomy<-as.data.frame(apply(X=reference_taxonomy,MARGIN=2,FUN=gsub,pattern="(_)\\1+",replacement="_"),stringsAsFactors=FALSE)
    # Name the taxonomic levels.
    colnames(reference_taxonomy)<-c("Domain","Phylum","Class","Order","Family","Genus","Species")
    # Collapse reference taxonomies by semi-colons.
    reference_taxonomy<-apply(X=reference_taxonomy,MARGIN=1,FUN=paste,collapse=";")
    # Add the formatted reference taxonomies back to the reference database.
    reference$Name<-reference_taxonomy
    
  } else { # If the reference database is from UNITE.
    
    # Get the taxonomy portion of the UNITE names.
    reference$Name<-sapply(X=strsplit(x=reference$Name,split="\\|"),FUN="[[",5)
    # Split up taxonomies by semi-colon.
    reference_taxonomy<-strsplit(x=reference$Name,split=";")
    # Create a data frame of the taxonomic levels.
    reference_taxonomy<-as.data.frame(do.call("rbind",reference_taxonomy),stringsAsFactors=FALSE)
    # Name the taxonomic levels.
    colnames(reference_taxonomy)<-c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
    # Remove leading taxonomic level indicators from taxon names.
    reference_taxonomy<-as.data.frame(apply(X=reference_taxonomy,MARGIN=2,FUN=sub,pattern="^[kpcofgs]__",replacement=""),stringsAsFactors=FALSE)
    # Collapse reference taxonomies by semi-colons.
    reference_taxonomy<-apply(X=reference_taxonomy,MARGIN=1,FUN=paste,collapse=";")
    # Add the formatted reference taxonomies back to the reference database.
    reference$Name<-reference_taxonomy
    # Replace bad characters (multiplication sign) which BLAST cannot handle.
    reference$Name<-gsub(pattern="\u00D7",replacement="x",x=reference$Name)
    
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
    
    # Loop through each taxonomy edit.
    for(i in 1:nrow(taxonomy_edits)){
      
      # Get the taxonomy edit.
      taxonomy_edit<-taxonomy_edits[i,]
      # Apply translation to reference database taxonomies.
      reference$Name<-sub(pattern=taxonomy_edit$Old_Taxonomy,
                          replacement=taxonomy_edit$New_Taxonomy,
                          x=reference$Name)
      
    }
    
  }
  
  # If a path to a sequence edit file is provided.
  if(!is.na(path_to_sequence_edits)){
    
    # Read in edits to reference sequences.
    sequence_edits<-utils::read.csv(file=path_to_sequence_edits,stringsAsFactors=FALSE)
    
    # Throw an error if the fields of the sequence edits file are not
    # Action, Common_Name, Domain, Phylum, Class, Order, Family, Genus, Species, Sequence, Notes.
    if(!identical(colnames(sequence_edits),c("Action","Common_Name","Domain","Phylum","Class","Order","Family","Genus","Species","Sequence","Notes"))) stop("The fields of the sequence edits file must be: 'Action', 'Common_Name', 'Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species', 'Sequence','Notes'.")
    
    # Throw an error if there are values in the Action field of the sequence
    # edits file which are not Add or Remove.
    if(!all(sequence_edits$Action %in% c("Add","Remove"))) stop("There are values in the 'Action' field of the sequence edits file which are not 'Add' or 'Remove'.")
    
    # Check that there are no NAs in the sequence edits fields.
    if(any(is.na(sequence_edits[,c("Domain","Phylum","Class","Order","Family","Genus","Species","Sequence")]) | sequence_edits[,c("Domain","Phylum","Class","Order","Family","Genus","Species","Sequence")]=="")) stop("There are NAs or blanks in the taxonomy or sequence fields of the sequence edits file. Please ensure that these fields have entries for all records.")
    
    # Check that there are no underscores in the taxonomy fields of the sequence edits file.
    if(any(t(apply(X=sequence_edits[,c("Domain","Phylum","Class","Order","Family","Genus","Species")],MARGIN=1,FUN=grepl,pattern="_")))) stop("There cannot be underscores in the taxonomy fields of the sequence edits file.")
    
    # Collapse sequence edit taxonomies by semi-colons.
    sequence_edits$Name<-apply(X=sequence_edits[,c("Domain","Phylum","Class","Order","Family","Genus","Species")],MARGIN=1,FUN=paste,collapse=";")
    
    # Replace spaces in sequence edit taxonomies with underscores.
    sequence_edits$Name<-gsub(pattern=" ",replacement="_",x=sequence_edits$Name)
    
    # Loop through each sequence edit action.
    for(i in unique(sequence_edits$Action)){
      
      # If the action is Add.
      if(i=="Add"){
        
        # Get sequences to add.
        seqs_add<-sequence_edits[sequence_edits$Action==i,]
        # Get just the name and sequence fields for the sequences to add.
        seqs_add<-seqs_add[,c("Name","Sequence")]
        # Add the sequences to the reference database.
        reference<-rbind(reference,seqs_add)
        
      } else { # If the action is Remove.
        
        # Get sequences to remove.
        seqs_remove<-sequence_edits[sequence_edits$Action==i,]
        
        # Loop through each sequence to remove.
        for(j in 1:nrow(seqs_remove)){
          
          # Get the sequence to remove.
          seq_remove<-seqs_remove[j,]
          
          # If the sequence to remove exists in the reference database.
          if(any(reference$Name==seq_remove$Name & reference$Sequence==seq_remove$Sequence)){
            
            # Remove the sequence from the reference database.
            reference<-reference[!(reference$Name==seq_remove$Name & reference$Sequence==seq_remove$Sequence),]
            
          } else { # If the sequence to remove does not exist in the reference database.
            
            # Throw a warning.
            warning(paste0("The following sequence to remove does not exist in the reference database: ",seq_remove$Species," (",seq_remove$Sequence,")"))
            
          }
          
        }
        
      }
      
    }
    
  }
  
  # Subset reference database to unique sequence-species combinations.
  reference<-unique(reference)
  
  # Throw an error if the reference database contains spaces in sequence names.
  if(any(grepl(pattern=" ",x=reference$Name))) stop("The reference database cannot contain spaces in sequence names.")
  
  # Get the number of taxonomic levels each sequence has.
  check_num_taxonomic_levels<-sapply(X=strsplit(x=reference$Name,split=";"),FUN=length)
  
  # If any sequences in the reference database lack a 7-level taxonomy.
  if(any(check_num_taxonomic_levels!=7)){
    
    # Get the sequences which lack a 7-level taxonomy.
    incorrect_num_taxonomic_levels<-reference[check_num_taxonomic_levels!=7,]
    # Create a field containing information on both the sequence name and nucleotide sequence.
    incorrect_num_taxonomic_levels$Message<-paste0(incorrect_num_taxonomic_levels$Name," (",incorrect_num_taxonomic_levels$Sequence,")")
    # Throw a error message mentioning the sequences which lack a 7-level taxonomy.
    stop(paste0("Not all sequences in the formatted database have 7 taxonomic levels. If applying taxonomy or sequence edits, please ensure that these edits preserve the 7-level taxonomy. The following sequences lack a 7-level taxonomy: ",paste(incorrect_num_taxonomic_levels$Message,collapse=", ")))
    
  }
  
  # If a local taxa list to subset the reference database to is provided.
  if(!is.na(path_to_list_of_local_taxa_to_subset)){
    
    # Read in local taxa list.
    local<-utils::read.csv(file=path_to_list_of_local_taxa_to_subset,stringsAsFactors=FALSE)
    
    # Check that the correct fields are present in the local taxa list.
    if(!identical(colnames(local),c("Common_Name","Domain","Phylum","Class","Order","Family","Genus","Species"))) stop('The field names in the local taxa list file should be: "Common_Name", "Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species".')
    
    # Check that there are no NAs in the taxonomies of the local taxa list.
    if(any(is.na(local[,c("Domain","Phylum","Class","Order","Family","Genus","Species")]) | local[,c("Domain","Phylum","Class","Order","Family","Genus","Species")]=="")) stop("There are NAs or blanks in the taxonomies of the local taxa list file. Please ensure that all taxonomy fields have entries for all records.")
    
    # Check that there are no underscores in the taxonomy fields of the local taxa list.
    if(any(t(apply(X=local[,c("Domain","Phylum","Class","Order","Family","Genus","Species")],MARGIN=1,FUN=grepl,pattern="_")))) stop("There cannot be underscores in the taxonomy fields of the local taxa list file.")
    
    # Collapse local taxa names by semi-colons.
    local$Name<-apply(X=local[,c("Domain","Phylum","Class","Order","Family","Genus","Species")],MARGIN=1,FUN=paste,collapse=";")
    
    # Check that all species are unique.
    if(any(duplicated(local$Name))) stop("There are duplicate species (i.e., there are multiple records with the same taxonomy) in the local taxa list file.")
    
    # Replace spaces in the local name field with underscores.
    local$Name<-gsub(pattern=" ",replacement="_",x=local$Name)
    
    # If any sequences of local species are present in the reference database.
    if(any(local$Name %in% reference$Name)){
      
      # Subset the reference database to sequences of local species.
      reference<-reference[reference$Name %in% local$Name,]
      
    } else { # If no sequences of local species are present in the reference database.
      
      # Throw an error.
      stop("No sequences of local species are present in the reference database.")
      
    }
    
  }
  
  # Write out the formatted reference database.
  write.fasta(names=reference$Name,sequences=reference$Sequence,
              file=path_to_output_BLAST_database)
  
  # Make a BLAST database out of the formatted reference database fasta file.
  system2(command=makeblastdb_command,args=c(paste0("-in ",path_to_output_BLAST_database),
                                             "-dbtype nucl"),
          stdout=FALSE)
  
}