#' Perform Geographically-Conscious Taxonomic Assignment
#' 
#' @description Performs taxonomic assignment of DNA metabarcoding sequences while considering geographic location.
#' @details Sequences are BLASTed against a global reference database, and the tool suggests locally occurring species which are most closely related (by taxonomy) to any of the best-matching BLAST hits (by bit score). Optionally, local sister taxonomic groups without reference sequences can be added to the local taxa suggestions by setting the `include_missing` argument to `TRUE`. If a local taxa list is not provided, then local taxa suggestions will be disabled, but all best-matching BLAST hits will still be returned. Alternatively, a reference database containing just the sequences of local species can be used, and local taxa suggestions can be disabled to return all best BLAST matches of local species. The reference database should be formatted with the [`format_reference_database`][format_reference_database()] function, and the local taxa lists can be prepared using the [`get_taxonomies.species_binomials`][get_taxonomies.species_binomials()] and [`get_taxonomies.IUCN`][get_taxonomies.IUCN()] functions. Output field definitions are:
#' * Sequence_name: The query sequence name.
#' * Sequence: The query sequence.
#' * Best_match_references: Species binomials of all best-matching BLAST hits (by bit score) from the reference database.
#' * Best_match_E_value: The E-value associated with the best-matching BLAST hits.
#' * Best_match_bit_score: The bit score associated with the best-matching BLAST hits.
#' * Best_match_query_cover.mean: The mean query cover of all best-matching BLAST hits.
#' * Best_match_query_cover.SD: The standard deviation of query cover of all best-matching BLAST hits.
#' * Best_match_PID.mean: The mean percent identity of all best-matching BLAST hits.
#' * Best_match_PID.SD: The standard deviation of percent identity of all best-matching BLAST hits.
#' * Local_taxa (Field only present if a path to a local taxa list is provided): The finest taxonomic unit(s) which include both any species of the best-matching BLAST hits and any local species. If the species of any of the best-matching BLAST hits are local, then the finest taxonomic unit(s) are at the species level.
#' * Local_species (Field only present if a path to a local taxa list is provided): Species binomials of all local species which belong to the taxonomic unit(s) in the Local_taxa field.
#' * Local_taxa.include_missing (Field only present if both a path to a local taxa list is provided and the `include_missing` argument is set to `TRUE`): Local sister taxonomic groups without reference sequences are added to the local taxa suggestions from the Local_taxa field.
#' * Local_species.include_missing (Field only present if both a path to a local taxa list is provided and `include_missing` argument is set to `TRUE`): Species binomials of all local species which belong to the taxonomic unit(s) in the Local_taxa.include_missing field.
#' @param path_to_sequences_to_classify String specifying path to FASTA file containing sequences to classify. File path cannot contain spaces.
#' @param path_to_BLAST_database String specifying path to BLAST reference database in FASTA format. File path cannot contain spaces.
#' @param path_to_output_file String specifying path to output file of classified sequences in CSV format.
#' @param path_to_list_of_local_taxa String specifying path to list of local species in CSV format. The file should contain the following fields: 'Common_Name', 'Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species'. There should be no 'NA's or blanks in the taxonomy fields. The species field should contain the binomial name without subspecies or other information below the species level. There should be no duplicate species (*i.e.*, multiple records with the same species binomial and taxonomy) in the local species list. If local taxa suggestions are not desired, set this variable to `NA` (the default).
#' @param include_missing Logical. If `TRUE`, then additional fields are included in the output CSV file in which local sister taxonomic groups without reference sequences are added to the local taxa suggestions. If `FALSE` (the default), then this is not performed.
#' @param blast_e_value Numeric. Maximum E-value of returned BLAST hits (lower E-values are associated with more 'significant' matches). The default is `1e-05`.
#' @param blast_max_target_seqs Numeric. Maximum number of BLAST target sequences returned per query sequence. Enough target sequences should be returned to ensure that all minimum E-value matches are returned for each query sequence. A warning will be produced if this value is not sufficient. The default is `2000`.
#' @param blast_task String specifying BLAST task specification. Use `'megablast'` (the default) to find very similar sequences (*e.g.*, intraspecies or closely related species). Use `'blastn-short'` for sequences shorter than 50 bases. See the blastn program help documentation for additional options and details.
#' @param full_names Logical. If `TRUE`, then full taxonomies are returned in the output CSV file. If `FALSE` (the default), then only the lowest taxonomic levels (*e.g.*, species binomials instead of the full species taxonomies) are returned in the output CSV file.
#' @param underscores Logical. If `TRUE`, then taxa names in the output CSV file use underscores instead of spaces. If `FALSE` (the default), then taxa names in the output CSV file use spaces.
#' @param separator String specifying the separator to use between taxa names in the output CSV file. The default is `', '`.
#' @param blastn_command String specifying path to the blastn program. The default (`'blastn'`) should work for standard BLAST installations. The user can provide a path to the blastn program for non-standard BLAST installations.
#' @returns No return value. Writes an output CSV file with fields defined in the details section.
#' @references A manuscript describing this taxonomic assignment method is in preparation.
#' @examplesIf blast_command_found(blast_command="blastn")
#' # Get path to example query sequences FASTA file.
#' path_to_query_sequences<-system.file("extdata",
#'                                      "example_query_sequences.fasta",
#'                                      package="LocaTT",
#'                                      mustWork=TRUE)
#' 
#' # Get path to example reference database FASTA file.
#' path_to_reference_database<-system.file("extdata",
#'                                         "example_blast_database.fasta",
#'                                         package="LocaTT",
#'                                         mustWork=TRUE)
#' 
#' # Get path to example local taxa list CSV file.
#' path_to_local_taxa_list<-system.file("extdata",
#'                                      "example_local_taxa_list.csv",
#'                                      package="LocaTT",
#'                                      mustWork=TRUE)
#' 
#' # Create a temporary file path for the output CSV file.
#' path_to_output_CSV_file<-tempfile(fileext=".csv")
#' 
#' # Run the local taxa tool.
#' local_taxa_tool(path_to_sequences_to_classify=path_to_query_sequences,
#'                 path_to_BLAST_database=path_to_reference_database,
#'                 path_to_output_file=path_to_output_CSV_file,
#'                 path_to_list_of_local_taxa=path_to_local_taxa_list,
#'                 include_missing=TRUE,
#'                 full_names=TRUE,
#'                 underscores=TRUE)
#' @export
local_taxa_tool<-function(path_to_sequences_to_classify,path_to_BLAST_database,path_to_output_file,path_to_list_of_local_taxa=NA,include_missing=FALSE,blast_e_value=1e-5,blast_max_target_seqs=2000,blast_task="megablast",full_names=FALSE,underscores=FALSE,separator=", ",blastn_command="blastn"){
  
  # Throw an error if the blastn command cannot not be found.
  if(!blast_command_found(blast_command=blastn_command)) stop("The blastn command could not be found. If using a non-standard installation of BLAST, set the path to the blastn command using the blastn_command argument.")
  
  # Throw an error if the path to the user-defined fasta file to classify contains spaces.
  if(grepl(pattern=" ",x=path_to_sequences_to_classify)) stop("There cannot be spaces in path_to_sequences_to_classify.")
  
  # Throw an error if the path to the user-defined BLAST database contains spaces.
  if(grepl(pattern=" ",x=path_to_BLAST_database)) stop("There cannot be spaces in path_to_BLAST_database.")
  
  # Throw an error if the include missing argument is not TRUE or FALSE.
  if(!is.logical(include_missing) | is.na(include_missing)) stop("The include_missing argument must be TRUE or FALSE.")
  
  # Throw an error if the full names argument is not TRUE or FALSE.
  if(!is.logical(full_names) | is.na(full_names)) stop("The full_names argument must be TRUE or FALSE.")
  
  # Throw an error if the underscores argument is not TRUE or FALSE.
  if(!is.logical(underscores) | is.na(underscores)) stop("The underscores argument must be TRUE or FALSE.")
  
  # Throw an error if the separator argument is not a character string.
  if(!is.character(separator)) stop("The separator argument must be a character string.")
  
  # Throw a warning if include missing is enabled but a local taxa list is not provided.
  if(include_missing & is.na(path_to_list_of_local_taxa)) warning("The include_missing argument is ignored if a local taxa list is not provided.")
  
  # Read in sequences to classify.
  sequences_to_classify<-read.fasta(file=path_to_sequences_to_classify)
  colnames(sequences_to_classify)[1]<-"Sequence_name"
  
  # Throw an error if there are NAs or blanks in the fasta file of sequences to classify.
  if(any(is.na(sequences_to_classify) | sequences_to_classify=="")) stop("NAs or blanks are present in the sequence names or sequences of the fasta file containing the sequences to classify.")
  
  # Throw an error if not all sequence names are unique in the fasta file of sequences
  # to classify.
  if(any(duplicated(sequences_to_classify$Sequence_name))) stop("Not all sequence names are unique in the fasta file of sequences to classify.")
  
  # Perform the command line BLAST algorithm between the query and reference sequences.
  blast_output<-system2(command=blastn_command,
                        args=c(paste0('-task ',blast_task),
                               paste0('-db ',path_to_BLAST_database),
                               paste0('-query ',path_to_sequences_to_classify),
                               '-outfmt "6 qseqid sseqid evalue bitscore qcovs pident"',
                               paste0('-max_target_seqs ',blast_max_target_seqs),
                               paste0('-evalue ',blast_e_value)),
                        stdout=TRUE)
  
  # Check that BLAST returned results for at least one sequence.
  if(length(blast_output)==0) stop("BLAST did not return results for any sequences. Try increasing the user-defined minimum E-value parameter and run this function again.")
  
  # Format the BLAST output.
  ## Get output as a data frame.
  blast_output<-as.data.frame(do.call("rbind",strsplit(x=blast_output,split="\t")),
                              stringsAsFactors=FALSE)
  ## Provide field names.
  colnames(blast_output)<-c("Sequence_name","Matches","E_value",
                            "Bit_score","Query_cover","PID")
  ## If spaces in taxa names are desired, replace underscores in matches with spaces.
  if(!underscores){
    blast_output$Matches<-gsub(pattern="_",replacement=" ",x=blast_output$Matches)
  }
  ## Format the fields for E-value and after as numeric.
  blast_output[,3:ncol(blast_output)]<-sapply(X=blast_output[,3:ncol(blast_output)],
                                              FUN=as.numeric)
  
  # Get best matches for each query sequence.
  # We'll use the maximum bit score because the E-values
  # can have underflow issues (i.e., for a given query sequence, there can
  # be multiple hits with very low E-values, which underflow to zero, and
  # these hits can have different bit scores). Using the maximum bit score
  # is equivalent to using the minimum E-value.
  ## Get the maximum bit score match for each query sequence.
  max_bit_score<-stats::aggregate(Bit_score~Sequence_name,data=blast_output,FUN=max)
  ## Rename field to maximum bit score.
  colnames(max_bit_score)[2]<-"max_bit_score"
  ## Merge maximum bit score with the BLAST output.
  blast_output<-merge(blast_output,max_bit_score)
  ## Keep just records which have the maximum bit score for each query sequence.
  blast_output<-blast_output[blast_output$Bit_score==blast_output$max_bit_score,]
  
  # Get number of best BLAST matches returned for each query sequence.
  num_best_BLAST_matches_by_sequence<-table(blast_output$Sequence_name)
  
  # Check whether each query sequence had the maximum number of returned best BLAST matches.
  max_best_BLAST_matches_by_sequence<-num_best_BLAST_matches_by_sequence==blast_max_target_seqs
  
  # If any of the unique sequences had no bit scores less than their maximum bit scores,
  # then we may not have all of the best-matching taxa for these sequences!
  if(any(max_best_BLAST_matches_by_sequence)){
    
    # Get the names of unique sequences which had no bit scores less than their maximum
    # bit scores.
    sequences_at_max_target_seq_matches<-data.frame(
      Sequence_name=names(which(max_best_BLAST_matches_by_sequence)),
      stringsAsFactors=FALSE)
    # Add the actual sequences to the data frame.
    sequences_at_max_target_seq_matches<-merge(x=sequences_at_max_target_seq_matches,
                                               sequences_to_classify,all.x=TRUE)
    # Create a single string for each sequence.
    sequences_at_max_target_seq_matches$Message<-paste0(
      sequences_at_max_target_seq_matches$Sequence_name,
      " (",sequences_at_max_target_seq_matches$Sequence,")")
    # Issue a warning that blast_max_target_seqs is not high enough for these sequences.
    warning(paste0("There are ",nrow(sequences_at_max_target_seq_matches),
                   " sequences for which a blast_max_target_seqs of ",blast_max_target_seqs,
                   " was insufficient. These sequences are the following: ",
                   paste(sequences_at_max_target_seq_matches$Message,collapse=", ")))
    
  }
  
  # Get average query cover and PID for the best-matching references for each sequence.
  ## Calculate means.
  best_match_means<-stats::aggregate(.~Sequence_name,
                                data=blast_output[,c("Sequence_name","Query_cover","PID")],
                                FUN=mean)
  ## Add a mean label to the column names.
  colnames(best_match_means)[2:3]<-paste0(colnames(best_match_means)[2:3],".mean")
  
  # Get standard deviation of query cover and PID for the best-matching references for
  # each sequence.
  ## Calculate standard deviation.
  best_match_SDs<-stats::aggregate(.~Sequence_name,
                                data=blast_output[,c("Sequence_name","Query_cover","PID")],
                                FUN=stats::sd)
  ## Add a standard deviation label to the column names.
  colnames(best_match_SDs)[2:3]<-paste0(colnames(best_match_SDs)[2:3],".SD")
  
  # Merge the averages and standard deviations of query cover and PID for the
  # best-matching references for each sequence.
  best_match_means_and_SDs<-merge(best_match_means,best_match_SDs)
  
  # Turn query to lower case.
  colnames(best_match_means_and_SDs)<-gsub(pattern="^Query_",replacement="query_",
                                           x=colnames(best_match_means_and_SDs))
  
  # Append best match to the numeric fields.
  colnames(best_match_means_and_SDs)[2:ncol(best_match_means_and_SDs)]<-paste0(
    "Best_match_",colnames(best_match_means_and_SDs)[2:ncol(best_match_means_and_SDs)])
  
  # Get just the sequence, matches, E-value, and bit-score fields in the BLAST output.
  blast_output<-blast_output[,c("Sequence_name","Matches","E_value","Bit_score")]
  
  # Define function for getting a given taxonomic level.
  getTaxonomicLevel<-function(character_vector,desired_level_numeric){
    # Split character vector names by semi-colon.
    split_character_vector_names<-strsplit(x=character_vector,split=";")
    # Define internal function for collapsing split names at the desired taxonomic level.
    collapse_split_name<-function(split_name) paste(split_name[1:desired_level_numeric],
                                                    collapse=";")
    # Apply the collapse names function across the list of split names.
    collapsed_names<-sapply(X=split_character_vector_names,FUN=collapse_split_name)
    # Return the collapsed names.
    return(collapsed_names)
  }
  
  # If a local taxa list is provided.
  if(!is.na(path_to_list_of_local_taxa)){
    
    # Read in local taxa list.
    local<-utils::read.csv(file=path_to_list_of_local_taxa,stringsAsFactors=FALSE)
    
    # Check that the correct fields are present in the local taxa list.
    if(!identical(colnames(local),c("Common_Name","Domain","Phylum","Class","Order",
                                    "Family","Genus","Species"))){
      stop('The field names in the local taxa list file should be: "Common_Name", "Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species".')
    }
    
    # Check that there are no NAs in the taxonomies of the local taxa list.
    if(any(is.na(local[,c("Domain","Phylum","Class","Order","Family","Genus","Species")]) | local[,c("Domain","Phylum","Class","Order","Family","Genus","Species")]=="")){
      stop("There are NAs or blanks in the taxonomies of the local taxa list file. Please ensure that all taxonomy fields have entries for all records.")
    }
    
    # Check that there are no underscores in the taxonomy fields of the local taxa list.
    if(any(t(apply(X=local[,c("Domain","Phylum","Class","Order","Family","Genus","Species")],
                   MARGIN=1,FUN=grepl,pattern="_")))){
      stop("There cannot be underscores in the taxonomy fields of the local taxa list file.")
    }
    
    # If underscores are desired, replace spaces with underscores in the
    # taxonomy fields of the local taxa list.
    if(underscores){
      local[,c("Domain","Phylum","Class","Order","Family","Genus","Species")]<-data.frame(
        lapply(X=local[,c("Domain","Phylum","Class","Order",
                          "Family","Genus","Species")],
               FUN=gsub,pattern=" ",replacement="_"))
    }
    
    # Collapse local taxa names by semi-colons.
    local$Name<-apply(X=local[,c("Domain","Phylum","Class","Order",
                                 "Family","Genus","Species")],
                      MARGIN=1,FUN=paste,collapse=";")
    
    # Check that all species are unique.
    if(any(duplicated(local$Name))) stop("There are duplicate species (i.e., there are multiple records with the same taxonomy) in the local taxa list file.")
    
    # If include missing is enabled.
    if(include_missing){
      
      # Read in the reference database.
      ref<-read.fasta(file=path_to_BLAST_database)
      
      # If underscores are not desired.
      if(!underscores){
        # Replace underscres in reference names with spaces.
        ref$Name<-gsub(pattern="_",replacement=" ",x=ref$Name)
      }
      
    }
    
    # Get unique sequence names from BLAST output. Sequences not included
    # in the BLAST output had no significant matches.
    matched_sequences_to_classify<-unique(blast_output$Sequence_name)
    
    # Get the number of unique sequences from BLAST output.
    num_seq<-length(matched_sequences_to_classify)
    
    # Create an empty data frame for storing sequence classifications.
    sequence_classifications<-data.frame(NULL)
    
    # Loop through each sequence to classify.
    for(i in 1:num_seq){
      
      # Get the sequence to classify.
      sequence_to_classify<-matched_sequences_to_classify[i]
      
      # Get the BLAST output for the sequence to classify.
      sequence_to_classify_df<-blast_output[blast_output$Sequence_name==sequence_to_classify,]
      
      # Get the maximum bit score and its associated minimum E-value.
      # Using just the first values since all of them will always be the same.
      ## Minimum E-value.
      min_E_value<-sequence_to_classify_df$E_value[1]
      ## Maximum bit-score.
      bit_score<-sequence_to_classify_df$Bit_score[1]
      
      # Get the unique taxa associated with the highest bit score.
      best_matches<-sort(unique(sequence_to_classify_df$Matches))
      
      # Get the species names for the best-matching reference sequences.
      if(full_names){ # If full names are requested, get the full taxonomies.
        best_match_references<-best_matches
      } else { # If full names are not requested, get the species binomials.
        best_match_references<-sapply(X=strsplit(x=best_matches,split=";"),FUN="[[",7)
      }
      
      # Set the while loop continuation variable.
      inclusive_taxonomic_group_found<-"No"
      
      # Set an initial current taxonomic level.
      current_taxonomic_level<-7
      
      # Simplify and check taxonomic names until matches are found between the taxonomies
      # of the maximum bit score reference sequences and the sequence to classify.
      while(inclusive_taxonomic_group_found=="No"){
        
        # Get the current taxonomic level of the best-matches.
        best_matches.current_level<-unique(getTaxonomicLevel(character_vector=best_matches,
                                           desired_level_numeric=current_taxonomic_level))
        
        # Get the current taxonomic level of the local species.
        local.current_level<-getTaxonomicLevel(character_vector=local$Name,
                                               desired_level_numeric=current_taxonomic_level)
        
        # If any best matches share the current taxonomic level with local species.
        if(any(best_matches.current_level %in% local.current_level)){
          
          # Change the while loop continuation variable.
          inclusive_taxonomic_group_found<-"Yes"
          
          # Get the inclusive taxa.
          inclusive_taxa<-sort(
            best_matches.current_level[best_matches.current_level %in% local.current_level])
          
          # Get potential local species.
          potential_spp<-sort(local$Name[local.current_level %in% inclusive_taxa])
          
          # If include missing is enabled.
          if(include_missing){
            
            # If the current taxonomic level is at phylum or below.
            if(current_taxonomic_level > 1){
              
              # Get the inclusive taxonomic level above the current one.
              inclusive_taxa.upper_level<-unique(
                getTaxonomicLevel(character_vector=inclusive_taxa,
                                  desired_level_numeric=current_taxonomic_level-1))
              
              # Get the taxonomic level above the current one for local species.
              local.upper_level<-getTaxonomicLevel(character_vector=local$Name,
                                    desired_level_numeric=current_taxonomic_level-1)
              
              # Create a logical vector specifying whether the local taxonomic
              # group is a sister group to the inclusive taxa.
              sister_group<-(local.upper_level %in% inclusive_taxa.upper_level)
              
            } else { # If the current taxonomic level is at domain.
              
              # Set the sister group logical vector to TRUE, since all domains
              # are sister groups.
              sister_group<-TRUE
              
            }
            
            # If not all sister taxonomic groups are the inclusive taxa
            # (i.e., if there are sister groups to the inclusive taxa).
            if(!all(local.current_level[sister_group] %in% inclusive_taxa)){
              
              # Get the current taxonomic level of the reference sequences.
              ref.current_level<-getTaxonomicLevel(character_vector=ref$Name,
                                          desired_level_numeric=current_taxonomic_level)
              
              # Create a logical vector specifying whether the local taxonomic
              # group is missing from the reference database.
              missing_ref_seqs<-(!(local.current_level %in% ref.current_level))
              
              # Create a logical vector specifying whether the local taxonomic
              # group is both a sister group to the inclusive taxa and
              # is missing from the reference database.
              missing_sister_group<-sister_group & missing_ref_seqs
              
              # If there are any local sister groups to the inclusive taxa
              # which are missing from the reference database.
              if(any(missing_sister_group)){
                
                # Get the local taxonomic groups sister to the local taxa
                # which are missing from the reference database.
                local_absent.current_level<-unique(local.current_level[missing_sister_group])
                
                # Get the local species belonging to the sister groups
                # which are missing from the reference database.
                local_absent.spp<-local$Name[missing_sister_group]
                
                # Add the missing taxa to the inclusive taxa.
                inclusive_taxa.include_missing<-sort(c(inclusive_taxa,
                                                       local_absent.current_level))
                
                # Add the missing species to the potential species.
                potential_spp.include_missing<-sort(c(potential_spp,local_absent.spp))
                
              } else {
                
                # If there are no local sister groups to the inclusive taxa
                # which are missing from the reference database.
                
                # Set the include missing taxa to the inclusive taxa.
                inclusive_taxa.include_missing<-inclusive_taxa
                
                # Set the include missing species to the potential species.
                potential_spp.include_missing<-potential_spp
                
              }
              
            } else {
              
              # If all sister taxonomic groups are the inclusive taxa
              # (i.e., if there are no sister groups to the inclusive taxa),
              # then there are no additional taxa to add.
              
              # Set the include missing taxa to the inclusive taxa.
              inclusive_taxa.include_missing<-inclusive_taxa
              
              # Set the include missing species to the potential species.
              potential_spp.include_missing<-potential_spp
              
            }
            
          }
          
          # If full names are not requested, report simplified taxonomies.
          if(!full_names){
            
            # Potential local species.
            potential_spp<-sapply(X=strsplit(x=potential_spp,split=";"),
                                  FUN="[[",7)
            # Inclusive taxa.
            inclusive_taxa<-sapply(X=strsplit(x=inclusive_taxa,split=";"),
                                   FUN="[[",current_taxonomic_level)
            
            # If include missing is enabled.
            if(include_missing){
              
              # Potential local species (including missing).
              potential_spp.include_missing<-sapply(X=strsplit(
                                                      x=potential_spp.include_missing,
                                                      split=";"),
                                                    FUN="[[",7)
              # Inclusive taxa (including missing).
              inclusive_taxa.include_missing<-sapply(X=strsplit(
                                                       x=inclusive_taxa.include_missing,
                                                       split=";"),
                                                     FUN="[[",current_taxonomic_level)
              
            }
            
          }
          
          # Collect classification information in a data frame.
          if(include_missing){
            # If include missing is enabled.
            sequence_information<-data.frame(Local_taxa=paste(inclusive_taxa,
                                               collapse=separator),
                                             Local_taxa.include_missing=paste(
                                               inclusive_taxa.include_missing,
                                               collapse=separator),
                                             Local_species=paste(potential_spp,
                                               collapse=separator),
                                             Local_species.include_missing=paste(
                                               potential_spp.include_missing,
                                               collapse=separator),
                                             Best_match_references=paste(
                                               best_match_references,
                                               collapse=separator),
                                             Best_match_E_value=min_E_value,
                                             Best_match_bit_score=bit_score,
                                             Sequence_name=sequence_to_classify,
                                             stringsAsFactors=FALSE)
          } else {
            # If include missing is not enabled.
            sequence_information<-data.frame(Local_taxa=paste(inclusive_taxa,
                                               collapse=separator),
                                             Local_species=paste(potential_spp,
                                               collapse=separator),
                                             Best_match_references=paste(
                                               best_match_references,collapse=separator),
                                             Best_match_E_value=min_E_value,
                                             Best_match_bit_score=bit_score,
                                             Sequence_name=sequence_to_classify,
                                             stringsAsFactors=FALSE)
          }
          
        } else {
          
          # If no best matches share the current taxonomic level with local species.
          
          # Move up a taxonomic level for the next iteration.
          current_taxonomic_level<-current_taxonomic_level-1
          
          # If the current taxonomic level counter reaches zero.
          # Then the inclusive taxonomic group is root.
          if(current_taxonomic_level==0){
            
            # Change the while loop continuation variable.
            inclusive_taxonomic_group_found<-"Yes"
            
            # Collect classification information in a data frame.
            if(include_missing){
              # If include missing is enabled.
              sequence_information<-data.frame(Local_taxa="root",
                                               Local_taxa.include_missing="root",
                                               Local_species=ifelse(full_names,
                                                 paste(sort(local$Name),collapse=separator),
                                                 paste(local$Species[order(local$Name)],
                                                       collapse=separator)),
                                               Local_species.include_missing=ifelse(full_names,
                                                 paste(sort(local$Name),collapse=separator),
                                                 paste(local$Species[order(local$Name)],
                                                       collapse=separator)),
                                               Best_match_references=paste(
                                                 best_match_references,collapse=separator),
                                               Best_match_E_value=min_E_value,
                                               Best_match_bit_score=bit_score,
                                               Sequence_name=sequence_to_classify,
                                               stringsAsFactors=FALSE)
            } else {
              # If include missing is not enabled.
              sequence_information<-data.frame(Local_taxa="root",
                                               Local_species=ifelse(full_names,
                                                 paste(sort(local$Name),collapse=separator),
                                                 paste(local$Species[order(local$Name)],
                                                       collapse=separator)),
                                               Best_match_references=paste(
                                                 best_match_references,collapse=separator),
                                               Best_match_E_value=min_E_value,
                                               Best_match_bit_score=bit_score,
                                               Sequence_name=sequence_to_classify,
                                               stringsAsFactors=FALSE)
            }
            
          }
          
        }
        
      }
      
      # Add sequence information to the storage data frame.
      sequence_classifications<-rbind(sequence_classifications,sequence_information)
      
    }
    
    # Merge taxa suggestions with means and SDs of query cover and PID.
    sequence_classifications<-merge(sequence_classifications,best_match_means_and_SDs)
    
    # Add sequences to the classification data frame.
    df<-merge(sequences_to_classify,sequence_classifications,all=TRUE)
    
    # Keep certain fields in the classification data frame.
    if(include_missing){
      # If include missing is enabled.
      df<-df[,c("Sequence_name","Sequence","Best_match_references","Best_match_E_value",
                "Best_match_bit_score","Best_match_query_cover.mean",
                "Best_match_query_cover.SD","Best_match_PID.mean","Best_match_PID.SD",
                "Local_taxa","Local_species","Local_taxa.include_missing",
                "Local_species.include_missing")]
    } else {
      # If include missing is not enabled.
      df<-df[,c("Sequence_name","Sequence","Best_match_references","Best_match_E_value",
                "Best_match_bit_score","Best_match_query_cover.mean",
                "Best_match_query_cover.SD","Best_match_PID.mean","Best_match_PID.SD",
                "Local_taxa","Local_species")]
    }
    
  } else { # If a local taxa list is not provided.
    
    # Get unique best-matches for each sequence.
    blast_output<-unique(blast_output)
    
    # Order the records in the best-matches by sequence name and taxon name.
    blast_output<-blast_output[order(blast_output$Sequence_name,blast_output$Matches),]
    
    # Get the species names for the best-matching reference sequences.
    if(full_names){ # If full names are requested, get the full taxonomies.
      blast_output$Species_names<-blast_output$Matches
    } else { # If full names are not requested, get the species binomials.
      blast_output$Species_names<-sapply(X=strsplit(x=blast_output$Matches,split=";"),
                                         FUN="[[",7)
    }
    
    # Get best-matching references.
    best_matches<-stats::aggregate(Species_names~Sequence_name,data=blast_output,
                                   FUN=paste,collapse=separator)
    
    # Get the E-values of best-matches.
    best_matches_E_value<-stats::aggregate(E_value~Sequence_name,data=blast_output,
                                           FUN="[",1)
    
    # Get the bit-scores of best-matches.
    best_matches_Bit_score<-stats::aggregate(Bit_score~Sequence_name,data=blast_output,
                                             FUN="[",1)
    
    # Combine information on best-matching references.
    best_matches<-merge(best_matches,best_matches_E_value)
    best_matches<-merge(best_matches,best_matches_Bit_score)
    best_matches<-merge(best_matches,best_match_means_and_SDs)
    
    # Add sequences to the classification data frame.
    df<-merge(sequences_to_classify,best_matches,all=TRUE)
    
    # Rename certain fields.
    colnames(df)[which(colnames(df)=="Species_names")]<-"Best_match_references"
    colnames(df)[which(colnames(df)=="E_value")]<-"Best_match_E_value"
    colnames(df)[which(colnames(df)=="Bit_score")]<-"Best_match_bit_score"
    
    # Keep certain fields in the classification data frame.
    df<-df[,c("Sequence_name","Sequence","Best_match_references","Best_match_E_value",
              "Best_match_bit_score","Best_match_query_cover.mean",
              "Best_match_query_cover.SD","Best_match_PID.mean","Best_match_PID.SD")]
    
  }
  
  # Set NA values in rows without an E-value to suggest that no significant similarity
  # was found.
  df[is.na(df$Best_match_E_value),3:ncol(df)]<-"No significant similarity found"
  
  # Write out classified sequences.
  utils::write.csv(x=df,file=path_to_output_file,row.names=FALSE)
  
}