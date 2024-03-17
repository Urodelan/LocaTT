#' Filter DNA Sequences by PCR Replicates
#'
#' @description Filters DNA sequences by minimum read count within a PCR replicate, minimum proportion within a PCR replicate, and number of detections across PCR replicates.
#' @details For each set of input polymerase chain reaction (PCR) replicate FASTA files associated with a sample, writes out DNA sequences which are detected across a minimum number of PCR replicates (`minimum_PCR_replicates` argument). Detection within a PCR replicate is defined as a sequence having at least a minimum read count *and* exceeding a minimum proportion of reads (`minimum_reads.sequence` and `minimum_proportion.sequence` arguments, respectively). When `binomial_test.enabled = TRUE`, a sequence must significantly exceed the minimum proportion within a PCR replicate at the provided alpha level (`binomial_test.alpha_level` argument) based on a one-sided binomial test (*i.e.*, [`binomial_test`][binomial_test()] with `alternative = "greater"`). Within a PCR replicate, p-values can be adjusted for multiple hypothesis testing by setting the `binomial_test.p.adjust.method` argument (see `stats::p.adjust.methods` and [`p.adjust`][stats::p.adjust()] in the [`stats`][stats::stats] package). PCR replicates which contain fewer than a minimum number of reads are discarded (`minimum_reads.PCR_replicate` argument) and do not contribute detections to any sequence. \cr \cr
#' DNA sequences in the input FASTA files are assumed to be summarized by frequency of occurrence, with each FASTA header line beginning with "Frequency: " and followed by the sequence's read count. Output FASTA files from [`truncate_and_merge_pairs`][truncate_and_merge_pairs()] have this format and can be used directly with this function. Each input FASTA file is assumed to contain the DNA sequence reads for a single PCR replicate for a single sample. \cr \cr
#' For pipeline calibration purposes, a data frame containing unfiltered DNA sequences with their read counts, proportions, and p-values in each PCR replicate is invisibly returned (see return value section). While the primary output of this function is the written CSV file of filtered sequences (described below), the invisibly returned data frame of unfiltered sequences can be helpful when calibrating or troubleshooting filtering parameters. To aid in troubleshooting filtering parameters, the data frame is invisibly returned even if the error "Filtering removed all sequences" is received. \cr \cr
#' For the primary output, this function writes a CSV file of filtered DNA sequences with the following field definitions:
#' * Sample: The sample name.
#' * Sequence: The DNA sequence.
#' * Detections_across_PCR_replicates: The number of PCR replicates the sequence was detected in.
#' * Read_count_by_PCR_replicate: The sequence's read count in each PCR replicate the sequence was detected in.
#' * Sequence_read_count: The sequence's total read count across the PCR replicates the sequence was detected in. Calculated as the sum of the read counts in the Read_count_by_PCR_replicate field.
#' * Sample_read_count: The sample's total read count across all sequences detected in the PCR replicates. Calculated as the sum of the read counts in Sequence_read_count field associated with the sample.
#' * Proportion_of_sample: The proportion of sample reads comprised by the sequence. Calculated by dividing the Sequence_read_count field by the Sample_read_count field. Equivalent to the weighted average of the sequence's proportion in each PCR replicate, with weights given by the proportion of the sample's total reads contained in each PCR replicate.
#' @param input_files A character vector of file paths to input FASTA files. DNA sequences in the input FASTA files are assumed to be summarized by frequency of occurrence, with each FASTA header line beginning with "Frequency: " and followed by the sequence's read count. Output FASTA files from [`truncate_and_merge_pairs`][truncate_and_merge_pairs()] have this format and can be used directly with this function. Each input FASTA file is assumed to contain the DNA sequence reads for a single PCR replicate for a single sample.
#' @param samples A character vector of sample identifiers, with one element for each element of `input_files`.
#' @param PCR_replicates A character vector of PCR replicate identifiers, with one element for each element of `input_files`.
#' @param output_file String specifying path to output file of filtered sequences in CSV format.
#' @param minimum_reads.PCR_replicate Numeric. PCR replicates which contain fewer reads than this value are discarded and do not contribute detections to any sequence. The default is `1` (*i.e.*, no PCR replicates discarded).
#' @param minimum_reads.sequence Numeric. For a sequence to be considered detected within a PCR replicate, the sequence's read count within the PCR replicate must match or exceed this value. The default is `1` (*i.e.*, no filtering by minimum read count within PCR replicates).
#' @param minimum_proportion.sequence Numeric. For a sequence to be considered detected within a PCR replicate, the proportion of reads in the PCR replicate comprised by the sequence must exceed this value. If `binomial_test.enabled = TRUE`, then this argument is used as the null hypothesis for a one-sided binomial test, and a significance test is used to determine whether the minimum proportion requirement for detection is satisfied instead. See the `binomial_test.enabled` argument below. The default is `0.005` (*i.e.*, 0.5%). To disable sequence filtering by minimum proportion within PCR replicates, set to `0`.
#' @param binomial_test.enabled Logical. If `TRUE` (the default), then for a sequence to be considered detected within a PCR replicate, the proportion of reads in the PCR replicate comprised by the sequence must significantly exceed the value of the `minimum_proportion.sequence` argument at the provided alpha level (`binomial_test.alpha_level` argument) based on a one-sided binomial test (*i.e.*, [`binomial_test`][binomial_test()] with `alternative = "greater"`). Optionally, p-values within a PCR replicate can be adjusted for multiple hypothesis testing by setting the `binomial_test.p.adjust.method` argument below. To disable significance testing, set to `FALSE` (minimum proportion filtering will still occur if `minimum_proportion.sequence > 0`, see above).
#' @param binomial_test.p.adjust.method String specifying the p-value adjustment method for multiple hypothesis testing. p-value adjustments are performed within each PCR replicate for each sample. Passed to the `method` argument of [`p.adjust`][stats::p.adjust()] in the [`stats`][stats::stats] package. Available methods are contained within the `stats::p.adjust.methods` vector. If `"none"` (the default), then p-value adjustments are not performed. Ignored if `binomial_test.enabled = FALSE`.
#' @param binomial_test.alpha_level Numeric. The alpha level used in deciding whether the proportion of reads in a PCR replicate comprised by a sequence significantly exceeds a minimum threshold required for detection. See the `binomial_test.enabled` argument. The default is `0.05`. Ignored if `binomial_test.enabled = FALSE`.
#' @param minimum_PCR_replicates Numeric. The minimum number of PCR replicates in which a sequence must be detected in order to be considered present in a sample (*i.e.*, not erroneous). The default is `2`.
#' @param delimiter.read_counts String specifying the delimiter between PCR replicate identifiers and sequence read counts in the Read_count_by_PCR_replicate field of the output CSV file (see details section). The default is `": "`.
#' @param delimiter.PCR_replicates String specifying the delimiter between PCR replicates in the Read_count_by_PCR_replicate field of the output CSV file (see details section). The default is `", "`.
#' @returns Invisibly returns a data frame containing unfiltered DNA sequences with their read counts, proportions, and p-values in each PCR replicate. While the primary output of this function is the written CSV file of filtered sequences described in the details section, the invisibly returned data frame of unfiltered sequences can be helpful when calibrating or troubleshooting filtering parameters. To aid in troubleshooting filtering parameters, the data frame is invisibly returned even if the error "Filtering removed all sequences" is received. Field definitions for the invisibly returned data frame of unfiltered sequences are:
#' * Sample: The sample name.
#' * PCR_replicate: The PCR replicate identifier.
#' * Sequence: The DNA sequence.
#' * Read_count.sequence: The sequence's read count within the PCR replicate.
#' * Read_count.PCR_replicate: The number of reads in the PCR replicate.
#' * Proportion_of_PCR_replicate.observed: The proportion of reads in the PCR replicate comprised by the sequence.
#' * Proportion_of_PCR_replicate.null (Field only present if `binomial_test.enabled = TRUE`): The null hypothesis for a one-sided binomial test (inherited from the `minimum_proportion.sequence` argument). See the p.value field below.
#' * p.value (Field only present if `binomial_test.enabled = TRUE`): The p-value from a one-sided binomial test of whether the proportion of reads in the PCR replicate comprised by the sequence exceeds the null hypothesis (_i.e._, [`binomial_test`][binomial_test()] with `alternative = "greater"`).
#' * p.value.adjusted (Field only present if `binomial_test.enabled = TRUE`): The p-value from the one-sided binomial test adjusted for multiple comparisons within each PCR replicate for each sample. See the p.value_adjustment_method field below.
#' * p.value_adjustment_method (Field only present if `binomial_test.enabled = TRUE`): The p-value adjustment method (inherited from the `binomial_test.p.adjust.method` argument).
#' @seealso
#' [`binomial_test`][binomial_test()] for performing vectorized one-sided binomial tests. \cr \cr
#' [`truncate_and_merge_pairs`][truncate_and_merge_pairs()] for truncating and merging read pairs prior to sequence filtering. \cr \cr
#' [`local_taxa_tool`][local_taxa_tool()] for performing geographically-conscious taxonomic assignment of filtered sequences.
#' @references A manuscript describing these methods is in preparation.
#' @examples
#' print("Insert example here.")
#' @export
filter_sequences<-function(input_files,samples,PCR_replicates,output_file,minimum_reads.PCR_replicate=1,minimum_reads.sequence=1,minimum_proportion.sequence=0.005,binomial_test.enabled=TRUE,binomial_test.p.adjust.method="none",binomial_test.alpha_level=0.05,minimum_PCR_replicates=2,delimiter.read_counts=": ",delimiter.PCR_replicates=", "){
  
  # Check user inputs.
  
  # Check input files argument.
  ## Throw an error if the input files are not a character vector.
  if(!is.character(input_files)) stop("The input files argument must be a character vector.")
  ## Throw an error if the input files contain any NAs.
  if(any(is.na(input_files))) stop("The input files argument cannot contain NAs.")
  
  # Check samples argument.
  ## Throw an error if the samples argument is not a character vector.
  if(!is.character(samples)) stop("The samples argument must be a character vector.")
  ## Throw an error if the samples argument contains any NAs.
  if(any(is.na(samples))) stop("The samples argument cannot contain NAs.")
  
  # Check PCR replicates argument.
  ## Throw an error if the PCR replicates argument is not a character vector.
  if(!is.character(PCR_replicates)) stop("The PCR replicates argument must be a character vector.")
  ## Throw an error if the samples argument contains any NAs.
  if(any(is.na(PCR_replicates))) stop("The PCR replicates argument cannot contain NAs.")
  
  # Throw an error if the length of the samples argument
  # does not match the number of input files.
  if(length(input_files)!=length(samples)) stop("The length of the samples argument must match the number of input files.")
  
  # Throw an error if the length of the PCR replicates argument
  # does not match the number of input files.
  if(length(input_files)!=length(PCR_replicates)) stop("The length of the PCR replicates argument must match the number of input files.")
  
  # Check output file argument.
  ## Throw an error if the output file argument is not a character vector.
  if(!is.character(output_file)) stop("The output file argument must be a character string.")
  ## Throw an error if the output file argument contains multiple elements.
  if(length(output_file) > 1) stop("The output file argument cannot have multiple elements.")
  ## Throw an error if the output file argument is NA.
  if(is.na(output_file)) stop("The output file argument cannot be NA.")
  
  # Check minimum_reads.PCR_replicate argument.
  ## Throw an error if the minimum_reads.PCR_replicate argument is not numeric.
  if(!is.numeric(minimum_reads.PCR_replicate)) stop("The minimum_reads.PCR_replicate argument must be numeric.")
  ## Throw an error if the minimum_reads.PCR_replicate argument contains multiple elements.
  if(length(minimum_reads.PCR_replicate) > 1) stop("The minimum_reads.PCR_replicate argument cannot have multiple elements.")
  ## Throw an error if the minimum_reads.PCR_replicate argument is NA.
  if(is.na(minimum_reads.PCR_replicate)) stop("The minimum_reads.PCR_replicate argument cannot be NA.")
  ## Throw an error if the minimum_reads.PCR_replicate argument is not an integer value.
  if(round(x=minimum_reads.PCR_replicate,digits=0)!=minimum_reads.PCR_replicate) stop("The minimum_reads.PCR_replicate argument must be an integer value.")
  ## Throw an error if the minimum_reads.PCR_replicate argument is less than one.
  if(minimum_reads.PCR_replicate < 1) stop("The minimum_reads.PCR_replicate argument must be greater than zero.")
  
  # Check minimum_reads.sequence argument.
  ## Throw an error if the minimum_reads.sequence argument is not numeric.
  if(!is.numeric(minimum_reads.sequence)) stop("The minimum_reads.sequence argument must be numeric.")
  ## Throw an error if the minimum_reads.sequence argument contains multiple elements.
  if(length(minimum_reads.sequence) > 1) stop("The minimum_reads.sequence argument cannot have multiple elements.")
  ## Throw an error if the minimum_reads.sequence argument is NA.
  if(is.na(minimum_reads.sequence)) stop("The minimum_reads.sequence argument cannot be NA.")
  ## Throw an error if the minimum_reads.sequence is not an integer value.
  if(round(x=minimum_reads.sequence,digits=0)!=minimum_reads.sequence) stop("The minimum_reads.sequence argument must be an integer value.")
  ## Throw an error if the minimum_reads.sequence argument is less than one.
  if(minimum_reads.sequence < 1) stop("The minimum_reads.sequence argument must be greater than zero.")
  
  # Check minimum_proportion.sequence argument.
  ## Throw an error if the minimum_proportion.sequence argument is not numeric.
  if(!is.numeric(minimum_proportion.sequence)) stop("The minimum_proportion.sequence argument must be numeric.")
  ## Throw an error if the minimum_proportion.sequence argument contains multiple elements.
  if(length(minimum_proportion.sequence) > 1) stop("The minimum_proportion.sequence argument cannot have multiple elements.")
  ## Throw an error if the minimum_proportion.sequence argument is NA.
  if(is.na(minimum_proportion.sequence)) stop("The minimum_proportion.sequence argument cannot be NA.")
  ## Throw an error if the minimum_proportion.sequence argument is less than zero.
  if(minimum_proportion.sequence < 0) stop("The minimum_proportion.sequence argument cannot be less than zero.")
  ## Throw an error if the minimum_proportion.sequence argument is greater than one.
  if(minimum_proportion.sequence > 1) stop("The minimum_proportion.sequence argument cannot be greater than one.")
  
  # Check binomial_test.enabled argument.
  ## Throw an error if the binomial_test.enabled argument is not logical.
  if(!is.logical(binomial_test.enabled)) stop("The binomial_test.enabled argument must be logical.")
  ## Throw an error if the binomial_test.enabled argument contains multiple elements.
  if(length(binomial_test.enabled) > 1) stop("The binomial_test.enabled argument cannot have multiple elements.")
  ## Throw an error if the binomial_test.enabled argument is NA.
  if(is.na(binomial_test.enabled)) stop("The binomial_test.enabled argument cannot be NA.")
  
  # Check p-value adjustment method.
  ## Throw an error if the p-value adjustment method is not a character string.
  if(!is.character(binomial_test.p.adjust.method)) stop("The p-value adjustment method must be a character string.")
  ## Throw an error if the p-value adjustment method has multiple elements.
  if(length(binomial_test.p.adjust.method) > 1) stop("The p-value adjustment method cannot have multiple elements.")
  ## Throw an error if the p-value adjustment method is not valid.
  if(!(binomial_test.p.adjust.method %in% stats::p.adjust.methods)) stop(paste0('The p-value adjustment method must be one of the following: "',paste(stats::p.adjust.methods,collapse='", "'),'"'))
  
  # Check binomial_test.alpha_level argument.
  ## Throw an error if the binomial_test.alpha_level argument is not numeric.
  if(!is.numeric(binomial_test.alpha_level)) stop("The binomial_test.alpha_level argument must be numeric.")
  ## Throw an error if the binomial_test.alpha_level argument contains multiple elements.
  if(length(binomial_test.alpha_level) > 1) stop("The binomial_test.alpha_level argument cannot have multiple elements.")
  ## Throw an error if the binomial_test.alpha_level argument is NA.
  if(is.na(binomial_test.alpha_level)) stop("The binomial_test.alpha_level argument cannot be NA.")
  ## Throw an error if the binomial_test.alpha_level argument is less than zero.
  if(binomial_test.alpha_level < 0) stop("The binomial_test.alpha_level argument cannot be less than zero.")
  ## Throw an error if the binomial_test.alpha_level argument is greater than one.
  if(binomial_test.alpha_level > 1) stop("The binomial_test.alpha_level argument cannot be greater than one.")
  
  # Check minimum_PCR_replicates argument.
  ## Throw an error if the minimum_PCR_replicates argument is not numeric.
  if(!is.numeric(minimum_PCR_replicates)) stop("The minimum_PCR_replicates argument must be numeric.")
  ## Throw an error if the minimum_PCR_replicates argument contains multiple elements.
  if(length(minimum_PCR_replicates) > 1) stop("The minimum_PCR_replicates argument cannot have multiple elements.")
  ## Throw an error if the minimum_PCR_replicates argument is NA.
  if(is.na(minimum_PCR_replicates)) stop("The minimum_PCR_replicates argument cannot be NA.")
  ## Throw an error if the minimum_PCR_replicates is not an integer value.
  if(round(x=minimum_PCR_replicates,digits=0)!=minimum_PCR_replicates) stop("The minimum_PCR_replicates argument must be an integer value.")
  ## Throw an error if the minimum_PCR_replicates argument is less than one.
  if(minimum_PCR_replicates < 1) stop("The minimum_PCR_replicates argument must be greater than zero.")
  
  # Check delimiter.read_counts argument.
  ## Throw an error if the delimiter.read_counts argument is not a character vector.
  if(!is.character(delimiter.read_counts)) stop("The delimiter.read_counts argument must be a character string.")
  ## Throw an error if the delimiter.read_counts argument contains multiple elements.
  if(length(delimiter.read_counts) > 1) stop("The delimiter.read_counts argument cannot have multiple elements.")
  ## Throw an error if the delimiter.read_counts argument is NA.
  if(is.na(delimiter.read_counts)) stop("The delimiter.read_counts argument cannot be NA.")
  
  # Check delimiter.PCR_replicates argument.
  ## Throw an error if the delimiter.PCR_replicates argument is not a character vector.
  if(!is.character(delimiter.PCR_replicates)) stop("The delimiter.PCR_replicates argument must be a character string.")
  ## Throw an error if the delimiter.PCR_replicates argument contains multiple elements.
  if(length(delimiter.PCR_replicates) > 1) stop("The delimiter.PCR_replicates argument cannot have multiple elements.")
  ## Throw an error if the delimiter.PCR_replicates argument is NA.
  if(is.na(delimiter.PCR_replicates)) stop("The delimiter.PCR_replicates argument cannot be NA.")
  
  # Begin operations.
  
  # Create a data frame containing the file names.
  files<-data.frame(File=input_files,
                    Sample=samples,
                    PCR=PCR_replicates,
                    stringsAsFactors=FALSE)
  
  # Throw an error of any duplicate combinations of sample and PCR replicate occur.
  if(any(duplicated(files[,c("Sample","PCR")]))) stop("There are duplicate combinations of sample and PCR replicate.")
  
  # Get unique sample names.
  sample.names<-unique(files$Sample)
  
  # Create an empty data frame for storing sequence information.
  df<-data.frame(NULL)
  
  # Loop through each sample.
  for(i in 1:length(sample.names)){
    
    # Get the sample name.
    sample.name<-sample.names[i]
    
    # Get the fasta files associated with the sample.
    sample.files<-files[files$Sample==sample.name,]
    
    # Loop through each fasta file associated with the sample.
    for(j in 1:nrow(sample.files)){
      
      # Get the record for the fasta file.
      sample.file<-sample.files[j,]
      
      # Read in the fasta file.
      sample<-read.fasta(file=sample.file$File)
      
      # Throw an error of not all header lines in the input file begin with "Frequency: ".
      if(!all(grepl(pattern="^Frequency: ",x=sample$Name))) stop(paste0('Not all header lines begin with "Frequency: " in input file: ',sample.file$File))
      
      # Get sequence counts from the sequence names.
      sample$Count<-sub(pattern="^Frequency: ",replacement="",x=sample$Name)
      
      # Throw an error if any counts cannot be converted to numeric.
      if(suppressWarnings(any(is.na(as.numeric(sample$Count))))) stop(paste0("Not all counts in the header lines are numeric in input file: ",sample.file$File))
      
      # Format counts as numeric.
      sample$Count<-as.numeric(sample$Count)
      
      # Throw an error if any counts are not integer values.
      if(!identical(round(x=sample$Count,digits=0),sample$Count)) stop(paste0("Not all counts are integer values in input file: ",sample.file$File))
      
      # Throw an error if any counts are not at least one.
      if(any(sample$Count < 1)) stop(paste0("Not all counts are greater than zero in input file: ",sample.file$File))
      
      # Create a field denoting the sample.
      sample$Sample<-sample.name
      
      # Create a field denoting the PCR replicate.
      sample$PCR<-sample.file$PCR
      
      # Get fields of interest.
      sample<-sample[,c("Sample","PCR","Sequence","Count")]
      
      # Add sample sequence information to the storage data frame.
      df<-rbind(df,sample)
      
    }
    
  }
  
  # Prepare for sequence filtering.
  
  # Calculate the total read count for each PCR replicate.
  ## Perform calculation.
  totals<-stats::aggregate(Count~Sample+PCR,data=df,FUN=sum)
  ## Rename the count field to total.
  colnames(totals)[colnames(totals)=="Count"]<-"Total"
  
  # Add total PCR read counts to the main data frame.
  df<-merge(x=df,y=totals)
  
  # Calculate proportion of reads each sequence comprises in each PCR replicate.
  df$Proportion<-df$Count/df$Total
  
  # If binomial testing is enabled.
  if(binomial_test.enabled){
    
    # Create field containing the null hypothesis.
    df$Proportion.null<-minimum_proportion.sequence
    
    # Perform binomial tests.
    df$p.value<-binomial_test(k=df$Count,n=df$Total,
                              p=minimum_proportion.sequence,
                              alternative="greater")
    
    # Apply p-value adjustment.
    ## Create empty storage field for adjusted p-values.
    df$p.value.adjusted<-NA
    ## Loop through each sample.
    for(i in unique(df$Sample)){
      ## Loop through each PCR replicate within the sample.
      for(j in unique(df$PCR[df$Sample==i])){
        ## Apply p-value adjustment for sample-PCR replicate combination.
        df$p.value.adjusted[df$Sample==i & df$PCR==j]<-stats::p.adjust(
          p=df$p.value[df$Sample==i & df$PCR==j],
          method=binomial_test.p.adjust.method)
      }
    }
    
    # Throw an error if p-value calculation failed for any sequences.
    if(any(is.na(df$p.value.adjusted))) stop("p-value calculation failed for some sequences.")
    
  }
  
  # Store unfiltered data set.
  
  # Store data frame prior to filtering.
  # This unfiltered data set will be invisibly returned by the function.
  df.unfiltered<-df
  
  # Re-order records in unfiltered data set.
  ## Order records by sample, PCR replicate, and sequence.
  df.unfiltered<-df.unfiltered[order(df.unfiltered$Sample,
                                     df.unfiltered$PCR,
                                     df.unfiltered$Sequence),]
  ## Reset row names.
  row.names(df.unfiltered)<-1:nrow(df.unfiltered)
  
  # Rename fields of the unfiltered data set.
  ## PCR replicate.
  colnames(df.unfiltered)[colnames(df.unfiltered)=="PCR"]<-"PCR_replicate"
  ## Sequence read count.
  colnames(df.unfiltered)[colnames(df.unfiltered)=="Count"]<-"Read_count.sequence"
  ## Sequence read count.
  colnames(df.unfiltered)[colnames(df.unfiltered)=="Total"]<-"Read_count.PCR_replicate"
  ## Sequence read count.
  colnames(df.unfiltered)[colnames(df.unfiltered)=="Proportion"]<-"Proportion_of_PCR_replicate.observed"
  ## Sequence read count.
  if(binomial_test.enabled) colnames(df.unfiltered)[colnames(df.unfiltered)=="Proportion.null"]<-"Proportion_of_PCR_replicate.null"
  
  # Add field for p-value adjustment method to the unfiltered data frame.
  if(binomial_test.enabled) df.unfiltered$p.value_adjustment_method<-binomial_test.p.adjust.method
  
  # Set expression to evaluate upon exiting the function (e.g., completion or error).
  # From this point forward in the function, the data frame of unfiltered sequences
  # will always be invisibly returned even if an error is encountered.
  on.exit(expr=return(invisible(df.unfiltered)),add=TRUE)
  
  # Perform sequence filtering.
  
  # Remove PCR replicates which have less than the minimum number of reads.
  df<-df[df$Total >= minimum_reads.PCR_replicate,]
  
  # Remove sequences which have fewer reads than the user-specified threshold.
  df<-df[df$Count >= minimum_reads.sequence,]
  
  # If binomial testing is enabled.
  if(binomial_test.enabled){
    
    # Subset to sequences whose observed proportions within PCR replicates significantly
    # exceed the user-defined threshold at the user-defined alpha level.
    df<-df[df$p.value.adjusted <= binomial_test.alpha_level,]
    
  } else { # If binomial testing is not enabled.
    
    # Subset to sequences whose observed proportions within PCR replicates
    # exceed the user-defined threshold.
    df<-df[df$Proportion > minimum_proportion.sequence,]
    
  }
  
  # Throw an error if filtering up to this point has removed all sequences.
  if(nrow(df)==0) stop("Filtering removed all sequences.")
  
  # Get sequence counts by PCR replicate.
  ## Get fields of interest.
  pcr.cnts<-df[,c("Sample","PCR","Sequence","Count")]
  ## Get total read count across PCR replicates for each sequence in each sample.
  pcr.cnts.totals<-stats::aggregate(Count~Sample+Sequence,data=pcr.cnts,FUN=sum)
  ## Rename totals field.
  colnames(pcr.cnts.totals)[colnames(pcr.cnts.totals)=="Count"]<-"Read_count.sequence"
  ## Create a joint PCR-count field.
  pcr.cnts$Label<-paste0(pcr.cnts$PCR,delimiter.read_counts,pcr.cnts$Count)
  ## Combine PCR-counts into a single field.
  pcr.cnts<-stats::aggregate(Label~Sample+Sequence,
                             data=pcr.cnts,
                             FUN=function(x) paste(sort(x),collapse=delimiter.PCR_replicates))
  ## Rename combined PCR-counts field.
  colnames(pcr.cnts)[colnames(pcr.cnts)=="Label"]<-"Read_count.PCR_replicates"
  ## Combine PCR-counts and sum of the PCR-counts into a single data frame.
  pcr.cnts<-merge(x=pcr.cnts,y=pcr.cnts.totals)
  
  # Create a contingency table of the number of PCR replicates each
  # sequence is detected in for each sample.
  ## Create table.
  df<-as.data.frame(table(df$Sample,df$Sequence),stringsAsFactors=FALSE)
  ## Adjust field names.
  colnames(df)<-c("Sample","Sequence","Number.PCR_replicates")
  
  # Subset to sequences detected in the minimum number of PCR replicates for each sample.
  df<-df[df$Number.PCR_replicates >= minimum_PCR_replicates,]
  
  # Throw an error if PCR filtering removed all sequences.
  if(nrow(df)==0) stop("Filtering removed all sequences.")
  
  # Add read counts to filtered data frame.
  df<-merge(x=df,y=pcr.cnts)
  
  # Get total number of reads for each sample.
  ## Calculate total number of reads for each sample.
  sample.cnts<-stats::aggregate(Read_count.sequence~Sample,data=df,FUN=sum)
  ## Update field name.
  colnames(sample.cnts)[colnames(sample.cnts)=="Read_count.sequence"]<-"Read_count.sample"
  
  # Add total number of reads for each sample to the filtered data frame.
  df<-merge(x=df,y=sample.cnts)
  
  # Calculate proportion of reads each sequence comprises within each sample.
  df$Proportion.sequence<-df$Read_count.sequence/df$Read_count.sample
  
  # Rename fields of filtered sequences data frame.
  colnames(df)<-c("Sample","Sequence","Detections_across_PCR_replicates","Read_count_by_PCR_replicate","Sequence_read_count","Sample_read_count","Proportion_of_sample")
  
  # Order records by sample and sequence.
  df<-df[order(df$Sample,df$Sequence),]
  
  # Write out filtered sequences.
  utils::write.csv(x=df,file=output_file,row.names=FALSE)
  
}