#' Merge Forward and Reverse DNA Sequence Reads
#'
#' @description Merges forward and reverse DNA sequence reads.
#' @details For each pair of forward and reverse DNA sequence reads, the reverse complement of the reverse read is internally derived using the [`reverse_complement`][reverse_complement()] function, and the read pair is merged into a single sequence if an overlap of at least the minimum length is found between the end of the forward read and the start of the reverse complement of the reverse read. If an overlap of the minimum length is not found, then an `NA` is returned for the merged read pair.
#' @param forward_reads A character vector of forward DNA sequence reads.
#' @param reverse_reads A character vector of reverse DNA sequence reads.
#' @param minimum_overlap Numeric. The minimum length of an overlap that must be found between the end of the forward read and the start of the reverse complement of the reverse read in order for a read pair to be merged. The default is `10`.
#' @returns A character vector of merged DNA sequence read pairs. `NA`s are returned for read pairs which could not be merged, which occurs when an overlap of at least the minimum length is not found between the end of the forward read and the start of the reverse complement of the reverse read.
#' @seealso
#' [`truncate_and_merge_pairs`][truncate_and_merge_pairs()] for truncating and merging forward and reverse DNA sequence reads.
#' @examples
#' merge_pairs(forward_reads=c("CCTTACGAATCCTGT","TTCTCCACCCGCGGATA","CGCCCGGAGTCCCTGTAGTA"),
#'             reverse_reads=c("GACAAACAGGATTCG","CAATATCCGCGGGTG","TACTACAGGGACTCC"))
#' @export
merge_pairs<-function(forward_reads,reverse_reads,minimum_overlap=10){
  
  # Throw an error if the forward reads are not a character vector.
  if(!is.character(forward_reads)) stop("The forward reads must be a character vector.")
  
  # Throw an error if the forward reads contain NAs.
  if(any(is.na(forward_reads))) stop("The forward reads cannot contain NAs.")
  
  # Throw an error if the reverse reads are not a character vector.
  if(!is.character(reverse_reads)) stop("The reverse reads must be a character vector.")
  
  # Throw an error if the reverse reads contain NAs.
  if(any(is.na(reverse_reads))) stop("The reverse reads cannot contain NAs.")
  
  # Throw an error if the forward and reverse read vectors are not the same length.
  if(length(forward_reads)!=length(reverse_reads)) stop("The forward and reverse read vectors must be the same length.")
  
  # Throw an error if the minimum overlap is not of class numeric.
  if(!is.numeric(minimum_overlap)) stop("The minimum overlap must be class numeric.")
  
  # Throw an error if the minimum overlap contains multiple elements.
  if(length(minimum_overlap) > 1) stop("The minimum overlap cannot have multiple elements.")
  
  # Throw an error if the minimum overlap is not an integer value.
  if(round(x=minimum_overlap,digits=0)!=minimum_overlap) stop("The minimum overlap must be an integer value.")
  
  # Throw an error if the minimum overlap is less than one.
  if(minimum_overlap < 1) stop("The minimum overlap must be greater than zero.")
  
  # Define function for merging a read pair.
  merge_pair<-function(forward_read,reverse_read,min_overlap){
    
    # Get the reverse complement of the reverse read.
    rev_comp<-reverse_complement(sequence=reverse_read)
    
    # Set the initial overlap length to the minimum overlap.
    overlap<-min_overlap
    
    # Set the finished flag to FALSE.
    finished<-FALSE
    
    # While the finished flag is set to FALSE.
    while(!finished){
      
      # If the overlap length does not exceed the length of the
      # reverse complement of the reverse sequence.
      if(overlap <= nchar(rev_comp)){
        
        # Get the first overlap-length-number of base pairs from the
        # reverse complement of the reverse read.
        pattern<-paste0(substr(x=rev_comp,start=1,stop=overlap),"$")
        
        # If the forward read ends in the first overlap-length-number of
        # base pairs from the reverse complement of the reverse read.
        if(grepl(pattern=pattern,x=forward_read)){
          
          # Get the portion of the reverse complement of the reverse read
          # to add to the forward read to merge the read pairs.
          sec_part<-substr(x=rev_comp,start=overlap+1,stop=nchar(rev_comp))
          
          # Add the portion of the reverse complement of the reverse read
          # to the forward read to merge the read pairs.
          merged<-paste0(forward_read,sec_part)
          
          # Set the finished flag to TRUE.
          finished<-TRUE
          
        } else {
          
          # If the forward read does not end in the first overlap-length-number
          # of base pairs from the reverse complement of the reverse read.
          
          # Increase the overlap length by one base pair.
          overlap<-overlap+1
          
        }
        
      } else {
        
        # If the overlap length exceeds the length of the
        # reverse complement of the reverse sequence.
        
        # The reads cannot be merged. Set the sequence to return to NA.
        merged<-NA
        
        # Set the finished flag to TRUE.
        finished<-TRUE
        
      }
      
    }
    
    # Return the merged read pair.
    return(merged)
    
  }
  
  # Apply the merge pair function to multiple read pairs.
  merged_pairs<-sapply(X=1:length(forward_reads),
                       FUN=function(i) merge_pair(forward_read=forward_reads[i],
                                                  reverse_read=reverse_reads[i],
                                                  min_overlap=minimum_overlap))
  
  # Return the merged read pairs.
  return(merged_pairs)
  
}