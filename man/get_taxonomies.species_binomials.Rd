% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_taxonomies.species_binomials.R
\name{get_taxonomies.species_binomials}
\alias{get_taxonomies.species_binomials}
\title{Get NCBI Taxonomies from Species Binomials}
\usage{
get_taxonomies.species_binomials(
  path_to_input_species_binomials,
  path_to_output_local_taxa_list,
  path_to_taxonomy_edits = NA,
  print_taxize_queries = TRUE
)
}
\arguments{
\item{path_to_input_species_binomials}{String specifying path to input species list with common and scientific names. The file should be in CSV format with a '.csv' extension. The file must contain the following fields: 'Common_Name', 'Scientific_Name'. Values in the 'Common_Name' field are optional. Values in the 'Scientific_Name' field are required.}

\item{path_to_output_local_taxa_list}{String specifying path to output species list with added NCBI taxonomies. The output file will be in CSV format, and the file path should have a '.csv' extension.}

\item{path_to_taxonomy_edits}{String specifying path to taxonomy edits file in CSV format with a '.csv' extension. The file must contain the following fields: 'Old_Taxonomy', 'New_Taxonomy', 'Notes'. Old taxonomies are replaced with new taxonomies in the order the records appear in the file. The taxonomic levels in the 'Old_Taxonomy' and 'New_Taxonomy' fields should be delimited by a semi-colon. If no taxonomy edits are desired, then set this variable to \code{NA} (the default).}

\item{print_taxize_queries}{Logical. Whether taxa queries should be printed. The default is \code{TRUE}.}
}
\description{
Remotely fetches taxonomies from the NCBI taxonomy database for a list of species binomials.
}
