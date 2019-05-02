# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:


#' @export
gene_info <- function(
	gene_names,
	output_fields=c("name", "description", "major_class_name", "sub_class_name"),
  ...
){
	dplyr::data_frame(gene_name=gene_names) %>% plyr::adply(1, function(row){
		tryCatch({
			zinc_REST(
					path=paste0(
						"genes/", row$gene_name[1], ".csv?",
						"output_fields=", paste(output_fields, collapse=" ")),
					...)
		}, error=function(e){
			print(paste0("ERROR: ", e))
			return(data.frame())
		})
	})
}

#' @export
target_info <- function(
	uniprot_entries,
	output_fields=c("name", "description", "major_class_name", "sub_class_name"),
	...
){
	dplyr::data_frame(uniprot_entry=uniprot_entries) %>% plyr::adply(1, function(row){
		zinc_REST(
			path=paste0(
				"genes/", row$uniprot_entry[1], ".csv?",
				"output_fields=", paste(output_fields, collapse=" ")),
			...)
	})
}
