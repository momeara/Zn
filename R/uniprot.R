# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:


#' @export
prepare_data <- function(data, col){
	if(!("data.frame" %in% class(data))){
		data <- data.frame(col=data)
		names(data) <- c(col)
	}
	if(nrow(data) == 0){
		cat("WARNING: No uniprot_entry values to convert\n")
	}
	data
}

uniprot_entry_web_lookup_debug <- function(
	request_bin_host,
	columns,
	user_agent_arg # e.g. user_agent("httr name@example.com")
){
	column_arg <- paste0(columns, collapse=",")
	r <- httr::GET(
		request_bin_host,
		user_agent_arg,
		query = list(
			query="DRD2_HUMAN",
			format = 'tab',
			columns = column_arg)) %>%
	content
}

#' See http://www.uniprot.org/help/uniprotkb_column_names for available columns
uniprot_entry_web_lookup <- function(
	uniprot_entries,
	columns,
	user_agent_arg,
	verbose=F
){
	column_arg <- paste0(columns, collapse=",")
	if(verbose){
		cat("Getting data for ", length(uniprot_entries), " uniprot entries ... \n", sep="")
	}
	r2 <- plyr::adply(
		dplyr::data_frame(uniprot_entry = uniprot_entries),
		1,
		function(df){
			if(verbose){
				cat("Looking up entry for '", df$uniprot_entry, "' ... ", sep="")
			}
			r <- httr::GET(
				"http://www.uniprot.org/uniprot/",
				user_agent_arg,
				query = list(
					query=df$uniprot_entry,
					format = 'tab',
					columns = column_arg)) %>%
				httr::content()
			if(r %>% is.null){
				if(verbose){
					cat(" MISSING\n")
				}
				return(dplyr::data_frame())
			} else {
				if(verbose){
					cat(" GOT IT\n")
				}
				r %>%
					readr::read_tsv() %>%
					dplyr::filter(`Entry name` == df$uniprot_entry) %>%
					head(1) %>%
					return
			}
		})
	if(nrow(r2) == 0){
		cat("WARNING: no entries were found\n")
	}
	r2
}


#' @export
uniprot_entry_full <- function(
	uniprot_entries,
	user_agent_arg
){
	uniprot_entries_fname <- tempfile()
	write.table(uniprot_entries, uniprot_entries_fname, quote=F, col.names=F, row.names=F, sep="")
	uniprot_entries_f <- httr::upload_file(uniprot_entries_fname, type="text/plain")
	z <- httr::POST(
		"http://www.uniprot.org/batch/",
		user_agent_arg,
		body=list(
			file = uniprot_entries_f,
			format='txt')) %>%
		httr::content %>%
		stringr::str_split("//\n")
	file.remove(uniprot_entries_fname)
	z
}

# retrieve taxa from uniprot descending from given ancestor
# as of Feb 2016, columns are
# Taxon	Mnemonic	Scientific name	Common name	Synonym	Other Names	Reviewed	Rank	Lineage	Parent	Virus hosts
# 89096	ABRANA	brothrix andinus	Andean mouse	Abrothrix andinus (Philippi, 1858); Akodon andinus; Andean akodont	reviewed	Species	Eukaryota; Metazoa; Chordata; Craniata; Vertebrata; Euteleostomi; Mammalia; Eutheria; Euarchontoglires; Glires; Rodentia; Sciurognathi; Muroidea; Cricetidae; Sigmodontinae; Abrothrix	156196
uniprot_taxa <- function(
	user_agent, # e.g. user_agent("httr name@example.com")
	ancestor=7711, # chordate
	reviewed=TRUE,
	format='tab',
	compress=TRUE,
	verbose=F
){
	library(httr)
	library(plyr)
	library(dplyr)
	library(readr)
	if(verbose){
		cat("Retriving taxa descending from '", ancestor, "' from uniprot...\n")
	}
	query_arg <- paste0(
		"rank:species ",
		ifelse(reviewed, "uniprot:(reviewed:yes) ", ""),
		"ancestor:", ancestor)

	r <- httr::GET(
		"http://www.uniprot.org/taxonomy/",
		user_agent,
		query = list(
			query=query_arg,
			format = 'tab',
			compress=ifelse(compress, 'yes', 'no'))) %>%
		httr::content()

	if(r %>% is.null){
		if(verbose){
			cat(" None retrieved\n")
		}
		return(data_frame())
	}

	r %>%
		rawConnection %>%
		gzcon %>%
		readr::read_tsv() %>%
		return
}


