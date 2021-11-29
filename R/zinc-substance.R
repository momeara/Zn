# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:

#' Standardize ZINC ID to have format ZINC[0-9]{12}
#' Specifically take all numeric values, on the right, pad with zeros to get to 12 numbers and then prepend with 'ZINC'
#' @export
standardize_zinc_ids <- function(zinc_ids){
	plyr::llply(zinc_ids, function(zinc_id){
		if(is.na(zinc_id)){
			return(NA)
		} else {
			paste0("ZINC", zinc_id %>% stringr::str_extract("[0-9]+$") %>% stringr::str_pad(
			  width = 12,
			  side = "left",
			  pad = "0"))
		}
	}) %>% unlist()
}


#' Process results of substance_info to be more user-friendly
#' @export
process_substance_info <- function(sub_info){
	fields <- names(sub_info)
	if("gene_names" %in% fields){
		# gene_names looks like "[u'GENEID1', u'GENEID2', ...]"
		# split this into an R list
		# add an n_genes column
		sub_info <- sub_info %>%
			dplyr::mutate(
				gene_names = gene_names %>%
					stringr::str_replace("^[\\[]", "") %>%
					stringr::str_replace("[\\]]$", "") %>%
					stringr::str_replace("^u'", "") %>%
					stringr::str_replace("'$", "") %>%
					strsplit("', u'", fixed=T)) %>%
			dplyr::mutate(
				n_genes = gene_names %>% vapply(length, 1L))
	}

	if("features" %in% fields){
		sub_info <- sub_info %>%
			dplyr::mutate(
				aggregator = !is.na(features) & stringr::str_detect(features, 'aggregator'),
				drug_code =
					ifelse(stringr::str_detect(features, 'fda'), 210,
					ifelse(stringr::str_detect(features, 'world'), 211,
					ifelse(stringr::str_detect(features, 'investigationa'), 212,
					ifelse(stringr::str_detect(features, 'in-man'), 213,
					ifelse(stringr::str_detect(features, 'in-vivo'), 213,
					ifelse(stringr::str_detect(features, 'in-cells'), 215,
					ifelse(stringr::str_detect(features, 'in_vitro'), 216,
					NA))))))),
				drug_level =
					ifelse(stringr::str_detect(features, 'fda'), 'fda',
					ifelse(stringr::str_detect(features, 'world'), 'world',
					ifelse(stringr::str_detect(features, 'investigationa'), 'investigational',
					ifelse(stringr::str_detect(features, 'in-man'), 'in-man',
					ifelse(stringr::str_detect(features, 'in-vivo'), 'in-vivo',
					ifelse(stringr::str_detect(features, 'in-cells'), 'in-cells',
					ifelse(stringr::str_detect(features, 'in_vitro'), 'in-vitro',
					NA))))))),
				biological_code =
					ifelse(stringr::str_detect(features, 'endogenous'), 203,
					ifelse(stringr::str_detect(features, 'metabolite'), 202,
					ifelse(stringr::str_detect(features, 'biogenic'), 201,
					NA))),
				biological_level =
					ifelse(stringr::str_detect(features, 'endogenous'), 'endogenous',
					ifelse(stringr::str_detect(features, 'metabolite'), 'metabolite',
					ifelse(stringr::str_detect(features, 'biogenic'), 'biogenenic',
					NA))))
	}
	if("purchasable" %in% fields){
		sub_info <- sub_info %>%
			dplyr::mutate(purchasable_code = purchasable)
	}

	if("purchasablility" %in% fields){
		sub_info <- sub_info %>%
			dplyr::mutate(purchasable_level = purchasability)
	}

	sub_info
}



#' Lookup info on substances by zinc_ids
#' @export
substance_info <- function(
	zinc_ids,
	output_fields = c("zinc_id", "preferred_name", "smiles", "purchasability", "features"),
	batch_size = length(zinc_ids),
	raw = FALSE,
	...) {
	raw_results <- tibble::tibble(
		zinc_id = zinc_ids,
		batch_id = rep_len(1:ceiling(length(zinc_ids) / batch_size), length(zinc_ids))) %>%
		plyr::ddply(c("batch_id"), function(df) {
			zinc_REST(
				path = "substances.csv",
				post_data = list(
					`zinc_id-in` = paste(df$zinc_id, collapse = " "),
					output_fields = paste(output_fields, collapse = " ")),
				...)}) %>%
		dplyr::select(-batch_id)
	if (!raw) {
		results <- process_substance_info(raw_results)
	} else {
	  results <- raw_results
	}
	results
}


#' Search for substances by search terms
#'
#'   http://zinc15.docking.org/substances/search/?q=N%23CC1%3DCC%3DC%28C%3DC1%29C%28N1C%3DNC%3DN1%29C1%3DCC%3DC%28C%3DC1%29C%23N
#'   http://zinc15.docking.org/substances/search/?count=all&output_format=cvs&q=N%23CC1%3DCC%3DC%28C%3DC1%29C%28N1C%3DNC%3DN1%29C1%3DCC%3DC%28C%3DC1%29C%23N&output_fields=zinc_id%20preferred_name%20smiles%20purchasability%20features]
#' @export
search_for_substances <- function(
	search_terms,
	output_fields = c("zinc_id", "preferred_name", "smiles", "purchasability", "features"),
	raw = FALSE,
	...) {
		raw_results <- tibble::tibble(search_term = search_terms) %>%
			plyr::adply(1, function(row) {
  		zinc_REST(
  			path = "substances/search",
  			query = list(
  				output_format = "csv",
  				q = row$search_term[1],
  				output_fields = paste(output_fields, collapse = " ")),
  			...)
  	})
	if (!raw) {
		results <- process_substance_info(raw_results)
	} else {
	  results <- raw_results
	}
	results
}

#' Search for substances using smiles and fine grained tolerance criteria
#' @export
resolve_substances <- function(
  smiles,
	output_fields = c("zinc_id", "smiles", "preferred_name", "purchasability", "features"),
	structures = TRUE,
	names = TRUE,
	suppliers = FALSE,
	analogs = FALSE,
	retired = FALSE,
	charges = FALSE,
	scaffolds = FALSE,
	fulltext = FALSE,
	multiple = FALSE,
	raw = FALSE,
	...
) {
	raw_results <- zinc_REST(
			path = "substances/resolved/",
			post_data = list(
				paste = paste(smiles, collapse = "\n"),
				output_fields = paste(output_fields, collapse = " "),
				output_format = 'csv',
				structures = structures,
				names = names,
				suppliers = suppliers,
				analogs = analogs,
				retired = retired,
				charges = charges,
				scaffolds = scaffolds,
				fulltext = fulltext,
				multiple = multiple),
			...)
	if (!raw) {
		results <- process_substance_info(raw_results)
	} else {
	  results <- raw_results
	}
	results
}

#' Search for analogs of the given substances
#' @export
substance_analogs <- function(
	zinc_ids,
	output_fields=c("zinc_id", "smiles", "preferred_name", "purchasability", "similarity", "features"),
	fingerprint = "ecfp4_fp",
	similarity = "tanimoto",
	threshold = 30,
	subsets = NULL,
	ref_batch_size = 20 * 10^6,
	ref_page = NULL,
	max_zinc_id = 1.2 * 10^9,
	raw = FALSE,
	verbose = FALSE,
	...
){

	if (is.null(subsets)) {
		path <- "substances.csv"
	} else {
		path <- paste0("substances/subsets/", subsets, ".csv")
	}

	plyr::ldply(zinc_ids, function(zinc_id) {
		post_data <- list(
			output_fields = paste(output_fields, collapse = " "))
		post_data[[paste0(fingerprint, "-", similarity, "-", threshold)]] <- zinc_id
		if(is.null(ref_batch_size)) {
			raw_results <- zinc_REST(
				path = path,
				post_data = post_data,
				verbose = verbose,
				...)
		} else {

			# search over slices of the zinc database to avoid time out errors
			if (is.null(ref_page)) {
				ref_page <- 1
			}
			done <- FALSE
			raw_results <- NULL
			while (!done) {
				post_data['zinc_id-between']= paste(
					format(ref_batch_size * (ref_page-1) + 1, scientific=FALSE),
					format(min(ref_batch_size * ref_page, max_zinc_id), scientific=FALSE))
				result_page <- zinc_REST(
					path=path,
					post_data=post_data,
					verbose=verbose,
					...)

				if (verbose) {
					cat("Found ", nrow(result_page), " analogs.\n", sep="")
				}

				if(is.null(result_page) || nrow(result_page) == 0){
					done <- TRUE
				} else {
					raw_results <- raw_results %>% rbind(result_page)
					if(ref_batch_size * ref_page >= max_zinc_id){
						done <- TRUE
					} else {
						ref_page <- ref_page + 1
					}
				}
			}
		}
		if(!raw){
			results <- process_substance_info(raw_results)
		} else{
			results <- raw_results
		}
		results <- results %>%
			dplyr::mutate(
				query_zinc_id=zinc_id)
	})
}
