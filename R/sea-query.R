# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:



#' @export
sea_base_url <- function(){"http://sea16.docking.org"}



#' @export
sea_score <- function(
	ref_set,
	ref_smi,
	query_set,
	query_smi,
	run_tag,
	verbose = F,
	url_path="custom/search/single"
){
	input_base <- paste(tempfile(), run_tag, sep = "_")
	url <- paste(sea_base_url(), url_path, sep = "/") %>%
		httr::parse_url() %>%
		httr::build_url()

	if (verbose) {
		cat("url: ", url, "\n")
	}

	if (is.character(ref_smi)) {
		if (!stringr::str_detect(ref_smi, ".smi$")) {
			cat("WARNING: referene .smi file '", ref_smi, "', does not end in '.smi'\n")
		}
		ref_smi_fname <- ref_smi
	} else {
		ref_smi_fname <- paste0(input_base, ".ref.smi")
		cat("writing ref_smi to -> '", ref_smi_fname, "' ... ", sep = "")
		write.table(
			ref_smi %>% dplyr::select(smiles, compound),
			ref_smi_fname, quote = F, sep = " ", row.names = F, col.names = F)
		cat("DONE\n")
	}

	if (is.character(query_smi)) {
		if (!stringr::str_detect(query_smi, ".smi$")) {
			cat("WARNING: query .smi file '", query_smi, "', does not end in '.smi'\n")
		}
		query_smi_fname <- query_smi
	} else {
		query_smi_fname <- paste0(input_base, ".query.smi")
		cat("writing query_smi to -> '", query_smi_fname, "' ... ", sep = "")
		write.table(
			query_smi %>% dplyr::select(smiles, compound),
			query_smi_fname, quote = F, sep = " ", row.names = F, col.names = F)
		cat("DONE\n")
	}

	ref_smi_f <- httr::upload_file(ref_smi_fname, type="application/smil")
	query_smi_f <- httr::upload_file(query_smi_fname, type="application/smil")

	tryCatch({
		request <- httr::POST(
			url,
			httr::add_headers(Accept = "application/json"),
			body = list(
				target_file = ref_smi_f,
				query_file = query_smi_f,
				as_set = "yes",
				wait = "yes")) %>%
			httr::content()
	}, error = function(e) {
		cat("ERROR parsing url='", url, "'\n", sep = "")
		print(e)
		return(data.frame())
	})
}

#source("~/work/sea/scripts/sea.R")
#source("~/work/sea/scripts/zinc_tools.R")
#z <- sea_target_vs_target(
#	ref_set = "target1",
#	ref_smi = data.frame(compound=c("ZINC4258247"), smiles=c("N[C@@H](Cc1cc(I)c(Oc2ccc(O)cc2)c(I)c1)C(=O)O")),
#	query_set = "target2",
#	query_smi = data.frame(compound=c("ZINC4258247"), smiles=c("N[C@@H](Cc1cc(I)c(Oc2ccc(O)cc2)c(I)c1)C(=O)O")),
#	run_tag = "test1",
#	verbose=T)

#' WARNING this may not work well
#' @export
tc_matrix <- function(
	ref_smiles,
	query_smiles,
	run_tag,
	verbose = F,
	url_path="tools/tc_matrix"
){
	input_base <- paste(tempfile(), run_tag, sep = "_")
	url <- paste(sea_base_url(), url_path, sep = "/") %>%
		httr::parse_url() %>%
		httr::build_url()

	if (verbose) {
		cat("url: ", url, "\n")
	}

	if (is.character(ref_smiles)) {
		if (!stringr::str_detect(ref_smiles, ".smi$")) {
			cat("WARNING: referene .smi file '", ref_smiles, "', does not end in '.smi'\n")
		}
		ref_smi_fname <- ref_smiles
	} else {
		ref_smi_fname <- paste0(input_base, ".ref.smi")
		cat("writing ref_smi to -> '", ref_smi_fname, "' ... ", sep = "")
		write.table(
			ref_smiles %>% dplyr::select(smiles, compound),
			ref_smi_fname, quote = F, sep = " ", row.names = F, col.names = F)
		cat("DONE\n")
	}

	if (is.character(query_smiles)) {
		if (!stringr::str_detect(query_smiles, ".smi$")) {
			cat("WARNING: query .smi file '", query_smi, "', does not end in '.smi'\n")
		}
		query_smi_fname <- query_smiles
	} else {
		query_smi_fname <- paste0(input_base, ".query.smi")
		cat("writing query_smi to -> '", query_smi_fname, "' ... ", sep = "")
		write.table(
			query_smiles %>% dplyr::select(smiles, compound),
			query_smi_fname, quote = F, sep = " ", row.names = F, col.names = F)
		cat("DONE\n")
	}

	ref_smi_f <- httr::upload_file(ref_smi_fname, type = "application/smil")
	query_smi_f <- httr::upload_file(query_smi_fname, type = "application/smil")

	tryCatch({
		request <- httr::POST(
			url,
			httr::add_headers(Accept = "application/json"),
			body = list(
				target_file = ref_smi_f,
				query_file = query_smi_f,
				wait = "yes")) %>%
			httr::content()
	}, error = function(e) {
		cat("ERROR parsing url='", url, "'\n", sep = "")
		print(e)
		return(data.frame())
	})
}

#' compute sea score given uniprot entries for targets
#' for the substance fields
#' @export
sea_target_vs_target <- function(
	ref_target,
	query_target,
	substance_fields = c("preferred_name", "smiles", "purchasability"),
	Tc_matrix = TRUE,
	verbose
){
	if (!(smiles %in% substance_fields)){
		substance_fields = c(substance_fields, "smiles")
	}

	substances <- target_to_substances(c(ref_target, query_target), fields, verbose)
	ref_substances <- substances %>% dplyr::filter(target == ref_target)
	query_substances <- substances %>% dplyr::filter(target == query_target)

	score <-
		sea_score(
			ref_set = ref_target,
			ref_smi = data.frame(ref_substances %>% dplyr::select(compound=zinc_id, smiles)),
			query_set = query_target,
			query_smi = data.frame(query_substances %>% dplyr::select(compound=zinc_id, smiles)),
			run_tag=paste0(ref_target, "_", query_target),
			verbose=verbose)

	if(tc_matrix){
		tcs <- tc_matrix(
			score$ref_substances$smiles,
			score$query_substances$smiles)
	} else {
		tcs <- NULL
	}
	list(
		ref_target = ref_target,
		query_target = query_target,
		ref_substances = ref_substances,
		query_substances = query_substances,
		tcs,
		score = score)
}

#' @export
sea_illustrate_smiles <- function(
	smiles,
	verbose = False,
	url_path="/custom/search/single"
){
	url <- paste(sea_base_url(), url_path, sep="/") %>%
		httr::parse_url() %>%
		httr::build_url()

	if(verbose){
		cat("url: ", url, "\n")
	}

	tryCatch({
		r <- httr::GET(url)
	}, error=function(e) {
		cat("ERROR getting url: url='", url, "'\n", sep="")
		print(e)
		return(data.frame())
	})

	if((httr::status_code(r) < 200) | httr::status_code(r) >= 300){
		cat("ERROR: url='", url, "'\n", sep="")
		cat("status_code='", r %>% httr::status_code(), "'\n", sep="")
		return(NULL)
	}

	r %>% content("raw")
}

#' @export
sea_illustrate_zinc_id <- function(
	zinc_ids,
	verbose
){
	zinc_ids %>%
		substance_info(fields=c("smiles"), verbose) %>%
		plyr::adply(1, function(substance){
			illustrate_smiles(substance$smiles[1], verbose)
		})
}
