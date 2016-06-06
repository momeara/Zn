# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:

#' @export
sea15_base_url <- function(){"http://sea15.docking.org"}

#' @export
zinc_base_url <- function(){"http://zinc15.docking.org"}

##TODO get post_data to work...
##this python works:
#import requests
#requests.post("http://zinc15.docking.org/substances.json", data={
#    "output_fields" : "zinc_id smiles"),
#    "zinc_id-in" : "ZINC000042835023 ZINC000003815700"}).json()
##but this R does does not:
#content(POST("http://zinc15.docking.org/substances.csv",body=list(
#    output_fields="zinc_id smiles",
#    `zinc_id-in`="ZINC000042835023 ZINC000003815700")), "text")

#' @export
zinc_REST <- function(
	path,
	post_data=NULL,
	count='all',
	result_batch_size=NULL,
	page=NULL,
	verbose=F
){
	url <- paste(zinc_base_url(), path, sep="/") %>%
		httr::parse_url() %>%
		httr::build_url()


	make_request <- function(url, post_data){
		if(verbose){
			cat("making request:\n")
			cat("  url: ", url, "\n")
			if(!is.null(post_data)){
				cat(
					"  post data:\n",
					paste(
						"    ", names(post_data), " having ", vapply(post_data, length, 1), " elements\n",
						sep=""),
					sep="")
			}
		}

		if(is.null(post_data)){
			tryCatch({
				r <- httr::GET(url)
			}, error=function(e) {
				cat("ERROR getting data from url: url='", url, "'\n", sep="")
				print(e)
				stop()
			})
		} else {
			tryCatch({
				r <- post(url, body=post_data)
			}, error=function(e) {
				cat("ERROR posting to url: url='", url, "'\n", sep="")
				print(e)
				stop()
			})
		}
		if((httr::status_code(r) < 200) | httr::status_code(r) >= 300){
			cat("ERROR: url='", url, "'\n", sep="")
			cat("status_code='", r %>% httr::status_code(), "'\n", sep="")
			stop()
		}

		if( r %>% httr::content("text") == "") {
			return(data.frame())
		}

		tryCatch({
			readr::read_delim(
				r %>% httr::content("text"),
				delim=",") %>%
			return
		}, error = function(e){
			cat("ERROR parsing url='", url, "'\n", sep="")
			print(e)
			return(data.frame())
		})
	}

	if(!is.null(result_batch_size)){
		if(is.null(page)){
			page <- 1
		}
		done <- FALSE
		data <- NULL
		data_page_tmp_fname <- paste0(tempfile(), "_zinc_query")
		if(verbose){
			cat("writing temporary results to '", data_page_tmp_fname, "_<page>.csv\n", sep="")
		}
		while(!done){
			url_page <- url %>% httr::parse_url()
			url_page$query$page <- page
			url_page$query$count <- result_batch_size
			url_page <- url_page %>% httr::build_url()
			data_page <- make_request(url_page, post_data)

			data_page %>% readr::write_csv(paste0(data_page_tmp_fname, "_", page, ".csv"))

			if(is.null(data_page) || nrow(data_page) == 0){
				# there is no new data
				done <- TRUE
			} else {
				cur_count <- ifelse(is.null(data), 0, nrow(data))
				new_count <- nrow(data_page)

				if(count != 'all' && (cur_count + new_count > count)){
					# the new data exceeds the requested number of entries
					data <- data %>% rbind(data_page[1:(count - cur_count),])
					done <- TRUE
				} else {
					# the new data does not exceed the requested number of entries
					data <- data %>% rbind(data_page)
					page <- page + 1
				}
			}
		}
	} else {
		url <- url %>% httr::parse_url()
		url$query$count <- count
		url <- url %>% httr::build_url()
		data <- make_request(url, post_data)
	}
	return(data)
}





#' @export
purchasable_compounds_for_ortholog <- function(orthologs, ...) {
	data_frame(ortholog=orthologs) %>% plyr::adply(1, function(row){
		zinc_REST(
			path=paste0("orthologs/", row$ortholog[1], "/substances.smi?purchasability=for-sale"),
			...)
	})
}

#' @export
substance_info <- function(
	zinc_ids,
	output_fields=c("preferred_name", "smiles", "purchasability"),
	...){
	zinc_REST(
		path="substances.json",
		post_data=list(
			`zinc_id-in`=paste(zinc_ids),
			output_fields=paste(output_fields, collapse=" ")))
}

#' @export
search_for_substances <- function(
	search_terms,
	output_fields=c("preferred_name", "smiles", "purchasability"),
	...){
	dplyr::data_frame(search_term=search_terms) %>% plyr::adply(1, function(row){
		zinc_REST(
			path=paste0(
				"substances/search?output_format=csv&q=", row$search_term[1], "&",
				"output_fields=", paste(output_fields, collapse=" ")),
			...)
	})
}

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

#' @export
gene_to_substances <- function(
	gene_names,
	output_fields=c("preferred_name", "smiles", "purchasability"),
	...
){
	dplyr::data_frame(gene_name=gene_names) %>% plyr::adply(1, function(row){
		zinc_REST(
			path=paste0(
				"genes/", row$gene_name[1], "/substances.csv?",
				"output_fields=", paste(output_fields, collapse=" ")),
			...)
	})
}

#' @export
ortholog_to_substances <- function(
	uniprot_entries,
	output_fields=c("preferred_name", "smiles", "purchasability"),
	...
){
	dplyr::data_frame(uniprot_entry=uniprot_entries) %>% plyr::adply(1, function(row){
		zinc_REST(
			path=paste0(
				"orthologs/", row$uniprot_entry[1], "/substances.csv?",
				"output_fields=", paste(output_fields, collapse=" ")),
			...)
	})
}

#' @export
substance_to_genes <- function(
	zinc_ids,
	output_fields=c("name", "description", "major_class_name", "sub_class_name"),
	...
){
	dplyr::data_frame(zinc_id=zinc_ids) %>% plyr::adply(1, function(row){
		zinc_REST(
			path=paste0(
				"substances/", row$zinc_id[1], "/genes.csv?",
				"output_fields=", paste(output_fields, collapse=" ")),
			...)
	})
}

#' @export
substance_to_orthologs <- function(
	zinc_ids,
	output_fields=c("preferred_name", "smiles", "purchasability"),
	...
){
	dplyr::data_frame(zinc_id=zinc_ids) %>% plyr::adply(1, function(row){
		zinc_REST(
			path=paste0(
				"substances/", row$zinc_id[1], "/orthologs.csv?",
				"output_fields=", paste(output_fields, collapse=" ")),
			...)
	})
}

substance_to_protomers <- function(
	zinc_ids,
	output_fields=c(
		"protomers.zinc_id",
		"protomers.prot_id",
		"protomers.net_charge",
		"protomers.desolv_apol",
		"protomers.desolv_pol",
		"protomers.ph_mod_fk",
		"protomers.true_logp",
		"protomers.true_mwt"),
	...){
	dplyr::data_frame(zinc_id=zinc_ids) %>% plyr::adply(1, function(row){
		zinc_REST(
			path=paste0(
				"substances/", row$zinc_id[1], "/protomers.csv?",
				"output_fields=", paste(output_fields, collapse=" ")),
			...)
	})
}


#' @export
catalog_info <- function(
	catalog_short_name,
	...
){
	zinc_REST(
		path = paste0("catalogs/", catalog_short_name, ".csv"),
		...)
}

#' @export
catalog_items <- function(
	catalog_short_name,
	output_fields=c(
		"zinc_id",
		"supplier_code",
		"substance.preferred_name",
		"substance.smiles",
		"substance.purchasable",
		"substance.purchasability",
		"substance.gene_names",
		"substance.rb",
		"substance.reactive",
		"substance.features"),
	...){
	zinc_REST(
		path=paste0(
			"catalogs/", catalog_short_name, "/items.csv?",
			"output_fields=", paste(output_fields, collapse=" ")),
		...)
}

#' @export
catalog_protomers <- function(
	catalog_short_name,
	output_fields=c(
		"zinc_id",
		"prot_id",
		"smiles",
		"net_charge",
		"desolv_apol",
		"desolv_pol",
		"ph_mod_fk",
		"true_logp",
		"true_mwt",
		"hba",
		"hbd",
		"num_aliphatic_rings",
		"num_aromatic_rings",
		"num_heavy_atoms",
		"num_rotatable_bonds",
		"chiral_centers",
		"reactive",
		"reactivity",
		"tpsa",
		"tranche_name",
		"tranche_prefix"
	),
	...){
	zinc_REST(
		path=paste0(
			"catalogs/", catalog_short_name, "/protomers.csv?",
			"output_fields=", paste(output_fields, collapse=" ")),
		...)
}


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
	input_base <- paste(tempfile(), run_tag, sep="_")
	url <- paste(sea_base_url, url_path, sep="/") %>%
		httr::parse_url() %>%
		httr::build_url()

	if(verbose){
		cat("url: ", url, "\n")
	}

	if(is.character(ref_smi)){
		if(!stringr::str_detect(ref_smi, ".smi$")){
			cat("WARNING: referene .smi file '", ref_smi, "', does not end in '.smi'\n")
		}
		ref_smi_fname <- ref_smi
	} else {
		ref_smi_fname <- paste0(input_base, ".ref.smi")
		cat("writing ref_smi to -> '", ref_smi_fname, "' ... ", sep="")
		write.table(
			ref_smi %>% dplyr::select(smiles, compound),
			ref_smi_fname, quote=F, sep=" ", row.names=F, col.names=F)
		cat("DONE\n")
	}

	if(is.character(query_smi)){
		if(!stringr::str_detect(query_smi, ".smi$")){
			cat("WARNING: query .smi file '", query_smi, "', does not end in '.smi'\n")
		}
		query_smi_fname <- query_smi
	} else {
		query_smi_fname <- paste0(input_base, ".query.smi")
		cat("writing query_smi to -> '", query_smi_fname, "' ... ", sep="")
		write.table(
			query_smi %>% dplyr::select(smiles, compound),
			query_smi_fname, quote=F, sep=" ", row.names=F, col.names=F)
		cat("DONE\n")
	}

	ref_smi_f <- httr::upload_file(ref_smi_fname, type="application/smil")
	query_smi_f <- httr::upload_file(query_smi_fname, type="application/smil")

	tryCatch({
		request <- httr::POST(
			url,
			httr::add_headers(Accept="application/json"),
			body=list(
				target_file = ref_smi_f,
				query_file = query_smi_f,
				as_set="yes",
				wait="yes")) %>%
			httr::content()
	}, error = function(e){
		cat("ERROR parsing url='", url, "'\n", sep="")
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
	verbose = F,
	url_path="tools/tc_matrix"
){
	input_base <- paste(tempfile(), run_tag, sep="_")
	url <- paste(sea_base_url, url_path, sep="/") %>%
		httr::parse_url() %>%
		httr::build_url()

	if(verbose){
		cat("url: ", url, "\n")
	}

	if(is.character(ref_smi)){
		if(!stringr::str_detect(ref_smi, ".smi$")){
			cat("WARNING: referene .smi file '", ref_smi, "', does not end in '.smi'\n")
		}
		ref_smi_fname <- ref_smi
	} else {
		ref_smi_fname <- paste0(input_base, ".ref.smi")
		cat("writing ref_smi to -> '", ref_smi_fname, "' ... ", sep="")
		write.table(
			ref_smi %>% dplyr::select(smiles, compound),
			ref_smi_fname, quote=F, sep=" ", row.names=F, col.names=F)
		cat("DONE\n")
	}

	if(is.character(query_smi)){
		if(!stringr::str_detect(query_smi, ".smi$")){
			cat("WARNING: query .smi file '", query_smi, "', does not end in '.smi'\n")
		}
		query_smi_fname <- query_smi
	} else {
		query_smi_fname <- paste0(input_base, ".query.smi")
		cat("writing query_smi to -> '", query_smi_fname, "' ... ", sep="")
		write.table(
			query_smi %>% dplyr::select(smiles, compound),
			query_smi_fname, quote=F, sep=" ", row.names=F, col.names=F)
		cat("DONE\n")
	}

	ref_smi_f <- httr::upload_file(ref_smi_fname, type="application/smil")
	query_smi_f <- httr::upload_file(query_smi_fname, type="application/smil")

	tryCatch({
		request <- httr::POST(
			url,
			httr::add_headers(Accept="application/json"),
			body=list(
				target_file = ref_smi_f,
				query_file = query_smi_f,
				wait="yes")) %>%
			httr::content()
	}, error = function(e){
		cat("ERROR parsing url='", url, "'\n", sep="")
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
	sea_base_url="http://sea15.docking.org",
	url_path="/custom/search/single"
){
	url <- paste(sea_base_url, url_path, sep="/") %>%
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

