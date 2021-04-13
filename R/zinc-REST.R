# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:



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
	post_data = NULL,
	count = "all",
	result_batch_size = NULL,
	temp_file_base = tempfile(),
	page = NULL,
	retry_attempts = 5,
	verbose = FALSE,
	...
) {
	url <- paste(zinc_base_url(), path, sep = "/") %>%
		httr::parse_url() %>%
		httr::build_url()
	url_args <- list(...)

	make_request <- function(url, post_data) {
		if (verbose) {
			cat("making request:\n")
			cat("  url: ", url, "\n")
			if (!is.null(post_data)) {
				post_args <- ""
				for (i in 1:length(post_data)) {
					elements <- stringr::str_split(post_data[i], "\\s+")[[1]]
					if (length(elements) <= 5) {
						value <- post_data[i]
					} else{
						value <- paste0(paste0(elements[1:5], collapse = " "), " ... <", length(elements) - 5, " more>")
					}
					post_args <- paste(post_args, "   ", names(post_data)[i], "=", value, "\n")
				}
				cat("  post data:\n", post_args, sep = "")
			}
		}

		succeeded <- FALSE
		for (i in 1:(retry_attempts + 1)) {
		  if (i > 1) {
				cat("retrying attempt ", i, " ... \n", sep = "")
			}
			r <- NULL
			if (is.null(post_data)) {
				r <- tryCatch(
					do.call(httr::GET, args = c(url = url, url_args)),
				  error = function(e) {
					  cat("ERROR getting data from url: url='", url, "'\n", sep = "")
					  print(e)
					  NULL
				  })
			} else {
  			r <- tryCatch(
					do.call(httr::POST, args = c(list(url = url, body = post_data), url_args)),
					error = function(e) {
						cat("ERROR posting to url: url='", url, "'\n", sep = "")
						print(e)
						NULL
					})
			}

			if (is.null(r) || (httr::status_code(r) < 200) || (httr::status_code(r) >= 300)) {
				cat("ERROR: url='", url, "'\n", sep = "")
				if (!is.null(post_data)) {
					cat("ERROR: post_data:\nERROR:  ",
						paste(
							paste0(names(post_data), "=", post_data),
							collapse = "\nERROR:  "),
						"\n", sep = "")
				}
				if (!is.null(r)) {
					cat("status_code='", r %>% httr::status_code(), "'\n", sep="")
				}
			} else{
				succeeded <- TRUE
				break
		  }
		}
		if(!succeeded){
			cat("Failed after", retry_attempts, "attempts\n")
			stop()
		}

		contents <- r %>% httr::content("text", encoding = "UTF-8")
		if( contents == "") {
			return(data.frame())
		}

		tryCatch({
			return(readr::read_delim(contents, delim = ","))
		}, error = function(e) {
			cat("ERROR parsing url='", url, "'\n", sep = "")
			print(e)
			return(data.frame())
		})
	}

	if (!is.null(result_batch_size)) {
		if (is.null(page)) {
			page <- 1
		}
		done <- FALSE
		data <- NULL
		data_page_tmp_fname <- paste0(temp_file_base, "_zinc_query")
		if (verbose) {
			cat("writing temporary results to '", data_page_tmp_fname, "_<page>.csv\n", sep="")
		}
		while (!done) {
			url_page <- url %>% httr::parse_url()
			url_page$query$page <- page
			url_page$query$count <- result_batch_size
			url_page <- url_page %>% httr::build_url()
			data_page <- make_request(url_page, post_data)

			data_page %>% readr::write_csv(paste0(data_page_tmp_fname, "_", page, ".csv"))

			if(is.null(data_page) || nrow(data_page) == 0) {
				# there is no new data
				done <- TRUE
			} else {
				cur_count <- ifelse(is.null(data), 0, nrow(data))
				new_count <- nrow(data_page)

				if(count != "all" && (cur_count + new_count >= count)) {
					# the new data exceeds the requested number of entries
					data <- data %>% rbind(data_page[1:(count - cur_count), ])
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
		if (!is.null(count)) {
			url$query$count <- count
		}
		url <- url %>% httr::build_url()
		data <- make_request(url, post_data)
	}
	return(data)
}
