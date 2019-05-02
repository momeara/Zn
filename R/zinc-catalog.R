# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:


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
		"tranche_prefix"),
	...){
	zinc_REST(
		path=paste0(
			"catalogs/", catalog_short_name, "/protomers.csv?",
			"output_fields=", paste(output_fields, collapse=" ")),
		...)
}
