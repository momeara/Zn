# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:


#' @export
purchasable_compounds_for_ortholog <- function(orthologs, ...) {
  tibble::tibble(ortholog = orthologs) %>%
    plyr::adply(1, function(row) {
      zinc_REST(
        path = paste0("orthologs/", row$ortholog[1], "/substances.smi"),
        query = list(purchasability = "for-sale"),
        ...)
    })
}



#' @export
gene_to_substances <- function(
  gene_names,
  output_fields=c("preferred_name", "smiles", "purchasability"),
  ...
) {
    tibble::tibble(gene_name = gene_names) %>%
        plyr::adply(1, function(row) {
        zinc_REST(
          path = paste0("genes/", row$gene_name[1], "/substances.csv"),
          query = list(output_fields =  paste(output_fields, collapse = " ")),
          ...)
      })
}

#' @export
ortholog_to_substances <- function(
  uniprot_entries,
  output_fields=c("preferred_name", "smiles", "purchasability"),
  ...
) {
    tibble::tibble(uniprot_entry = uniprot_entries) %>%
      plyr::adply(1, function(row) {
        zinc_REST(
          path = paste0("orthologs/", row$uniprot_entry[1], "/substances.csv"),
          query = list(output_fields = paste(output_fields, collapse = " ")),
          ...)
      })
}

#' @export
substance_to_genes <- function(
  zinc_ids,
  output_fields=c("name", "description", "major_class_name", "sub_class_name"),
  ...
) {
  tibble::tibble(zinc_id = zinc_ids) %>%
    plyr::adply(1, function(row) {
      zinc_REST(
        path = paste0("substances/", row$zinc_id[1], "/genes.csv"),
        query = list(output_fields = paste(output_fields, collapse = " ")),
        ...)
      })
}

#' @export
substance_to_orthologs <- function(
  zinc_ids,
  output_fields=c("preferred_name", "smiles", "purchasability"),
  ...
) {
  tibble::tibble(zinc_id = zinc_ids) %>%
    plyr::adply(1, function(row) {
      zinc_REST(
        path = paste0("substances/", row$zinc_id[1], "/orthologs.csv"),
        query = list(output_fields = paste(output_fields, collapse = " ")),
        ...)
    })
}

#' @export
substance_to_protomers <- function(
  zinc_ids,
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
  ...) {
  tibble::tibble(zinc_id = zinc_ids) %>%
    plyr::adply(1, function(row) {
      zinc_REST(
        path = paste0("substances/", row$zinc_id[1], "/protomers.csv"),
        query = list(output_fields = paste(output_fields, collapse = " ")),
        ...)
    })
}
