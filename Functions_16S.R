# Libraries -----------------------------------------------------------------

library(tidyverse)
library(readxl)
library(dplyr)
library(glue)

# Functions ---------------------------------------------------------------

## Format and Combine Tables -----------------------------------------------


import_counts <- function(.shared) {
  read_tsv(.shared) %>%
    select(-label, starts_with("Otu")) %>%
    rename(sample_id = Group) %>%
    pivot_longer(-sample_id, names_to="otu", values_to = "count")
}

import_taxonomy <- function(.taxonomy) {
  read_tsv(.taxonomy) %>%
    select("OTU", "Taxonomy") %>%
    rename_all(tolower) %>%
    mutate(taxonomy = str_replace_all(taxonomy, "\\(\\d+\\)", ""),
           taxonomy = str_replace(taxonomy, ";$", "")) %>%
    separate(taxonomy,
             into=c("kingdom", "phylum", "class", "order", "family", "genus", "species"),
             sep=";") %>%
    mutate(pretty_otu = str_replace(string=otu,
                                    pattern="tu0*",
                                    replacement = "TU "),
           species = str_replace(string=species,
                                 pattern=".+_unclassified",
                                 replacement=""),
           OTU = if_else(species == "", glue("*{genus}* ({pretty_otu})"),glue("*{genus} {species}* ({pretty_otu})")),
           OTU = str_replace(string=OTU,
                             pattern="\\*(.*)_unclassified\\*",
                             replacement="Unclassified *\\1*"), 
           genus = str_replace(string=genus,
                               pattern="(.*)",
                               replacement="*\\1*"),
           genus = str_replace(string=genus,
                               pattern="\\*(.*)_unclassified\\*",
                               replacement="Unclassified *\\1*")) 
}

Combine_tables <- function(.shared, .taxonomy, metadata.xlsx) {
  otu_counts <- import_counts(.shared)
  taxonomy <- import_taxonomy(.taxonomy)
  full_table <- read_excel(metadata.xlsx) %>%
    inner_join(., otu_counts, by="sample_id") %>%
    inner_join(., taxonomy, by="otu")
}

create_rel_abund <- function(.shared, .taxonomy, metadata.xlsx) {
  Combine_tables(.shared, .taxonomy, metadata.xlsx) %>%
    group_by(sample_id)%>%
    mutate(rel_abund = count / sum(count))
}

### Stacked Bar Chart -------------------------------------------------------

make_mean_rel_abund <- function(my_otu_rel_abund, my_group_list) {
  otu_rel_abund %>%  
    group_by(otu_rel_abund[,my_group_list]) %>%
    group_by(phylum,class,order,family,genus,OTU, .add = TRUE) %>%
    summarize(mean_rel_abund = 100*mean(rel_abund), .groups="drop")
}

taxon_filter <- function(mean_rel_abund, taxon_cutoff=3, my_filter=TRUE,taxon=phylum){
  mean_rel_abund %>%
    filter({{my_filter}}) %>% 
    group_by({{taxon}}) %>%
    summarise(pool = max(mean_rel_abund) < taxon_cutoff, 
              mean = mean(mean_rel_abund),
              .groups = "drop") %>%
    arrange(pool) %>% print(n=10)
}

genus_filter <- function(mean_rel_abund,taxon_pool,genus_cutoff=3,my_filter=TRUE){
  inner_join(mean_rel_abund, taxon_pool, by = "phylum") %>%
    mutate(taxon_phylum = if_else(pool,"Other", phylum)) %>%
    filter({{my_filter}}) %>% 
    group_by(taxon_phylum,family,genus) %>%
    summarise(pool = max(mean_rel_abund) < genus_cutoff, 
              mean = mean(mean_rel_abund),
              .groups = "drop") %>%
    arrange(pool) %>% print(n=20)  
}

apply_filters <- function(mean_rel_abund, taxon_pool2){
  inner_join(mean_rel_abund, taxon_pool2, by = c("OTU")) %>%
    mutate(taxon = if_else(pool, paste("Other",taxon_phylum), genus)) %>%
    mutate(filtered_taxon = if_else(taxon == "Other Other","Other Bacteria", taxon))
}

filtered_summary <- function(final, my_group_list){
  final %>%
    group_by(final[,my_group_list]) %>%
    group_by(phylum, family, filtered_taxon, .add = TRUE) %>%
    summarise(mean_rel_abund = sum(mean_rel_abund)) %>%
    arrange(desc(str_detect(filtered_taxon, "Other Bacteria")),phylum,desc(str_detect(filtered_taxon,"Other")),str_detect(family, "Lachnospiraceae"),mean_rel_abund)
}

process_rel_abund <- function(otu_rel_abund, my_group_list, phylum_cutoff, genus_cutoff){
  mean_rel_abund <- make_mean_rel_abund(otu_rel_abund, my_group_list)
  taxon_pool <- phylum_filter(mean_rel_abund, phylum_cutoff)
  taxon_pool2 <- genus_filter(mean_rel_abund, taxon_pool, genus_cutoff)
  final <- apply_filters(mean_rel_abund, taxon_pool2)
  message("You have ",  taxon_pool %>%filter(pool == FALSE)%>%nrow(), " phyla and ", taxon_pool2 %>%filter(pool == FALSE)%>%nrow(), " genera after filtering.")
  message("You will need ", length(unique(final$filtered_taxon))," different colors.")
  summary <- filtered_summary(final,my_group_list)
  
}

### OTU Heatmap -------------------------------------------------------------

calc_fold_change <- function(mean_rel_abund,taxon_pool,my_group_list,baseline){
  min_val <- min(mean_rel_abund[mean_rel_abund$mean_rel_abund > 0, "mean_rel_abund"])
  combined <- inner_join(mean_rel_abund, taxon_pool, by = "OTU") %>%
    mutate(taxon = if_else(pool,"Other Bacteria",OTU)) %>%
    filter(taxon != "Other Bacteria")
  combined %>%
    group_by(combined[,my_group_list]) %>%
    group_by(OTU, .add = TRUE) %>%
    mutate(log_delta = log2(mean_rel_abund/(case_when(
      mean_rel_abund[{{baseline}}] >= min_val ~ mean_rel_abund[{{baseline}}],
      TRUE ~  min_val)))) %>%
    mutate(log_delta = ifelse(log_delta == -Inf,NA,log_delta)) %>%
    arrange(phylum,str_detect(family,"Lachnospiraceae"),mean_rel_abund)
}

make_heatmap_ready <- function(mean_rel_abund, taxon_cutoff=1, my_filter = TRUE,my_group_list,baseline,taxon=OTU){
  taxon_pool <- taxon_filter(mean_rel_abund, taxon_cutoff, {{my_filter}},{{taxon}})
  min_val <- min(mean_rel_abund[mean_rel_abund$mean_rel_abund > 0, "mean_rel_abund"])
  combined <- calc_fold_change(mean_rel_abund,taxon_pool,my_group_list,{{baseline}})
}
