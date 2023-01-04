library(stringr)
library(taxize)
library(rlang)

#name_input that will be searched in this script (it substitutes text input in the Apply form in TeloBase)
name_input<- "Arabidopsis thaliana"

#save telomere_sequences.csv as reference that is addressed later when searching if genus or family is not already in the database (it can be downloaded from the TeloBase or GitHub repository)
reference<- read.csv2("telomere_sequences.csv", stringsAsFactors = FALSE)

#remove cf., af., aff. from the name that will be searched
searched_item<- str_replace_all(name_input, c("cf. " = "", "af. " = "", "aff. " = ""))

#cut more from the searched item than first 2 words (usually should be giving only the genus and species for better search in the db)
if (sapply(strsplit(searched_item, " "), length) > 1) {
  searched_item<- word(searched_item,1,2, sep=" ")
}

#create empty table with name rank and id that will be taken up instead of a query that was empty after 3 attempts by classification() (creates problem with atomic vector warning)
taxa_NA<- data.frame(NA,NA,NA)
colnames(taxa_NA)<- c("name", "rank", "id")

#search query in GBIF, NCBI and TOL db using classification() function from taxize package
taxa_gbif<- classification(as.character(searched_item), db = "gbif", rows = 1)
taxa_ncbi<- classification(as.character(searched_item), db = "ncbi", rows = 1)
taxa_tol<- classification(as.character(searched_item), db = "tol", rows = 1)


#filter classification results to only kingdom, phylum, class, order, family, genus and species taxonomy ranking and if empty, use taxa_NA
taxa_gbif_filtered<- if (is.atomic(taxa_gbif[[1]]) == FALSE) {
  taxa_gbif[[1]][taxa_gbif[[1]]$rank %in% c("kingdom", "phylum", "class", "order", "family", "genus", "species"),]
} else {
  taxa_NA
}

taxa_ncbi_filtered<- if (is.atomic(taxa_ncbi[[1]]) == FALSE) {
  taxa_ncbi[[1]][taxa_ncbi[[1]]$rank %in% c("family", "genus", "species"),]
} else {
  taxa_NA
}

taxa_tol_filtered<- if (is.atomic(taxa_tol[[1]]) == FALSE) {
  taxa_tol[[1]][taxa_tol[[1]]$rank %in% c("family", "genus", "species"),]
} else {
  taxa_NA
} 

#fill taxa by GBIF, NCBI or TOL in this order, otherwise NA

family<- if (!is_empty(taxa_gbif_filtered$name[taxa_gbif_filtered$rank == "family"]) && !is.na(taxa_gbif_filtered$name[taxa_gbif_filtered$rank == "family"])) {
  taxa_gbif_filtered$name[taxa_gbif_filtered$rank == "family"]
} else {
  if (!is_empty(taxa_ncbi_filtered$name[taxa_ncbi_filtered$rank == "family"]) && !is.na(taxa_ncbi_filtered$name[taxa_ncbi_filtered$rank == "family"])) {
    taxa_ncbi_filtered$name[taxa_ncbi_filtered$rank == "family"]
  } else {
    if (!is_empty(taxa_tol_filtered$name[taxa_tol_filtered$rank == "family"]) && !is.na(taxa_tol_filtered$name[taxa_tol_filtered$rank == "family"])) {
      taxa_tol_filtered$name[taxa_tol_filtered$rank == "family"]
    } else {
      NA
    }
  }
}

genus<- if (!is_empty(taxa_gbif_filtered$name[taxa_gbif_filtered$rank == "genus"]) && !is.na(taxa_gbif_filtered$name[taxa_gbif_filtered$rank == "genus"])) {
  taxa_gbif_filtered$name[taxa_gbif_filtered$rank == "genus"]
} else {
  if (!is_empty(taxa_ncbi_filtered$name[taxa_ncbi_filtered$rank == "genus"]) && !is.na(taxa_ncbi_filtered$name[taxa_ncbi_filtered$rank == "genus"])) {
    taxa_ncbi_filtered$name[taxa_ncbi_filtered$rank == "genus"]
  } else {
    if (!is_empty(taxa_tol_filtered$name[taxa_tol_filtered$rank == "genus"]) && !is.na(taxa_tol_filtered$name[taxa_tol_filtered$rank == "genus"])) {
      taxa_tol_filtered$name[taxa_tol_filtered$rank == "genus"]
    } else {
      NA
    }
  }
}

species<- if (!is_empty(taxa_gbif_filtered$name[taxa_gbif_filtered$rank == "species"]) && !is.na(taxa_gbif_filtered$name[taxa_gbif_filtered$rank == "species"])) {
  taxa_gbif_filtered$name[taxa_gbif_filtered$rank == "species"]
} else {
  if (!is_empty(taxa_ncbi_filtered$name[taxa_ncbi_filtered$rank == "species"]) && !is.na(taxa_ncbi_filtered$name[taxa_ncbi_filtered$rank == "species"])) {
    taxa_ncbi_filtered$name[taxa_ncbi_filtered$rank == "species"]
  } else {
    if (!is_empty(taxa_tol_filtered$name[taxa_tol_filtered$rank == "species"]) && !is.na(taxa_tol_filtered$name[taxa_tol_filtered$rank == "species"])) {
      taxa_tol_filtered$name[taxa_tol_filtered$rank == "species"]
    } else {
      NA
    }
  }
}

#unification of genus and species if they do not conform to each other (due to the use of different db, etc.)
if (!is.na(genus) && !is.na(species) & genus != word(species, 1)) {
  genus<- word(species, 1)
}

#check if genus is already present in the reference to take the higher order of taxa from the reference
if (!is.na(genus) & genus %in% reference$genus) {
  #save first encounter of the reference to the separate dataframe
  reference_to_higher_taxa<- reference[reference$genus == genus & !is.na(reference$genus),][1,]
  
  #save higher taxa information from the reference
  family<- as.character(reference_to_higher_taxa$family)
  order<- as.character(reference_to_higher_taxa$order)
  class<- as.character(reference_to_higher_taxa$class)
  phylum<- as.character(reference_to_higher_taxa$phylum)
  kingdom<- as.character(reference_to_higher_taxa$kingdom)
  
} else {
  
  #check if family is already present in the reference to take the higher order of taxa from the reference
  if (!is.na(family) & family %in% reference$family) {
    
    #save first encounter of the reference to the separate dataframe
    reference_to_higher_taxa<- reference[reference$family == family & !is.na(reference$family),][1,]
    
    #save higher taxa information from the reference
    order<- as.character(reference_to_higher_taxa$order)
    class<- as.character(reference_to_higher_taxa$class)
    phylum<- as.character(reference_to_higher_taxa$phylum)
    kingdom<- as.character(reference_to_higher_taxa$kingdom)
    
  } else {
    if (!is.na(family)) {
      #do taxa search from the family again in GBIF (allows more unified taxonomy as differences are expected in higher order taxa)
      taxa_gbif_second<- classification(as.character(family), db = "gbif", rows = 1)
      
      
      #filter classification results using family to only kingdom, phylum, class, order and if empty, use taxa_NA
      taxa_gbif_filtered_second<- if (is.atomic(taxa_gbif_second[[1]]) == FALSE) {
        taxa_gbif_second[[1]][taxa_gbif_second[[1]]$rank %in% c("kingdom", "phylum", "class", "order"),]
      } else {
        taxa_NA
      }
      
      #fill taxa by second GBIF search results, otherwise NA
      
      kingdom<- if (!is_empty(taxa_gbif_filtered_second$name[taxa_gbif_filtered_second$rank == "kingdom"]) && 
                    !is.na(taxa_gbif_filtered_second$name[taxa_gbif_filtered_second$rank == "kingdom"])) {
        taxa_gbif_filtered_second$name[taxa_gbif_filtered_second$rank == "kingdom"]
      } else {
        NA
      }
      
      phylum<- if (!is_empty(taxa_gbif_filtered_second$name[taxa_gbif_filtered_second$rank == "phylum"]) && 
                   !is.na(taxa_gbif_filtered_second$name[taxa_gbif_filtered_second$rank == "phylum"])) {
        taxa_gbif_filtered_second$name[taxa_gbif_filtered_second$rank == "phylum"]
      } else {
        NA
      }
      
      class<- if (!is_empty(taxa_gbif_filtered_second$name[taxa_gbif_filtered_second$rank == "class"]) && 
                  !is.na(taxa_gbif_filtered_second$name[taxa_gbif_filtered_second$rank == "class"])) {
        taxa_gbif_filtered_second$name[taxa_gbif_filtered_second$rank == "class"]
      } else {
        NA
      }
      
      order<- if (!is_empty(taxa_gbif_filtered_second$name[taxa_gbif_filtered_second$rank == "order"]) && 
                  !is.na(taxa_gbif_filtered_second$name[taxa_gbif_filtered_second$rank == "order"])) {
        taxa_gbif_filtered_second$name[taxa_gbif_filtered_second$rank == "order"]
      } else {
        NA
      }
    } else {
      
      #use first GBIF search to fill the higher taxa (if family had NA) 
      
      kingdom<- if (!is_empty(taxa_gbif_filtered$name[taxa_gbif_filtered$rank == "kingdom"]) && 
                    !is.na(taxa_gbif_filtered$name[taxa_gbif_filtered$rank == "kingdom"])) {
        taxa_gbif_filtered$name[taxa_gbif_filtered$rank == "kingdom"]
      } else {
        NA
      }
      
      phylum<- if (!is_empty(taxa_gbif_filtered$name[taxa_gbif_filtered$rank == "phylum"]) && 
                   !is.na(taxa_gbif_filtered$name[taxa_gbif_filtered$rank == "phylum"])) {
        taxa_gbif_filtered$name[taxa_gbif_filtered$rank == "phylum"]
      } else {
        NA
      }
      
      class<- if (!is_empty(taxa_gbif_filtered$name[taxa_gbif_filtered$rank == "class"]) && 
                  !is.na(taxa_gbif_filtered$name[taxa_gbif_filtered$rank == "class"])) {
        taxa_gbif_filtered$name[taxa_gbif_filtered$rank == "class"]
      } else {
        NA
      }
      
      order<- if (!is_empty(taxa_gbif_filtered$name[taxa_gbif_filtered$rank == "order"]) && 
                  !is.na(taxa_gbif_filtered$name[taxa_gbif_filtered$rank == "order"])) {
        taxa_gbif_filtered$name[taxa_gbif_filtered$rank == "order"]
      } else {
        NA
      }
    }
  }
}
