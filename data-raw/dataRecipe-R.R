## code to prepare `dataRecipe.R` dataset goes here
IUPAC_CODE_MAP_extended <- c(Biostrings::IUPAC_CODE_MAP,
                             structure("-" ,names =  "-"),
                             structure(paste0(Biostrings::IUPAC_CODE_MAP, "-"),
                                       names = names(Biostrings::IUPAC_CODE_MAP)))
use_data(IUPAC_CODE_MAP_extended)
IUPAC_code_2nucs <- c("A", "C", "G", "T", "M", "R", "W", "S", "Y", "K", "N", "-", ".")
use_data(IUPAC_code_2nucs)
IUPAC_code_simpl <- c("A", "C", "G", "T", "N", "-", ".")
use_data(IUPAC_code_simpl)
