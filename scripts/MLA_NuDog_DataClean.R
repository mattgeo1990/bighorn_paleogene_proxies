# 01_clean_NuDog_session22.R
# Purpose: Standardize Type.1 and Type.2 classifications for the
# NuDog clumped-isotope dataset.

# ---- Load packages ----
library(tidyverse)
library(here)


# ---- File paths ----
input_file <- here(
  "data", "raw",
  "Nu_Dog_Clump_Session22_Oct 2025-July2026_Matlab-2.csv"
)

output_file <- here(
  "data", "processed",
  "Nu_Dog_Clump_Session22_Oct 2025-July2026_Matlab-2_MLAcleaned.csv"
)


# ---- Classification lists ----
standard_names <- c(
  "102-GC-AZ01",
  "102-GC-AZ01 HEATED",
  "ETH-1",
  "ETH-1 HEATED",
  "ETH-2",
  "ETH-2 HEATED",
  "ETH-3",
  "ETH-3 HEATED",
  "ETH-4",
  "ETH-4 HEATED",
  "IAEA-C1",
  "IAEA-C1 HEATED"
)

heated_gas_names <- c(
  "ONLINE DEPLETED HEATED GAS",
  "ONLINE ENRICHED HEATED GAS"
)


# ---- Load data ----
NuDog <- read_csv(
  input_file,
  show_col_types = FALSE
)


# ---- Clean and classify data ----
NuDog <- NuDog %>%
  mutate(
    # Remove accidental leading, trailing, or repeated spaces
    Sample_Name = str_squish(Sample_Name),
    
    # Uppercase copy used only for consistent matching
    Sample_Name_match = str_to_upper(Sample_Name),
    
    # Identify Bighorn Basin carbonate samples
    is_BHB_sample = coalesce(
      str_detect(Sample, regex("PK95|PB", ignore_case = TRUE)),
      FALSE
    ),
    
    # Classify the broad analysis type
    Type.1 = case_when(
      Sample_Name_match %in% heated_gas_names ~ "HeatedGas",
      
      str_detect(
        Sample_Name,
        regex("Eq Gas", ignore_case = TRUE)
      ) ~ "EquilibratedGas",
      
      Sample_Name_match %in% standard_names ~ "Standard",
      
      TRUE ~ "Sample"
    ),
    
    # Classify the specific analysis subtype
    Type.2 = case_when(
      # Standards are labeled by standard name
      Type.1 == "Standard" ~ Sample_Name,
      
      # Equilibrated and heated gases are enriched or depleted
      Type.1 %in% c("EquilibratedGas", "HeatedGas") &
        str_detect(
          Sample_Name,
          regex("Enriched", ignore_case = TRUE)
        ) ~ "Enriched",
      
      Type.1 %in% c("EquilibratedGas", "HeatedGas") &
        str_detect(
          Sample_Name,
          regex("Depleted", ignore_case = TRUE)
        ) ~ "Depleted",
      
      # Identify Matthew's BHB carbonate samples
      Type.1 == "Sample" & is_BHB_sample ~ "Matthew's BHB Carbs",
      
      # Remove incorrect BHB labels from unrelated samples
      Type.1 == "Sample" &
        Type.2 == "Matthew's BHB Carbs" &
        !is_BHB_sample ~ NA_character_,
      
      # Remove enriched/depleted labels from ordinary samples
      Type.1 == "Sample" &
        Type.2 %in% c("Enriched", "Depleted") ~ NA_character_,
      
      # Preserve all remaining existing values
      TRUE ~ Type.2
    )
  ) %>%
  select(
    -Sample_Name_match,
    -is_BHB_sample
  )


# ---- Check classifications ----
NuDog %>%
  count(Type.1, Type.2, sort = TRUE) %>%
  print(n = Inf)

NuDog %>%
  filter(str_detect(
    Sample_Name,
    regex("Heated Gas", ignore_case = TRUE)
  )) %>%
  count(Sample_Name, Type.1, Type.2)


# ---- Export cleaned dataset ----
write_csv(
  NuDog,
  output_file,
  na = ""
)