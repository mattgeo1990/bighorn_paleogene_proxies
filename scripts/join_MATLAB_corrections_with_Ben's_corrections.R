library(tidyverse)
library(here)

# ---- Load data ----

ben <- read_csv(
  here(
    "data",
    "junkyard data",
    "IPL_D47_Bens_corrections.csv"
  )
)

will <- read_csv(
  here(
    "data",
    "junkyard data",
    "IPL_D47_Wills_corrections.csv"
  )
)

summary <- read_csv(
  here(
    "data",
    "raw",
    "IPL_D47_BHB_Pg_Summary_June2026.csv"
  )
)

# ---- Normalize sample names ----

normalize_sample_name <- function(x) {
  x %>%
    str_replace_all("'", "") %>%
    str_squish() %>%
    str_to_upper() %>%
    str_replace("^PKQS-", "PK95-") %>%
    str_replace("^PK95-(?=[0-9])", "PK95-SC-") %>%
    str_replace("-80-?N-?([12])", "-80-N\\1") %>%
    str_remove("-SPAR$")
}

# ---- Build accepted horizon lookup ----

horizon_lookup <- summary %>%
  filter(!is.na(MLA_horizon_id)) %>%
  distinct(MLA_horizon_id) %>%
  mutate(
    normalized_id = normalize_sample_name(MLA_horizon_id)
  )

# ---- Guess horizon ID from sample name ----

guess_horizon_id <- function(sample_name, lookup) {
  x <- normalize_sample_name(sample_name)
  
  exact <- lookup %>%
    filter(normalized_id == x)
  
  if (nrow(exact) == 1) {
    return(exact$MLA_horizon_id)
  }
  
  replicate_match <- lookup %>%
    filter(
      map_lgl(
        normalized_id,
        function(id) {
          prefix <- paste0(id, "-")
          
          startsWith(x, prefix) &&
            str_detect(
              str_remove(x, fixed(prefix)),
              "^[0-9]{1,2}$"
            )
        }
      )
    )
  
  if (nrow(replicate_match) == 1) {
    return(replicate_match$MLA_horizon_id)
  }
  
  NA_character_
}

# ---- Prepare Ben's preferred T47 values ----

ben_keys <- ben %>%
  mutate(
    Sample_Name_normalized = normalize_sample_name(Sample_Name),
    
    MLA_horizon_id_ben = map_chr(
      Sample_Name,
      guess_horizon_id,
      lookup = horizon_lookup
    )
  ) %>%
  transmute(
    IPLnum = as.integer(IPLnum),
    Ben_Sample_Name = Sample_Name,
    MLA_horizon_id_ben,
    T47_ben = T47_preferred
  )

# ---- Prepare Will's corrected mineral isotope values ----

will_mineral <- will %>%
  rename(
    Sample_Name = `'SampleID'`,
    d13Ccarb_VPDB = `'d13CPDB'`,
    d18Ocarb_VPDB = `'d18OPDBmin'`
  ) %>%
  mutate(
    Sample_Name = str_replace_all(Sample_Name, "'", ""),
    Sample_Name_normalized = normalize_sample_name(Sample_Name),
    
    MLA_horizon_id_will = map_chr(
      Sample_Name,
      guess_horizon_id,
      lookup = horizon_lookup
    )
  ) %>%
  transmute(
    IPLnum = as.integer(IPLnum),
    Sample_Name,
    MLA_horizon_id_will,
    T47_will = T47_preferred,
    
    source_d18O_VPDB = d18Ocarb_VPDB,
    source_d18O_VSMOW =
      1.03092 * d18Ocarb_VPDB + 30.92,
    source_d13C_VPDB = d13Ccarb_VPDB
  )

# ---- Combine Ben's T47 with Will's mineral isotopes ----

correction_lookup <- will_mineral %>%
  left_join(
    ben_keys,
    by = "IPLnum",
    relationship = "one-to-one"
  ) %>%
  transmute(
    IPLnum,
    Sample_Name,
    
    MLA_horizon_id = coalesce(
      MLA_horizon_id_ben,
      MLA_horizon_id_will
    ),
    
    T47_match = round(
      coalesce(T47_ben, T47_will),
      1
    ),
    
    source_d18O_VSMOW,
    source_d18O_VPDB,
    source_d13C_VPDB,
    
    correction_source = if_else(
      !is.na(T47_ben),
      "Ben T47; Will mineral isotopes",
      "Will T47; Will mineral isotopes"
    )
  )

# ---- Review guessed horizon IDs ----

horizon_id_review <- correction_lookup %>%
  select(
    IPLnum,
    Sample_Name,
    MLA_horizon_id,
    T47_match,
    correction_source
  ) %>%
  mutate(
    horizon_match_status = if_else(
      is.na(MLA_horizon_id),
      "CHECK",
      "matched"
    )
  ) %>%
  arrange(
    desc(horizon_match_status == "CHECK"),
    MLA_horizon_id,
    T47_match
  )

print(horizon_id_review, n = Inf)

horizon_id_review %>%
  count(horizon_match_status) %>%
  print()

# ---- Check for duplicate matching keys ----

duplicate_keys <- correction_lookup %>%
  filter(
    !is.na(MLA_horizon_id),
    !is.na(T47_match)
  ) %>%
  count(
    MLA_horizon_id,
    T47_match
  ) %>%
  filter(n > 1)

print(duplicate_keys, n = Inf)

if (nrow(duplicate_keys) > 0) {
  stop(
    "Correction lookup contains duplicate horizon-temperature keys."
  )
}

# ---- Join IPL numbers and mineral isotopes to summary ----

summary <- summary %>%
  select(
    -any_of(
      c(
        "IPLnum",
        "Sample_Name",
        "correction_source"
      )
    )
  ) %>%
  mutate(
    T47_match = round(T47_preferred, 1)
  ) %>%
  left_join(
    correction_lookup,
    by = c(
      "MLA_horizon_id",
      "T47_match"
    ),
    relationship = "many-to-one"
  ) %>%
  mutate(
    IPL_NuDog_d18Ocarb_VSMOW = coalesce(
      source_d18O_VSMOW,
      IPL_NuDog_d18Ocarb_VSMOW
    ),
    
    IPL_NuDog_d18Ocarb_VPDB = coalesce(
      source_d18O_VPDB,
      IPL_NuDog_d18Ocarb_VPDB
    ),
    
    IPL_NuDog_d13Ccarb_VPDB = coalesce(
      source_d13C_VPDB,
      IPL_NuDog_d13Ccarb_VPDB
    )
  )

# ---- Review unmatched summary rows ----

unmatched_summary <- summary %>%
  filter(is.na(IPLnum)) %>%
  select(
    MLA_sample_id,
    MLA_horizon_id,
    T47_preferred,
    Session
  )

print(unmatched_summary, n = Inf)

if (nrow(unmatched_summary) > 0) {
  stop(
    "Some summary rows did not match a correction record."
  )
}

summary %>%
  count(correction_source) %>%
  print()

# ---- Clean temporary columns ----

IPL_D47_BHB_Pg_Summary <- summary %>%
  select(
    -T47_match,
    -source_d18O_VSMOW,
    -source_d18O_VPDB,
    -source_d13C_VPDB
  )

# ---- Overwrite original summary ----

write_csv(
  IPL_D47_BHB_Pg_Summary,
  here(
    "data",
    "raw",
    "IPL_D47_BHB_Pg_Summary_June2026.csv"
  ),
  na = ""
)


plot(IPL_D47_BHB_Pg_Summary$T47_preferred, IPL_D47_BHB_Pg_Summary$IPL_NuDog_d18Ocarb_VSMOW)

