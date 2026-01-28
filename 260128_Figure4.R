library(tidyverse)
library(patchwork)

# ------------------------------------------------------------------
# 1. Path 
# ------------------------------------------------------------------
path_theoretical <- "./InSilico_TheoreticalPeptides"
path_results     <- "./data/1.8.1"

# ------------------------------------------------------------------
# 2. Functions
# ------------------------------------------------------------------

# 2.1 Load Theoretical Peptides
theo_list <- lapply(1:12, function(i) {
  read.csv(file.path(path_theoretical, paste0("TheoreticalPep_Mix-", i, ".txt")), sep = "\t")
})

# 2.2 Read diann
read_data <- function(path) {
  lapply(1:12, function(i) {
    read.csv(file.path(path, sprintf("Mix-%d.pr_matrix.tsv", i)), sep = "\t") %>%
      filter(Proteotypic == 1) %>%
      select(Protein.Group, Stripped.Sequence, Precursor.Id) %>%
      rename(Sequence = Stripped.Sequence) %>%
      mutate(Mix = paste0("Mix", i))
  })
}

# 2.3 Process Matrix 
process_matrix <- function(data_list) {
  bind_rows(data_list) %>%
    pivot_wider(
      names_from  = Mix,
      values_from = Mix,
      values_fill = NA_character_
    ) %>%
    group_by(Precursor.Id) %>%
    summarise(
      across(starts_with("Mix"), ~ {
        cnt <- sum(!is.na(.x))
        if(cnt == 0) NA_integer_ else cnt
      }),
      Sequence = first(Sequence),
      Protein.Group = first(Protein.Group),
      .groups = "drop"
    ) %>%
    mutate(
      Count = rowSums(!is.na(select(., starts_with("Mix"))))
    )
}

# ------------------------------------------------------------------
# 3. Without Carryover Lists
# ------------------------------------------------------------------

# Load and process
pr_list_raw <- read_data(path_results)
pr_matrix   <- process_matrix(pr_list_raw)

# Create base lists
list_wo_base <- lapply(1:12, function(i) {
  curr <- paste0("Mix", i)
  prev <- paste0("Mix", i - 1)
  
  if (i == 1) {
    target <- pr_matrix %>% filter(!is.na(!!sym(curr)))
  } else {
    target <- pr_matrix %>% filter(!is.na(!!sym(curr)) & is.na(!!sym(prev)))
  }
  
  target %>% arrange(desc(Count))
})

# ------------------------------------------------------------------
# 4. Standardize Protein Groups
# ------------------------------------------------------------------

list_wo <- list()

for (i in 1:12) {
  # 1. Reread raw files to get a unique Protein Group per Sequence
  raw_map <- read.csv(file.path(path_results, sprintf("Mix-%d.pr_matrix.tsv", i)), sep = "\t") %>%
    select(Protein.Group, Stripped.Sequence) %>%
    rename(Sequence = Stripped.Sequence) %>%
    distinct(Sequence, .keep_all = TRUE) 
  
  # 2. Assign standardized Protein Groups to the list
  list_wo[[i]] <- list_wo_base[[i]] %>%
    select(-Protein.Group) %>% 
    left_join(raw_map, by = "Sequence")
}

# ------------------------------------------------------------------
# 5. Figure 4A: Precursor FDR
# ------------------------------------------------------------------

# Function to calculate FDR per Mix at a fixed threshold (Count <= 3)
calculate_precursor_fdr_mixwise <- function(data_list, ref_list, threshold = 3) {
  stats <- map2_dfr(data_list, ref_list, function(df, theo) {
    
    # Filter by detection frequency
    target <- df %>% filter(Count <= threshold)
    if(nrow(target) == 0) return(tibble(All=0, False=0))
    
    # Check against theoretical
    n_false <- sum(!target$Sequence %in% theo$Sequence)
    tibble(All = nrow(target), False = n_false)
  })
  
  stats %>%
    mutate(
      FDR = ifelse(All == 0, 0, (False/All)*100),
      Mix = factor(paste0("Mix-", 1:12), levels = paste0("Mix-", 1:12))
    )
}

# Calculate stats
fdr_precursor_plot <- calculate_precursor_fdr_mixwise(list_wo, theo_list, threshold = 3)

# Plot Figure 4A
fig4a <- ggplot(fdr_precursor_plot, aes(x = Mix, y = FDR)) +
  geom_col(fill = "#4C7AF1", width = 0.7) +
  labs(title = "Precursor FDR (Count <= 3)", x = NULL, y = "Precursor FDR (%)") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1.5))

# ------------------------------------------------------------------
# 6. Figure 4B: Protein FDR
# ------------------------------------------------------------------

pg_list_clean <- list()

for (i in 1:12) {
  df <- list_wo[[i]]
  
  # Classify Protein Groups: Simple vs. Complex
  type1 <- df %>% filter(!grepl("/", Protein.Group))
  
  type2 <- df %>% 
    filter(grepl("/", Protein.Group)) %>%
    filter(grepl(paste0("Mix", i, "_"), Protein.Group)) %>%
    separate_rows(Protein.Group, sep = "/") %>%
    filter(grepl(paste0("Mix", i, "_"), Protein.Group))
  
  type3 <- df %>% 
    filter(grepl("/", Protein.Group)) %>%
    filter(!grepl(paste0("Mix", i, "_"), Protein.Group))
  
  # Clean strings and remove duplicates
  pg_list_clean[[i]] <- bind_rows(type1, type2, type3) %>%
    mutate(Protein.Group = str_replace_all(Protein.Group, pattern = "Mix[0-9]+_", replacement = "")) %>%
    distinct(Protein.Group, Sequence, .keep_all = TRUE) %>%
    mutate(MatchKey = paste(Protein.Group, Sequence, sep = ",")) 
}

# Validation against Theoretical
pg_list_checked <- lapply(1:12, function(i) {
  theo_keys <- paste(theo_list[[i]]$Protein.Group, theo_list[[i]]$Sequence, sep = ",")
  
  pg_list_clean[[i]] %>%
    mutate(Exist = ifelse(MatchKey %in% theo_keys, 1, 0))
})

# Calculate Stats (Count <= 3)
fdr_protein_stats <- map_dfr(pg_list_checked, function(df) {
  target <- df %>% filter(Count <= 3)
  
   #A protein is False if it contains any false peptide
  prot_status <- target %>%
    group_by(Protein.Group) %>%
    arrange(Exist) %>% 
    distinct(Protein.Group, .keep_all = TRUE)
  
  tibble(
    All = nrow(prot_status),
    False = sum(prot_status$Exist == 0)
  )
})

# Prepare Data for Plotting
fdr_protein_plot <- fdr_protein_stats %>%
  mutate(
    FDR = (False / All) * 100,
    Mix = factor(paste0("Mix-", 1:12), levels = paste0("Mix-", 1:12))
  )

fig4b <- ggplot(fdr_protein_plot, aes(x = Mix, y = FDR)) +
  geom_col(fill = "#4C7AF1", width = 0.7) +
  labs(title = "Protein FDR (Count <= 3)", x = NULL, y = "Protein FDR (%)") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 5))

# ------------------------------------------------------------------

# ------------------------------------------------------------------
# 7. Save Final Figure
# ------------------------------------------------------------------
final_plot <- fig4a / fig4b
final_plot
ggsave("Figure4_Combined_Final.pdf", plot = final_plot, width = 10, height = 8)
