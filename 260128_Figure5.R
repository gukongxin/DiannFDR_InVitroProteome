library(tidyverse)
library(patchwork)

# ------------------------------------------------------------------
# 1. Path 
# ------------------------------------------------------------------
path_theoretical <- "./InSilico_TheoreticalPeptides"
path_results     <- "./data/1.8.1"

# ------------------------------------------------------------------
# 2. Replicating Figure 4 
# ------------------------------------------------------------------

# 2.1 Load Theoretical Peptides
theo_list <- lapply(1:12, function(i) {
  read.csv(file.path(path_theoretical, paste0("TheoreticalPep_Mix-", i, ".txt")), sep = "\t")
})

# 2.2 Read & Summarize Precursors
pr_list_raw <- lapply(1:12, function(i) {
  read.csv(file.path(path_results, sprintf("Mix-%d.pr_matrix.tsv", i)), sep = "\t") %>%
    filter(Proteotypic == 1) %>%
    select(Protein.Group, Stripped.Sequence, Precursor.Id) %>%
    rename(Sequence = Stripped.Sequence) %>%
    mutate(Mix = paste0("Mix", i))
})

pr_matrix <- bind_rows(pr_list_raw) %>%
  pivot_wider(names_from = Mix, values_from = Mix, values_fill = NA_character_) %>%
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
  mutate(Count = rowSums(!is.na(select(., starts_with("Mix")))))

# 2.3 Without Carryover Lists
list_wo_base <- lapply(1:12, function(i) {
  curr <- paste0("Mix", i)
  prev <- paste0("Mix", i - 1)
  if (i == 1) target <- pr_matrix %>% filter(!is.na(!!sym(curr)))
  else        target <- pr_matrix %>% filter(!is.na(!!sym(curr)) & is.na(!!sym(prev)))
  target %>% arrange(desc(Count))
})

# ------------------------------------------------------------------
# 3. Type 1/2/3
# ------------------------------------------------------------------

list_true_hits <- list()

for (i in 1:12) {
  #  Standardize Protein Groups
  raw_map <- read.csv(file.path(path_results, sprintf("Mix-%d.pr_matrix.tsv", i)), sep = "\t") %>%
    select(Protein.Group, Stripped.Sequence) %>%
    rename(Sequence = Stripped.Sequence) %>%
    distinct(Sequence, .keep_all = TRUE) 
  
  df_std <- list_wo_base[[i]] %>%
    select(-Protein.Group) %>% 
    left_join(raw_map, by = "Sequence")
  
  # Type 1/2/3 Classification & Splitting
  type1 <- df_std %>% filter(!grepl("/", Protein.Group))
  
  type2 <- df_std %>% 
    filter(grepl("/", Protein.Group)) %>%
    filter(grepl(paste0("Mix", i, "_"), Protein.Group)) %>%
    separate_rows(Protein.Group, sep = "/") %>%  
    filter(grepl(paste0("Mix", i, "_"), Protein.Group))
  
  type3 <- df_std %>% 
    filter(grepl("/", Protein.Group)) %>%
    filter(!grepl(paste0("Mix", i, "_"), Protein.Group))
  
  # Combine & Clean
  df_clean <- bind_rows(type1, type2, type3) %>%
    mutate(Protein.Group = str_replace_all(Protein.Group, pattern = "Mix[0-9]+_", replacement = "")) %>%
    distinct(Protein.Group, Sequence, .keep_all = TRUE) %>% 
    mutate(MatchKey = paste(Protein.Group, Sequence, sep = ","))
  
  #  Validate against Theoretical 
  theo_keys <- paste(theo_list[[i]]$Protein.Group, theo_list[[i]]$Sequence, sep = ",")
  
  df_final <- df_clean %>%
    mutate(Exist = ifelse(MatchKey %in% theo_keys, 1, 0)) %>%
    filter(Count <= 3, Exist == 1) 
  
  list_true_hits[[i]] <- df_final
}


# ------------------------------------------------------------------
# 4. Figure 5: Theoretical Filtering
# ------------------------------------------------------------------

list_theo_filtered <- lapply(theo_list, function(df) {
  df %>%
    mutate(length = nchar(Sequence)) %>%
    filter(length >= 7 & length <= 30) %>%
    distinct(Sequence, .keep_all = TRUE)
})

# ------------------------------------------------------------------
# 5. Figure 5A:  Identification
# ------------------------------------------------------------------

# Experimental Counts
stats_pep_exp <- map_dfr(list_true_hits, ~ data.frame(Value = nrow(.))) %>%
  mutate(Type = "Experimental (1.8.1)", Mix = paste0("Mix-", 1:12))

# Theoretical Counts
stats_pep_theo <- map_dfr(list_theo_filtered, ~ data.frame(Value = nrow(.))) %>%
  mutate(Type = "Theoretical", Mix = paste0("Mix-", 1:12))

# Combine
plot_data_5a <- bind_rows(stats_pep_theo, stats_pep_exp) %>%
  mutate(Mix = factor(Mix, levels = paste0("Mix-", 1:12)))

# Verify Count 
print(head(plot_data_5a %>% filter(Type == "Experimental (1.8.1)")))

# Plot 5A Left
fig5a_left <- ggplot() +
  geom_col(data = filter(plot_data_5a, Type == "Theoretical"), 
           aes(x = Mix, y = Value), fill = "#E8E9EB", width = 0.9) +
  geom_col(data = filter(plot_data_5a, Type == "Experimental (1.8.1)"), 
           aes(x = Mix, y = Value), fill = "#4C7AF1", width = 0.6) +
  labs(y = "Number of Peptides", x = NULL) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, color = "black"), panel.grid = element_blank()) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 42000))


# ------------------------------------------------------------------
# 6. Figure 5B: Protein Identification
# ------------------------------------------------------------------

# Experimental Proteins
stats_pg_exp <- map_dfr(list_true_hits, function(df) {
  data.frame(Value = n_distinct(df$Protein.Group))
}) %>% mutate(Type = "Experimental (1.8.1)", Mix = paste0("Mix-", 1:12))

# Theoretical Proteins
stats_pg_theo <- map_dfr(list_theo_filtered, function(df) {
  data.frame(Value = n_distinct(df$Protein.Group))
}) %>% mutate(Type = "Theoretical", Mix = paste0("Mix-", 1:12))

# Combine
plot_data_5b <- bind_rows(stats_pg_theo, stats_pg_exp) %>%
  mutate(Mix = factor(Mix, levels = paste0("Mix-", 1:12)))

# Plot 5B Left
fig5b_left <- ggplot() +
  geom_col(data = filter(plot_data_5b, Type == "Theoretical"), 
           aes(x = Mix, y = Value), fill = "#E8E9EB", width = 0.9) +
  geom_col(data = filter(plot_data_5b, Type == "Experimental (1.8.1)"), 
           aes(x = Mix, y = Value), fill = "#4C7AF1", width = 0.6) +
  labs(y = "Number of Proteins", x = NULL) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, color = "black"), panel.grid = element_blank()) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 2000))





