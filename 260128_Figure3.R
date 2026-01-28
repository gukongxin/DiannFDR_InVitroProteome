library(tidyverse)

# ------------------------------------------------------------------
# 1. Path 
# ------------------------------------------------------------------
# Path to the theoretical peptide
path_theoretical <- "./InSilico_TheoreticalPeptides"
path_results     <- "./data/1.8.1"

# ------------------------------------------------------------------
# 2. Read
# ------------------------------------------------------------------

# Load theoretical peptide
theo_pep_list <- lapply(1:12, function(i) {
  read.csv(file.path(path_theoretical, paste0("TheoreticalPep_Mix-", i, ".txt")), sep = "\t") %>% 
    select(Sequence)
})

# Load diann search results
pr_list <- lapply(1:12, function(i) {
  # Read the precursor matrix
  read.delim(file.path(path_results, sprintf("Mix-%d.pr_matrix.tsv", i))) %>%
    filter(Proteotypic == 1) %>%               # Filter for proteotypic
    select(Precursor.Id, Stripped.Sequence) %>%
    mutate(Mix = paste0("Mix", i))
})

# ------------------------------------------------------------------
# 3. Data Processing
# ------------------------------------------------------------------
# Convert to a wide format
pr_all <- bind_rows(pr_list) %>%
  pivot_wider(names_from = Mix, values_from = Mix, values_fill = NA) %>%
  mutate(
    # Calculate detection frequency (Count) across 12 mixtures
    Count = rowSums(!is.na(select(., starts_with("Mix")))),
    Sequence = Stripped.Sequence
  )

# ------------------------------------------------------------------
# 4. Define Peptide Lists (With / Without Carryover)
# ------------------------------------------------------------------

# With Carryovers
# Include all peptides identified in the current mixture
list_with <- lapply(1:12, function(i) {
  pr_all %>% filter(!is.na(!!sym(paste0("Mix", i))))
})

# Without Carryovers
# Include peptides identified in the current mixture but not in the previous one
list_wo <- lapply(1:12, function(i) {
  curr <- paste0("Mix", i)
  prev <- paste0("Mix", i - 1)
  
  if (i == 1) {
    pr_all %>% filter(!is.na(!!sym(curr)))
  } else {
    pr_all %>% filter(!is.na(!!sym(curr)) & is.na(!!sym(prev)))
  }
})

# ------------------------------------------------------------------
# 5. FDR Calculation Function
# ------------------------------------------------------------------
calculate_fdr <- function(mix_list, theo_list) {
  map_dfr(1:12, function(threshold) {
    
      # Calculate FDR for each individual mixture
    daily_stats <- map2_dfr(mix_list, theo_list, function(df, theo) {
      # Filter peptides based on the frequency threshold
      target <- df %>% filter(Count <= threshold)
      

      if(nrow(target) == 0) return(tibble(FDR = 0))
      
      # Identify false discoveries
      n_false <- sum(!target$Sequence %in% theo$Sequence)
      n_total <- nrow(target)
      
      # Calculate FDR 
      tibble(FDR = (n_false / n_total) * 100)
    })
    
     # Calculate the average FDR 
    daily_stats %>%
      summarise(Average_FDR = mean(FDR, na.rm = TRUE)) %>%
      mutate(Threshold = threshold)
  })
}

# ------------------------------------------------------------------
# 6. Resuts and Plotting
# ------------------------------------------------------------------

# Calculate FDR for both strategies
res_w  <- calculate_fdr(list_with, theo_pep_list) %>% mutate(Type = "With Carryovers")
res_wo <- calculate_fdr(list_wo, theo_pep_list)   %>% mutate(Type = "Without Carryovers")

# Combine results
plot_data <- bind_rows(res_w, res_wo)

# Generate Plot
ggplot(plot_data, aes(x = factor(Threshold, levels = 12:1), y = Average_FDR, group = Type)) +
  geom_line(aes(linetype = Type), color = "#4C7AF1", linewidth = 0.8) +
  geom_point(aes(fill = Type), shape = 21, size = 3) +
  scale_fill_manual(values = c("With Carryovers" = "#4C7AF1", "Without Carryovers" = "white")) +
  scale_linetype_manual(values = c("With Carryovers" = "dashed", "Without Carryovers" = "solid")) +
  labs(
    x = "Peptides Detection count",
    y = "Average FDR (%)"
  ) +
  theme_bw() +
  theme(
    legend.position = c(0.7, 0.8), 
    legend.title = element_blank(),
    axis.text = element_text(color = "black")
  )

# Save the figure
#ggsave("Figure3_FDR.pdf", width = 5, height = 4)