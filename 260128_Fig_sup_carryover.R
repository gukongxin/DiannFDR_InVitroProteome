library(tidyverse)

# ------------------------------------------------------------------
# 1.  Directly 
# ------------------------------------------------------------------
setwd("/data/1.8.1")

# ------------------------------------------------------------------
# 2. Read
# ------------------------------------------------------------------

# Mix-1 to Mix-12
pr_list <- lapply(1:12, function(i) {
  
  fname <- sprintf("Mix-%d.pr_matrix.tsv", i)
  
  read.delim(fname) %>%
    filter(Proteotypic == 1) %>%
    select(Precursor.Id) %>%
    mutate(Mix = paste0("Mix", i), Present = 1)
})

# Convert to matrix format
pr_matrix <- bind_rows(pr_list) %>%
  pivot_wider(names_from = Mix, values_from = Present)

# ------------------------------------------------------------------
# 3. Calculate Carryover
# ------------------------------------------------------------------

res <- tibble(Sample_ID = 2:12) %>%
  mutate(
    Sample = paste0("Mix-", Sample_ID),
    
    # Numerator: Precursors identified in BOTH the current and the previous Mix (Intersection)
    Carried = map_int(Sample_ID, ~{
      prev <- paste0("Mix", .x - 1)
      curr <- paste0("Mix", .x)
      sum(pr_matrix[[prev]] == 1 & pr_matrix[[curr]] == 1, na.rm = TRUE)
    }),
    
    # Denominator: Total number of identifications in the current Mix
    Total = map_int(Sample_ID, ~sum(!is.na(pr_matrix[[paste0("Mix", .x)]]))),
    
    # Calculate Percentage
    Ratio = (Carried / Total) * 100,
    Version = "1.8.1"
  )

# ------------------------------------------------------------------
# 4. Visualization
# ------------------------------------------------------------------
ggplot(res, aes(x = factor(Sample, levels = Sample), y = Carried)) +
  geom_col(fill = "#4C7AF1", width = 0.7) +
  geom_text(aes(label = sprintf("%.2f%%", Ratio)), vjust = -0.5, size = 3.5) +
  labs(title = "Precursor Carryover Analysis",
       x = NULL, y = "Carried Precursors") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank()
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))

# Save plot to disk
#ggsave("carryovers.pdf", width = 8, height = 4, dpi = 300)