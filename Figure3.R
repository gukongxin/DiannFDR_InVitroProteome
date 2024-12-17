library(tidyverse)
library(patchwork)

# -----------------------------------------------------------------------------
### figure 3a

# read TheoreticalPeptides

# Read files
TheoreticalPepList <- lapply(1:12, function(i) {
  read.csv(paste0('TheoreticalPep_Mix-', i, ".txt"), sep = "\t")
})

# Combine list and process
Theoretical_all <- bind_rows(TheoreticalPepList) %>%
  group_by(Sequence) %>%
  summarise(across(
    starts_with("Mix"), 
    ~ ifelse(sum(., na.rm = TRUE) == 0, NA, sum(., na.rm = TRUE))
    ),
    Protein.Group = first(Protein.Group)) %>%
  mutate(Count = rowSums(!is.na(select(., -Sequence, -Protein.Group)))) %>%
  select(Protein.Group, Sequence, starts_with("Mix"), Count) %>%  
  mutate(PepLength = nchar(Sequence)) %>%  
  filter(PepLength <= 30 & PepLength >= 7)  

# Calculate frequency/count
SharePep_Theor <- Theoretical_all %>%
  count(Count) %>%
  mutate(Ratio = round(n / sum(n) * 100, 2),
         AccPercentage = cumsum(Ratio),
         log10 = log10(n))

# figure
figure_1a <- ggplot(SharePep_Theor) +
  geom_bar(
    aes(x = Count, y = log10), 
    stat = "identity", 
    width = 0.85, 
    fill = "#CCCCCC"
  ) +
  labs(
    x = "Peptides presence frequency", 
    y = "log10, peptides"
  ) +
  theme_bw() +
  theme(
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 10, color = "black"),
    axis.ticks = element_line(linewidth = 0.25),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  scale_y_continuous(
    expand = c(0, 0.1), 
    limits = c(0, 6)
  ) +
  geom_text(
    aes(x = Count, y = log10, label = n), 
    vjust = 0.5, 
    hjust = 1.2, 
    size = 3.5, 
    color = "black", 
    angle = 90
  )

figure_1a

# ------------------------------------------------------------------------------

# process peptide list
process_pep <- function(prlist, mix_prefix = "Mix") {
  # Combine and process data
  prall <- bind_rows(prlist) %>%
    pivot_wider(
      names_from = Mix, 
      values_from = Mix
    ) %>%
    group_by(Precursor.Id) %>%
    summarise(
      across(
        starts_with(mix_prefix), 
        ~ ifelse(sum(!is.na(.), na.rm = TRUE) == 0, NA, sum(!is.na(.), na.rm = TRUE))
      ),
      Stripped.Sequence = first(Stripped.Sequence)
    ) %>%
    mutate(
      Count = rowSums(!is.na(select(., -Stripped.Sequence, -Precursor.Id)))
    ) %>%
    rename(
      Sequence = Stripped.Sequence
    ) %>%
    select(
      Precursor.Id, Sequence, starts_with(mix_prefix), Count
    )
  
  # Calculate frequency/counts
  SharedPep <- prall %>%
    distinct(Sequence, .keep_all = TRUE) %>%
    count(Count) %>%
    mutate(
      Ratio = n / sum(n) * 100,
      AccPercentage = cumsum(Ratio),
      log10 = log10(n)
    )
  
  # Return the processed data
  list(
    prall,
    SharedPep
  )
}

# peptide frequency figure
p_shared_pep <- function(df, fill_color) {
  ggplot(df) +
    geom_bar(
      aes(x = Count, y = log10), 
      stat = "identity", 
      width = 0.85, 
      fill = fill_color
    ) +
    theme_bw() +
    theme(
      axis.title = element_text(size = 10),
      axis.text = element_text(size = 10, color = "black"),
      axis.ticks = element_line(linewidth = 0.25),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    scale_y_continuous(
      expand = c(0, 0.1), 
      limits = c(0, 6)
    ) +
    geom_text(
      aes(x = Count, y = log10, label = n), 
      vjust = 0.5, 
      hjust = 1.2, 
      size = 3.5, 
      color = "black", 
      angle = 90
    )
}

#-------------------------------------------------------------------------------
### figure 3b

# read diann 181 

# Read pr.matrix report
prlist_181 <- lapply(1:12, function(i) {
  read.csv(paste0("Mix-", i, ".pr_matrix.tsv"), sep = "\t") %>%
    filter(Proteotypic == 1) %>%
    select(Stripped.Sequence, Precursor.Id) %>%
    mutate(Mix = paste0("Mix", i))
})


# Calculate 
Pep_181 <- process_pep(prlist_181)

# Create figure
figure_1b <- p_shared_pep(Pep_181[[2]], "#FFCD92")
figure_1b
# ------------------------------------------------------------------------------
### figure 3c

# read diann 192

# Read pr.matrix report
prlist_19 <- lapply(1:12, function(i) {
  read.csv(paste0("Mix-", i, ".pr_matrix.tsv"), sep = "\t") %>%
    filter(Proteotypic == 1) %>%
    select(Stripped.Sequence, Precursor.Id) %>%
    mutate(Mix = paste0("Mix", i))
})

# Calculate 
Pep_19 <- process_pep(prlist_19)

# Create figure
figure_1c <- p_shared_pep(Pep_19[[2]], "#91DA91")
figure_1c

# ------------------------------------------------------------------------------
### figure 4d

# Function to create peptide lists for each mix (with carryovers)
create_mixlist <- function(prall) {
  lapply(1:12, function(i) {
    Mixname <- paste0("Mix", i)
    prall %>%
      filter(!is.na(!!sym(Mixname))) %>%
      select(Sequence, all_of(Mixname), Count) %>%
      arrange(desc(Count)) %>%
      distinct(Sequence, .keep_all = TRUE)
  })
}

# Function to calculate FDR
calculateFDR <- function(mixList, theoreticalPepList) {
  FDR_list <- lapply(1:12, function(countLimit) {
    temp_list <- lapply(1:length(mixList), function(i) {
      df <- mixList[[i]]
      fasta <- theoreticalPepList[[i]]
      df %>%
        filter(Count <= countLimit) %>%
        mutate(Exist = ifelse(Sequence %in% fasta$Sequence, 1, 0)) %>%
        summarise(Alldiscovery = n(), Fdiscovery = sum(Exist == 0)) %>%
        mutate(SharedNumber = countLimit)
    })
    bind_rows(temp_list)
  })
  list(FDR_list = FDR_list,
       FDR_summary = bind_rows(FDR_list) %>%
         mutate(FDR = round(Fdiscovery / Alldiscovery * 100, 2)) %>%
         group_by(SharedNumber) %>%
         summarise(FDR = round(mean(FDR), 2)))
}

# calculate FDR
Mix_181_w <- create_mixlist(Pep_181[[1]])
Mix_19_w <- create_mixlist(Pep_19[[1]])

FDR_181_w <- calculateFDR(Mix_181_w, TheoreticalPepList)
FDR_19_w <- calculateFDR(Mix_19_w, TheoreticalPepList)

# FDR summary
FDR_sumy_181_w <- FDR_181_w$FDR_summary
FDR_list_181_w <- FDR_181_w$FDR_list

FDR_sumy_19_w <- FDR_19_w$FDR_summary
FDR_list_19_w <- FDR_19_w$FDR_list

# Combine FDR results (with carryovers)
FDR_sumy_w <- bind_cols(FDR_sumy_181_w, FDR_sumy_19_w) %>%
  select(1, 2, 4) %>%
  setNames(c("Count", "diann181", "diann19"))



# Function to create peptide lists without carryovers
create_mixlist_wo <- function(prall) {
  Mix1 <- prall %>%
    filter(!is.na(Mix1)) %>%
    select(Sequence, Mix1, Count) %>%
    arrange(desc(Count)) %>%
    distinct(Sequence, .keep_all = TRUE)
  
  Mix2to12 <- lapply(2:12, function(i) {
    Mixname <- paste0("Mix", i)
    Mixname_1 <- paste0("Mix", i - 1)
    prall %>%
      filter(!is.na(!!sym(Mixname))) %>%
      select(Sequence, !!sym(Mixname_1), !!sym(Mixname), Count) %>%
      arrange(desc(Count)) %>%
      distinct(Sequence, .keep_all = TRUE) %>%
      filter(is.na(!!sym(Mixname_1))) %>%
      select(Sequence, !!sym(Mixname), Count)
  })
  c(list(Mix1), Mix2to12)
}

#  FDR without carryovers 
Mix_181_wo <- create_mixlist_wo(Pep_181[[1]])
Mix_19_wo  <- create_mixlist_wo(Pep_19[[1]])

FDR_181_wo <- calculateFDR(Mix_181_wo, TheoreticalPepList)
FDR_19_wo <- calculateFDR(Mix_19_wo, TheoreticalPepList)

# Extract FDR summary and FDR list
FDR_sumy_181_wo <- FDR_181_wo $FDR_summary
FDR_list_181_wo <- FDR_181_wo $FDR_list

FDR_sumy_19_wo <- FDR_19_wo$FDR_summary
FDR_list_19_wo <- FDR_19_wo$FDR_list

# Combine FDR results (without carryovers)
FDR_sumy_wo <- bind_cols(FDR_sumy_181_wo, FDR_sumy_19_wo) %>%
  select(1, 2, 4) %>%
  setNames(c("Count", "diann181", "diann19"))

# Combine FDR results for plotting
FDR_figure <- bind_cols(FDR_sumy_w, FDR_sumy_wo) %>%
  select(1, 2, 3, 5, 6) %>%
  setNames(c("Count",
             "DIA-NN 1.8.1 with carryovers",
             "DIA-NN 1.9 with carryovers",
             "DIA-NN 1.8.1 without carryovers",
             "DIA-NN 1.9 without carryovers"))

# Convert data to long format for plotting
FDR_figure_long <- pivot_longer(FDR_figure, cols = -Count,
                                names_to = "version",
                                values_to = "FDR")



# ==============================================================================
# Set factor levels for Count
FDR_figure_long$Count <- factor(as.character(FDR_figure_long$Count), levels = as.character(12:1))

# Create figure
figure_1d <- ggplot(FDR_figure_long, aes(x = Count, y = FDR, group = version, color = version)) +
  geom_line(aes(linetype = ifelse(grepl("with carryovers", version), "solid", "dashed")),
            linewidth = 1) +  
  geom_point(size = 2, shape = 21, aes(fill = version), color = "black", stroke = 0.8) +
  theme_bw() +
  scale_x_discrete() +
  scale_color_manual(values = c(
    "DIA-NN 1.8.1 without carryovers" = "#FFCD92",
    "DIA-NN 1.8.1 with carryovers" = "#FFCD92",
    "DIA-NN 1.9 without carryovers" = "#91DA91",
    "DIA-NN 1.9 with carryovers" = "#91DA91"
  )) +
  scale_fill_manual(values = c(
    "DIA-NN 1.8.1 without carryovers" = "#FFCD92",
    "DIA-NN 1.8.1 with carryovers" = "#FFCD92",
    "DIA-NN 1.9 without carryovers" = "#91DA91",
    "DIA-NN 1.9 with carryovers" = "#91DA91"
  )) +  
  theme(
    legend.position = c(0.7, 0.8),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 10.1),
    axis.text.y = element_text(size = 10.1),
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    axis.text = element_text(size = 5, color = "black"),
    axis.ticks = element_line(linewidth = 0.25)
  ) +
  labs(x = "Peptides detection frequency", y = "Average FDR (%)")
figure_1d

