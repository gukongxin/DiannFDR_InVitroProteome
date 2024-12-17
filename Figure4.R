source("Figure3.R")

#-------------------------------------------------------------------------------
### figure 4a

# Peptides FDR count <= 3
pepFDR_181 <- FDR_list_181_wo[[3]]
pepFDR_19 <- FDR_list_19_wo[[3]]

# calculate FDR
pepFDR <- bind_cols(pepFDR_181,pepFDR_19) %>% 
  mutate(FDR_181 = .[[2]]/.[[1]]*100,
         FDR_19 = .[[5]]/.[[4]]*100,
         Mixtures = paste0("Mix-",1:12)) %>% 
  select(7,8,9)

#signif difference
t.test(pepFDR$FDR_181, pepFDR$FDR_19, paired = TRUE) 

#covert to long
pepFDR_long <- pivot_longer(pepFDR,cols = -Mixtures,names_to = "version", values_to = "FDR")

#factor levels
pepFDR_long$Mixtures <- factor(pepFDR_long$Mixtures,levels = c(paste0("Mix-",1:12)))

figure_4a_left <- ggplot(pepFDR_long,aes(x=Mixtures,y=FDR,fill=version))+
  geom_bar(stat = "identity", position = position_dodge(), width = 0.7,color= "black",linewidth=0.3)+
  ylim(0,1)+
  scale_fill_manual(values = c("#FFCD92","#91DA91"))+
  theme_bw()+
  labs(y="Peptides FDR(%)",x=NULL)+
  theme(legend.position = c(0.8,0.9),
        axis.text.x = element_text(size=12,angle=45,color = "black",hjust = 1),  
        axis.text.y = element_text(size = 10),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

figure_4a_right <- ggplot(pepFDR_long, aes(x = version, y = FDR, fill = version)) +
  stat_summary(fun = mean, geom = "bar", width = 0.6, color = "black",linewidth=0.3) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2,linewidth=0.3) +
  geom_point(position = position_jitter(width = 0.1), size = 1, color = "black") + 
  theme_bw()+
  theme(legend.position = "none",
        axis.text.x = element_text(size=12,angle=45,color = "black",hjust = 1),  
        axis.text.y = element_text(size = 10,color = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+
  scale_fill_manual(values = c("#FFCD92","#91DA91"))+
  labs(x=NULL,y="Average peptides FDR,%")+
  ylim(0,1)

figure_4a <- figure_4a_left + figure_4a_right + plot_layout(widths = c(3, 1))
figure_4a

# ------------------------------------------------------------------------------

### Figure 4b

# versionb 1.8.1

# Read output for DIANN181
pglist_181 <- list()
for (i in 1:12) {
  filename <- paste0("Mix-", i, ".pr_matrix.tsv")
  pglist_181[[i]] <- read.csv(filename, sep = "\t") %>%
    select(Protein.Group, Stripped.Sequence) %>%
    rename(Sequence = Stripped.Sequence) %>%
    distinct(Sequence, .keep_all = TRUE)
}

# Match peptides and protein groups for DIANN181
for (i in 1:12) {
  Mix_181_wo[[i]] <- Mix_181_wo[[i]] %>%
    left_join(pglist_181[[i]], by = "Sequence")
}

# Group by type and process Protein.Group for DIANN181
type1_181 <- list()
type2_181 <- list()
type3_181 <- list()
Mix_pg_all_181 <- list()

for (i in 1:12) {
  type1_181[[i]] <- Mix_181_wo[[i]] %>%
    filter(!grepl("/", Protein.Group))
  
  type2_181[[i]] <- Mix_181_wo[[i]] %>%
    filter(grepl("/", Protein.Group)) %>%
    filter(grepl(paste0("Mix", i, "_"), Protein.Group)) %>%
    separate_rows(Protein.Group, sep = "/") %>%
    filter(grepl(paste0("Mix", i, "_"), Protein.Group))
  
  type3_181[[i]] <- Mix_181_wo[[i]] %>%
    filter(grepl("/", Protein.Group)) %>%
    filter(!grepl(paste0("Mix", i, "_"), Protein.Group))
  
  Mix_181_wo[[i]] <- bind_rows(type1_181[[i]], type2_181[[i]], type3_181[[i]]) %>%
    mutate(Protein.Group = str_replace_all(Protein.Group, pattern = "Mix[0-9]+_", replacement = "")) %>%
    distinct(Protein.Group, Sequence, .keep_all = TRUE)
  
  Mix_pg_all_181[[i]] <- Mix_181_wo[[i]] %>%
    mutate(Name = paste(Protein.Group, Sequence, sep = ",")) 
}

# Prepare TheoreticalPeptideList for DIANN181
TheoreticalPepList_for181 <- list()
for (i in 1:12) {
  TheoreticalPepList_for181[[i]] <- TheoreticalPepList[[i]] %>%
    mutate(Name = paste(Protein.Group, Sequence, sep = ",")) %>%
    select(Name)
}

# Search and update existence in DIANN181
for (i in 1:12) {
  Mix_pg_all_181[[i]] <- Mix_pg_all_181[[i]] %>%
    mutate(Exist = ifelse(Name %in% TheoreticalPepList_for181[[i]]$Name, 1, 0)) %>%
    separate(Name, c("Protein.Group", "Sequence"), sep = ",")
}

# Calculate Protein FDR for DIANN181 with count <= 3
pgFDR_count3_181 <- lapply(Mix_pg_all_181, function(x) {
  filter(x, Count <= 3) %>%
    group_by(Protein.Group) %>%
    arrange(Exist) %>%
    distinct(Protein.Group, .keep_all = TRUE)
})

PgFDR_181 <- map_dfr(pgFDR_count3_181, ~ data.frame(
  All = nrow(.x),
  False = sum(.x$Exist == 0)
)) %>% mutate(FDR = False / All * 100)

# version 1.9

# Read output for DIANN19
pglist_19 <- list()
for (i in 1:12) {
  filename <- paste0("Mix-", i, ".pr_matrix.tsv")
  pglist_19[[i]] <- read.csv(filename, sep = "\t") %>%
    select(Protein.Group, Stripped.Sequence) %>%
    rename(Sequence = Stripped.Sequence) %>%
    distinct(Sequence, .keep_all = TRUE)
}

# Match peptides and protein groups for DIANN19
for (i in 1:12) {
  Mix_19_wo[[i]] <- Mix_19_wo[[i]] %>%
    left_join(pglist_19[[i]], by = "Sequence")
}

# Group by type and process Protein.Group for DIANN19
type1_19 <- list()
type2_19 <- list()
type3_19 <- list()
Mix_pg_all_19 <- list()

for (i in 1:12) {
  type1_19[[i]] <- Mix_19_wo[[i]] %>%
    filter(!grepl("/", Protein.Group))
  
  type2_19[[i]] <- Mix_19_wo[[i]] %>%
    filter(grepl("/", Protein.Group)) %>%
    filter(grepl(paste0("Mix", i, "_"), Protein.Group)) %>%
    separate_rows(Protein.Group, sep = "/") %>%
    filter(grepl(paste0("Mix", i, "_"), Protein.Group))
  
  type3_19[[i]] <- Mix_19_wo[[i]] %>%
    filter(grepl("/", Protein.Group)) %>%
    filter(!grepl(paste0("Mix", i, "_"), Protein.Group))
  
  Mix_19_wo[[i]] <- bind_rows(type1_19[[i]], type2_19[[i]], type3_19[[i]]) %>%
    mutate(Protein.Group = str_replace_all(Protein.Group, pattern = "Mix[0-9]+_", replacement = "")) %>%
    distinct(Protein.Group, Sequence, .keep_all = TRUE)
  
  Mix_pg_all_19[[i]] <- Mix_19_wo[[i]] %>%
    mutate(Name = paste(Protein.Group, Sequence, sep = ",")) %>%
    select(Name, Count)
}

# Prepare TheoreticalPeptideList for DIANN19
TheoreticalPep_for19 <- list()
for (i in 1:12) {
  TheoreticalPep_for19[[i]] <- TheoreticalPepList[[i]] %>%
    mutate(Name = paste(Protein.Group, Sequence, sep = ",")) %>%
    select(Name)
}

# Search and update existence in DIANN19
for (i in 1:12) {
  Mix_pg_all_19[[i]] <- Mix_pg_all_19[[i]] %>%
    mutate(Exist = ifelse(Name %in% TheoreticalPep_for19[[i]]$Name, 1, 0)) %>%
    separate(Name, c("Protein.Group", "Sequence"), sep = ",")
}

# Calculate Protein FDR for DIANN19 with count <= 3
pgFDR_count3_19 <- lapply(Mix_pg_all_19, function(x) {
  filter(x, Count <= 3) %>%
    group_by(Protein.Group) %>%
    arrange(Exist) %>%
    distinct(Protein.Group, .keep_all = TRUE)
})

PgFDR_19 <- map_dfr(pgFDR_count3_19, ~ data.frame(
  All = nrow(.x),
  False = sum(.x$Exist == 0)
)) %>% mutate(FDR = False / All * 100)

# combine and make figure

#combine cols
pgFDR <- bind_cols(PgFDR_181,PgFDR_19) %>% 
  mutate(Mixtures = paste0("Mix-",1:12)) %>% 
  select(3,6,7) %>% 
  setNames(c("FDR_181","FDR19","Mixtures"))

#signif difference
t.test(pgFDR$FDR_181,pgFDR$FDR19,paired = T)

#change to long
pgFDR_long <- pivot_longer(pgFDR,cols = -Mixtures,names_to = "version",values_to = "FDR")

#factor levels
pgFDR_long$Mixtures <- factor(pgFDR_long$Mixtures,levels = c(paste0("Mix-",1:12)))

figure_4b_right <- ggplot(pgFDR_long,aes(x=Mixtures,y=FDR,fill=version))+
  geom_bar(stat = "identity", position = position_dodge(), width = 0.7,color= "black",linewidth=0.3)+
  ylim(0,4)+
  scale_fill_manual(values = c("#FFCD92","#91DA91"))+
  theme_bw()+
  labs(y="Proteins FDR(%)",x=NULL)+
  theme(legend.position = c(0.8,0.9),
        axis.text.x = element_text(size=12,angle=45,color = "black",hjust = 1),  
        axis.text.y = element_text(size = 10),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

figure_4b_left <- ggplot(pgFDR_long, aes(x = version, y = FDR, fill = version)) +
  stat_summary(fun = mean, geom = "bar", width = 0.6, color = "black",linewidth=0.3) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2,linewidth=0.3) +
  geom_point(position = position_jitter(width = 0.1), size = 1, color = "black") +
  theme_bw()+
  theme(legend.position = "none",
        axis.text.x = element_text(size=12,angle=45,color = "black",hjust = 1),  
        axis.text.y = element_text(size = 10,color = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+
  scale_fill_manual(values = c("#FFCD92","#91DA91"))+
  labs(x=NULL,y="Average proteins FDR,%")+
  ylim(0,4)

figure_4b <- figure_4b_right + figure_4b_left + plot_layout(widths = c(3, 1))
figure_4b

#-------------------------------------------------------------------------------
###  Figure 4C
#clean peptides list without contaminnat

cleanPep_181 <- lapply(Mix_pg_all_181, function(x){
  x %>% filter(Count <= 3)
})

cleanPep_19 <- lapply(Mix_pg_all_19, function(x){
  x %>% filter(Count <= 3)
})


#only true
cleanPep_181_True <- lapply(cleanPep_181, function(x){
  x %>% filter(Exist==1)
})

cleanPep_19_True <- lapply(cleanPep_19, function(x){
  x %>% filter(Exist==1)
})

# frequency
#181
Truelist181 <- lapply(cleanPep_181_True , function(x) {
  
  freq_data <- table(x$Protein.Group) %>% as.data.frame()
  
})

Truelist181 <- bind_rows(Truelist181)

#19
Truelist19 <- lapply(cleanPep_19_True , function(x) {
  
  freq_data <- table(x$Protein.Group) %>% as.data.frame()
  
})

Truelist19 <- bind_rows(Truelist19)

#figure
figure_4c <- ggplot()+
  geom_histogram(data=Truelist181,aes(Freq),binwidth = 1, fill = "#FFCD92",color = "black",linewidth=0.3,alpha=0.7)+
  geom_histogram(data=Truelist19,aes(Freq),binwidth = 1, fill = "#91DA91",color = "black",linewidth=0.3,alpha=0.7)+
  theme_bw()+
  labs(y="Proteins",x= "Peptides per protein")+
  theme(legend.position="none",
        axis.text.x = element_text(size=12,color = "black",hjust = 1),
        axis.text.y = element_text(size = 10,color = "black"),
        axis.title.y = element_text(color = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

#  -----------------------------------------------------------------------------
###  Figure 4D 

#only false
cleanPep_181_False <- lapply(cleanPep_181, function(x){
  x %>% filter(Exist==0)
})

cleanPep_19_False <- lapply(cleanPep_19, function(x){
  x %>% filter(Exist==0)
})

# frequency
#181
Falselist181 <- lapply(cleanPep_181_False , function(x) {
  
  freq_data <- table(x$Protein.Group) %>% as.data.frame()
  
})

Falselist181 <- bind_rows(Falselist181)

#19
Falselist19 <- lapply(cleanPep_19_False , function(x) {
  
  freq_data <- table(x$Protein.Group) %>% as.data.frame()
  
})

Falselist19 <- bind_rows(Falselist19)

#figure
figure_4d <- ggplot()+
  geom_histogram(data=Falselist181,aes(Freq),binwidth = 1, fill = "#FFCD92",color = "black",linewidth=0.3,alpha=0.7)+
  geom_histogram(data=Falselist19,aes(Freq),binwidth = 1, fill = "#91DA91",color = "black",linewidth=0.3,alpha=0.7)+
  theme_bw()+
  labs(y="Proteins",x= "Peptides per protein")+
  theme(legend.position="none",
        axis.text.x = element_text(size=12,color = "black",hjust = 1),
        axis.text.y = element_text(size = 10,color = "black"),
        axis.title.y = element_text(color = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

#-------------------------------------------------------------------------------
###  Figure 4E

# diann 181

# True 
True_onepep_181 <- Truelist181 %>% filter(Freq==1) %>% nrow() %>% as.data.frame()
# false
False_onepep_181 <- Falselist181 %>% filter(Freq==1) %>% nrow() %>% as.data.frame()
OnePep_pg_181 <- bind_rows(True_onepep_181,False_onepep_181)%>% mutate(version=c("True","False"))

# diann 19
# True 
True_onepep_19 <- Truelist19 %>% filter(Freq==1) %>% nrow() %>% as.data.frame()
# false
False_onepep_19 <- Falselist19 %>% filter(Freq==1) %>% nrow() %>% as.data.frame()
OnePep_pg_19 <- bind_rows(True_onepep_19,False_onepep_19)%>% mutate(version=c("True","False"))

# combine
onepepprotein <- bind_rows(OnePep_pg_181,OnePep_pg_19) %>% 
  setNames(c("number","type")) %>% 
  mutate(version= c(rep("DIA-NN 1.8.1",2),rep("DIA-NN 1.9",2)))

# Figure 
figure_4e <- ggplot(onepepprotein) +
  aes(x = version, y = number, fill = type) +
  geom_col(width = 0.75) +
  scale_fill_hue(direction = 1) +
  theme_bw()+
  scale_fill_manual(values = c("#AD483C","#F3EED1"))+
  geom_text(aes(label=number))+
  theme(legend.position = "none",
        axis.text.x = element_text(size=12,angle=45,color = "black",hjust = 1),
        axis.text.y = element_text(size = 10,color = "black"),
        axis.title.y = element_text(color = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+
  labs(x=NULL,y="Proteins")

figure_4cde <- figure_4c+figure_4d+figure_4e+plot_layout(widths = c(1.5,1.5,1))

#-------------------------------patchwork-----------------------------------
figure_4 <- figure_4a/figure_4b/figure_4cde
figure_4

