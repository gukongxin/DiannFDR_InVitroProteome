library(ggvenn)
source("Figure4.R")

# ------------------------------------------------------------------------------
###  Figure 5A

pep_theor <- list()

for (i in 1:12) {
  FastaName <- paste0('TheoreticalPep_Mix-',i,".txt")
  pep_theor[[i]] <- read.csv(FastaName,sep = "\t") %>% 
    select(Protein.Group,Sequence) %>% 
    mutate(length = nchar(Sequence)) %>% 
    filter(length >=7 & length <=30)  %>% 
    select(Sequence) %>% 
    distinct(Sequence)
}


# identification
PepId_181 <- cleanPep_181_True %>% lapply(nrow) %>% unlist() %>% data.frame()
PepId_19 <- cleanPep_19_True %>% lapply(nrow) %>% unlist() %>% data.frame()
PepId_theor <- pep_theor %>% lapply(nrow) %>% unlist() %>% data.frame() %>% 
  mutate(Sample = paste0("Mix-",1:12)) %>% 
  setNames(c("values","Sample"))

PepId <- bind_cols(PepId_181,PepId_19) %>% 
  setNames(c("diann181","diann19")) %>% 
  mutate(Sample = paste0("Mix-",1:12)) 
  

#signif difference
t.test(PepId$diann181,PepId$diann19,paired = T)

#recovery
pepreovery <- bind_cols(PepId,PepId_theor) %>% 
  mutate(recovery181 = diann181/values*100,
         recovery19 = diann19/values*100)
  
#convert to long
PepId_long <- pivot_longer(PepId,cols = -Sample, names_to = "version", values_to = "values")

#factor
PepId_theor$Sample <- factor(PepId_theor$Sample,levels = c(paste0("Mix-",1:12)))

#figure
figure_5a_left <- ggplot() +
  geom_col(data = PepId_theor, aes(x = Sample, y = values), fill = "#E8E9EB",width = 0.9) +
  geom_col(data = PepId_long, aes(x = Sample, y = values, fill = version), position = "dodge") +
  scale_fill_manual(values = c("#FFCD92", "#91DA91")) +
  theme_bw() +
  labs(x = NULL, y = "Number of peptides")+
  theme(legend.position="none",
        axis.text.x = element_text(size=12,angle=45,color = "black",hjust = 1),
        axis.text.y = element_text(size = 10,color = "black"),
        axis.title.y = element_text(color = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())


figure_5a_right <- ggplot(PepId_long, aes(x = version, y = values, fill = version)) +
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
  labs(x=NULL,y="Average number of peptides")+
  ylim(0,15000)

figure_5a <- figure_5a_left+figure_5a_right+plot_layout(widths = c(3,1))
figure_5a
#-------------------------------------------------------------------------------
###  Figure 5B

#total number of proteins, read theoretical proteins and filter

Pg_theor <- list()

for (i in 1:12) {
  FastaName <- paste0('TheoreticalPep_Mix-',i,".txt")
  Pg_theor[[i]] <- read.csv(FastaName,sep = "\t") %>% 
    select(Protein.Group,Sequence) %>% 
    mutate(length = nchar(Sequence)) %>% 
    filter(length >=7 & length <=30)  %>% 
    select(Protein.Group) %>% 
    distinct(Protein.Group)
}

# identification
PgId_181 <- cleanPep_181_True %>% lapply(function(x){
  x %>% select(Protein.Group) %>% 
    distinct() %>% 
    nrow()
}) %>% 
  unlist() %>% 
  data.frame()

PgId_19 <- cleanPep_19_True %>% lapply(function(x){
  x %>% select(Protein.Group) %>% 
    distinct() %>% 
    nrow()
}) %>% 
  unlist() %>% 
  data.frame()

PgId_theor <- Pg_theor %>% lapply(nrow) %>% unlist() %>% data.frame() %>% 
  mutate(Sample = paste0("Mix-",1:12)) %>% 
  setNames(c("values","Sample"))

PgId <- bind_cols(PgId_181,PgId_19) %>% 
  setNames(c("diann181","diann19")) %>% 
  mutate(Sample = paste0("Mix-",1:12)) 


#signif difference
t.test(PgId$diann181,PgId$diann19,paired = T)

#recovery
pgreovery <- bind_cols(PgId,PgId_theor) %>% 
  mutate(recovery181 = diann181/values*100,
         recovery19 = diann19/values*100)

#convert to long
PgId_long <- pivot_longer(PgId,cols = -Sample, names_to = "version", values_to = "values")

#factor
PgId_theor$Sample <- factor(PgId_theor$Sample,levels = c(paste0("Mix-",1:12)))

#figure
figure_5b_left <- ggplot() +
  geom_col(data = PgId_theor, aes(x = Sample, y = values), fill = "#E8E9EB",width = 0.9) +
  geom_col(data = PgId_long, aes(x = Sample, y = values, fill = version), position = "dodge") +
  scale_fill_manual(values = c("#FFCD92", "#91DA91")) +
  theme_bw() +
  labs(x = NULL, y = "Number of peptides")+
  theme(legend.position="none",
        axis.text.x = element_text(size=12,angle=45,color = "black",hjust = 1),
        axis.text.y = element_text(size = 10,color = "black"),
        axis.title.y = element_text(color = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())


figure_5b_right <- ggplot(PgId_long, aes(x = version, y = values, fill = version)) +
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
  labs(x=NULL,y="Average number of proteins")+
  ylim(0,2000)

figure_5b <- figure_5b_left+figure_5b_right+plot_layout(widths = c(3,1))
figure_5b

# ------------------------------------------------------------------------------
### Figure 5C

# peptides
venn_pep_181 <- bind_rows(cleanPep_181_True) %>% 
  select(Sequence) %>% 
  distinct() %>% 
  pull()

venn_pep_19 <- bind_rows(cleanPep_19_True) %>% 
  select(Sequence) %>% 
  distinct() %>% 
  pull()

# combine
venn_pep <- list(
  Pep_181 = venn_pep_181,
  Pep_19 = venn_pep_19
)

p_venn_pep <- ggvenn(venn_pep, fill_color = c("#FFCD92","#91DA91"))

# proteins
venn_pg_181 <- bind_rows(cleanPep_181_True) %>% 
  select(Protein.Group) %>% 
  distinct()%>% 
  pull()

venn_pg_19 <- bind_rows(cleanPep_19_True) %>% 
  select(Protein.Group) %>% 
  distinct()%>% 
  pull()

# combine
venn_pg <- list(
  Pg_181 = venn_pg_181,
  Pg_19 = venn_pg_19
)

p_venn_pg <-ggvenn(venn_pg, fill_color = c("#FFCD92","#91DA91"))
p_venn <- p_venn_pep + p_venn_pg

# ------------------------------------------------------------------------------
###  peptide length

# diann 1.8.1
length181 <- lapply(cleanPep_181_True, function(x){
  x %>%  mutate(peplength = nchar(Sequence))
}) %>% bind_rows() 

length181_freq <- table(length181$peplength) %>% 
  as.data.frame() %>% 
  mutate(ratio = Freq/ 130187*100)

# diann 1.9
length19 <- lapply(cleanPep_19_True, function(x){
  x %>% mutate(peplength = nchar(Sequence))
}) %>% bind_rows()

length19_freq <- table(length19$peplength) %>% 
  as.data.frame() %>% 
  mutate(ratio = Freq/ 134873*100)

# theoretical
lengthFasta <- lapply(pep_theor, function(x){
  x %>% mutate(peplength = nchar(Sequence))
}) %>% 
  bind_rows() %>% 
  distinct(Sequence,.keep_all = T)

lengthFasta_freq <- table(lengthFasta$peplength) %>% 
  as.data.frame() %>% 
  mutate(ratio = Freq/ 339463*100)

# figure
figure_5d <- ggplot() +
  geom_smooth(data = length181_freq, aes(x = Var1, y = ratio,group=1), method = "loess", se = FALSE,color="#FFCD92")+
  geom_smooth(data = length19_freq, aes(x = Var1, y = ratio,group=1), method = "loess", se = FALSE,color="#91DA91")+
  geom_smooth(data = lengthFasta_freq, aes(x = Var1, y = ratio,group=1), method = "loess", se = FALSE,color="gray")+
  theme_bw()+
  theme(
        axis.text.x = element_text(size=12,color = "black",hjust = 1),
        axis.text.y = element_text(size = 10,color = "black"),
        axis.title.y = element_text(color = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+
  labs(x= "Peptides length",y="Percentage (%)")

#------------------------pathchwork------------------------------------------
figure_5cde <- p_venn + figure_5d

figure_5 <- figure_5a/figure_5b/figure_5cde
figure_5

