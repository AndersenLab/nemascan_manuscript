require(tidyverse)
require(RColorBrewer)
require(magrittr)
require(ggtree)
require(ggplotify)
setwd(paste0(dirname(rstudioapi::getActiveDocumentContext()$path),"/.."))

#########################################################################################################
### Comparisons of Local Performance (Chromosome Regions, Divergent Regions, Swept/Divergent Strains) ###
#########################################################################################################
load(file = "data/File_S5.RData")
load(file = "data/File_S6.RData")

designations <- dat.genomic.locations %>%
  dplyr::filter(algorithm == "MIXED") %>%
  droplevels() %>%
  dplyr::mutate(designation = case_when(Simulated == TRUE & Detected == TRUE & aboveBF == TRUE ~ "Detected.CV",
                                        Simulated == TRUE & Detected == FALSE & aboveBF == FALSE ~ "Missed.CV",
                                        Simulated == TRUE & Detected == TRUE & aboveBF == FALSE ~ "CV.Not.Significant.In.Interval",
                                        Simulated == FALSE & Detected == TRUE & aboveBF == TRUE ~ "False.Discovery"),
                sim.loc = as.factor(sim.loc),
                strains = case_when(strains == "swept" ~ "Swept",
                                    strains == "unswept" ~ "Divergent")) %>%
  tidyr::separate(col = detected.peak,
                  into = c("peak.CHROM","peak.POS"), 
                  sep = ":", remove = F) %>%
  dplyr::group_by(h2, strains, sim.loc, Simulated.Divergence, Simulated.Region, batch, designation, Rep) %>%
  dplyr::summarise(n = n()) %>%
  tidyr::pivot_wider(names_from = designation, values_from = n)
designations[is.na(designations)] <- 0
levels(designations$sim.loc) <- c("I","II","III","IV","V","X")

power <- designations %>%
  dplyr::mutate(Detected = Detected.CV + CV.Not.Significant.In.Interval + False.Discovery, 
                Simulated = Detected.CV + CV.Not.Significant.In.Interval + Missed.CV,
                Power = Detected.CV/Simulated, 
                FDR = False.Discovery/Detected) %>%
  dplyr::group_by(h2, strains, sim.loc, Simulated.Divergence, Simulated.Region) %>% # ditch replicate and batch; individual mappings have been evaluated
  dplyr::summarise(mean.Power = mean(Power),
                   sd.Power = sd(Power)) %>%
  dplyr::mutate(a = 1 - mean.Power) %>%
  dplyr::mutate(sd.top.Power = if_else(a < sd.Power, true = a, false = sd.Power)) %>%
  dplyr::mutate(strainset.chrom = paste(strains, sim.loc, sep = "_"))
levels(power$Simulated.Region) <- c("Chromosome Arm","Chromosome Center")
levels(power$Simulated.Divergence) <- c("Hyperdivergent Regions","All Other Loci")


FDR <- designations %>%
  dplyr::mutate(Detected = Detected.CV + CV.Not.Significant.In.Interval + False.Discovery, 
                Simulated = Detected.CV + CV.Not.Significant.In.Interval + Missed.CV,
                Power = Detected.CV/Simulated, 
                FDR = False.Discovery/Detected) %>%
  dplyr::filter(!is.na(FDR)) %>%
  dplyr::group_by(h2, strains, sim.loc, Simulated.Divergence, Simulated.Region) %>% # ditch replicate and batch; individual mappings have been evaluated
  dplyr::summarise(mean.FDR = mean(FDR),
                   sd.FDR = sd(FDR)) %>%
  dplyr::mutate(strainset.chrom = paste(strains, sim.loc, sep = "_"))
levels(FDR$Simulated.Region) <- c("Chromosome Arm","Chromosome Center")
levels(FDR$Simulated.Divergence) <- c("Hyperdivergent Regions","All Other Loci")


################
### Figure 5 ###
################
combined.metrics <- power %>%
  dplyr::full_join(., FDR)
combined.performance.figure <- ggplot(combined.metrics) + 
  theme_bw() +
  # FDR
  geom_line(mapping = aes(x = h2, group = interaction(strains,sim.loc), y = mean.FDR, colour = strainset.chrom),
            position = position_dodge(width = 0.3), linetype = 2) +
  geom_pointrange(mapping = aes(x = h2, group = interaction(strains,sim.loc), y = mean.FDR, 
                                ymin = mean.FDR - sd.FDR, ymax = mean.FDR + sd.FDR, 
                                fill = strainset.chrom, colour = strainset.chrom), 
                  position = position_dodge(width = 0.3),shape = 23, size = 0.3) +
  # power
  geom_line(mapping = aes(x = h2, group = interaction(strains,sim.loc), y = mean.Power, colour = strainset.chrom), 
            position = position_dodge(width = 0.3), linetype = 1) +
  geom_pointrange(mapping = aes(x = h2, group = interaction(strains,sim.loc), y = mean.Power, 
                                ymin = mean.Power - sd.Power, ymax = mean.Power + sd.top.Power, 
                                fill = strainset.chrom, colour = strainset.chrom), 
                  position = position_dodge(width = 0.3), shape = 21) +
  facet_grid(Simulated.Divergence ~ Simulated.Region) +
  scale_colour_manual(values = c(brewer.pal(9, "Blues")[4:9], brewer.pal(9, name = "Oranges")[4:9]), 
                      name = "Population") + 
  scale_fill_manual(values = c(brewer.pal(9, "Blues")[4:9], brewer.pal(9, name = "Oranges")[4:9]), 
                    name = "Population") +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.title.y = element_blank()) + 
  labs(x = expression(italic(h^2)))


gg.divergent.tree <- as.ggplot(function() plot(divergent.tree, 
                                               type = "u", 
                                               edge.color = "blue", 
                                               show.tip.label = F, edge.width = 0.5, no.margin = TRUE))

gg.swept.tree <- as.ggplot(function() plot(swept.tree, 
                                           type = "u", 
                                           edge.color = "orange", 
                                           show.tip.label = F, edge.width = 0.5, no.margin = TRUE))

stacked.trees <- cowplot::plot_grid(gg.divergent.tree + 
                                      theme(plot.background = element_rect(fill = "white",colour = NA)), 
                                    gg.swept.tree + 
                                      theme(plot.background = element_rect(fill = "white",colour = NA)), 
                                    ncol = 1, labels = c("B","C"))
fig5 <- cowplot::plot_grid(combined.performance.figure,stacked.trees, ncol = 2, labels = c("A",""), rel_widths = c(1.8,1))
ggsave(fig5 + theme(plot.background = element_rect(fill = "white",colour = NA)), filename = "plots/figure.5.png", width = 7.5, height = 4)


##############
### Stats  ###
##############
overall <- designations %>%
  dplyr::mutate(Detected = Detected.CV + CV.Not.Significant.In.Interval + False.Discovery, 
                Simulated = Detected.CV + CV.Not.Significant.In.Interval + Missed.CV,
                Power = Detected.CV/Simulated, 
                FDR = False.Discovery/Detected)
levels(overall$Simulated.Region) <- c("Arm","Center")
levels(overall$Simulated.Divergence) <- c("Hyperdivergent","Non-Hyperdivergent")

# Comparing divergent and swept strain sets among all other factors
tests <- overall %>%
  dplyr::mutate(all.factors = strains) %>%
  dplyr::group_by(h2, Simulated.Divergence, Simulated.Region) %>%
  tidyr::nest()
for(i in 1:length(tests$data)){
  print(paste(as.character(tests$h2[[i]]), 
              as.character(tests$Simulated.Divergence[[i]]),
              as.character(tests$Simulated.Region[[i]]), sep = "; "))
  KW <- tests$data[[i]] %>%
    data.frame()
  
  if(length(unique(KW$Power)) == 1){
    print("No variance in power among mappings")
  } else {
    print(kruskal.test(Power ~ strains, data = KW))
    
  }
  
}
for(i in 1:length(tests$data)){
  print(paste(as.character(tests$h2[[i]]), 
              as.character(tests$Simulated.Divergence[[i]]),
              as.character(tests$Simulated.Region[[i]]), sep = "; "))
  KW <- tests$data[[i]] %>%
    data.frame()
  
  if(length(unique(KW$FDR)) == 1){
    print("No variance in power among mappings")
  } else {
    print(kruskal.test(FDR ~ strains, data = KW))
    
  }
  
}


#############################
### Supplemental Table 6  ###
#############################
tests <- overall %>%
  dplyr::mutate(all.factors = strains) %>%
  dplyr::group_by(sim.loc, h2, Simulated.Divergence, Simulated.Region) %>%
  tidyr::nest()

KW.power.results <- list()
for(i in 1:length(tests$data)){
  # print(paste(as.character(tests$h2[[i]]),
  #             as.character(tests$Simulated.Divergence[[i]]),
  #             as.character(tests$Simulated.Region[[i]]),
  #             as.character(tests$sim.loc[[i]]), sep = "; "))
  KW <- tests$data[[i]] %>%
    data.frame()
  
  if(length(unique(KW$strains)) == 1){
    print("Only one strain set among compared groups")
    
  } else if(length(unique(KW$Power)) == 1){
    print("No variance in power among mappings")
    jw <- kruskal.test(Power ~ strains, data = KW)
    
    jw.2 <- data.frame(as.character(tests$h2[[i]]), 
                       as.character(tests$Simulated.Divergence[[i]]),
                       as.character(tests$Simulated.Region[[i]]),
                       as.character(tests$sim.loc[[i]]), jw$p.value)
    colnames(jw.2) <- c("h2","Simulated.Divergence","Simulated.Region","sim.loc","KW.p.POWER")
    KW.power.results[[i]] <- jw.2
    
    
  } else {
    
    jw <- kruskal.test(Power ~ strains, data = KW)
    
    jw.2 <- data.frame(as.character(tests$h2[[i]]), 
                       as.character(tests$Simulated.Divergence[[i]]),
                       as.character(tests$Simulated.Region[[i]]),
                       as.character(tests$sim.loc[[i]]), jw$p.value)
    colnames(jw.2) <- c("h2","Simulated.Divergence","Simulated.Region","sim.loc","KW.p.POWER")
    print(dunn_test(data = KW, Power ~ strains))
    KW.power.results[[i]] <- jw.2
  }
  
}
KW.FDR.results <- list()
for(i in 1:length(tests$data)){
  # print(paste(as.character(tests$h2[[i]]),
  #             as.character(tests$Simulated.Divergence[[i]]),
  #             as.character(tests$Simulated.Region[[i]]),
  #             as.character(tests$sim.loc[[i]]), sep = "; "))
  KW <- tests$data[[i]] %>%
    data.frame()
  
  if(length(unique(KW$strains)) == 1){
    print("Only one strain set among compared groups")
    
  } else if(length(unique(KW$FDR)) == 1){
    print("No variance in FDR among mappings")
    jw <- kruskal.test(FDR ~ strains, data = KW)
    
    jw.2 <- data.frame(as.character(tests$h2[[i]]), 
                       as.character(tests$Simulated.Divergence[[i]]),
                       as.character(tests$Simulated.Region[[i]]),
                       as.character(tests$sim.loc[[i]]), jw$p.value)
    colnames(jw.2) <- c("h2","Simulated.Divergence","Simulated.Region","sim.loc","KW.p.FDR")
    KW.FDR.results[[i]] <- jw.2
    
    
  } else {
    
    jw <- kruskal.test(FDR ~ strains, data = KW)
    
    jw.2 <- data.frame(as.character(tests$h2[[i]]), 
                       as.character(tests$Simulated.Divergence[[i]]),
                       as.character(tests$Simulated.Region[[i]]),
                       as.character(tests$sim.loc[[i]]), jw$p.value)
    colnames(jw.2) <- c("h2","Simulated.Divergence","Simulated.Region","sim.loc","KW.p.FDR")
    print(dunn_test(data = KW, FDR ~ strains))
    KW.FDR.results[[i]] <- jw.2
  }
  
}

strain.effect.padj.power <- KW.power.results %>%
  dplyr::bind_rows() %>%
  dplyr::mutate(sig = if_else(KW.p.POWER < 0.05, true = "SIG" , false = "not-SIG"),
                KW.p.POWER = round(KW.p.POWER, digits = 5)) %>%
  dplyr::select(-sig) %>%
  tidyr::pivot_wider(names_from = c(sim.loc), values_from = KW.p.POWER) %>%
  dplyr::arrange(Simulated.Region, Simulated.Divergence) %>%
  data.frame()
strain.effect.padj.power[is.na(strain.effect.padj.power)] <- ""
colnames(strain.effect.padj.power)[1:3] <- c("Trait Heritability","Simulated Divergence","Simulated Chromosome Region")
write.csv(strain.effect.padj.power, file = "tables/supplemental.table.6.csv")


#############################
### Supplemental Table 7  ###
#############################
tests <- overall %>%
  dplyr::mutate(all.factors = paste(Simulated.Divergence, Simulated.Region, sep = "_")) %>%
  dplyr::group_by(h2, strains) %>%
  tidyr::nest()

KW.power.results <- list()
KW.power.dunn.tests <- list()
for(i in 1:length(tests$data)){
  # print(paste(as.character(tests$h2[[i]]),
  #             as.character(tests$strains[[i]]), sep = "; "))
  KW <- tests$data[[i]] %>%
    data.frame()
  if(length(unique(KW$Power)) == 1){
    print("No variance in power among mappings")
    jw <- kruskal.test(Power ~ all.factors, data = KW)
    jw.2 <- data.frame(as.character(tests$h2[[i]]), 
                       as.character(tests$strains[[i]]),
                       jw$p.value)
    colnames(jw.2) <- c("h2","strains","KW.p.POWER")
    
    KW.power.results[[i]] <- jw.2
  } else {
    jw <- kruskal.test(Power ~ all.factors, data = KW)
    jw.2 <- data.frame(as.character(tests$h2[[i]]), 
                       as.character(tests$strains[[i]]),
                       jw$p.value)
    colnames(jw.2) <- c("h2","strains","KW.p.POWER")
    d <-dunn_test(data = KW, Power ~ all.factors)
    dunn <- d %>%
      # dplyr::filter(p.adj < 0.05) %>%
      dplyr::select(-c(.y.,statistic,p)) %>%
      dplyr::mutate(strain.set = tests$strains[i],
                    h2 = tests$h2[i])
    print(dunn)
    
    KW.power.dunn.tests[[i]] <- dunn
    
    KW.power.results[[i]] <- jw.2
  }
}

KW.FDR.results <- list()
KW.FDR.dunn.tests <- list()
for(i in 1:length(tests$data)){
  # print(paste(as.character(tests$h2[[i]]),
  #             as.character(tests$strains[[i]]), sep = "; "))
  KW <- tests$data[[i]] %>%
    data.frame()
  if(length(unique(KW$FDR)) == 1){
    print("No variance in FDR among mappings")
    jw <- kruskal.test(FDR ~ all.factors, data = KW)
    jw.2 <- data.frame(as.character(tests$h2[[i]]), 
                       as.character(tests$strains[[i]]),
                       jw$p.value)
    colnames(jw.2) <- c("h2","strains","KW.p.FDR")
    
    KW.FDR.results[[i]] <- jw.2
  } else {
    jw <- kruskal.test(FDR ~ all.factors, data = KW)
    jw.2 <- data.frame(as.character(tests$h2[[i]]), 
                       as.character(tests$strains[[i]]),
                       jw$p.value)
    colnames(jw.2) <- c("h2","strains","KW.p.FDR")
    d <-dunn_test(data = KW, FDR ~ all.factors)
    dunn <- d %>%
      # dplyr::filter(p.adj < 0.05) %>%
      dplyr::select(-c(.y.,statistic,p)) %>%
      dplyr::mutate(strain.set = tests$strains[i],
                    h2 = tests$h2[i])
    print(dunn)
    
    KW.FDR.dunn.tests[[i]] <- dunn
    KW.FDR.results[[i]] <- jw.2
  }
}

options(scipen = 999999)
KW.chr_region.divergence.effect.power <- KW.power.results %>%
  dplyr::bind_rows() %>%
  dplyr::mutate(sig = if_else(KW.p.POWER < 0.05, true = "SIG" , false = "not-SIG")) %>%
  dplyr::select(-sig) %>%
  tidyr::pivot_wider(names_from = strains, values_from = KW.p.POWER) %>%
  data.frame()
KW.chr_region.divergence.effect.power

chr_regions.divergence.posthoc.padj.power <- KW.power.dunn.tests %>%
  dplyr::bind_rows() %>%
  dplyr::select(h2, strain.set, everything()) %>%
  dplyr::mutate(p.adj = round(p.adj, digits = 5))
colnames(chr_regions.divergence.posthoc.padj.power) <- c("Trait Heritability","Strain Set","Group 1","Group 2","n1","n2","Adjusted p-value","Significance")
write.csv(chr_regions.divergence.posthoc.padj.power, file = "tables/supplemental.table.7.csv")


#############################
### Supplemental Table 8  ###
#############################
KW.chr_region.divergence.effect.FDR <- KW.FDR.results %>%
  dplyr::bind_rows() %>%
  dplyr::mutate(sig = if_else(KW.p.FDR < 0.05, true = "SIG" , false = "not-SIG")) %>%
  dplyr::select(-sig) %>%
  tidyr::pivot_wider(names_from = strains, values_from = KW.p.FDR) %>%
  data.frame()
KW.chr_region.divergence.effect.FDR

chr_regions.divergence.posthoc.padj.FDR <- KW.FDR.dunn.tests %>%
  dplyr::bind_rows() %>%
  dplyr::select(h2, strain.set, everything()) %>%
  dplyr::mutate(p.adj = round(p.adj, digits = 5))
colnames(chr_regions.divergence.posthoc.padj.FDR) <- c("Trait Heritability","Strain Set","Group 1","Group 2","n1","n2","Adjusted p-value","Significance")
write.csv(chr_regions.divergence.posthoc.padj.FDR, file = "tables/supplemental.table.8.csv")


#############################
### Supplemental Table 9  ###
#############################
tests <- overall %>%
  dplyr::mutate(all.factors = sim.loc) %>%
  dplyr::group_by(h2,strains, Simulated.Divergence, Simulated.Region) %>%
  tidyr::nest()

KW.power.results <- list()
KW.power.dunn.tests <- list()
for(i in 1:length(tests$data)){
  KW <- tests$data[[i]] %>%
    data.frame()
  if(length(unique(KW$Power)) == 1){
    print("No variance in power among mappings")
    jw <- kruskal.test(Power ~ all.factors, data = KW)
    jw.2 <- data.frame(as.character(tests$h2[[i]]), 
                       as.character(tests$strains[[i]]), 
                       as.character(tests$Simulated.Divergence[[i]]),
                       as.character(tests$Simulated.Region[[i]]),
                       jw$p.value)
    colnames(jw.2) <- c("h2","strains","Simulated.Divergence","Simulated.Chromosomal.Region","KW.p.POWER")
    
    KW.power.results[[i]] <- jw.2
  } else {
    jw <- kruskal.test(Power ~ all.factors, data = KW)
    jw.2 <- data.frame(as.character(tests$h2[[i]]), 
                       as.character(tests$strains[[i]]), 
                       as.character(tests$Simulated.Divergence[[i]]),
                       as.character(tests$Simulated.Region[[i]]),
                       jw$p.value)
    colnames(jw.2) <- c("h2","strains","Simulated.Divergence","Simulated.Chromosomal.Region","KW.p.POWER")
    d <-dunn_test(data = KW, Power ~ all.factors)
    dunn <- d %>%
      dplyr::filter(p.adj < 0.05) %>%
      dplyr::select(-c(.y.,statistic,p)) %>%
      dplyr::mutate(strain.set = tests$strains[i],
                    h2 = tests$h2[i],
                    Simulated.Divergence = tests$Simulated.Divergence[i],
                    Simulated.Region = tests$Simulated.Region[i])
    print(dunn)
    
    KW.power.dunn.tests[[i]] <- dunn
    
    KW.power.results[[i]] <- jw.2
  }
}

KW.FDR.results <- list()
KW.FDR.dunn.tests <- list()
for(i in 1:length(tests$data)){
  KW <- tests$data[[i]] %>%
    data.frame()
  if(length(unique(KW$FDR)) == 1){
    print("No variance in FDR among mappings")
    jw <- kruskal.test(FDR ~ all.factors, data = KW)
    jw.2 <- data.frame(as.character(tests$h2[[i]]), 
                       as.character(tests$strains[[i]]), 
                       as.character(tests$Simulated.Divergence[[i]]),
                       as.character(tests$Simulated.Region[[i]]),
                       jw$p.value)
    colnames(jw.2) <- c("h2","strains","Simulated.Divergence","Simulated.Chromosomal.Region","KW.p.FDR")
    
    KW.FDR.results[[i]] <- jw.2
  } else {
    jw <- kruskal.test(FDR ~ all.factors, data = KW)
    jw.2 <- data.frame(as.character(tests$h2[[i]]), 
                       as.character(tests$strains[[i]]), 
                       as.character(tests$Simulated.Divergence[[i]]),
                       as.character(tests$Simulated.Region[[i]]),
                       jw$p.value)
    colnames(jw.2) <- c("h2","strains","Simulated.Divergence","Simulated.Chromosomal.Region","KW.p.FDR")
    d <- dunn_test(data = KW, FDR ~ all.factors)
    dunn <- d %>%
      dplyr::filter(p.adj < 0.05) %>%
      dplyr::select(-c(.y.,statistic,p)) %>%
      dplyr::mutate(strain.set = tests$strains[i],
                    h2 = tests$h2[i],
                    Simulated.Divergence = tests$Simulated.Divergence[i],
                    Simulated.Region = tests$Simulated.Region[i])
    print(dunn)
    
    KW.FDR.dunn.tests[[i]] <- dunn
    
    KW.FDR.results[[i]] <- jw.2
  }
}

KW.chromosome.effect.power <- KW.power.results %>%
  dplyr::bind_rows() %>%
  dplyr::mutate(sig = if_else(KW.p.POWER < 0.05, true = "SIG" , false = "not-SIG")) %>%
  dplyr::select(-sig) %>%
  tidyr::pivot_wider(names_from = strains, values_from = KW.p.POWER) %>%
  dplyr::arrange(h2, Simulated.Divergence,Simulated.Chromosomal.Region) %>%
  data.frame()
KW.chromosome.effect.power

chromosome.effect.posthoc.padj.power <- KW.power.dunn.tests %>%
  dplyr::bind_rows() %>%
  dplyr::select(h2, strain.set, Simulated.Divergence, Simulated.Region, everything()) %>%
  dplyr::mutate(p.adj = round(p.adj, digits = 5)) %>%
  dplyr::arrange(h2, Simulated.Divergence,Simulated.Region)
colnames(chromosome.effect.posthoc.padj.power) <- c("Trait Heritability","Strain Set","Simulated Divergence","Simulated Chromosome Region",
                                                    "Group 1","Group 2","n1","n2","Adjusted p-value","Significance")
write.csv(chromosome.effect.posthoc.padj.power, file = "tables/supplemental.table.9.csv")


##############################
### Supplemental Table 10  ###
##############################
KW.chromosome.effect.FDR <- KW.FDR.results %>%
  dplyr::bind_rows() %>%
  dplyr::mutate(sig = if_else(KW.p.FDR < 0.05, true = "SIG" , false = "not-SIG")) %>%
  dplyr::select(-sig) %>%
  tidyr::pivot_wider(names_from = strains, values_from = KW.p.FDR) %>%
  dplyr::arrange(h2, Simulated.Divergence,Simulated.Chromosomal.Region) %>%
  data.frame()

chromosome.effect.posthoc.padj.FDR <- KW.FDR.dunn.tests %>%
  dplyr::bind_rows() %>%
  dplyr::select(h2, strain.set, Simulated.Divergence, Simulated.Region, everything())%>%
  dplyr::mutate(p.adj = round(p.adj, digits = 5)) %>%
  dplyr::arrange(h2, Simulated.Divergence, Simulated.Region)
colnames(chromosome.effect.posthoc.padj.FDR) <- c("Trait Heritability","Strain Set","Simulated Divergence","Simulated Chromosome Region",
                                                    "Group 1","Group 2","n1","n2","Adjusted p-value","Significance")
write.csv(chromosome.effect.posthoc.padj.FDR, file = "tables/supplemental.table.10.csv")
