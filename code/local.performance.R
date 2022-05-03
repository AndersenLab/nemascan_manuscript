require(tidyverse)
require(valr)
require(RColorBrewer)
require(magrittr)
require(ggtree)
require(ggplotify)
require(ggh4x)
require(rstatix)
setwd(paste0(dirname(rstudioapi::getActiveDocumentContext()$path),"/.."))
arms.centers <- data.table::fread("data/ARMS_CENTERS.tsv") %>%
  `colnames<-`(c("chrom","start","end","region"))
#########################################################################################################
### Comparisons of Local Performance (Chromosome Regions, Divergent Regions, Swept/Divergent Strains) ###
#########################################################################################################
load(file = "data/File_S8.RData")
load(file = "data/File_S9.RData")
sim.nested <- all.performance %>%
  dplyr::group_by(strain_set, divergence, sim) %>%
  tidyr::nest()

assign.location.metadata <- function(data, x, y){
  
  simulated.chrom <- data %>%
    dplyr::filter(Simulated == TRUE) %>%
    dplyr::select(CHROM)
  
  if(nrow(simulated.chrom) == 0){
  print("Simulated QTL filtered from arms/centers")
  } else {
    for.arms.centers <- data %>%
      dplyr::mutate(simulated.CHROM = simulated.chrom$CHROM,
                    strain_set = x,
                    divergence = y) %>%
      dplyr::mutate(chrom = as.factor(CHROM),
                    start = POS,
                    end = POS)
    levels(for.arms.centers$chrom) <- c("I","II","III","IV","V","X")
    simulated.chrom.region <- for.arms.centers %>%
      dplyr::filter(Simulated == TRUE) %>%
      valr::bed_intersect(., arms.centers) %>% 
      dplyr::select(region.y) %>% as.character()
    
    for.arms.centers %>%
      dplyr::mutate(simulated.region = simulated.chrom.region,
                    simulated.region = if_else(simulated.region != "character(0)", 
                                               true = simulated.region, 
                                               false = "TIP"))
    
  }
}

dat.genomic.locations.chrom.recoded <- purrr::pmap(.l = list(sim.nested$data, 
                                                             sim.nested$strain_set, 
                                                             sim.nested$divergence),
                                                   .f = assign.location.metadata) %>%
  purrr::keep(., is.tibble) %>%
  Reduce(rbind,.)
dat.genomic.locations.chrom.recoded %<>%
  dplyr::mutate(simulated.CHROM = as.factor(simulated.CHROM),
                simulated.region = as.factor(simulated.region)) %>% 
  dplyr::select(-start, -end)
levels(dat.genomic.locations.chrom.recoded$simulated.CHROM) <- c("I","II","III","IV","V","X")
designations <- dat.genomic.locations.chrom.recoded %>%
  dplyr::mutate(designation = case_when(Simulated == TRUE & Detected == TRUE & aboveBF == TRUE ~ "Detected.CV",
                                        Simulated == TRUE & Detected == FALSE & aboveBF == FALSE ~ "Missed.CV",
                                        Simulated == TRUE & Detected == TRUE & aboveBF == FALSE ~ "CV.Not.Significant.In.Interval",
                                        Simulated == FALSE & Detected == TRUE & aboveBF == TRUE ~ "False.Discovery")) %>%
  tidyr::separate(col = detected.peak,
                  into = c("peak.CHROM","peak.POS"), 
                  sep = ":", remove = F) %>%
  dplyr::group_by(h2, strain_set, simulated.region, divergence, designation, Rep, simulated.CHROM) %>%
  dplyr::summarise(n = n()) %>%
  tidyr::pivot_wider(names_from = designation, values_from = n)
designations[is.na(designations)] <- 0





power <- designations %>%
  dplyr::mutate(Detected = Detected.CV + CV.Not.Significant.In.Interval + False.Discovery, 
                Simulated = Detected.CV + CV.Not.Significant.In.Interval + Missed.CV,
                Power = Detected.CV/Simulated, 
                FDR = False.Discovery/Detected) %>%
  dplyr::group_by(strain_set, simulated.CHROM, divergence, simulated.region) %>% # ditch replicate and batch; individual mappings have been evaluated
  dplyr::summarise(mean.Power = mean(Power),
                   sd.Power = sd(Power)) %>%
  dplyr::mutate(a = 1 - mean.Power) %>%
  dplyr::mutate(sd.top.Power = if_else(mean.Power + sd.Power > 1, true = a, false = sd.Power)) %>%
  dplyr::mutate(strainset.chrom = paste(strain_set, simulated.CHROM, sep = "_"))
levels(power$simulated.region) <- c("Chromosome Arm","Chromosome Center","Chromosome Tip")
levels(power$divergence) <- c("Hyper-divergent Regions","All Other Loci")


FDR <- designations %>%
  dplyr::mutate(Detected = Detected.CV + CV.Not.Significant.In.Interval + False.Discovery, 
                Simulated = Detected.CV + CV.Not.Significant.In.Interval + Missed.CV,
                Power = Detected.CV/Simulated, 
                FDR = False.Discovery/Detected) %>%
  dplyr::filter(!is.na(FDR)) %>%
  dplyr::group_by(strain_set, simulated.CHROM, divergence, simulated.region) %>%
  dplyr::summarise(mean.FDR = mean(FDR),
                   sd.FDR = sd(FDR)) %>%
  dplyr::mutate(b = 1 - mean.FDR) %>%
  dplyr::mutate(sd.top.FDR = if_else(mean.FDR + sd.FDR > 1, true = b, false = sd.FDR)) %>%
  dplyr::mutate(strainset.chrom = paste(strain_set, simulated.CHROM, sep = "_"))
levels(FDR$simulated.region) <- c("Chromosome Arm","Chromosome Center","Chromosome Tip")
levels(FDR$divergence) <- c("Hyper-divergent Regions","All Other Loci")


################
### Figure 5 ###
################
combined.metrics <- power %>%
  dplyr::full_join(., FDR) %>%
  dplyr::mutate(divergence = if_else(divergence == "Hyperdivergent", true = "Hyper-divergent Regions", false = divergence))

levels(combined.metrics$simulated.region) <- c("Arm","Center","Tip")
combined.metrics$simulated.region <- factor(combined.metrics$simulated.region, levels = c("Tip","Arm","Center"))

# V1 Figure
# combined.performance.figure <- ggplot(combined.metrics) + 
#   theme_bw() +
#   # FDR
#   geom_line(mapping = aes(x = h2, group = interaction(strain_set,simulated.CHROM), y = mean.FDR, colour = strainset.chrom),
#             position = position_dodge(width = 0.3), linetype = 2) +
#   geom_pointrange(mapping = aes(x = h2, group = interaction(strain_set,simulated.CHROM), y = mean.FDR, 
#                                 ymin = mean.FDR - sd.FDR, ymax = mean.FDR + sd.top.FDR, 
#                                 fill = strainset.chrom, colour = strainset.chrom), 
#                   position = position_dodge(width = 0.3),shape = 23, size = 0.3) +
#   # power
#   geom_line(mapping = aes(x = h2, group = interaction(strain_set,simulated.CHROM), y = mean.Power, colour = strainset.chrom), 
#             position = position_dodge(width = 0.3), linetype = 1) +
#   geom_pointrange(mapping = aes(x = h2, group = interaction(strain_set,simulated.CHROM), y = mean.Power, 
#                                 ymin = mean.Power - sd.Power, ymax = mean.Power + sd.top.Power, 
#                                 fill = strainset.chrom, colour = strainset.chrom), 
#                   position = position_dodge(width = 0.3), shape = 21) +
#   facet_grid(divergence ~ simulated.region) +
#   scale_colour_manual(values = c(brewer.pal(9, name = "Oranges")[4:9], brewer.pal(9, "Blues")[4:9]),
#                       name = "Population") + 
#   scale_fill_manual(values = c(brewer.pal(9, name = "Oranges")[4:9], brewer.pal(9, "Blues")[4:9]), 
#                     name = "Population") +
#   theme(legend.position = "none",
#         panel.grid = element_blank(),
#         axis.title.y = element_blank()) + 
#   labs(x = expression(italic(h^2)))


# combined.performance.figure.2 <- ggplot(combined.metrics)  + 
#   theme_bw() +
#   geom_pointrange(mapping = aes(x = h2, group = interaction(strain_set,simulated.CHROM), y = mean.FDR, 
#                                 ymin = mean.FDR - sd.FDR, ymax = mean.FDR + sd.top.FDR, 
#                                 colour = strain_set), 
#                   position = position_dodge(width = 0.3),shape = 23, size = 0.3) + 
#   geom_pointrange(mapping = aes(x = h2, group = interaction(strain_set,simulated.CHROM), y = mean.Power, 
#                                 ymin = mean.Power - sd.Power, ymax = mean.Power + sd.Power, 
#                                 fill = strain_set, colour = strain_set), 
#                   position = position_dodge(width = 0.5), shape = 21, size = 0.3) + 
#   scale_colour_manual(values = c("orange","blue"),
#                       name = "Population") + 
#   scale_fill_manual(values = c("orange","blue"), 
#                     name = "Population") + 
#   scale_y_continuous(breaks = seq(0,1,0.5)) +
#   facet_nested(divergence + simulated.region ~ simulated.CHROM) + 
#   theme(legend.position = "none",
#         panel.grid = element_blank(),
#         axis.title.y = element_blank(),
#         axis.text = element_text(size = 8),
#         strip.text = element_text(size = 8)) + 
#   labs(x = expression(italic(h^2)))

combined.performance.figure.3 <- ggplot(combined.metrics)  + 
  theme_bw() +
  geom_pointrange(mapping = aes(x = simulated.region, group = interaction(strain_set,simulated.CHROM), y = mean.FDR, 
                                ymin = mean.FDR - sd.FDR, ymax = mean.FDR + sd.top.FDR, 
                                colour = strain_set), 
                  position = position_dodge(width = 0.5),shape = 23, size = 0.3) + 
  geom_pointrange(mapping = aes(x = simulated.region, group = interaction(strain_set,simulated.CHROM), y = mean.Power, 
                                ymin = mean.Power - sd.Power, ymax = mean.Power + sd.top.Power, 
                                fill = strain_set, colour = strain_set), 
                  position = position_dodge(width = 0.7), shape = 21, size = 0.3) + 
  scale_colour_manual(values = c("orange","blue"),
                      name = "Population") + 
  scale_fill_manual(values = c("orange","blue"), 
                    name = "Population") + 
  scale_y_continuous(breaks = seq(0,1,0.5)) +
  facet_nested(divergence ~ simulated.CHROM) + 
  theme(legend.position = "none",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 8, angle = 45, hjust = 1, vjust = 1),
        strip.text = element_text(size = 8)) + 
  labs(x = "Chromosome Region")
combined.performance.figure.3




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
fig5 <- cowplot::plot_grid(combined.performance.figure.3,stacked.trees, ncol = 2, labels = c("A",""), rel_widths = c(2,1))
ggsave(fig5 + theme(plot.background = element_rect(fill = "white",colour = NA)), filename = "plots/figure.5.png", width = 7, height = 4)
ggsave(fig5 + theme(plot.background = element_rect(fill = "white",colour = NA)), filename = "plots/figure.5.pdf", width = 7, height = 4)


##############
### Stats  ###
##############
overall <- designations %>%
  dplyr::mutate(Detected = Detected.CV + CV.Not.Significant.In.Interval + False.Discovery, 
                Simulated = Detected.CV + CV.Not.Significant.In.Interval + Missed.CV,
                Power = Detected.CV/Simulated, 
                FDR = False.Discovery/Detected) %>%
  dplyr::mutate(divergence = as.factor(divergence))
levels(overall$simulated.region) <- c("Arm","Center","Tip")
levels(overall$divergence) <- c("Hyperdivergent","Non-Hyperdivergent")

# Comparing divergent and swept strain sets among all other factors
tests <- overall %>%
  dplyr::mutate(all.factors = strain_set) %>%
  dplyr::group_by(divergence, simulated.region) %>%
  tidyr::nest()
for(i in 1:length(tests$data)){
  print(paste(as.character(tests$divergence[[i]]),
              as.character(tests$simulated.region[[i]]), sep = "; "))
  KW <- tests$data[[i]] %>%
    data.frame()
  
  if(length(unique(KW$Power)) == 1){
    print("No variance in power among mappings")
  } else {
    print(kruskal.test(Power ~ strain_set, data = KW))
    
  }
  
}
for(i in 1:length(tests$data)){
  print(paste(as.character(tests$divergence[[i]]),
              as.character(tests$simulated.region[[i]]), sep = "; "))
  KW <- tests$data[[i]] %>%
    data.frame()
  
  if(length(unique(KW$FDR)) == 1){
    print("No variance in FDR among mappings")
  } else {
    print(kruskal.test(FDR ~ strain_set, data = KW))
    
  }
  
}


#############################
### Supplemental Table 6  ###
#############################
tests <- overall %>%
  dplyr::mutate(all.factors = strain_set) %>%
  dplyr::group_by(simulated.CHROM, divergence, simulated.region) %>%
  tidyr::nest()

KW.power.results <- list()
for(i in 1:length(tests$data)){
  KW <- tests$data[[i]] %>%
    data.frame()
  
  if(length(unique(KW$strain_set)) == 1){
    print("Only one strain set among compared groups")
    
  } else if(length(unique(KW$Power)) == 1){
    print("No variance in power among mappings")
    jw <- kruskal.test(Power ~ strain_set, data = KW)
    
    jw.2 <- data.frame(as.character(tests$divergence[[i]]),
                       as.character(tests$simulated.region[[i]]),
                       as.character(tests$simulated.CHROM[[i]]), jw$p.value)
    colnames(jw.2) <- c("Simulated.Divergence","Simulated.Region","Simulated.Chrom","KW.p.POWER")
    KW.power.results[[i]] <- jw.2
    
    
  } else {
    
    jw <- kruskal.test(Power ~ strain_set, data = KW)
    
    jw.2 <- data.frame(as.character(tests$divergence[[i]]),
                       as.character(tests$simulated.region[[i]]),
                       as.character(tests$simulated.CHROM[[i]]), jw$p.value)
    colnames(jw.2) <- c("Simulated.Divergence","Simulated.Region","Simulated.Chrom","KW.p.POWER")
    print(dunn_test(data = KW, Power ~ strain_set))
    KW.power.results[[i]] <- jw.2
  }
  
}

KW.FDR.results <- list()
for(i in 1:length(tests$data)){
  KW <- tests$data[[i]] %>%
    data.frame()
  
  if(length(unique(KW$strains)) == 1){
    print("Only one strain set among compared groups")
    
  } else if(length(unique(KW$FDR)) == 1){
    print("No variance in FDR among mappings")
    jw <- kruskal.test(FDR ~ strain_set, data = KW)
    
    jw.2 <- data.frame(as.character(tests$divergence[[i]]),
                       as.character(tests$simulated.region[[i]]),
                       as.character(tests$simulated.CHROM[[i]]), jw$p.value)
    colnames(jw.2) <- c("Simulated.Divergence","Simulated.Region","simulated.CHROM","KW.p.FDR")
    KW.FDR.results[[i]] <- jw.2
    
    
  } else {
    
    jw <- kruskal.test(FDR ~ strain_set, data = KW)
    
    jw.2 <- data.frame(as.character(tests$divergence[[i]]),
                       as.character(tests$simulated.region[[i]]),
                       as.character(tests$simulated.CHROM[[i]]), jw$p.value)
    colnames(jw.2) <- c("Simulated.Divergence","Simulated.Region","simulated.CHROM","KW.p.FDR")
    print(dunn_test(data = KW, FDR ~ strain_set))
    KW.FDR.results[[i]] <- jw.2
  }
  
}

strain.effect.padj.power <- KW.power.results %>%
  dplyr::bind_rows() %>%
  dplyr::mutate(sig = if_else(KW.p.POWER < 0.05, true = "SIG" , false = "not-SIG"),
                KW.p.POWER = round(KW.p.POWER, digits = 5)) %>%
  dplyr::select(-sig) %>%
  tidyr::pivot_wider(names_from = c(Simulated.Chrom), values_from = KW.p.POWER) %>%
  dplyr::select(Simulated.Divergence, Simulated.Region, I, II, III, IV, V, X) %>%
  dplyr::arrange(Simulated.Region, Simulated.Divergence) %>%
  data.frame()
strain.effect.padj.power[is.na(strain.effect.padj.power)] <- ""
colnames(strain.effect.padj.power)[1:2] <- c("Simulated Divergence","Simulated Chromosome Region")
write.csv(strain.effect.padj.power, file = "tables/supplemental.table.6.csv")


#############################
### Supplemental Table 7  ###
#############################
tests <- overall %>%
  dplyr::mutate(all.factors = paste(divergence, simulated.region, sep = "_")) %>%
  dplyr::group_by(strain_set) %>%
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
    jw.2 <- data.frame(as.character(tests$strain_set[[i]]),
                       jw$p.value)
    colnames(jw.2) <- c("strains","KW.p.POWER")
    
    KW.power.results[[i]] <- jw.2
  } else {
    jw <- kruskal.test(Power ~ all.factors, data = KW)
    jw.2 <- data.frame(as.character(tests$strain_set[[i]]),
                       jw$p.value)
    colnames(jw.2) <- c("strains","KW.p.POWER")
    d <-dunn_test(data = KW, Power ~ all.factors)
    dunn <- d %>%
      # dplyr::filter(p.adj < 0.05) %>%
      dplyr::select(-c(.y.,statistic,p)) %>%
      dplyr::mutate(strain.set = tests$strain_set[i])
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
    jw.2 <- data.frame(as.character(tests$strain_set[[i]]),
                       jw$p.value)
    colnames(jw.2) <- c("strains","KW.p.FDR")
    
    KW.FDR.results[[i]] <- jw.2
  } else {
    jw <- kruskal.test(FDR ~ all.factors, data = KW)
    jw.2 <- data.frame(as.character(tests$strain_set[[i]]),
                       jw$p.value)
    colnames(jw.2) <- c("strains","KW.p.FDR")
    d <-dunn_test(data = KW, FDR ~ all.factors)
    dunn <- d %>%
      # dplyr::filter(p.adj < 0.05) %>%
      dplyr::select(-c(.y.,statistic,p)) %>%
      dplyr::mutate(strain.set = tests$strain_set[i])
    print(dunn)
    
    KW.FDR.dunn.tests[[i]] <- dunn
    KW.FDR.results[[i]] <- jw.2
  }
}

options(scipen = 999999)
KW.chr_region.divergence.effect.power <- KW.power.results %>%
  dplyr::bind_rows() %>%
  dplyr::filter(!is.na(KW.p.POWER)) %>%
  dplyr::mutate(sig = if_else(KW.p.POWER < 0.05, true = "SIG" , false = "not-SIG")) %>%
  dplyr::select(-sig) %>%
  tidyr::pivot_wider(names_from = strains, values_from = KW.p.POWER) %>%
  data.frame()
KW.chr_region.divergence.effect.power

chr_regions.divergence.posthoc.padj.power <- KW.power.dunn.tests %>%
  dplyr::bind_rows() %>%
  dplyr::select(strain.set, everything()) %>%
  dplyr::mutate(p.adj = round(p.adj, digits = 5))
colnames(chr_regions.divergence.posthoc.padj.power) <- c("Strain Set","Group 1","Group 2","n1","n2","Adjusted p-value","Significance")
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
  dplyr::mutate(all.factors = simulated.CHROM) %>%
  dplyr::group_by(strain_set, divergence, simulated.region) %>%
  tidyr::nest()

KW.power.results <- list()
KW.power.dunn.tests <- list()
for(i in 1:length(tests$data)){
  KW <- tests$data[[i]] %>%
    data.frame()
  if(length(unique(KW$Power)) == 1){
    print("No variance in power among mappings")
    jw <- kruskal.test(Power ~ all.factors, data = KW)
    jw.2 <- data.frame(as.character(tests$strain_set[[i]]), 
                       as.character(tests$divergence[[i]]),
                       as.character(tests$simulated.region[[i]]),
                       jw$p.value)
    colnames(jw.2) <- c("strains","Simulated.Divergence","Simulated.Chromosomal.Region","KW.p.POWER")
    
    KW.power.results[[i]] <- jw.2
  } else {
    jw <- kruskal.test(Power ~ all.factors, data = KW)
    jw.2 <- data.frame(as.character(tests$strain_set[[i]]), 
                       as.character(tests$divergence[[i]]),
                       as.character(tests$simulated.region[[i]]),
                       jw$p.value)
    colnames(jw.2) <- c("strains","Simulated.Divergence","Simulated.Chromosomal.Region","KW.p.POWER")
    d <-dunn_test(data = KW, Power ~ all.factors)
    dunn <- d %>%
      dplyr::filter(p.adj < 0.05) %>%
      dplyr::select(-c(.y.,statistic,p)) %>%
      dplyr::mutate(strain.set = tests$strain_set[i],
                    divergence = tests$divergence[i],
                    Simulated.Region = tests$simulated.region[i])
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
    jw.2 <- data.frame(as.character(tests$strain_set[[i]]), 
                       as.character(tests$divergence[[i]]),
                       as.character(tests$simulated.region[[i]]),
                       jw$p.value)
    colnames(jw.2) <- c("strains","Simulated.Divergence","Simulated.Chromosomal.Region","KW.p.FDR")
    
    KW.FDR.results[[i]] <- jw.2
  } else {
    jw <- kruskal.test(FDR ~ all.factors, data = KW)
    jw.2 <- data.frame(as.character(tests$strain_set[[i]]), 
                       as.character(tests$divergence[[i]]),
                       as.character(tests$simulated.region[[i]]),
                       jw$p.value)
    colnames(jw.2) <- c("strains","Simulated.Divergence","Simulated.Chromosomal.Region","KW.p.FDR")
    d <- dunn_test(data = KW, FDR ~ all.factors)
    dunn <- d %>%
      dplyr::filter(p.adj < 0.05) %>%
      dplyr::select(-c(.y.,statistic,p)) %>%
      dplyr::mutate(strain.set = tests$strain_set[i],
                    divergence = tests$divergence[i],
                    Simulated.Region = tests$simulated.region[i])
    print(dunn)
    
    KW.FDR.dunn.tests[[i]] <- dunn
    
    KW.FDR.results[[i]] <- jw.2
  }
}

KW.chromosome.effect.power <- KW.power.results %>%
  dplyr::bind_rows() %>%
  dplyr::filter(!is.na(KW.p.POWER)) %>%
  dplyr::mutate(sig = if_else(KW.p.POWER < 0.05, true = "SIG" , false = "not-SIG")) %>%
  dplyr::select(-sig) %>%
  tidyr::pivot_wider(names_from = strains, values_from = KW.p.POWER) %>%
  dplyr::arrange(Simulated.Divergence, Simulated.Chromosomal.Region) %>%
  data.frame()
KW.chromosome.effect.power

chromosome.effect.posthoc.padj.power <- KW.power.dunn.tests %>%
  dplyr::bind_rows() %>%
  dplyr::select(strain.set, divergence, Simulated.Region, everything()) %>%
  dplyr::mutate(p.adj = round(p.adj, digits = 5)) %>%
  dplyr::arrange(divergence,Simulated.Region)
colnames(chromosome.effect.posthoc.padj.power) <- c("Strain Set","Simulated Divergence","Simulated Chromosome Region",
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
  dplyr::arrange(Simulated.Divergence,Simulated.Chromosomal.Region) %>%
  data.frame()

chromosome.effect.posthoc.padj.FDR <- KW.FDR.dunn.tests %>%
  dplyr::bind_rows() %>%
  dplyr::select(strain.set, divergence, Simulated.Region, everything())%>%
  dplyr::mutate(p.adj = round(p.adj, digits = 5)) %>%
  dplyr::arrange(divergence, Simulated.Region)
colnames(chromosome.effect.posthoc.padj.FDR) <- c("Strain Set","Simulated Divergence","Simulated Chromosome Region",
                                                  "Group 1","Group 2","n1","n2","Adjusted p-value","Significance")
write.csv(chromosome.effect.posthoc.padj.FDR, file = "tables/supplemental.table.10.csv")
