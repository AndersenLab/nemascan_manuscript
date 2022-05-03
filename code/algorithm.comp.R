require(tidyverse)
require(cowplot)
require(rstatix)
setwd(paste0(dirname(rstudioapi::getActiveDocumentContext()$path),"/.."))

#################################
### Algorithm Simulation Data ###
#################################

load(file = "data/File_S1.RData") # EMMA, Algorithm Sims
load(file = "data/File_S2.RData") # Algorithm Genomic Inflation
dat <- simulation.metrics.df %>%
  tidyr::separate(col = sim,
                  into = c("nQTL","Rep","h2","MAF","effect_range","strain_set"), 
                  sep = "_", remove = F) %>%
  tidyr::separate(col = QTL,
                  into = c("CHROM","POS"), 
                  sep = ":", remove = F) %>%
  dplyr::mutate(h2 = as.factor(h2),
                nQTL = as.factor(nQTL),
                POS = as.numeric(POS),
                startPOS = as.numeric(startPOS),
                peakPOS = as.numeric(peakPOS),
                endPOS = as.numeric(endPOS),
                interval.var.exp  = as.numeric(interval.var.exp),
                Simulated.QTL.VarExp = as.numeric(Simulated.QTL.VarExp), 
                peak_id = as.numeric(peak_id),
                BETA = as.numeric(BETA),
                Effect = as.numeric(Effect),
                Frequency = as.numeric(Frequency),
                log10p = dplyr::if_else(Simulated == FALSE, 
                                        true = interval.log10p, 
                                        false = log10p), # false discoveries inherit the log10p value of the peak marker for the interval
                log10p = as.numeric(log10p),
                interval_size = as.numeric(interval_size),
                aboveBF = dplyr::case_when(aboveBF == 1 ~ TRUE, 
                                           aboveBF == 0 ~ FALSE,
                                           is.na(aboveBF) ~ TRUE), # false discoveries by definition exceed significance threshold
                aboveBF = as.factor(aboveBF)) %>%
  dplyr::filter(!c(Detected == FALSE & aboveBF == TRUE),
                CHROM != 7, 
                algorithm != "LMM-EXACT-INBRED-LOCO")
dat$algorithm <- as.factor(dat$algorithm)
levels(dat$algorithm) <- c("EMMA",
                           "fastGWA-lmm-exact",
                           "fastGWA-lmm-exact-INBRED",
                           "fastGWA-lmm-exact-LOCO")



#############################
### Supplemental Figure 1 ###
#############################
supp.figure.2 <- dat %>%
  dplyr::filter(Simulated == TRUE) %>%
  dplyr::mutate(nQTL = paste0(nQTL," QTL")) %>%
  ggplot(., mapping = aes(x = abs(Effect))) +
  theme_bw() + 
  geom_histogram() + 
  facet_grid(nQTL~., scales = "free") + 
  theme(legend.position = "none")+
  labs(x = "Assigned QTL Effect Magnitude",
       y = "Frequency of Causal Variants")
ggsave(supp.figure.2, filename = "plots/supp.fig.2.png", height = 6, width = 6)


################
### Figure 1 ###
################
# Power and FDR calculation
designations <- dat %>%
  dplyr::mutate(designation = case_when(Simulated == TRUE & Detected == TRUE & aboveBF == TRUE ~ "Detected.CV",
                                        Simulated == TRUE & Detected == FALSE & aboveBF == FALSE ~ "Missed.CV",
                                        Simulated == TRUE & Detected == TRUE & aboveBF == FALSE ~ "CV.Not.Significant.In.Interval",
                                        Simulated == FALSE & Detected == TRUE & aboveBF == TRUE ~ "False.Discovery")) %>%
  tidyr::separate(col = detected.peak,
                  into = c("peak.CHROM","peak.POS"), 
                  sep = ":", remove = F) %>%
  dplyr::mutate(QTL.v.peak = abs(as.numeric(POS)-as.numeric(peak.POS))) %>%
  dplyr::group_by(designation, algorithm, nQTL, h2, Rep) %>%
  dplyr::summarise(n = n()) %>%
  tidyr::pivot_wider(names_from = designation, values_from = n)
designations[is.na(designations)] <- 0
designations$Detected <- rowSums(designations[,5:7])


Power <- designations %>%
  dplyr::mutate(Simulated = as.numeric(as.character(nQTL)),
                Power = Detected.CV/Simulated,
                Artefact.Rate = False.Discovery/Detected,
                Detected.CV.NS.Rate = CV.Not.Significant.In.Interval/Detected) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(h2, nQTL, algorithm) %>%
  dplyr::summarise(mean.Power = mean(Power),
                   sd.Power = sd(Power),
                   n = n()) %>%
  dplyr::mutate(a = 1 - mean.Power) %>%
  dplyr::mutate(sd.top.Power = if_else(a < sd.Power, true = a, false = sd.Power)) %>%
  dplyr::mutate(sd.bottom.Power = if_else(mean.Power-sd.Power < 0, true = (0-mean.Power)*-1, false = sd.Power))
Power$h2 <- as.numeric(as.character(Power$h2))


algorithm.palette <- c("#F94144", 
                       "#F9C74F",
                       "#43AA8B",
                       "#277DA1", 
                       "grey50", # EMMA
                       "#000000") 
# names(algorithm.palette) <- unique(summarizd.all.gifs$ALGORITHM)
names(algorithm.palette) <- c("fastGWA-lmm-exact-LOCO","fastGWA-lmm-exact-LOCO-PCA",
                              "fastGWA-lmm-exact-INBRED","fastGWA-lmm-exact-INBRED-PCA",
                              "EMMA","fastGWA-lmm-exact")
algorithm.A <- Power %>%
  dplyr::mutate(nQTL = paste0(nQTL, " QTL")) %>%
  ggplot(., mapping = aes(x = h2, y = mean.Power, colour = algorithm, group = algorithm )) + 
  theme_bw(base_size = 10) + 
  geom_line(position = position_dodge(width = 0.03), size = 0.25) +
  geom_pointrange(aes(ymin=mean.Power-sd.bottom.Power, ymax=mean.Power+sd.top.Power), position = position_dodge(width = 0.03), size = 0.25) +
  # geom_point(position = position_dodge(width = 0.03)) +
  
  facet_grid(~nQTL, scales = "free", space = "free") + 
  scale_colour_manual(values = algorithm.palette, name = "Algorithm") + 
  guides(colour = guide_legend(ncol = 1, title.position = "top", title.hjust = 0.5)) + 
  # scale_y_continuous(breaks = seq(0,1,0.25), labels = seq(0,1,0.25), limits = c(0,1)) + 
  theme(strip.text = element_text(size = 8),
        legend.position = "right", 
        panel.grid = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8)) +
  labs(y = "Power",
       x = expression(italic(h^2)))


AF <- designations %>%
  dplyr::mutate(Simulated = as.numeric(as.character(nQTL)),
                Power = Detected.CV/Simulated,
                Artefact.Rate = if_else(is.na(False.Discovery/Detected),
                                        true = 0,
                                        false = (False.Discovery/Detected)),
                Detected.CV.NS.Rate = if_else(is.na(CV.Not.Significant.In.Interval/Detected),
                                              true = 0,
                                              false = (CV.Not.Significant.In.Interval/Detected))) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(h2, nQTL, algorithm) %>%
  dplyr::summarise(mean.AR = mean(Artefact.Rate),
                   sd.AR = sd(Artefact.Rate),
                   n = n()) %>%
  dplyr::mutate(a = 1 - mean.AR) %>%
  dplyr::mutate(sd.top.AR = if_else(a < sd.AR, true = a, false = sd.AR)) %>%
  dplyr::mutate(sd.bottom.AR = if_else(mean.AR-sd.AR < 0, true = (0-mean.AR)*-1, false = sd.AR))
AF$h2 <- as.numeric(as.character(AF$h2))
algorithm.B <- AF %>%
  dplyr::mutate(nQTL = paste0(nQTL, " QTL")) %>%
  ggplot(., mapping = aes(x = h2, y = mean.AR, colour = algorithm, group = algorithm )) + 
  theme_bw(base_size = 10) + 
  geom_line(position = position_dodge(width = 0.03), size = 0.25) +
  geom_pointrange(aes(ymin=mean.AR-sd.bottom.AR, ymax=mean.AR+sd.top.AR), position = position_dodge(width = 0.03), size = 0.25) +
  # geom_point(position = position_dodge(width = 0.03)) +
  
  facet_grid(~nQTL, scales = "free", space = "free") + 
  scale_colour_manual(values = algorithm.palette, name = "Algorithm") + 
  guides(colour = guide_legend(nrow = 2)) + 
  # scale_y_continuous(breaks = seq(0,1,0.25), labels = seq(0,1,0.25), limits = c(0,1)) +
  theme(strip.text = element_text(size = 8),
        legend.position = "top", 
        panel.grid = element_blank()) +
  labs(y = "FDR",
       x = expression(italic(h^2)))


all.gifs.df <- all.gifs %>%
  tidyr::separate(col = SIM,
                  into = c("nQTL","Rep","h2","MAF","effect_range","strain_set"), 
                  sep = "_", remove = F)
summarizd.all.gifs <- all.gifs.df %>%
  dplyr::group_by(nQTL, h2, ALGORITHM) %>%
  dplyr::summarise(mean.gif = mean(GIF),
                   sd.gif = sd(GIF),
                   n = n()) %>%
  dplyr::mutate(ALGORITHM = gsub(ALGORITHM, pattern = "_", replacement = "-"),
                ALGORITHM = paste0("LMM-EXACT-", ALGORITHM))
summarizd.all.gifs$ALGORITHM <- as.factor(summarizd.all.gifs$ALGORITHM)
levels(summarizd.all.gifs$ALGORITHM) <- c("fastGWA-lmm-exact-INBRED",
                                          "fastGWA-lmm-exact-INBRED-PCA",
                                          "fastGWA-lmm-exact-LOCO",
                                          "fastGWA-lmm-exact-LOCO-PCA")

algorithm.C <- ggplot() + 
  theme_bw(base_size = 10) + 
  geom_hline(yintercept = 1, linetype = 3) + 
  geom_pointrange(data = summarizd.all.gifs, mapping = aes(x = h2, y = mean.gif,
                                                           ymin = mean.gif-sd.gif, 
                                                           ymax = mean.gif+sd.gif,
                                                           colour = ALGORITHM, 
                                                           fill = ALGORITHM),
                  position = position_dodge(width = 0.3), shape = 21) + 
  scale_color_manual(values = algorithm.palette, name = "Algorithm") +
  scale_fill_manual(values = algorithm.palette, name = "Algorithm") + 
  theme(panel.grid = element_blank()) + 
  labs(x = expression(italic(h^2)),
       y = expression(λ[GC]))
pre.plots <- cowplot::plot_grid(algorithm.A + theme(legend.position = "none",
                                                plot.background = element_rect(fill = "white",colour = NA)),
                            algorithm.B + theme(legend.position = "none",
                                                plot.background = element_rect(fill = "white",colour = NA)),
                            labels = "AUTO",
                            ncol = 1, rel_heights = c(0.9,1))
legends <- cowplot::get_legend(algorithm.A)
pre.plots.2 <- cowplot::plot_grid(algorithm.C + theme(legend.position = "none",
                                                plot.background = element_rect(fill = "white",colour = NA)),
                                  legends,
                                  labels = c("C",""),
                                  ncol = 1, rel_heights = c(1,0.8))

figure.1 <- cowplot::plot_grid(pre.plots, pre.plots.2, ncol = 2)
figure.1
ggsave(figure.1 + theme(plot.background = element_rect(fill = "white",colour = NA)), filename = "plots/figure.1.png", height = 4, width = 7)


#############################
### Supplemental Table 2  ###
#############################
power.FDR.stats <- designations %>%
  dplyr::mutate(Detected = Detected.CV + CV.Not.Significant.In.Interval + False.Discovery, 
                Simulated = Detected.CV + CV.Not.Significant.In.Interval + Missed.CV,
                Power = Detected.CV/Simulated, 
                FDR = False.Discovery/Detected) %>%
  dplyr::group_by(h2, nQTL) %>%
  tidyr::nest()

KW.power.results <- list()
KW.power.dunn.tests <- list()
for(i in 1:length(power.FDR.stats$data)){
  
  KW <- power.FDR.stats$data[[i]] %>%
    data.frame()
  if(length(unique(KW$Power)) == 1){
    print("No variance in power among mappings")
    jw <- kruskal.test(Power ~ algorithm, data = KW)
    jw.2 <- data.frame(as.character(power.FDR.stats$h2[[i]]),
                       as.character(power.FDR.stats$nQTL[[i]]),
                       jw$p.value)
    colnames(jw.2) <- c("h2","nQTL","KW.p.POWER")
    
    KW.power.results[[i]] <- jw.2
  } else {
    jw <- kruskal.test(Power ~ algorithm, data = KW)
    jw.2 <- data.frame(as.character(power.FDR.stats$h2[[i]]),
                       as.character(power.FDR.stats$nQTL[[i]]),
                       jw$p.value)
    colnames(jw.2) <- c("h2","nQTL","KW.p.POWER")
    d <-rstatix::dunn_test(data = KW, Power ~ algorithm)
    dunn <- d %>%
      # dplyr::filter(p.adj < 0.05) %>%
      dplyr::select(-c(.y.,statistic,p)) %>%
      dplyr::mutate(h2 = power.FDR.stats$h2[i],
                    nQTL = power.FDR.stats$nQTL[i])
    print(data.frame(dunn))
    
    KW.power.dunn.tests[[i]] <- dunn
    KW.power.results[[i]] <- jw.2
  }
}

KW.power.results %>%
  dplyr::bind_rows() %>%
  dplyr::mutate(sig = if_else(KW.p.POWER < 0.05, true = "SIG" , false = "not-SIG"),
                KW.p.POWER = round(KW.p.POWER, digits = 5)) %>%
  dplyr::select(-sig) %>%
  dplyr::arrange(nQTL, h2)

algorithm.posthoc.tests.power <- KW.power.dunn.tests %>%
  dplyr::bind_rows() %>%
  dplyr::mutate(p.adj = round(p.adj, digits = 5)) %>%
  dplyr::arrange(nQTL, h2) %>%
  dplyr::select(h2, nQTL, everything())
colnames(algorithm.posthoc.tests.power) <- c("Trait Heritability","Number of Simulated QTL","Group 1","Group 2","n1","n2","Adjusted p-value","Significance")
write.csv(algorithm.posthoc.tests.power, "tables/supplemental.table.2.csv", row.names = F, quote = F)


#############################
### Supplemental Table 3  ###
#############################
KW.FDR.results <- list()
KW.FDR.dunn.tests <- list()
for(i in 1:length(power.FDR.stats$data)){
  
  KW <- power.FDR.stats$data[[i]] %>%
    data.frame()
  if(length(unique(KW$FDR)) == 1){
    print("No variance in power among mappings")
    jw <- kruskal.test(FDR ~ algorithm, data = KW)
    jw.2 <- data.frame(as.character(power.FDR.stats$h2[[i]]),
                       as.character(power.FDR.stats$nQTL[[i]]),
                       jw$p.value)
    colnames(jw.2) <- c("h2","nQTL","KW.p.FDR")
    
    KW.FDR.results[[i]] <- jw.2
  } else {
    jw <- kruskal.test(FDR ~ algorithm, data = KW)
    jw.2 <- data.frame(as.character(power.FDR.stats$h2[[i]]),
                       as.character(power.FDR.stats$nQTL[[i]]),
                       jw$p.value)
    colnames(jw.2) <- c("h2","nQTL","KW.p.FDR")
    d <-rstatix::dunn_test(data = KW, FDR ~ algorithm)
    dunn <- d %>%
      # dplyr::filter(p.adj < 0.05) %>%
      dplyr::select(-c(.y.,statistic,p)) %>%
      dplyr::mutate(h2 = power.FDR.stats$h2[i],
                    nQTL = power.FDR.stats$nQTL[i])
    print(data.frame(dunn))
    
    KW.FDR.dunn.tests[[i]] <- dunn
    KW.FDR.results[[i]] <- jw.2
  }
}


KW.FDR.results %>%
  dplyr::bind_rows() %>%
  dplyr::mutate(sig = if_else(KW.p.FDR < 0.05, true = "SIG" , false = "not-SIG"),
                KW.p.FDR = round(KW.p.FDR, digits = 5)) %>%
  dplyr::select(-sig)%>%
  dplyr::arrange(nQTL, h2)

algorithm.posthoc.tests.FDR <-KW.FDR.dunn.tests %>%
  dplyr::bind_rows() %>%
  dplyr::mutate(p.adj = round(p.adj, digits = 5)) %>%
  dplyr::arrange(nQTL, h2) %>%
  dplyr::select(h2, nQTL, everything())
colnames(algorithm.posthoc.tests.FDR) <- c("Trait Heritability","Number of Simulated QTL","Group 1","Group 2","n1","n2","Adjusted p-value","Significance")
write.csv(algorithm.posthoc.tests.FDR, "tables/supplemental.table.3.csv", row.names = F, quote = F)

## Genomic Inflation Stats ##
gifs.aov <- aov(all.gifs.df, formula = GIF ~ ALGORITHM)
TukeyHSD(gifs.aov)

