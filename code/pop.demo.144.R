require(tidyverse)
require(RColorBrewer)
require(magrittr)
setwd(paste0(dirname(rstudioapi::getActiveDocumentContext()$path),"/.."))

#####################################################################
### Subsampled, Swept, and Divergent Strain Comparisons @ n = 144 ###
#####################################################################
load(file = "data/File_S6.RData")
load(file = "data/File_S7.RData")

n144.df <- simulation.metrics.df %>%
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
                # log10p = dplyr::if_else(Simulated == FALSE,
                #                         true = interval.log10p,
                #                         false = log10p), # false discoveries inherit the log10p value of the peak marker for the interval
                # log10p = as.numeric(log10p),
                interval_size = as.numeric(interval_size),
                aboveBF = dplyr::case_when(aboveBF == 1 ~ TRUE, 
                                           aboveBF == 0 ~ FALSE,
                                           is.na(aboveBF) ~ TRUE), # false discoveries by definition exceed significance threshold
                aboveBF = as.factor(aboveBF)) %>%
  dplyr::filter(!c(Detected == FALSE & aboveBF == TRUE),
                CHROM != 7)

population.metrics.df[is.na(population.metrics.df)] <- 0
curated.population.metrics.df <- population.metrics.df %>%
  dplyr::select(-pct.population.swept) %>%
  dplyr::mutate(sweptI = I > 0.4,
                sweptIV = IV > 0.4,
                sweptV = V > 0.4,
                sweptX = X > 0.4) %>%
  tidyr::unite("sweptchroms",c(sweptI,sweptIV,sweptV,sweptX)) %>%
  dplyr::mutate(population.group = case_when(sweptchroms == "TRUE_TRUE_TRUE_TRUE" ~ "Swept",
                                             sweptchroms == "FALSE_FALSE_FALSE_FALSE" ~ "Divergent",
                                             TRUE ~ "Subsampled"))




designations.demo <- n144.df %>%
  dplyr::full_join(.,curated.population.metrics.df) %>%
  dplyr::mutate(designation = case_when(Simulated == TRUE & Detected == TRUE & aboveBF == TRUE ~ "Detected.CV",
                                        Simulated == TRUE & Detected == FALSE & aboveBF == FALSE ~ "Missed.CV",
                                        Simulated == TRUE & Detected == TRUE & aboveBF == FALSE ~ "CV.Not.Significant.In.Interval",
                                        Simulated == FALSE & Detected == TRUE & aboveBF == TRUE ~ "False.Discovery")) %>%
  dplyr::filter(designation != "False.Discovery") %>%
  dplyr::mutate(Sim.Var.Exp.Bin = cut(Simulated.QTL.VarExp, breaks = unique(c(seq(0,0.1,by = 0.025), 
                                                                              seq(0.1,0.9,by = 0.1))))) %>%
  tidyr::separate(col = detected.peak,
                  into = c("peak.CHROM","peak.POS"), 
                  sep = ":", remove = F) %>%
  dplyr::mutate(QTL.v.peak = abs(as.numeric(POS)-as.numeric(peak.POS))) %>%
  dplyr::group_by(Sim.Var.Exp.Bin, strain_set, population.group, designation) %>%
  dplyr::summarise(n = n()) %>%
  tidyr::pivot_wider(names_from = designation, values_from = n)
designations.demo[is.na(designations.demo)] <- 0

power.var.exp.demo <- designations.demo %>%
  dplyr::mutate(Detected = Detected.CV + CV.Not.Significant.In.Interval, 
                # normally would also add in False.Discovery, but figure is focusing on simulated QTL and leaving out doesn't change Power calculation
                Simulated = Detected.CV + CV.Not.Significant.In.Interval + Missed.CV,
                Power = Detected.CV/Simulated) %>%
  dplyr::group_by(Sim.Var.Exp.Bin, population.group) %>%
  dplyr::summarise(mean.Power = mean(Power),
                   sd.Power = sd(Power)) %>%
  dplyr::mutate(Sim.Var.Exp.Bin = gsub(as.character(Sim.Var.Exp.Bin), pattern = "\\(", replacement = "") %>%
                  gsub(., pattern = "\\]", replacement = "") %>%
                  gsub(., pattern = ",", replacement = "-")) %>%
  tidyr::separate(col = Sim.Var.Exp.Bin, into = c("top","bottom"), sep = "-", remove = FALSE) %>%
  dplyr::mutate(top = as.numeric(top)*100,
                bottom = as.numeric(bottom)*100) %>%
  tidyr::unite("Sim.Var.Exp.Bin", top:bottom, sep = "-", remove = FALSE)

#############################
### Supplemental Table 5  ###
#############################
supp.table.5 <- power.var.exp.demo %>%
  dplyr::ungroup() %>%
  dplyr::select(population.group, Sim.Var.Exp.Bin, mean.Power, sd.Power) %>%
  dplyr::mutate(power = if_else(is.na(sd.Power), 
                                true = as.character(round(mean.Power,2)), 
                                false = paste(round(mean.Power,2), round(sd.Power,2), sep = " ± "))) %>%
  dplyr::select(-c(mean.Power, sd.Power)) %>%
  tidyr::pivot_wider(names_from = Sim.Var.Exp.Bin, values_from = power)
colnames(supp.table.5) <- c("Population Type",colnames(supp.table.5)[2:length(colnames(supp.table.5))])
write.csv(supp.table.5, "tables/supplemental.table.5.csv", quote = F, row.names = F)


################
### Figure 4 ###
################
power.var.exp.demo %<>%
  dplyr::mutate(max = mean.Power + sd.Power,
                min = mean.Power - sd.Power,
                sd.Power.max = if_else(max > 1, 
                                       true = 1 - mean.Power,
                                       false = sd.Power),
                sd.Power.min = if_else(min < 0, 
                                       true = mean.Power,
                                       false = sd.Power)) %>%
  dplyr::mutate(detail = paste0(population.group," Strains; n = 144"))
A <- ggplot(power.var.exp.demo, mapping = aes(x = reorder(Sim.Var.Exp.Bin, bottom), y = mean.Power, 
                                                fill = detail,
                                                group = detail)) + 
  theme_bw(base_size = 11) + 
  geom_hline(yintercept = 0.8, linetype = 4, alpha = 0.2) + 
  geom_line(aes(colour = detail),position = position_dodge(width = 0.5)) +
  scale_size_manual(values = c(1,0.5,0.5), guide = "none") + 
  scale_colour_manual(values = c("blue","darkgreen","orange"), name = "Population Composition") +
  geom_errorbar(data = power.var.exp.demo,
                width = 1, alpha = 0.8,
                mapping = aes(y = mean.Power, ymax = mean.Power+sd.Power.max, ymin = mean.Power-sd.Power.min, colour = detail),
                position=position_dodge(width=0.5)) +
  geom_point(position = position_dodge(width = 0.5), shape = 21) + 
  scale_fill_manual(values = c("blue","darkgreen","orange"), name = "Population Composition") +
  ylim(c(0,1)) + 
  theme(legend.position = "right",
        legend.title = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.x = element_text()) +
  labs(x = "Variance Explained by Simulated QTL (%)",
       y = "Power")


MAF.by.VarExp.Bin.demo <- n144.df %>%
  dplyr::full_join(.,curated.population.metrics.df) %>%
  dplyr::mutate(designation = case_when(Simulated == TRUE & Detected == TRUE & aboveBF == TRUE ~ "Detected.CV",
                                        Simulated == TRUE & Detected == FALSE & aboveBF == FALSE ~ "Missed.CV",
                                        Simulated == TRUE & Detected == TRUE & aboveBF == FALSE ~ "CV.Not.Significant.In.Interval",
                                        Simulated == FALSE & Detected == TRUE & aboveBF == TRUE ~ "False.Discovery")) %>%
  dplyr::filter(designation != "False.Discovery") %>%
  dplyr::mutate(Sim.Var.Exp.Bin = cut(Simulated.QTL.VarExp, breaks = unique(c(seq(0,0.1,by = 0.025),
                                                                              seq(0.1,0.9,by = 0.1)))),
                Sim.Var.Exp.Bin = gsub(as.character(Sim.Var.Exp.Bin), pattern = "\\(", replacement = "") %>%
                  gsub(., pattern = "\\]", replacement = "") %>%
                  gsub(., pattern = ",", replacement = "-")) %>%
  tidyr::separate(col = Sim.Var.Exp.Bin, into = c("top","bottom"), sep = "-", remove = FALSE) %>%
  dplyr::mutate(top = as.numeric(top)*100,
                bottom = as.numeric(bottom)*100) %>%
  tidyr::unite("Sim.Var.Exp.Bin", top:bottom, sep = "-", remove = FALSE)
B <- ggplot(MAF.by.VarExp.Bin.demo, mapping = aes(x = reorder(Sim.Var.Exp.Bin,bottom), 
                                                    y = Frequency,
                                                    fill = population.group)) + 
  theme_bw(base_size = 11) + 
  geom_boxplot(outlier.alpha = 0.05, position=position_dodge(0.92)) +
  # geom_jitter(shape = 21 ,position=position_dodge()) +
  scale_fill_manual(values = c("blue","darkgreen","orange"), name = "Population Composition") + 
  theme(legend.position = "none",
        legend.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  labs(x = "Variance Explained by Simulated QTL (%)",
       y = "MAF")


AF.dists.demo <- n144.df %>%
  dplyr::full_join(.,curated.population.metrics.df) %>%
  dplyr::filter(Simulated == TRUE) %>%
  dplyr::mutate(detail = paste0(population.group," Strains; n = 144")) %>%
  dplyr::mutate(detail = as.factor(detail)) %>%
  dplyr::rename(`Minor Allele Frequency` = Frequency, 
                `Variance Explained by Simulated QTL (%)` = Simulated.QTL.VarExp) %>%
  dplyr::select(`Minor Allele Frequency`, 
                `Variance Explained by Simulated QTL (%)`, detail) %>%
  tidyr::pivot_longer(cols = c(`Minor Allele Frequency`, `Variance Explained by Simulated QTL (%)`), names_to = "metric", values_to = "value")
C <- AF.dists.demo %>%
  dplyr::filter(metric == "Minor Allele Frequency",
                ) %>%
  ggplot(., mapping = aes(x = value, colour = detail)) +
  geom_density(size = 0.4, adjust = 0.5) +
  theme_bw(base_size = 11) +
  scale_colour_manual(values = c("blue","darkgreen","orange"), name = "Population Composition") +
  # facet_grid(.~metric, scales = "free_x") +
  theme(legend.position = "none",
        legend.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank()) +
  labs(y = "Number of QTL",
       x = "Minor Allele Frequency")


A.legend <- cowplot::get_legend(A)
pre <- cowplot::plot_grid(A + theme(axis.text.x = element_blank(),
                                        axis.title.x = element_blank(),
                                        axis.ticks.x = element_blank(),
                                        legend.position = "none"), 
                            B, 
                            C, ncol = 1, align = 'v', axis = 'l', labels = "AUTO")
fig4 <- cowplot::ggdraw(pre + cowplot::draw_plot(A.legend, .5, .8, .6, 0))
ggsave(fig4, filename = "plots/figure.4.png", width = 7.5, height = 7)
