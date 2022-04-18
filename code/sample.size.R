require(tidyverse)
require(RColorBrewer)
require(magrittr)
setwd(paste0(dirname(rstudioapi::getActiveDocumentContext()$path),"/.."))

###################################
### Sample Size Simulation Data ###
###################################
load(file = "data/File_S4.RData")
load(file = "data/File_S5.RData")
  
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
                # log10p = dplyr::if_else(Simulated == FALSE,
                #                         true = interval.log10p,
                #                         false = log10p), # false discoveries inherit the log10p value of the peak marker for the interval
                # log10p = as.numeric(log10p),
                interval_size = as.numeric(interval_size),
                aboveBF = dplyr::case_when(aboveBF == 1 ~ TRUE, 
                                           aboveBF == 0 ~ FALSE,
                                           is.na(aboveBF) ~ TRUE), # false discoveries by definition exceed significance threshold
                aboveBF = as.factor(aboveBF)) %>%
  dplyr::filter(CHROM != 7)


dat.population.features <- dat %>%
  dplyr::full_join(., population.metrics.df)



dat.population.features$pop.size <- factor(dat.population.features$pop.size, 
                                           levels = c("100","200","300","400","500"))
pop.size.pal <- c("#8CDBCF","#BDDBD0","#CBC9AD","#656839","#71861D")
names(pop.size.pal) <- levels(dat.population.features$pop.size)



designations <- dat.population.features %>%
  droplevels() %>%
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
  dplyr::group_by(Sim.Var.Exp.Bin, nQTL, strain_set, designation, pop.size) %>%
  dplyr::summarise(n = n()) %>%
  tidyr::pivot_wider(names_from = designation, values_from = n)
designations[is.na(designations)] <- 0

power.var.exp <- designations %>%
  dplyr::mutate(Detected = Detected.CV + CV.Not.Significant.In.Interval, 
                # normally would also add in False.Discovery, but figure is focusing on simulated QTL and leaving out doesn't change Power calculation
                Simulated = Detected.CV + CV.Not.Significant.In.Interval + Missed.CV,
                Power = Detected.CV/Simulated) %>%
  dplyr::group_by(Sim.Var.Exp.Bin, pop.size) %>%
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
### Supplemental Figure 3 ###
#############################
AF.dists <- dat.population.features %>%
  dplyr::filter(Simulated == TRUE,
                !duplicated(QTL)) %>%
  tidyr::unite("detail", pop.size, remove = F, sep = " Strains; n = ") %>%
  dplyr::mutate(detail = as.factor(detail))
pop.size.pal.B <- pop.size.pal
names(pop.size.pal.B) <- levels(AF.dists$detail)

AF.dist.plot <- AF.dists %>%
  tidyr::unite("detail", pop.size, remove = F, sep = " Strains; n = ") %>%
  ggplot(., mapping = aes(x = Frequency, fill = detail)) +
  geom_density(alpha = 0.6, adjust = 0.8, position = "stack") + 
  theme_bw(base_size = 12) +
  scale_fill_manual(values = pop.size.pal.B, name = "Sample Size") + 
  theme(legend.position = "right",
        legend.title = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) + 
  scale_x_continuous(breaks = seq(0,0.5,0.1)) +
  labs(x = "Minor Allele Frequency",
       y = "Smoothed Density")
ggsave(AF.dist.plot, filename = "plots/supp.fig.4.png", width = 7.5, height = 3)


#############################
### Supplemental Table 4  ###
#############################
supp.table.4 <- power.var.exp %>%
  dplyr::ungroup() %>%
  dplyr::select(pop.size, Sim.Var.Exp.Bin, mean.Power, sd.Power) %>%
  dplyr::mutate(power = if_else(is.na(sd.Power), 
                                true = as.character(round(mean.Power,2)), 
                                false = paste(round(mean.Power,2), round(sd.Power,2), sep = " ± "))) %>%
  dplyr::select(-c(mean.Power, sd.Power)) %>%
  tidyr::pivot_wider(names_from = Sim.Var.Exp.Bin, values_from = power)
colnames(supp.table.4) <- c("Sample Size",colnames(supp.table.4)[2:length(colnames(supp.table.4))])
write.csv(supp.table.4, "tables/supplemental.table.4.csv", quote = F, row.names = F)


#################
### Figure 3  ###
#################
power.var.exp %<>%
  dplyr::mutate(max = mean.Power + sd.Power,
                min = mean.Power - sd.Power,
                sd.Power.max = if_else(max > 1, 
                                       true = 1 - mean.Power,
                                       false = sd.Power),
                sd.Power.min = if_else(min < 0, 
                                       true = mean.Power,
                                       false = sd.Power)) %>%
  dplyr::mutate(pop.size = paste0("n = ",pop.size))

pop.size.pal.B <- pop.size.pal
names(pop.size.pal.B) <- unique(power.var.exp$pop.size)

A <- power.var.exp %>%
  ggplot(., mapping = aes(x = reorder(Sim.Var.Exp.Bin, bottom), y = mean.Power, 
                          fill = pop.size,
                          group = pop.size)) + 
  theme_bw(base_size = 12) + 
  geom_line(position=position_dodge(width=0.5), mapping = aes(colour = pop.size)) +
  geom_hline(yintercept = 0.8, linetype = 4, alpha = 0.2) + 
  # geom_smooth(se = F, span = 0.5, aes(size = pop.type, colour = pop.size)) + 
  scale_size_manual(values = c(1,0.5,0.5), guide = "none") + 
  scale_colour_manual(values = pop.size.pal.B, name = "Sample Size") +
  geom_errorbar(data = power.var.exp,
                width = 1, alpha = 0.8,
                mapping = aes(y = mean.Power, ymax = mean.Power+sd.Power.max, ymin = mean.Power-sd.Power.min, colour = pop.size),
                position=position_dodge(width=0.5)) +
  geom_point(position = position_dodge(width = 0.5), shape = 21, size = 2) + 
  ylim(c(0,1)) + 
  scale_fill_manual(values = pop.size.pal.B, name = "Sample Size") +
  theme(legend.position = "right",
        legend.title = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  labs(x = "Variance Explained by Simulated QTL (%)",
       y = "Power")
A.legend <- cowplot::get_legend(A)
A.plain <- cowplot::plot_grid(A + theme(legend.position = "none"))
fig2 <- cowplot::ggdraw(A.plain + cowplot::draw_plot(A.legend, .5, .5, .6, 0))
ggsave(fig2, filename = "plots/figure.3.png", width = 7.5, height = 4)


################
### Table 1  ###
################
overall.metrics <- dat.population.features %>%
  droplevels() %>%
  dplyr::mutate(designation = case_when(Simulated == TRUE & Detected == TRUE & aboveBF == TRUE ~ "Detected.CV",
                                        Simulated == TRUE & Detected == FALSE & aboveBF == FALSE ~ "Missed.CV",
                                        Simulated == TRUE & Detected == TRUE & aboveBF == FALSE ~ "CV.Not.Significant.In.Interval",
                                        Simulated == FALSE & Detected == TRUE & aboveBF == TRUE ~ "False.Discovery")) %>%
  dplyr::group_by(strain_set, designation, pop.size, Rep) %>%
  dplyr::summarise(n = n()) %>%
  tidyr::pivot_wider(names_from = designation, values_from = n)
overall.metrics[is.na(overall.metrics)] <- 0

overall.metrics %<>%
  dplyr::mutate(Detected = Detected.CV + CV.Not.Significant.In.Interval + False.Discovery, 
                Simulated = Detected.CV + CV.Not.Significant.In.Interval + Missed.CV,
                Power = Detected.CV/Simulated,
                FDR = False.Discovery/Detected)
overall.metrics[is.na(overall.metrics)] <- 0

sample.size.performance.table <- overall.metrics %>% 
  dplyr::group_by(pop.size) %>%
  dplyr::summarise(mean.Power = mean(Power),
                   sd.Power = sd(Power),
                   mean.FDR = mean(FDR),
                   sd.FDR = sd(FDR)) %>% 
  dplyr::mutate(Power = paste(round(mean.Power,2), round(sd.Power,2), sep = " ± "),
                FDR = paste(round(mean.FDR,2), round(sd.FDR,2), sep = " ± ")) %>%
  dplyr::select(pop.size, Power, FDR)
colnames(sample.size.performance.table) <- c("Sample Size","Power","FDR")
write.csv(sample.size.performance.table, "tables/table.1.csv", quote = F, row.names = F)

