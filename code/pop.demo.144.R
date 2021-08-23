require(tidyverse)
require(RColorBrewer)
require(magrittr)
setwd(paste0(dirname(rstudioapi::getActiveDocumentContext()$path),"/.."))

#####################################################################
### Subsampled, Swept, and Divergent Strain Comparisons @ n = 144 ###
#####################################################################
load(file = "data/File_S4.RData")

designations.demo <- n144.df %>%
  dplyr::filter(algorithm == "MIXED",
                h2 == 0.8) %>%
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
