require(data.table)
require(tidyverse)
require(RColorBrewer)
require(cowplot)
setwd(paste0(dirname(rstudioapi::getActiveDocumentContext()$path),"/.."))
####################################
### Recapitulating Previous QTVs ###
####################################
load(file = "data/File_S7.RData")
load(file = "data/File_S8.RData")
expected.P.val.dists <- function(x,y){
  my.ps <- x %>%
    dplyr::select(marker, P)
  my.ps$exp.P <- (rank(my.ps$P, ties.method="first"))/(length(my.ps$P))
  
  x %>%
    dplyr::left_join(.,my.ps) %>%
    dplyr::mutate(trait = y)
}


##################
### Figure 5A  ###
##################
genome <- data.table::fread("data/genome.bed.tsv") %>%
  dplyr::rename(CHROM = chr, startPOS = start, endPOS = stop) %>%
  tidyr::pivot_longer(cols = c(startPOS, endPOS), names_to = "x", values_to = "POS") %>%
  dplyr::select(-x)

combined.significant.hits <- processed.nemascan.mappings %>%
  dplyr::full_join(., processed.cegwas.mappings) %>%
  dplyr::mutate(status = "real") %>%
  dplyr::select(CHROM, POS, log10p, BF, aboveBF,  trait, platform, startPOS, peakPOS, endPOS, status) %>%
  dplyr::distinct() %>%
  dplyr::filter(CHROM != "MtDNA",
                aboveBF == 1) %>%
  dplyr::arrange(log10p)

all.combos <- combined.significant.hits %>%
  dplyr::group_by(trait, platform, CHROM) %>%
  tidyr::nest()

lims <- list()
for(i in 1:length(all.combos$data)){
  lims[[i]] <- genome %>%
    dplyr::filter(genome$CHROM %in% all.combos$CHROM[i]) %>%
    dplyr::mutate(platform = all.combos$platform[i],
                  trait = all.combos$trait[i]) 
}
lims.df <- lims %>%
  dplyr::bind_rows(.) %>%
  dplyr::mutate(log10p = NA,
                aboveBF = 1,
                status = "fake")

mapping.summary <- combined.significant.hits %>%
  dplyr::full_join(., lims.df) %>%
  ggplot(., mapping = aes(x = POS/1000000, y = platform, colour = log10p, alpha = status)) + 
  theme_bw(base_size = 11) + 
  geom_jitter(height = 0.01, width = 0, size = 2) + 
  # scale_colour_distiller(palette = "BuPu", direction = 1, name = expression(-log[10](italic(p)))) +
  scale_color_gradient(high = "#D7263D", low = "#C6EEFA", name = expression(-log[10](italic(p)))) + 
  scale_alpha_manual(values = c(0,1), guide = "none") + 
  theme(panel.grid = element_blank(),
        axis.title.y = element_blank(),
        strip.text.y = element_text(angle=0, size = 7),
        legend.position = "top") + 
  facet_grid(rows = trait~CHROM, scales = "free", drop = F) + 
  labs(x = "Genomic position (Mb)")


##################
### Figure 5B  ###
##################
unique.hits <- processed.nemascan.mappings %>%
  dplyr::select(CHROM, marker, POS, AF1, BETA, SE, P, algorithm, log10p, trait, BF, aboveBF, platform) %>%
  dplyr::distinct()
nested.unique.hits <- unique.hits %>%
  dplyr::group_by(trait) %>% 
  tidyr::nest()
nested.unique.hits.qq <- purrr::map2(nested.unique.hits$data,
                                     nested.unique.hits$trait,
                                     expected.P.val.dists) %>%
  dplyr::bind_rows(.)

unique.previous.hits <- processed.cegwas.mappings %>%
  dplyr::select(CHROM, marker, POS, log10p, P, trait, BF, aboveBF, platform) %>%
  dplyr::distinct()

nested.unique.previous.hits <- unique.previous.hits %>%
  dplyr::group_by(trait) %>% 
  tidyr::nest()

nested.unique.previous.hits.qq <- purrr::map2(nested.unique.previous.hits$data,
                                              nested.unique.previous.hits$trait,
                                              expected.P.val.dists) %>%
  dplyr::bind_rows(.)

combined.exp.log10ps <- nested.unique.hits.qq %>%
  dplyr::full_join(., nested.unique.previous.hits.qq) %>%
  dplyr::mutate(marker = gsub(marker, pattern = "_", replacement = ":"))

thresh <- combined.exp.log10ps %>%
  dplyr::group_by(trait, platform) %>%
  dplyr::summarise(thresh = unique(BF))

qq_platform <- ggplot(combined.exp.log10ps, mapping = aes(x = -log10(exp.P), y = -log10(P), colour = platform)) + 
  theme_bw(base_size = 11) + 
  geom_point(alpha = 0.2) + 
  scale_colour_manual(values = c("blue","red"), name = "Platform") + 
  geom_abline(slope = 1) +
  geom_hline(data = thresh, aes(yintercept = thresh)) + 
  facet_grid(.~trait) + 
  theme(legend.position = "top",
        panel.grid = element_blank(),
        strip.text = element_text(size = 7)) + 
  labs(y = expression(Observed -log[10](italic(p))),
       x = expression(Expected -log[10](italic(p))))


#################
### Figure 5  ###
#################
fig5 <- cowplot::plot_grid(mapping.summary + theme(legend.position = "right"),
                           qq_platform + theme(legend.position = "bottom"), 
                           ncol = 1,rel_heights = c(1,1), labels = "AUTO")
ggsave(fig5, filename = "plots/figure.5.png", width = 7.5, height = 5.5)


##############################
### Supplemental Figure 4  ###
##############################
combined.processed.mappings <- processed.cegwas.mappings %>%
  dplyr::full_join(processed.nemascan.mappings) %>%
  dplyr::select(CHROM, POS, log10p, trait, BF, aboveBF, platform) %>%
  dplyr::filter(CHROM != "MtDNA") %>%
  dplyr::mutate(significance = if_else(aboveBF == 1, true = trait, false = "NA"),
                platform = as.factor(platform),
                significance = as.factor(significance)) %>%
  dplyr::distinct() %>%
  dplyr::arrange(aboveBF)
combined.processed.mappings$platform <- factor(combined.processed.mappings$platform, levels = c("Previous Mappings","NemaScan"))
combined.processed.mappings$significance <- factor(combined.processed.mappings$significance, 
                                                   levels = c("Abamectin Resistance","Arsenic Resistance","Dauer Pheromone Response","Etoposide Resistance","Telomere Length","NA"))
multitrait.manplot <- ggplot() + 
  geom_point(combined.processed.mappings[which(combined.processed.mappings$significance == "NA"),], 
             mapping = aes(x = POS/1000000, y = log10p)) + 
  geom_point(combined.processed.mappings[which(combined.processed.mappings$significance != "NA"),], 
             mapping = aes(x = POS/1000000, y = log10p, colour = significance)) + 
  theme_bw(base_size = 12) + 
  facet_grid(platform ~ CHROM, space = "free", scales = "free_x") + 
  scale_colour_manual(values = c(brewer.pal(n = 5, name = "Set2")), name = "Trait") + 
  theme(legend.position = "top",
        panel.grid = element_blank(),
        legend.text = element_text(size = 8)) + 
  labs(x = "Genomic position (Mb)",
       y = expression(-log[10](italic(p))))
ggsave(multitrait.manplot + guides(colour = guide_legend(nrow = 2)), filename = "plots/supp.fig.4.png", width = 7.5, height = 4.75)


##############################
### Supplemental Figure 5  ###
##############################
raw.abamectin.loco.qq <- expected.P.val.dists(raw.abamectin.loco, "Abamectin Resistance") %>%
  dplyr::select(marker, P, algorithm, exp.P, trait)
raw.abamectin.inbred.qq <- expected.P.val.dists(raw.abamectin.inbred, "Abamectin Resistance")  %>%
  dplyr::select(marker, P, algorithm, exp.P, trait)


raw.arsenic.loco.qq <- expected.P.val.dists(raw.arsenic.loco, "Arsenic Resistance") %>%
  dplyr::select(marker, P, algorithm, exp.P, trait)
raw.arsenic.inbred.qq <- expected.P.val.dists(raw.arsenic.inbred, "Arsenic Resistance")  %>%
  dplyr::select(marker, P, algorithm, exp.P, trait)


raw.dauer.loco.qq <- expected.P.val.dists(raw.dauer.loco, "Dauer Pheromone Response") %>%
  dplyr::select(marker, P, algorithm, exp.P, trait)
raw.dauer.inbred.qq <- expected.P.val.dists(raw.dauer.inbred, "Dauer Pheromone Response")  %>%
  dplyr::select(marker, P, algorithm, exp.P, trait)


raw.etoposide.loco.qq <- expected.P.val.dists(raw.etoposide.loco, "Etoposide Resistance") %>%
  dplyr::select(marker, P, algorithm, exp.P, trait)
raw.etoposide.inbred.qq <- expected.P.val.dists(raw.etoposide.inbred, "Etoposide Resistance")  %>%
  dplyr::select(marker, P, algorithm, exp.P, trait)


raw.telomere.loco.qq <- expected.P.val.dists(raw.telomere.loco, "Telomere Length") %>%
  dplyr::select(marker, P, algorithm, exp.P, trait)
raw.telomere.inbred.qq <- expected.P.val.dists(raw.telomere.inbred, "Telomere Length")  %>%
  dplyr::select(marker, P, algorithm, exp.P, trait)


raw.qqs <- raw.abamectin.inbred.qq %>%
  dplyr::full_join(., raw.abamectin.loco.qq) %>%
  
  dplyr::full_join(., raw.arsenic.inbred.qq) %>%
  dplyr::full_join(., raw.arsenic.loco.qq) %>%
  
  dplyr::full_join(., raw.dauer.inbred.qq) %>%
  dplyr::full_join(., raw.dauer.loco.qq) %>%
  
  dplyr::full_join(., raw.etoposide.inbred.qq) %>%
  dplyr::full_join(., raw.etoposide.loco.qq) %>%
  
  dplyr::full_join(., raw.telomere.inbred.qq) %>%
  dplyr::full_join(., raw.telomere.loco.qq)

threshes <- processed.nemascan.mappings %>%
  dplyr::select(trait, BF) %>%
  dplyr::distinct()

loco_inbred_trait.comps <- raw.qqs %>%
  dplyr::full_join(., threshes) %>%
  dplyr::mutate(aboveBF = if_else(-log10(P) > BF, true = "Significant Association", false = "Below Significance Threshold")) %>%
  ggplot(., mapping = aes(x = -log10(exp.P), y = -log10(P), colour = as.factor(aboveBF))) + 
  theme_bw(base_size = 12) +
  geom_point() + 
  scale_colour_manual(values = c("black","red"), name = "Above Significance Threshold") +
  geom_abline(slope = 1) +
  facet_grid(trait~algorithm) +
  theme(legend.position = "top",
        legend.justification = "center",
        panel.grid = element_blank(),
        strip.text.y = element_text(size = 7.5)) + 
  labs(y = expression(Observed -log[10](italic(p))),
       x = expression(Expected -log[10](italic(p))))

ggsave(loco_inbred_trait.comps, filename = "plots/supp.fig.5.png", width = 7.5, height = 9)
