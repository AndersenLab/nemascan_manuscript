require(data.table)
require(tidyverse)
require(RColorBrewer)
require(cowplot)
setwd(paste0(dirname(rstudioapi::getActiveDocumentContext()$path),"/.."))
####################################
### Recapitulating Previous QTVs ###
####################################
load(file = "data/File_S10.RData")
expected.P.val.dists <- function(x,y){
  my.ps <- x %>%
    dplyr::select(marker, P)
  my.ps$exp.P <- (rank(my.ps$P, ties.method="first"))/(length(my.ps$P))
  
  x %>%
    dplyr::left_join(.,my.ps) %>%
    dplyr::mutate(trait = y)
}
calc.genomic.inflation <- function(x,y,z){
  gif.chisq <- qchisq(1-x$P,1)
  gif <- median(gif.chisq)/qchisq(0.5,1)
  trait <- y
  gif.df <- data.frame(gif, y, z)
  colnames(gif.df) <- c("GIF","trait","platform")
  return(gif.df)
}

##################
### Figure 6A  ###
##################
genome <- data.table::fread("data/genome.bed.tsv") %>%
  dplyr::rename(CHROM = chr, startPOS = start, endPOS = stop) %>%
  tidyr::pivot_longer(cols = c(startPOS, endPOS), names_to = "x", values_to = "POS") %>%
  dplyr::select(-x)

combined.significant.hits <- processed.nemascan.mappings.inbred %>%
  dplyr::full_join(., processed.nemascan.mappings.loco) %>%
  dplyr::full_join(., processed.cegwas.mappings) %>%
  dplyr::mutate(trait = if_else(trait == "Dauer Pheromone Response",
                                true = "Pheromone Response", 
                                false = trait)) %>%
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
  facet_grid(rows = trait~CHROM, scales = "free_x", space = "free", drop = F) + 
  labs(x = "Genomic position (Mb)")


#################
### Figure 6B ###
#################
unique.hits <- processed.nemascan.mappings.inbred %>%
  dplyr::full_join(.,processed.nemascan.mappings.loco) %>%
  dplyr::mutate(trait = if_else(trait == "Dauer Pheromone Response",
                                true = "Pheromone Response", 
                                false = trait)) %>%
  dplyr::select(CHROM, marker, POS, AF1, BETA, SE, P, log10p, trait, BF, aboveBF, platform) %>%
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
  dplyr::mutate(trait = if_else(trait == "Dauer Pheromone Response",
                                true = "Pheromone Response", 
                                false = trait)) %>%
  dplyr::mutate(marker = gsub(marker, pattern = "_", replacement = ":"),
                aboveBF = if_else(log10p > BF, true = 1, false = 0)) %>%
  dplyr::arrange(aboveBF)

thresh <- combined.exp.log10ps %>%
  dplyr::group_by(trait, platform) %>%
  dplyr::summarise(thresh = unique(BF))

nested.for.gif <- combined.exp.log10ps %>%
  dplyr::group_by(trait, platform) %>%
  tidyr::nest()

gif.plotting <- purrr::pmap(.l = list(nested.for.gif$data, 
                 nested.for.gif$trait, 
                 nested.for.gif$platform), 
            
            .f = calc.genomic.inflation) %>%
  Reduce(rbind,.) %>%
  dplyr::mutate(lambda = "λ[GC]",
                GIF = paste(" =", round(GIF,digits = 2)))



qq_platform <- ggplot() + 
  theme_bw(base_size = 11) + 
  geom_hline(data = thresh, aes(yintercept = thresh)) + 
  geom_text(data = gif.plotting, 
            mapping = aes(x = 1, y = 13, label = lambda),
            parse = T, 
            size = 3.3) + 
  geom_text(data = gif.plotting, 
            mapping = aes(x = 2.3, y = 13, label = GIF),
            size = 3.3) + 
  geom_point(data = combined.exp.log10ps,
             mapping = aes(x = -log10(exp.P), y = -log10(P), colour = as.factor(aboveBF)),
             size = 2) + 
  scale_colour_manual(values = c("black","red"), name = "Platform") +
  geom_abline(slope = 1) +
  facet_grid(platform~trait) + 
  theme(legend.position = "top",
        panel.grid = element_blank(),
        strip.text = element_text(size = 7)) + 
  labs(y = expression(Observed -log[10](italic(p))),
       x = expression(Expected -log[10](italic(p))))

################
### Figure 6 ###
################
fig6 <- cowplot::plot_grid(mapping.summary + theme(legend.position = "right"),
                           qq_platform + theme(legend.position = "none"), 
                           ncol = 1,rel_heights = c(1,1), labels = "AUTO")
ggsave(fig6, filename = "plots/figure.6.png", width = 7.5, height = 7)
ggsave(fig6, filename = "plots/figure.6.pdf", width = 7.5, height = 7)


#############################
### Supplemental Figure 4 ###
#############################
combined.processed.mappings <- processed.cegwas.mappings %>%
  dplyr::full_join(processed.nemascan.mappings.inbred) %>%
  dplyr::full_join(processed.nemascan.mappings.loco) %>%
  dplyr::mutate(trait = if_else(trait == "Dauer Pheromone Response",
                                true = "Pheromone Response", 
                                false = trait)) %>%
  dplyr::select(CHROM, POS, log10p, trait, BF, aboveBF, platform) %>%
  dplyr::filter(CHROM != "MtDNA") %>%
  dplyr::mutate(significance = if_else(aboveBF == 1, true = trait, false = "NA"),
                platform = as.factor(platform),
                significance = as.factor(significance)) %>%
  dplyr::distinct() %>%
  dplyr::arrange(aboveBF)
combined.processed.mappings$platform <- factor(combined.processed.mappings$platform, levels = c("Previous Mappings","NemaScan: Inbred","NemaScan: LOCO"))
combined.processed.mappings$significance <- factor(combined.processed.mappings$significance, 
                                                   levels = c("Abamectin Resistance","Arsenic Resistance","Pheromone Response","Etoposide Resistance","Telomere Length","NA"))
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
ggsave(multitrait.manplot + guides(colour = guide_legend(nrow = 2)), filename = "plots/supp.fig.5.png", width = 7.5, height = 6.5)


