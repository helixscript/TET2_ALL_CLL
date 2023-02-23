library(dplyr)
library(ggplot2)
library(grDevices)
library(RColorBrewer)
library(scales)
source('lib.R')

d <- bind_rows(lapply(list.files('.', pattern = 'intSites.csv', full.names = TRUE, recursive = TRUE), readr::read_csv))

# TET2 including UTRs.
w <- 50000
o <- subset(d, seqnames == 'chr4' & start >= (105146293-w) & start <= (105279816+w))

o$timePointMonths2 <- paste0('M', o$timePointMonths)

o[o$timePoint == 'D14',]$timePointMonths2 <- 'M0.5'
o[o$timePoint == 'D21',]$timePointMonths2 <- 'M0.5'

o <- o[o$posid != 'chr4+105272458',]

o$timePointMonths2 <- factor(o$timePointMonths2, levels = arrange(o, timePointDays) %>% pull(timePointMonths2) %>% unique())
o[grepl('TCELLS:CAR\\+CD8\\-', o$cellType),]$cellType <- 'T CELLS'
o[grepl('TCELLS:CAR\\+CD8\\+', o$cellType),]$cellType <- 'T CELLS'


o <- select(o, patient, cellType, timePointMonths2, posid, estAbund, relAbund)
o$clone <- paste(o$patient, o$posid)

openxlsx::write.xlsx(select(o, -clone) %>% rename(timePoint = timePointMonths2) %>% arrange(timePoint), file = 'o.xlsx')

o$relAbund <- o$relAbund / 100

# Expanded: chr4+105272458

ggplot(arrange(o, clone), aes(timePointMonths2, relAbund, shape =  cellType, group = clone, fill = patient)) +
  theme_bw() +
  scale_shape_manual(values = c(21:25))+
  scale_fill_manual(values = colorRampPalette(brewer.pal(12, "Paired"))(n_distinct(o$patient))) +
  geom_point(alpha = 0.8, size = 4, position = position_jitter(width = 0.2,height = 0.0002, seed = 2)) +
  geom_line(position = position_jitter(width = 0.2, height = 0.0002, seed = 2)) +
  scale_y_continuous(labels = scales::percent_format()) +
  guides(fill = guide_legend(override.aes=list(shape=21))) +
  labs(x = 'Timepoint', y = 'Sample Relative Abundance') +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),  
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16))


o$position <- as.integer(stringr::str_extract(o$posid, '\\d+$'))
o$strand <- stringr::str_extract(o$posid, '[\\+\\-]')
o$start <- o$position
o$end <- o$position
o$seqnames <- stringr::str_extract(o$posid, 'chr\\d+')
o$siteLabel <- o$patient

#   chr4:105,146,293-105,279,816
o2 <- o[! duplicated(o$clone),]

o2$orientation <- stringr::str_extract(o2$posid, '[\\+\\-]')


ggplot(o2, aes(position, 0, fill = patient, shape = orientation)) +
  theme_bw() +
  geom_hline(yintercept = 0) +
  geom_rect(data=tibble(x1 =105146293, x2 = 105279816, y1 = 0.25, y2 = -0.25), 
            mapping=aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2), color="black", alpha=0.10, inherit.aes = FALSE) +
  scale_fill_manual(values = colorRampPalette(brewer.pal(12, "Paired"))(n_distinct(o2$patient))) +
  scale_shape_manual(values = c(21, 24))+
  geom_point(alpha = 0.8, size = 3, position = position_jitter(width = 0,height = 0.05, seed = 2)) +
  ylim(c(-0.75, 0.75)) +
  labs(x = 'Genomic position', y = '') +
  guides(fill = guide_legend(override.aes=list(shape=21))) +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(size = 12),
        axis.title.x = element_text(size = 16))

createIntUCSCTrack(o, siteLabel = 'siteLabel')


