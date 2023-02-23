library(dplyr)
library(ggplot2)
library(grDevices)
library(RColorBrewer)
library(scales)
source('lib.R')

d <- bind_rows(lapply(list.files('UPENN_CART19_ALL/', pattern = 'intSites.csv', full.names = TRUE, recursive = TRUE), readr::read_csv))

d[grepl('D28', d$timePoint),]$timePoint <- 'M1'
d[grepl('D13', d$timePoint),]$timePoint <- 'M0.5'
d[grepl('D15', d$timePoint),]$timePoint <- 'M0.5'
d[grepl('D10', d$timePoint),]$timePoint <- 'M0.5'
d[grepl('D14', d$timePoint),]$timePoint <- 'M0.5'
d[grepl('D6', d$timePoint),]$timePoint <- 'M0.5'
d[grepl('D7', d$timePoint),]$timePoint <- 'M0.5'
d[grepl('D17', d$timePoint),]$timePoint <- 'M0.5'
d[grepl('D21', d$timePoint),]$timePoint <- 'M1'
d[grepl('D25', d$timePoint),]$timePoint <- 'M1'


# TET2 including UTRs.
w <- 50000
o <- subset(d, seqnames == 'chr4' & start >= (105146293-w) & start <= (105279816+w))

o$timePoint <- factor(o$timePoint, levels = arrange(o, timePointDays) %>% pull(timePoint) %>% unique())
o[grepl('T\\-CELLS', o$cellType),]$cellType <- 'T CELLS'
o[grepl('TCELLS:CAR\\+', o$cellType),]$cellType <- 'T CELLS'
o[grepl('ALPHA\\-BETA\\sT\\sCELLS', o$cellType),]$cellType <- 'T CELLS'
o[grepl('GAMMA\\-DELTA\\sT\\sCELLS', o$cellType),]$cellType <- 'T CELLS'
o$patient <- sub('pCHP959-158', 'p959-158', o$patient)

o <- select(o, patient, cellType, timePoint, posid, estAbund, relAbund)
o$clone <- paste(o$patient, o$posid)

openxlsx::write.xlsx(select(o, -clone) %>% arrange(timePoint), file = 'ALL.xlsx')

o$relAbund <- o$relAbund / 100

# Expanded: chr4+105272458

p1 <- ggplot(arrange(o, clone), aes(timePoint, relAbund, shape =  cellType, group = clone, fill = patient)) +
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

ggsave(p1, file="ALL_plot1.png", units = 'in', width = 10)

o$position <- as.integer(stringr::str_extract(o$posid, '\\d+$'))
o$strand <- stringr::str_extract(o$posid, '[\\+\\-]')
o$start <- o$position
o$end <- o$position
o$seqnames <- stringr::str_extract(o$posid, 'chr\\d+')
o$siteLabel <- o$patient

createIntUCSCTrack(o, siteLabel = 'siteLabel', outputFile = 'ALL_TET2.ucsc')

o2 <- o[! duplicated(o$clone),]
o2$orientation <- stringr::str_extract(o2$posid, '[\\+\\-]')
o2$orientation <- ifelse(o2$orientation == '+', 'positive', 'negative')

p2 <- ggplot(o2, aes(position, 0, fill = patient, shape = orientation)) +
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

ggsave(p2, file="ALL_plot2.png", units = 'in', width = 10)
