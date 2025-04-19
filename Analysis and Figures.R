setwd('/Users/ahmedelewa/Google Drive/My Drive/01_projects/074_augsburg/05_BIO474manuscript/homo_elegans/')

##################################################################################################################
###################              Code associated with Augsburg BIO474 Manuscript 2025          ###################
###################                 For questions: Ahmed Elewa, elewa@email.com                ###################
##################################################################################################################

###################
#### Areflexia ####
###################

ejle <- read.csv('EJLE.csv',header=TRUE)
# Reshape to long format
long_data <- reshape(ejle,
                     varying = c('N2', 'RB'),
                     v.names = 'value',
                     timevar = 'strain',
                     times = c('N2', 'RB'),
                     direction = 'long')

# Set up plotting area: adjust rows/cols depending on number of react types
pdf('fig_EJLE.pdf')
par(mfrow = c(2,4),      # One plot per col
    mar = c(2, 4, 2, 1),    # Smaller margins: bottom, left, top, right
    oma = c(2, 0, 2, 0))    # Optional: outer margins
subset_data <- subset(ejle, react == 'bends')
boxplot(subset_data$N2, subset_data$RB,
        names = c('N2', expression(italic('ndrr-2'))), las =1,
        main = 'Bends', cex.lab=1.5, cex.main = 1.5,
        ylab = 'mean', ylim = c(0,5),cex.axis=1.5,
        col = c('#C0C0C0', '#F0F0F0'))

points(jitter(rep(1, nrow(subset_data))), subset_data$N2, pch = 16, col = 'black',cex=1.5)
points(jitter(rep(2, nrow(subset_data))), subset_data$RB, pch = 16, col = 'darkgray',cex=1.5)
t.test(subset_data[,2:3])
text(x = 1.5, y = 4.5, labels = '***', cex = 2)

subset_data <- subset(ejle, react == 'omega')
boxplot(subset_data$N2, subset_data$RB,
        names = c('N2', expression(italic('ndrr-2'))), las =1,
        main = 'Omega',cex.lab=1.5, cex.main = 1.5,
        ylab = '', ylim=c(0,9), cex.axis=1.5, 
        col = c('#C0C0C0', '#F0F0F0'))
points(jitter(rep(1, nrow(subset_data))), subset_data$N2, pch = 16, col = 'black',cex=1.5)
points(jitter(rep(2, nrow(subset_data))), subset_data$RB, pch = 16, col = 'darkgray',cex=1.5)
t.test(subset_data[,2:3])
text(x = 1.5, y = 8.5, labels = '***', cex = 2)

subset_data <- subset(ejle, react == 'pause')
boxplot(subset_data$N2, subset_data$RB,
        names = c('N2', expression(italic('ndrr-2'))), las =1,
        main = 'Pauses', cex.lab=1.5, cex.main = 1.5,
        ylab = '', ylim=c(0,6), cex.axis=1.5,
        col = c('#C0C0C0', '#F0F0F0'))
points(jitter(rep(1, nrow(subset_data))), subset_data$N2, pch = 16, col = 'black',cex=1.5)
points(jitter(rep(2, nrow(subset_data))), subset_data$RB, pch = 16, col = 'darkgray',cex=1.5)
t.test(subset_data[,2:3])
text(x = 1.5, y = 5.5, labels = '*', cex = 2)
dev.off()

#############################
#### Hip flexor weakness ####
#############################

ly <- read.csv('LY.csv',header=TRUE)
# Reshape to long format
long_data <- reshape(ly,
                     varying = c('N2', 'JD'),
                     v.names = 'value',
                     timevar = 'strain',
                     times = c('N2', 'JD'),
                     direction = 'long')
pdf('fig_LY.pdf')
par(mfrow = c(2,4),      # One plot per col
    mar = c(2, 4, 2, 1),    # Smaller margins: bottom, left, top, right
    oma = c(2, 0, 2, 0))    # Optional: outer margins
subset_data <- subset(ly, react == 'bends')
boxplot(subset_data$N2, subset_data$JD,
        names = c('N2', expression(italic('cca-1'))), las =1,
        main = 'Bends in 7 min', cex.lab=1.5, cex.main = 1.5,
        ylab = 'count',ylim=c(50,165),cex.axis=1.5,
        col = c('#C0C0C0', '#F0F0F0'))
points(jitter(rep(1, nrow(subset_data))), subset_data$N2, pch = 16, col = 'black',cex=1.5)
points(jitter(rep(2, nrow(subset_data))), subset_data$JD, pch = 16, col = 'darkgray',cex=1.5)
t.test(subset_data[,2:3])
text(x = 1.5, y = 160, labels = '***', cex = 2)
dev.off()

#########################
#### Aterial rupture ####
#########################
library(survival)
dkyd <- read.csv('DKYD.csv',header=TRUE)
pdf('fig_DKYD.pdf')
par(mfrow = c(2, 3), mar = c(4, 4, 2, 1), oma = c(0, 0, 0, 0))
# Unique weeks
weeks <- unique(dkyd$week)
# Plot for each week
for (w in weeks) {
  subset_data <- subset(dkyd, week == w)
  
  # Get cumulative deaths
  burst_N2 <- cumsum(subset_data$N2)
  burst_l15 <- cumsum(subset_data$l15)
  burst_l20 <- cumsum(subset_data$l20)
  burst_l25 <- cumsum(subset_data$l25)
  
  # Total at t=0
  start_N2 <- max(burst_N2)
  start_l15 <- max(burst_l15)
  start_l20 <- max(burst_l20)
  start_l25 <- max(burst_l25)
  
  # Calculate burst
  time <- subset_data$time
  burst_N2 <- burst_N2 / start_N2
  burst_l15 <-  burst_l15 / start_l15
  burst_l20 <-  burst_l20 / start_l20
  burst_l25 <- burst_l25 / start_l25
  
  # Plot survival
  plot(time, burst_N2, type = 'l', col = 'black', ylim = c(0, 1), las =1,
       xlab = 'Time (min)', ylab = 'Fraction Burst',
       main = paste('Burst Assay - Week', w), lwd = 3)
  lines(time, burst_l15, type = 'l', col = '#4575B4', lwd = 3)
  lines(time, burst_l20, type = 'l', col = '#FDAE61', lwd = 3)
  lines(time, burst_l25, type = 'l', col = '#D73027', lwd = 3)
  legend('bottomright', legend = c(expression(N[2]),
                                  expression(italic(let-2(ts)) ~ '15' * degree * 'C'),
                                  expression(italic(let-2(ts)) ~ '20' * degree * 'C'),
                                  expression(italic(let-2(ts)) ~ '25' * degree * 'C')),
         col = c('black', '#4575B4', '#FDAE61', '#D73026'),
         lty = 1, lwd = 3, bty = 'n', cex = 0.8)
}
dev.off()

# Count data
counts <- data.frame(
  time = 11:15,
  l15 = c(1, 0, 0, 10, 2),
  l25 = c(11, 2, 0, 0, 0)
)

# Expand counts to individual rows
expand_group <- function(time, count, group) {
  if (count == 0) return(NULL)
  data.frame(time = rep(time, count), group = group, status = 1)
}

df_long <- do.call(rbind, c(
  lapply(1:nrow(counts), function(i) expand_group(counts$time[i], counts$l15[i], 'l15')),
  lapply(1:nrow(counts), function(i) expand_group(counts$time[i], counts$l25[i], 'l25'))
))

# Build survival object
surv_obj <- Surv(df_long$time, df_long$status)

# Compare groups
survdiff(surv_obj ~ group, data = df_long)



######################
#### Clinodactyly ####
######################
wbpheno <- read.csv('WBPhenotypes.csv',header=TRUE,row.names = 1)
wbpheno <- wbpheno[wbpheno$RNAi.Phenotype.Observed!='N.A.'|wbpheno$Allele.Phenotype.Observed!='N.A.',]
dim(wbpheno)
#[1] 10411     2

# Convert the data frame into the desired format
wpo <- do.call(rbind, lapply(rownames(wbpheno), function(id) {
  # Extract features from columns A and B for each row
  features <- unlist(strsplit(c(wbpheno[id, 'RNAi.Phenotype.Observed'], wbpheno[id, 'Allele.Phenotype.Observed']), ',\\s*'))  # Split by comma and trim spaces
  features <- features[features != 'N.A.']  # Remove N.A.
  # Create a new data frame for this ID
  data.frame(ID = id, Feature = features)
}))

colnames(wpo) <- c('ID','phenotype')
wpo <- unique(wpo)

ortho <- read.csv('ortholist_master.csv',header=TRUE)
ortho <- ortho[ortho$WormBase.ID%in%rownames(wbpheno),]
dim(ortho)
#[1] 17334    11
w2h <- ortho[,c(1,6)]
colnames(w2h) <- c('worm','human')
w2h <- unique(w2h)

hpo <- read.csv('phenotype_to_genes.csv',header=TRUE)
hpo <- hpo[hpo$gene_symbol%in%ortho$HGNC.Symbol,]
dim(hpo)
#[1] 627480      5
hpo <- hpo[,c(4,2)]
colnames(hpo) <- c('ID','phenotype')
hpo <- unique(hpo)

hpoi <- c('Abnormal upper limb bone morphology') # Humann phenotype of interest
wpoi <- c('short') # Worm phenotype of interest

# Get human genes associated with hpoi phenotype
human_genes <- hpo$ID[hpo$phenotype == hpoi]
# Map human genes to worm orthologs
worm_genes_from_human <- unique(w2h$worm[w2h$human %in% human_genes])
# Number of orthologs linked to the human phenotype
n <- length(worm_genes_from_human)
# Get worm genes associated with wpoi
worm_genes_in_pheno <- unique(wpo$ID[wpo$phenotype == wpoi])
# Number of orthologs linked to the worm phenotype
m <- length(intersect(w2h$worm, worm_genes_in_pheno))
# Number of common orthologs between human and worm phenotypes
c <- length(intersect(worm_genes_from_human, worm_genes_in_pheno))
# Total number of orthologs shared between the two species
N <- nrow(w2h)
# Perform the hypergeometric test
p_value <- phyper(c - 1, m, N - m, n, lower.tail = FALSE)
p_value
genes_in_common = intersect(worm_genes_from_human, worm_genes_in_pheno)
unique(ortho[ortho$WormBase.ID%in%genes_in_common,'Common.Name'])
#[1] 'nca-2'   'sma-6'   'unc-73'  'unc-63'  'unc-77'  'ocr-4'   'pdl-1'   'unc-104' 'pezo-1'  'unc-38' 
#[11] 'egl-4'   'egl-27'  'osm-9'   'pqn-85'  'rnt-1'   'sem-4'   'unc-30'  'unc-89'  'lig-4'   'ssl-1'  
#[21] 'lev-1'   'egas-2' 

###########################
#### Cleft soft palate ####
###########################

mc <- read.csv('MC.csv')
mean(rowSums(mc[mc$strain=='LT',c(2,3)]))
sd(rowSums(mc[mc$strain=='LT',c(2,3)]))
mean(rowSums(mc[mc$strain=='N2',c(2,3)]))
sd(rowSums(mc[mc$strain=='N2',c(2,3)]))
t.test(c(rowSums(mc[mc$strain=='N2',c(2,3)]),rowSums(mc[mc$strain=='LT',c(2,3)])))

#################################
#### Abnormal cell phenotype ####
#################################

satw <- read.csv('SATW.csv')

# Convert MM:SS to total seconds
convert_to_seconds <- function(x) {
  parts <- strsplit(as.character(x), ':')[[1]]
  as.numeric(parts[1]) * 60 + as.numeric(parts[2])
}

# Apply conversion
satw$duration_sec <- sapply(satw$duration, convert_to_seconds)

mean(satw[satw$strain == 'N2', 'duration_sec' ])/60
sd(satw[satw$strain == 'N2', 'duration_sec' ])/60
mean(satw[satw$strain == 'polg', 'duration_sec' ])/60
sd(satw[satw$strain == 'polg', 'duration_sec' ])/60
wilcox.test(duration_sec ~ strain, data = satw[satw$student == 'S', ])
wilcox.test(duration_sec ~ strain, data = satw[satw$student == 'T', ])
wilcox.test(duration_sec ~ strain, data = satw[satw$student %in% c('S','T'), ])

pdf('fig_SATW.pdf')
par(mfrow = c(2,4), mar = c(2, 4, 2, 1), oma = c(2, 0, 2, 0)) 
# Filter for N2 and polg
filtered <- subset(satw, strain %in% c('N2', 'polg'))

# Map strain to numeric x-axis positions
x_pos <- ifelse(filtered$strain == 'N2', 1, 2)

# Jitter the x positions to prevent overlap
x_jittered <- jitter(x_pos, amount = 0.15)

# Draw boxplot
boxplot(duration_sec/60 ~ strain, data = filtered,
        main = 'Swimming Assay', xaxt = 'n',
        names = c('N2', expression(italic('polg-1/+'))), las =1,
        ylab = 'Swimming duration (min)', xlab = 'Strain', ylim = c(0,15),
         col = c('#C0C0C0', '#F0F0F0'))
axis(1, at = c(1, 2), labels = FALSE, tick = TRUE)
axis(1, at = c(1, 2), labels = c(expression(atop('N2')), expression(atop(italic('polg-1/+')))), tick = FALSE, line = 1.5, cex.axis = 1.2)

# Overlay individual points
points(x_jittered, filtered$duration_sec/60,
       pch = 16, col = c('black','darkgray'))
text(x = 1.5, y = 13, labels = '***', cex = 2)
legend('topleft', legend = c('Student 1', 'Student 2'),
       col = c('black','darkgray'),
       pch=16, bty = 'n')
dev.off()

#######################
#### Egg retention ####
#######################
pxtv <- read.csv('PXTV.csv')
pxtv <- pxtv[pxtv$total>3,]
pxtv$few = rowSums(pxtv[,4:7])

# Normalize stage counts by row total
pxtv[c('few_n','many_n')] <- pxtv[c('few','many')] / rowSums(pxtv[c('few','many')])

# Reshape to long format
df_long <- reshape(pxtv,
                   varying = list(c('few_n','many_n')),
                   v.names = 'proportion',
                   timevar = 'stage',
                   times = c('few_n','many_n'),
                   direction = 'long')

df_long$stage <- factor(df_long$stage, levels = c('few_n' ,'many_n'))
pdf('fig_PXTV.pdf')
par(mfrow = c(2,3), mar = c(2, 4, 2, 1), oma = c(2, 0, 2, 0)) 
boxplot(proportion ~ strain * stage,
        data = df_long,
        col = c('#C0C0C0', '#F0F0F0'),
        las = 2, ylim = c(0,1.15), xaxt='n', 
        xlab = '',
        ylab = 'Proportion of Embryos',
        main = 'Embryos retained in uterus')

# Overlay jittered data points, colored by student
xpos <- as.numeric(interaction(df_long$strain, df_long$stage))
points(jitter(xpos, amount = 0.15), df_long$proportion,
       pch = 16, col = c('black','darkgray'))

legend('topleft', legend = c('N2', expression(italic('Y64G10A.7'))),
       fill = c('#C0C0C0', '#F0F0F0'), bty = 'n')
legend('topright', legend = c('Student 1', 'Student 2'),
       col = c('black','darkgray'), pch = 16, bty = 'n')

axis(1, at = c(1.5, 3.5), labels = FALSE, tick = TRUE)
axis(1, at = c(1.5, 3.5), labels = c('1-4 cell', '>4 cell'), tick = FALSE, line = -0.3, cex.axis = 1.2)
text(x = 1.5, y = 0.95, labels = '**', cex = 2)
text(x = 3.5, y = 0.95, labels = '**', cex = 2)

dev.off()
wilcox.test(pxtv[pxtv$strain=='N2','few_n'],pxtv[pxtv$strain=='RB','few_n'])
wilcox.test(pxtv[pxtv$strain=='N2','many_n'],pxtv[pxtv$strain=='RB','many_n'])



#############################################################################
#### Modeling lethal short-limbed short stature C. elegans morphometrics ####
#############################################################################

blds <- read.csv('BLDS.csv')
# reorder strain factor
blds$strain <- factor(blds$strain, levels = c('N2', 'MT'))
# Normalize width to length of animal
blds$pTOl = blds$width_p/blds$length
blds$vTOl = blds$width_v/blds$length
# Get ratio of width at pharynx to vulva
blds$pTOv = blds$width_p/blds$width_v

# mean body length
mean(blds$length[blds$strain == 'N2' & !is.na(blds$length)])
sd(blds$length[blds$strain == 'N2' & !is.na(blds$length)])
length(blds$length[blds$strain == 'N2' & !is.na(blds$length)])
mean(blds$length[blds$strain == 'MT' & !is.na(blds$length)])
sd(blds$length[blds$strain == 'MT' & !is.na(blds$length)])
length(blds$length[blds$strain == 'MT' & !is.na(blds$length)])

# mean width at pharynx
mean(blds$width_p[blds$strain == 'N2' & !is.na(blds$width_p)])
sd(blds$width_p[blds$strain == 'N2' & !is.na(blds$width_p)])
length(blds$width_p[blds$strain == 'N2' & !is.na(blds$width_p)])
mean(blds$width_p[blds$strain == 'MT' & !is.na(blds$width_p)])
sd(blds$width_p[blds$strain == 'MT' & !is.na(blds$width_p)])
length(blds$width_p[blds$strain == 'MT' & !is.na(blds$width_p)])

# mean width at vulva
mean(blds$width_v[blds$strain == 'N2' & !is.na(blds$width_v)])
sd(blds$width_v[blds$strain == 'N2' & !is.na(blds$width_v)])
length(blds$width_v[blds$strain == 'N2' & !is.na(blds$width_v)])
mean(blds$width_v[blds$strain == 'MT' & !is.na(blds$width_v)])
sd(blds$width_v[blds$strain == 'MT' & !is.na(blds$width_v)])
length(blds$width_v[blds$strain == 'MT' & !is.na(blds$width_v)])

# mean width at normalized pharynx
mean(blds$pTOl[blds$strain == 'N2' & !is.na(blds$pTOl)])
sd(blds$pTOl[blds$strain == 'N2' & !is.na(blds$pTOl)])
length(blds$pTOl[blds$strain == 'N2' & !is.na(blds$pTOl)])
mean(blds$pTOl[blds$strain == 'MT' & !is.na(blds$pTOl)])
sd(blds$pTOl[blds$strain == 'MT' & !is.na(blds$pTOl)])
length(blds$pTOl[blds$strain == 'MT' & !is.na(blds$pTOl)])

# mean width at normalized vulva
mean(blds$vTOl[blds$strain == 'N2' & !is.na(blds$vTOl)])
sd(blds$vTOl[blds$strain == 'N2' & !is.na(blds$vTOl)])
length(blds$vTOl[blds$strain == 'N2' & !is.na(blds$vTOl)])
mean(blds$vTOl[blds$strain == 'MT' & !is.na(blds$vTOl)])
sd(blds$vTOl[blds$strain == 'MT' & !is.na(blds$vTOl)])
length(blds$vTOl[blds$strain == 'MT' & !is.na(blds$vTOl)])

wilcox.test(length ~ strain, data = blds[!is.na(blds$length),])
wilcox.test(width_p ~ strain, data = blds[!is.na(blds$width_p),])
wilcox.test(width_v ~ strain, data = blds[!is.na(blds$width_v),])
wilcox.test(pTOl ~ strain, data = blds[!is.na(blds$pTOl),])
wilcox.test(vTOl ~ strain, data = blds[!is.na(blds$vTOl),])
wilcox.test(pTOv ~ strain, data = blds[!is.na(blds$pTOv),])


pdf('fig_BLDS.pdf')
# Define layout matrix
par(mfrow = c(1, 1), mar = c(2, 4, 2, 1), oma = c(2, 0, 2, 0)) 

# Define a 2-row layout matrix: top row has plots, bottom row is all 0s (empty)
layout(matrix(c(
  1, 1, 2, 2, 2, 3, 3, 3,  # top row: plots
  0, 0, 0, 0, 0, 0, 0, 0   # bottom row: empty
), nrow = 2, byrow = TRUE),
widths = rep(1, 8),
heights = c(1, 1))  # equal height rows → canvas is split in half

# Define color mapping for trials
trial_colors_base <- c('1' = '#4575B4', '2' = '#FDAE61', '3' = '#D73027')
trial_colors <- sapply(trial_colors_base, adjustcolor, alpha.f = 0.5)

## Plot length
boxplot(blds[blds$strain=='N2','length'], blds[blds$strain=='MT','length'],
        xaxt = 'n',las=1,
        main = 'Length',cex.lab=1.2,
        ylab = 'millimeters', ylim=c(0,1.4), cex.axis=1.2, 
        col = c('#C0C0C0', '#F0F0F0'))
axis(1, at = c(1, 2), labels = FALSE, tick = TRUE)
axis(1, at = c(1, 2), labels = c(expression(atop('N2')), expression(atop(italic('egl-15')))), tick = FALSE, line = 1.5, cex.axis = 1.2)
subset <- blds$strain == 'N2' & !is.na(blds$length)
points(jitter(rep(1, sum(subset)), amount = 0.1), blds$length[subset], pch = 16,col = trial_colors[as.character(blds$trial[subset])])
subset <- blds$strain == 'MT' & !is.na(blds$length)
points(jitter(rep(2, sum(subset)), amount = 0.1), blds$length[subset], pch = 16,col = trial_colors[as.character(blds$trial[subset])])
legend('bottomleft', legend = c('Trial 1', 'Trial 2', 'Trial 3'), 
       col = trial_colors, pch = 16, bty = 'n')

## Plot width at pharynx and vulva
boxplot(blds[blds$strain == 'N2', 'width_p'],blds[blds$strain == 'MT', 'width_p'],
        blds[blds$strain == 'N2', 'width_v'],blds[blds$strain == 'MT', 'width_v'],
        main = 'Width', cex.lab=1.2,las = 1, xaxt = 'n',  # suppress x-axis
        ylab = 'millimeters', ylim = c(0, 0.125),
        cex.axis = 1.2,
        col = c('#C0C0C0', '#F0F0F0', '#C0C0C0', '#F0F0F0'))

# Plot data points 
subset <- blds$strain == 'N2' & !is.na(blds$width_p)
points(jitter(rep(1, sum(subset)), amount = 0.1), blds$width_p[subset], pch = 16,col = trial_colors[as.character(blds$trial[subset])])
subset <- blds$strain == 'MT' & !is.na(blds$width_p)
points(jitter(rep(2, sum(subset)), amount = 0.1), blds$width_p[subset], pch = 16,col = trial_colors[as.character(blds$trial[subset])])
subset <- blds$strain == 'N2' & !is.na(blds$width_v)
points(jitter(rep(3, sum(subset)), amount = 0.1), blds$width_v[subset], pch = 16,col = trial_colors[as.character(blds$trial[subset])])
subset <- blds$strain == 'MT' & !is.na(blds$width_v)
points(jitter(rep(4, sum(subset)), amount = 0.1), blds$width_v[subset], pch = 16,col = trial_colors[as.character(blds$trial[subset])])

# Add axis labels and legend
axis(1, at = c(1.5, 3.5), labels = FALSE, tick = TRUE)
axis(1, at = c(1.5, 3.5), labels = c('pharynx', 'vulva'), tick = FALSE, line = -0.3, cex.axis = 1.2)
text(x = 3.5, y = 0.095, labels = '*', cex = 2)
legend('topleft', legend = c('N2', expression(italic('egl-15'))),
       fill = c('#C0C0C0', '#F0F0F0'), bty = 'n', cex = 1)
legend('topright', legend = c('Trial 1', 'Trial 2', 'Trial 3'), 
       col = trial_colors, pch = 16, bty = 'n')


## Plot NORMALIZED width at pharynx and vulva (normalized to length)
boxplot(blds[blds$strain == 'N2', 'pTOl'],blds[blds$strain == 'MT', 'pTOl'],
        blds[blds$strain == 'N2', 'vTOl'],blds[blds$strain == 'MT', 'vTOl'],
        main = 'Width normalized to length', cex.lab=1.2,las = 1, xaxt = 'n',  # suppress x-axis
        ylab = 'millimeters', ylim = c(0, 0.125),
        cex.axis = 1.2,
        col = c('#C0C0C0', '#F0F0F0', '#C0C0C0', '#F0F0F0'))

# Plot data points 
subset <- blds$strain == 'N2' & !is.na(blds$pTOl)
points(jitter(rep(1, sum(subset)), amount = 0.1), blds$pTOl[subset], pch = 16,col = trial_colors[as.character(blds$trial[subset])])
subset <- blds$strain == 'MT' & !is.na(blds$pTOl)
points(jitter(rep(2, sum(subset)), amount = 0.1), blds$pTOl[subset], pch = 16,col = trial_colors[as.character(blds$trial[subset])])
subset <- blds$strain == 'N2' & !is.na(blds$vTOl)
points(jitter(rep(3, sum(subset)), amount = 0.1), blds$vTOl[subset], pch = 16,col = trial_colors[as.character(blds$trial[subset])])
subset <- blds$strain == 'MT' & !is.na(blds$vTOl)
points(jitter(rep(4, sum(subset)), amount = 0.1), blds$vTOl[subset], pch = 16,col = trial_colors[as.character(blds$trial[subset])])

# Add axis labels and legend
axis(1, at = c(1.5, 3.5), labels = FALSE, tick = TRUE)
axis(1, at = c(1.5, 3.5), labels = c('pharynx', 'vulva'), tick = FALSE, line = -0.3, cex.axis = 1.2)
text(x = 1.5, y = 0.095, labels = '*', cex = 2)
text(x = 3.5, y = 0.095, labels = '***', cex = 2)
legend('topleft', legend = c('N2', expression(italic('egl-15'))),
       fill = c('#C0C0C0', '#F0F0F0'), bty = 'n', cex = 1)
legend('topright', legend = c('Trial 1', 'Trial 2', 'Trial 3'), 
       col = trial_colors, pch = 16, bty = 'n')
dev.off()

