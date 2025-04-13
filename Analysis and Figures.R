setwd('/Users/ahmedelewa/Google Drive/My Drive/01_projects/074_augsburg/05_BIO474manuscript/homo_elegans/')

## Areflexia
ejle <- read.csv('EJLE.csv',header=TRUE)
# Reshape to long format
long_data <- reshape(ejle,
                     varying = c("N2", "RB"),
                     v.names = "value",
                     timevar = "strain",
                     times = c("N2", "RB"),
                     direction = "long")

# Set up plotting area: adjust rows/cols depending on number of react types
pdf('fig_EJLE.pdf')
par(mfrow = c(2,4),      # One plot per col
    mar = c(2, 4, 2, 1),    # Smaller margins: bottom, left, top, right
    oma = c(2, 0, 2, 0))    # Optional: outer margins
subset_data <- subset(ejle, react == 'bends')
boxplot(subset_data$N2, subset_data$RB,
        names = c("N2", expression(italic("ndrr-2"))), las =1,
        main = 'Bends', cex.lab=1.5, cex.main = 1.5,
        ylab = "mean", ylim = c(0,5),cex.axis=1.5,
        col = c("lightblue", "lightgreen"))
points(jitter(rep(1, nrow(subset_data))), subset_data$N2, pch = 16, col = "blue",cex=1.5)
points(jitter(rep(2, nrow(subset_data))), subset_data$RB, pch = 16, col = "darkgreen",cex=1.5)
t.test(subset_data[,2:3])
text(x = 1.5, y = 4.5, labels = "***", cex = 2)

subset_data <- subset(ejle, react == 'omega')
boxplot(subset_data$N2, subset_data$RB,
        names = c("N2", expression(italic("ndrr-2"))), las =1,
        main = 'Omega',cex.lab=1.5, cex.main = 1.5,
        ylab = "", ylim=c(0,9), cex.axis=1.5, 
        col = c("lightblue", "lightgreen"))
points(jitter(rep(1, nrow(subset_data))), subset_data$N2, pch = 16, col = "blue",cex=1.5)
points(jitter(rep(2, nrow(subset_data))), subset_data$RB, pch = 16, col = "darkgreen",cex=1.5)
t.test(subset_data[,2:3])
text(x = 1.5, y = 8.5, labels = "***", cex = 2)

subset_data <- subset(ejle, react == 'pause')
boxplot(subset_data$N2, subset_data$RB,
        names = c("N2", expression(italic("ndrr-2"))), las =1,
        main = 'Pauses', cex.lab=1.5, cex.main = 1.5,
        ylab = "", ylim=c(0,6), cex.axis=1.5,
        col = c("lightblue", "lightgreen"))
points(jitter(rep(1, nrow(subset_data))), subset_data$N2, pch = 16, col = "blue",cex=1.5)
points(jitter(rep(2, nrow(subset_data))), subset_data$RB, pch = 16, col = "darkgreen",cex=1.5)
t.test(subset_data[,2:3])
text(x = 1.5, y = 5.5, labels = "*", cex = 2)
dev.off()


## Hip flexor weakness
ly <- read.csv('LY.csv',header=TRUE)
# Reshape to long format
long_data <- reshape(ly,
                     varying = c("N2", "JD"),
                     v.names = "value",
                     timevar = "strain",
                     times = c("N2", "JD"),
                     direction = "long")
pdf('fig_LY.pdf')
par(mfrow = c(2,4),      # One plot per col
    mar = c(2, 4, 2, 1),    # Smaller margins: bottom, left, top, right
    oma = c(2, 0, 2, 0))    # Optional: outer margins
subset_data <- subset(ly, react == 'bends')
boxplot(subset_data$N2, subset_data$JD,
        names = c("N2", expression(italic("cca-1"))), las =1,
        main = 'Bends in 7 min', cex.lab=1.5, cex.main = 1.5,
        ylab = "count",ylim=c(50,165),cex.axis=1.5,
        col = c("lightblue", "lightgreen"))
points(jitter(rep(1, nrow(subset_data))), subset_data$N2, pch = 16, col = "blue",cex=1.5)
points(jitter(rep(2, nrow(subset_data))), subset_data$JD, pch = 16, col = "darkgreen",cex=1.5)
t.test(subset_data[,2:3])
text(x = 1.5, y = 160, labels = "***", cex = 2)
dev.off()
dim(ly)

## Aterial rupture
dkyd <- read.csv('DKYD.csv',header=TRUE)
par(mfrow = c(2, 3), mar = c(4, 4, 2, 1), oma = c(2, 1, 2, 1))
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
  plot(time, burst_N2, type = "l", col = "black", ylim = c(0, 1),
       xlab = "Time (min)", ylab = "Fraction Burst",
       main = paste("Burst Assay - Week", w), lwd = 3)
  lines(time, burst_l15, type = "l", col = "#4575B4", lwd = 3)
  lines(time, burst_l20, type = "l", col = "#FDAE61", lwd = 3)
  lines(time, burst_l25, type = "l", col = "#D73027", lwd = 3)
  
  legend("bottomright", legend = c(expression(N[2]),
                                  expression(italic(let-2(ts)) ~ "15" * degree * "C"),
                                  expression(italic(let-2(ts)) ~ "20" * degree * "C"),
                                  expression(italic(let-2(ts)) ~ "25" * degree * "C")),
         col = c("black", "#4575B4", "#FDAE61", "#D73027"),
         lty = 1, lwd = 3, bty = "n", cex = 0.8)
}

