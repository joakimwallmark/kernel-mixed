library(gridExtra)
library(ggplot2)
library(ggpubr)

xScorePlot <- function(df, ylab, y_min = NA, y_max = NA) {
  res_plot <- ggplot(df, aes(x = score, y = value, linetype = method)) +
    geom_line() +
    theme_bw() +
    xlab("X score") +
    ylab(ylab) +
    scale_linetype_manual(values = c(1, 2, 3), labels = c("EG", "CE", "PSE")) +
    theme(
      legend.position = "bottom",
      legend.title = element_blank(),
      axis.title.x = element_blank()
    ) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(limits = c(y_min, y_max), expand = c(0, 0))
}


# Splines ---------------------------------------------------------------------
load(file = "fitted_splines.RData")
n <- 1500
R <- 1000
bin_items <- c(50)
poly_items <- c(10)
bin_items_A <- c(25)
poly_items_A <- c(5)

harder_y <- c(F, T) # for filenames
diff_pop <- c(F, T) # for filenames

maxs <- 80
no_scen <- length(harder_y) * length(diff_pop) * length(bin_items) * length(bin_items_A)
biasplots <- vector(mode = "list", length = 0)
seeplots <- vector(mode = "list", length = 0)
# loop through scenarios
for (item_scen_A in 1:length(bin_items_A)) {
  for (item_scen in 1:length(bin_items)) {
    for (pop_scen in 1:length(diff_pop)) {
      for (test_scen in harder_y) {
        no_bin <- bin_items[item_scen]
        no_poly <- poly_items[item_scen]
        no_bin_A <- bin_items_A[item_scen_A]
        no_poly_A <- poly_items_A[item_scen_A]
        filename_err <- paste("data/spline/spline Error R", R, " n", n, " b", no_bin, "p", no_poly,
          " bA", no_bin_A, "pA", no_poly_A,
          " diffy ", test_scen, " diffpop ", diff_pop[pop_scen], ".RData",
          sep = ""
        )
        load(filename_err)

        score <- rep(0:maxs, 6)
        method <- c(
          rep("IRTKE EG", maxs + 1),
          rep("LLKE EG", maxs + 1),
          rep("IRTKE CE", maxs + 1),
          rep("LLKE CE", maxs + 1),
          rep("IRTKE PSE", maxs + 1),
          rep("LLKE PSE", maxs + 1)
        )
        bias <- c(
          res$spline$local$IRTKE2$bias,
          res$spline$local$IRTKE$bias,
          res$spline$local$KE$bias,
          res$spline$local$IRTKE2CE$bias,
          res$spline$local$IRTKECE$bias,
          res$spline$local$KECE$bias,
          res$spline$local$IRTKE2PSE$bias,
          res$spline$local$IRTKEPSE$bias,
          res$spline$local$KEPSE$bias
        )
        bias_df <- data.frame(score, method, value = bias)
        sees <- c(
          res$spline$local$IRTKE2$SEE,
          res$spline$local$IRTKE$SEE,
          res$spline$local$KE$SEE,
          res$spline$local$IRTKE2CE$SEE,
          res$spline$local$IRTKECE$SEE,
          res$spline$local$KECE$SEE,
          res$spline$local$IRTKE2PSE$SEE,
          res$spline$local$IRTKEPSE$SEE,
          res$spline$local$KEPSE$SEE
        )
        see_df <- data.frame(score, method, value = sees)
        bias_df$method <- factor(bias_df$method, levels = c("IRTKE EG", "LLKE EG", "IRTKE CE", "LLKE CE", "IRTKE PSE", "LLKE PSE"))
        see_df$method <- factor(see_df$method, levels = c("IRTKE EG", "LLKE EG", "IRTKE CE", "LLKE CE", "IRTKE PSE", "LLKE PSE"))


        y_max <- max(bias_df[!grepl("EG", bias_df$method), ]$value) + 0.05 * max(bias_df[!grepl("EG", bias_df$method), ]$value)
        y_min <- min(bias_df[!grepl("EG", bias_df$method), ]$value) + 0.05 * min(bias_df[!grepl("EG", bias_df$method), ]$value)
        y_max <- max(bias_df$value) + 0.05 * max(bias_df$value)
        y_min <- min(bias_df$value) + 0.05 * min(bias_df$value)

        pbiasirtke <- xScorePlot(bias_df[grepl("IRTKE ", bias_df$method), ], y_min = y_min, y_max = y_max, "Bias")
        pbiasllke <- xScorePlot(bias_df[grepl("LLKE", bias_df$method), ], y_min = y_min, y_max = y_max, "Bias")
        
        y_max <- max(see_df$value) + 0.05 * max(see_df$value)
        pseeirtke <- xScorePlot(see_df[grepl("IRTKE ", see_df$method), ], y_min = 0, y_max = y_max, "SEE")
        pseellke <- xScorePlot(see_df[grepl("LLKE", see_df$method), ], y_min = 0, y_max = y_max, "SEE")
        
        scenario <- paste(no_bin, "p", no_poly, " bA", no_bin_A, "pA", no_poly_A, " diffy ", test_scen, " diffpop ", diff_pop[pop_scen])
        biasplots[[scenario]] <- list(biasirtke = pbiasirtke, biasllke = pbiasllke)
        seeplots[[scenario]] <- list(seeirtke = pseeirtke, seellke = pseellke)
      }
    }
  }
}

b1 <- biasplots$`50 p 10  bA 25 pA 5  diffy  FALSE  diffpop  FALSE`$biasirtke +
  labs(title = expression(paste("IRTKE  ", X == Y, "  ", P == Q))) +
  theme(plot.title = element_text(size = 9, hjust = 0.5))
b2 <- biasplots$`50 p 10  bA 25 pA 5  diffy  FALSE  diffpop  FALSE`$biasllke +
  labs(title = expression(paste("LLKE  ", X == Y, "  ", P == Q))) +
  theme(plot.title = element_text(size = 9, hjust = 0.5))

b3 <- biasplots$`50 p 10  bA 25 pA 5  diffy  TRUE  diffpop  FALSE`$biasirtke +
  labs(title = expression(paste("IRTKE  ", X != Y, "  ", P == Q))) +
  theme(plot.title = element_text(size = 9, hjust = 0.5))
b4 <- biasplots$`50 p 10  bA 25 pA 5  diffy  TRUE  diffpop  FALSE`$biasllke +
  labs(title = expression(paste("LLKE  ", X != Y, "  ", P == Q))) +
  theme(plot.title = element_text(size = 9, hjust = 0.5))

b5 <- biasplots$`50 p 10  bA 25 pA 5  diffy  FALSE  diffpop  TRUE`$biasirtke +
  labs(title = expression(paste("IRTKE  ", X == Y, "  ", P != Q))) +
  theme(plot.title = element_text(size = 9, hjust = 0.5))
b6 <- biasplots$`50 p 10  bA 25 pA 5  diffy  FALSE  diffpop  TRUE`$biasllke +
  labs(title = expression(paste("LLKE  ", X == Y, "  ", P != Q))) +
  theme(plot.title = element_text(size = 9, hjust = 0.5))

b7 <- biasplots$`50 p 10  bA 25 pA 5  diffy  TRUE  diffpop  TRUE`$biasirtke +
  labs(title = expression(paste("IRTKE  ", X != Y, "  ", P != Q))) +
  theme(plot.title = element_text(size = 9, hjust = 0.5))
b8 <- biasplots$`50 p 10  bA 25 pA 5  diffy  TRUE  diffpop  TRUE`$biasllke +
  labs(title = expression(paste("LLKE  ", X != Y, "  ", P != Q))) +
  theme(plot.title = element_text(size = 9, hjust = 0.5))

ggarrange(b1, b2, b3, b4, b5, b6, b7, b8,
  nrow = 4, ncol = 2, common.legend = T, legend = "bottom"
)

s1 <- seeplots$`50 p 10  bA 25 pA 5  diffy  FALSE  diffpop  FALSE`$seeirtke +
  labs(title = expression(paste("IRTKE  ", X == Y, "  ", P == Q))) +
  theme(plot.title = element_text(size = 9, hjust = 0.5))
s2 <- seeplots$`50 p 10  bA 25 pA 5  diffy  FALSE  diffpop  FALSE`$seellke +
  labs(title = expression(paste("LLKE  ", X == Y, "  ", P == Q))) +
  theme(plot.title = element_text(size = 9, hjust = 0.5))

s3 <- seeplots$`50 p 10  bA 25 pA 5  diffy  TRUE  diffpop  FALSE`$seeirtke +
  labs(title = expression(paste("IRTKE  ", X != Y, "  ", P == Q))) +
  theme(plot.title = element_text(size = 9, hjust = 0.5))
s4 <- seeplots$`50 p 10  bA 25 pA 5  diffy  TRUE  diffpop  FALSE`$seellke +
  labs(title = expression(paste("LLKE  ", X != Y, "  ", P == Q))) +
  theme(plot.title = element_text(size = 9, hjust = 0.5))

s5 <- seeplots$`50 p 10  bA 25 pA 5  diffy  FALSE  diffpop  TRUE`$seeirtke +
  labs(title = expression(paste("IRTKE  ", X == Y, "  ", P != Q))) +
  theme(plot.title = element_text(size = 9, hjust = 0.5))
s6 <- seeplots$`50 p 10  bA 25 pA 5  diffy  FALSE  diffpop  TRUE`$seellke +
  labs(title = expression(paste("LLKE  ", X == Y, "  ", P != Q))) +
  theme(plot.title = element_text(size = 9, hjust = 0.5))

s7 <- seeplots$`50 p 10  bA 25 pA 5  diffy  TRUE  diffpop  TRUE`$seeirtke +
  labs(title = expression(paste("IRTKE  ", X != Y, "  ", P != Q))) +
  theme(plot.title = element_text(size = 9, hjust = 0.5))
s8 <- seeplots$`50 p 10  bA 25 pA 5  diffy  TRUE  diffpop  TRUE`$seellke +
  labs(title = expression(paste("LLKE  ", X != Y, "  ", P != Q))) +
  theme(plot.title = element_text(size = 9, hjust = 0.5))

ggarrange(s1, s2, s3, s4, s5, s6, s7, s8,
  nrow = 4, ncol = 2, common.legend = T, legend = "bottom"
)

ggarrange(b1, b2, b5, b6, s1, s2, s5, s6, 
          nrow = 4, ncol = 2, common.legend = T, legend = "bottom",
          labels = c("a", "b", "c", "d", "e", "f", "g", "h")
)


# IRT ---------------------------------------------------------------------
n <- 1500
R <- 1000
bin_items <- c(50)
poly_items <- c(10)
bin_items_A <- c(25)
poly_items_A <- c(5)
harder_y <- c(F, T)
diff_pop <- c(F, T)

maxs <- 80
no_scen <- length(harder_y) * length(diff_pop) * length(bin_items) * length(bin_items_A)
plots <- vector(mode = "list", length = 0)
rowi <- 1
biasplots <- vector(mode = "list", length = 0)
seeplots <- vector(mode = "list", length = 0)
# loop through scenarios
for (item_scen_A in 1:length(bin_items_A)) {
  for (item_scen in 1:length(bin_items)) {
    for (pop_scen in 1:length(diff_pop)) {
      for (test_scen in harder_y) {
        no_bin <- bin_items[item_scen]
        no_poly <- poly_items[item_scen]
        no_bin_A <- bin_items_A[item_scen_A]
        no_poly_A <- poly_items_A[item_scen_A]
        filename_err <- paste("Data/irt/irt Error R", R, " n", n, " b", no_bin, "p", no_poly,
          " bA", no_bin_A, "pA", no_poly_A,
          " diffy ", test_scen, " diffpop ", diff_pop[pop_scen], ".RData",
          sep = ""
        )
        load(filename_err)

        score <- rep(0:maxs, 6)
        method <- c(
          rep("IRTKE EG", maxs + 1),
          rep("LLKE EG", maxs + 1),
          rep("IRTKE CE", maxs + 1),
          rep("LLKE CE", maxs + 1),
          rep("IRTKE PSE", maxs + 1),
          rep("LLKE PSE", maxs + 1)
        )
        bias <- c(
          res$EEtrue$local$IRTKE$bias,
          res$EEtrue$local$KE$bias,
          res$EEtrue$local$IRTKECE$bias,
          res$EEtrue$local$KECE$bias,
          res$EEtrue$local$IRTKEPSE$bias,
          res$EEtrue$local$KEPSE$bias
        )
        bias_df <- data.frame(score, method, value = bias)
        SEEs <- c(
          res$EEtrue$local$IRTKE$SEE,
          res$EEtrue$local$KE$SEE,
          res$EEtrue$local$IRTKECE$SEE,
          res$EEtrue$local$KECE$SEE,
          res$EEtrue$local$IRTKEPSE$SEE,
          res$EEtrue$local$KEPSE$SEE
        )
        SEE_df <- data.frame(score, method, value = SEEs)

        SEE_df <- data.frame(score, method, value = SEEs)
        bias_df$method <- factor(bias_df$method, levels = c("IRTKE EG", "LLKE EG", "IRTKE CE", "LLKE CE", "IRTKE PSE", "LLKE PSE"))
        SEE_df$method <- factor(SEE_df$method, levels = c("IRTKE EG", "LLKE EG", "IRTKE CE", "LLKE CE", "IRTKE PSE", "LLKE PSE"))

        y_max <- max(bias_df[!grepl("EG", bias_df$method), ]$value) + 0.05 * max(bias_df[!grepl("EG", bias_df$method), ]$value)
        y_min <- min(bias_df[!grepl("EG", bias_df$method), ]$value) + 0.05 * min(bias_df[!grepl("EG", bias_df$method), ]$value)
        y_max <- max(bias_df$value) + 0.05 * max(bias_df$value)
        y_min <- min(bias_df$value) + 0.05 * min(bias_df$value)

        pbiasirtke <- xScorePlot(bias_df[!grepl("LLKE", bias_df$method), ], y_min = y_min, y_max = y_max, "Bias")
        pbiasllke <- xScorePlot(bias_df[!grepl("IRTKE", bias_df$method), ], y_min = y_min, y_max = y_max, "Bias")

        y_max <- max(SEE_df$value) + 0.05 * max(SEE_df$value)
        pseeirtke <- xScorePlot(SEE_df[!grepl("LLKE", SEE_df$method), ], y_min = 0, y_max = y_max, "SEE")
        pseellke <- xScorePlot(SEE_df[!grepl("IRTKE", SEE_df$method), ], y_min = 0, y_max = y_max, "SEE")

        scenario <- paste(no_bin, "p", no_poly, " bA", no_bin_A, "pA", no_poly_A, " diffy ", test_scen, " diffpop ", diff_pop[pop_scen])
        biasplots[[scenario]] <- list(biasirtke = pbiasirtke, biasllke = pbiasllke)
        seeplots[[scenario]] <- list(seeirtke = pseeirtke, seellke = pseellke)
      }
    }
  }
}

b1 <- biasplots$`50 p 10  bA 25 pA 5  diffy  FALSE  diffpop  FALSE`$biasirtke +
  labs(title = expression(paste("IRTKE  ", X == Y, "  ", P == Q))) +
  theme(plot.title = element_text(size = 9, hjust = 0.5))
b2 <- biasplots$`50 p 10  bA 25 pA 5  diffy  FALSE  diffpop  FALSE`$biasllke +
  labs(title = expression(paste("LLKE  ", X == Y, "  ", P == Q))) +
  theme(plot.title = element_text(size = 9, hjust = 0.5))
b3 <- biasplots$`50 p 10  bA 25 pA 5  diffy  TRUE  diffpop  FALSE`$biasirtke +
  labs(title = expression(paste("IRTKE  ", X != Y, "  ", P == Q))) +
  theme(plot.title = element_text(size = 9, hjust = 0.5))
b4 <- biasplots$`50 p 10  bA 25 pA 5  diffy  TRUE  diffpop  FALSE`$biasllke +
  labs(title = expression(paste("LLKE  ", X != Y, "  ", P == Q))) +
  theme(plot.title = element_text(size = 9, hjust = 0.5))
b5 <- biasplots$`50 p 10  bA 25 pA 5  diffy  FALSE  diffpop  TRUE`$biasirtke +
  labs(title = expression(paste("IRTKE  ", X == Y, "  ", P != Q))) +
  theme(plot.title = element_text(size = 9, hjust = 0.5))
b6 <- biasplots$`50 p 10  bA 25 pA 5  diffy  FALSE  diffpop  TRUE`$biasllke +
  labs(title = expression(paste("LLKE  ", X == Y, "  ", P != Q))) +
  theme(plot.title = element_text(size = 9, hjust = 0.5))
b7 <- biasplots$`50 p 10  bA 25 pA 5  diffy  TRUE  diffpop  TRUE`$biasirtke +
  labs(title = expression(paste("IRTKE  ", X != Y, "  ", P != Q))) +
  theme(plot.title = element_text(size = 9, hjust = 0.5))
b8 <- biasplots$`50 p 10  bA 25 pA 5  diffy  TRUE  diffpop  TRUE`$biasllke +
  labs(title = expression(paste("LLKE  ", X != Y, "  ", P != Q))) +
  theme(plot.title = element_text(size = 9, hjust = 0.5))

ggarrange(b1, b2, b3, b4, b5, b6, b7, b8,
  nrow = 4, ncol = 2, common.legend = T, legend = "bottom"
)

s1 <- seeplots$`50 p 10  bA 25 pA 5  diffy  FALSE  diffpop  FALSE`$seeirtke +
  labs(title = expression(paste("IRTKE  ", X == Y, "  ", P == Q))) +
  theme(plot.title = element_text(size = 9, hjust = 0.5))
s2 <- seeplots$`50 p 10  bA 25 pA 5  diffy  FALSE  diffpop  FALSE`$seellke +
  labs(title = expression(paste("LLKE  ", X == Y, "  ", P == Q))) +
  theme(plot.title = element_text(size = 9, hjust = 0.5))
s3 <- seeplots$`50 p 10  bA 25 pA 5  diffy  TRUE  diffpop  FALSE`$seeirtke +
  labs(title = expression(paste("IRTKE  ", X != Y, "  ", P == Q))) +
  theme(plot.title = element_text(size = 9, hjust = 0.5))
s4 <- seeplots$`50 p 10  bA 25 pA 5  diffy  TRUE  diffpop  FALSE`$seellke +
  labs(title = expression(paste("LLKE  ", X != Y, "  ", P == Q))) +
  theme(plot.title = element_text(size = 9, hjust = 0.5))
s5 <- seeplots$`50 p 10  bA 25 pA 5  diffy  FALSE  diffpop  TRUE`$seeirtke +
  labs(title = expression(paste("IRTKE  ", X == Y, "  ", P != Q))) +
  theme(plot.title = element_text(size = 9, hjust = 0.5))
s6 <- seeplots$`50 p 10  bA 25 pA 5  diffy  FALSE  diffpop  TRUE`$seellke +
  labs(title = expression(paste("LLKE  ", X == Y, "  ", P != Q))) +
  theme(plot.title = element_text(size = 9, hjust = 0.5))
s7 <- seeplots$`50 p 10  bA 25 pA 5  diffy  TRUE  diffpop  TRUE`$seeirtke +
  labs(title = expression(paste("IRTKE  ", X != Y, "  ", P != Q))) +
  theme(plot.title = element_text(size = 9, hjust = 0.5))
s8 <- seeplots$`50 p 10  bA 25 pA 5  diffy  TRUE  diffpop  TRUE`$seellke +
  labs(title = expression(paste("LLKE  ", X != Y, "  ", P != Q))) +
  theme(plot.title = element_text(size = 9, hjust = 0.5))

ggarrange(s1, s2, s3, s4, s5, s6, s7, s8,
  nrow = 4, ncol = 2, common.legend = T, legend = "bottom"
)

