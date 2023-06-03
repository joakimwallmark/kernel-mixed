library(gridExtra)
library(ggplot2)
library(ggpubr)

x_score_plot <- function(df, ylab, y_min = NA, y_max = NA, labels = c("EG", "CE", "PSE")) {
  res_plot <- ggplot(df, aes(x = score, y = value, linetype = method)) +
    geom_line() +
    theme_bw() +
    xlab("X score") +
    ylab(ylab) +
    scale_linetype_manual(values = c(1, 2, 3), labels = labels) +
    theme(
      legend.position = "bottom",
      legend.title = element_blank(),
      axis.title.x = element_blank()
    ) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(limits = c(y_min, y_max), expand = c(0, 0))
  return(res_plot)
}


# Splines ---------------------------------------------------------------------
load(file = "fitted_splines.RData")
n <- 1500
R <- 1000
bin_items <- c(50)
poly_items <- c(10)
bin_items_a <- c(25)
poly_items_a <- c(5)

harder_y <- c(FALSE, TRUE) # for filenames
diff_pop <- c(FALSE, TRUE) # for filenames

maxs <- 80
no_scen <- length(harder_y) * length(diff_pop) * length(bin_items) * length(bin_items_a)
biasplots <- vector(mode = "list", length = 0)
seeplots <- vector(mode = "list", length = 0)
# loop through scenarios
for (item_scen_A in seq_along(bin_items_a)) {
  for (item_scen in seq_along(bin_items)) {
    for (pop_scen in seq_along(diff_pop)) {
      for (test_scen in harder_y) {
        no_bin <- bin_items[item_scen]
        no_poly <- poly_items[item_scen]
        no_bin_a <- bin_items_a[item_scen_A]
        no_poly_a <- poly_items_a[item_scen_A]
        filename_err <- paste("data/spline/spline Error R", R, " n", n, " b", no_bin, "p", no_poly,
          " bA", no_bin_a, "pA", no_poly_a,
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
        bias_df$method <- factor(bias_df$method,
          levels = c("IRTKE EG", "LLKE EG", "IRTKE CE", "LLKE CE", "IRTKE PSE", "LLKE PSE")
        )
        see_df$method <- factor(see_df$method,
          levels = c("IRTKE EG", "LLKE EG", "IRTKE CE", "LLKE CE", "IRTKE PSE", "LLKE PSE")
        )


        y_max <- max(bias_df[!grepl("EG", bias_df$method), ]$value) +
          0.05 * max(bias_df[!grepl("EG", bias_df$method), ]$value)
        y_min <- min(bias_df[!grepl("EG", bias_df$method), ]$value) +
          0.05 * min(bias_df[!grepl("EG", bias_df$method), ]$value)
        y_max <- max(bias_df$value) + 0.05 * max(bias_df$value)
        y_min <- min(bias_df$value) + 0.05 * min(bias_df$value)

        pbiasirtke <- x_score_plot(bias_df[grepl("IRTKE ", bias_df$method), ], y_min = y_min, y_max = y_max, "Bias")
        pbiasllke <- x_score_plot(bias_df[grepl("LLKE", bias_df$method), ], y_min = y_min, y_max = y_max, "Bias")

        y_max <- max(see_df$value) + 0.05 * max(see_df$value)
        pseeirtke <- x_score_plot(see_df[grepl("IRTKE ", see_df$method), ], y_min = 0, y_max = y_max, "SEE")
        pseellke <- x_score_plot(see_df[grepl("LLKE", see_df$method), ], y_min = 0, y_max = y_max, "SEE")

        scenario <- paste(
          no_bin, "p", no_poly,
          " bA", no_bin_a, "pA", no_poly_a, " diffy ",
          test_scen, " diffpop ", diff_pop[pop_scen]
        )
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
  nrow = 4, ncol = 2, common.legend = TRUE, legend = "bottom"
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
  nrow = 4, ncol = 2, common.legend = TRUE, legend = "bottom"
)

ggarrange(b1, b2, b5, b6, s1, s2, s5, s6,
  nrow = 4, ncol = 2, common.legend = TRUE, legend = "bottom",
  labels = c("a", "b", "c", "d", "e", "f", "g", "h")
)


# IRT with IRTOSE incldued ---------------------------------------------------------------------
n <- 1500
R <- 1000
bin_items <- c(50)
poly_items <- c(10)
bin_items_a <- c(25)
poly_items_a <- c(5)
harder_y <- c(FALSE, TRUE)
diff_pop <- c(FALSE, TRUE)

maxs <- 80
no_scen <- length(harder_y) * length(diff_pop) * length(bin_items) * length(bin_items_a)
plots <- vector(mode = "list", length = 0)
rowi <- 1
biasplots <- vector(mode = "list", length = 0)
seeplots <- vector(mode = "list", length = 0)
# loop through scenarios
for (item_scen_a in seq_along(bin_items_a)) {
  for (item_scen in seq_along(bin_items)) {
    for (pop_scen in seq_along(diff_pop)) {
      for (test_scen in harder_y) {
        no_bin <- bin_items[item_scen]
        no_poly <- poly_items[item_scen]
        no_bin_a <- bin_items_a[item_scen_a]
        no_poly_a <- poly_items_a[item_scen_a]
        filename_err <- paste("Data/irt/irt Error R", R, " n", n, " b", no_bin, "p", no_poly,
                              " bA", no_bin_a, "pA", no_poly_a,
                              " diffy ", test_scen, " diffpop ", diff_pop[pop_scen], ".RData",
                              sep = ""
        )
        load(filename_err)
        
        score <- rep(0:maxs, 8)
        method <- c(
          rep("IRTOSE EG", maxs + 1),
          rep("IRTOSE NEAT", maxs + 1),
          rep("IRTKE EG", maxs + 1),
          rep("LLKE EG", maxs + 1),
          rep("IRTKE CE", maxs + 1),
          rep("LLKE CE", maxs + 1),
          rep("IRTKE PSE", maxs + 1),
          rep("LLKE PSE", maxs + 1)
        )
        bias <- c(
          res$EEtrue$local$IRTOSE$bias,
          res$EEtrue$local$IRTOSENEAT$bias,
          res$EEtrue$local$IRTKE$bias,
          res$EEtrue$local$KE$bias,
          res$EEtrue$local$IRTKECE$bias,
          res$EEtrue$local$KECE$bias,
          res$EEtrue$local$IRTKEPSE$bias,
          res$EEtrue$local$KEPSE$bias
        )
        bias_df <- data.frame(score, method, value = bias)
        sees <- c(
          res$EEtrue$local$IRTOSE$SEE,
          res$EEtrue$local$IRTOSENEAT$SEE,
          res$EEtrue$local$IRTKE$SEE,
          res$EEtrue$local$KE$SEE,
          res$EEtrue$local$IRTKECE$SEE,
          res$EEtrue$local$KECE$SEE,
          res$EEtrue$local$IRTKEPSE$SEE,
          res$EEtrue$local$KEPSE$SEE
        )
        see_df <- data.frame(score, method, value = sees)
        bias_df$method <- factor(
          bias_df$method,
          levels = c("IRTOSE EG", "IRTOSE NEAT", "IRTKE EG", "LLKE EG", "IRTKE CE", "LLKE CE", "IRTKE PSE", "LLKE PSE")
        )
        see_df$method <- factor(
          see_df$method,
          levels = c("IRTOSE EG", "IRTOSE NEAT", "IRTKE EG", "LLKE EG", "IRTKE CE", "LLKE CE", "IRTKE PSE", "LLKE PSE")
        )
        
        y_max <- max(bias_df[!grepl("EG", bias_df$method), ]$value) +
          0.05 * max(bias_df[!grepl("EG", bias_df$method), ]$value)
        y_min <- min(bias_df[!grepl("EG", bias_df$method), ]$value) +
          0.05 * min(bias_df[!grepl("EG", bias_df$method), ]$value)
        y_max <- max(bias_df$value) + 0.05 * max(bias_df$value)
        y_min <- min(bias_df$value) + 0.05 * min(bias_df$value)
        
        pbiasirtose <- x_score_plot(bias_df[grepl("IRTOSE", bias_df$method), ], y_min = y_min, y_max = y_max, "Bias", labels = c("EG", "NEAT"))
        pbiasirtke <- x_score_plot(bias_df[grepl("IRTKE", bias_df$method), ], y_min = y_min, y_max = y_max, "Bias")
        pbiasllke <- x_score_plot(bias_df[grepl("LLKE", bias_df$method), ], y_min = y_min, y_max = y_max, "Bias")
        
        y_max <- max(see_df$value) + 0.05 * max(see_df$value)
        pseeirtose <- x_score_plot(see_df[grepl("IRTOSE", see_df$method), ], y_min = 0, y_max = y_max, "SEE", labels = c("EG", "NEAT"))
        pseeirtke <- x_score_plot(see_df[grepl("IRTKE", see_df$method), ], y_min = 0, y_max = y_max, "SEE")
        pseellke <- x_score_plot(see_df[grepl("LLKE", see_df$method), ], y_min = 0, y_max = y_max, "SEE")
        
        scenario <- paste(
          no_bin, "p",
          no_poly, " bA", no_bin_a, "pA",
          no_poly_a, " diffy ",
          test_scen, " diffpop ", diff_pop[pop_scen]
        )
        biasplots[[scenario]] <- list(biasirtke = pbiasirtke, biasllke = pbiasllke, biasirtose = pbiasirtose)
        seeplots[[scenario]] <- list(seeirtke = pseeirtke, seellke = pseellke, seeirtose = pseeirtose)
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

b9 <- biasplots$`50 p 10  bA 25 pA 5  diffy  FALSE  diffpop  FALSE`$biasirtose +
  labs(title = expression(paste("IRTOSE  ", X == Y, "  ", P == Q))) +
  theme(plot.title = element_text(size = 9, hjust = 0.5))
b10 <- biasplots$`50 p 10  bA 25 pA 5  diffy  TRUE  diffpop  FALSE`$biasirtose +
  labs(title = expression(paste("IRTOSE  ", X != Y, "  ", P == Q))) +
  theme(plot.title = element_text(size = 9, hjust = 0.5))
b11 <- biasplots$`50 p 10  bA 25 pA 5  diffy  FALSE  diffpop  TRUE`$biasirtose +
  labs(title = expression(paste("IRTOSE  ", X == Y, "  ", P != Q))) +
  theme(plot.title = element_text(size = 9, hjust = 0.5))
b12 <- biasplots$`50 p 10  bA 25 pA 5  diffy  TRUE  diffpop  TRUE`$biasirtose +
  labs(title = expression(paste("IRTOSE  ", X != Y, "  ", P != Q))) +
  theme(plot.title = element_text(size = 9, hjust = 0.5))

ggarrange(b1, b2, b9, b3, b4, b10, b5, b6, b11, b7, b8, b12,
          nrow = 4, ncol = 3, common.legend = TRUE, legend = "bottom"
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

s9 <- seeplots$`50 p 10  bA 25 pA 5  diffy  FALSE  diffpop  FALSE`$seeirtose +
  labs(title = expression(paste("IRTOSE  ", X == Y, "  ", P == Q))) +
  theme(plot.title = element_text(size = 9, hjust = 0.5))
s10 <- seeplots$`50 p 10  bA 25 pA 5  diffy  TRUE  diffpop  FALSE`$seeirtose +
  labs(title = expression(paste("IRTOSE  ", X != Y, "  ", P == Q))) +
  theme(plot.title = element_text(size = 9, hjust = 0.5))
s11 <- seeplots$`50 p 10  bA 25 pA 5  diffy  FALSE  diffpop  TRUE`$seeirtose +
  labs(title = expression(paste("IRTOSE  ", X == Y, "  ", P != Q))) +
  theme(plot.title = element_text(size = 9, hjust = 0.5))
s12 <- seeplots$`50 p 10  bA 25 pA 5  diffy  TRUE  diffpop  TRUE`$seeirtose +
  labs(title = expression(paste("IRTOSE  ", X != Y, "  ", P != Q))) +
  theme(plot.title = element_text(size = 9, hjust = 0.5))

ggarrange(s1, s2, s9, s3, s4, s10, s5, s6, s11, s7, s8, s12,
          nrow = 4, ncol = 3, common.legend = TRUE, legend = "bottom"
)



# IRT ---------------------------------------------------------------------
n <- 1500
R <- 1000
bin_items <- c(50)
poly_items <- c(10)
bin_items_a <- c(25)
poly_items_a <- c(5)
harder_y <- c(FALSE, TRUE)
diff_pop <- c(FALSE, TRUE)

maxs <- 80
no_scen <- length(harder_y) * length(diff_pop) * length(bin_items) * length(bin_items_a)
plots <- vector(mode = "list", length = 0)
rowi <- 1
biasplots <- vector(mode = "list", length = 0)
seeplots <- vector(mode = "list", length = 0)
# loop through scenarios
for (item_scen_a in seq_along(bin_items_a)) {
  for (item_scen in seq_along(bin_items)) {
    for (pop_scen in seq_along(diff_pop)) {
      for (test_scen in harder_y) {
        no_bin <- bin_items[item_scen]
        no_poly <- poly_items[item_scen]
        no_bin_a <- bin_items_a[item_scen_a]
        no_poly_a <- poly_items_a[item_scen_a]
        filename_err <- paste("Data/irt/irt Error R", R, " n", n, " b", no_bin, "p", no_poly,
          " bA", no_bin_a, "pA", no_poly_a,
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
        sees <- c(
          res$EEtrue$local$IRTKE$SEE,
          res$EEtrue$local$KE$SEE,
          res$EEtrue$local$IRTKECE$SEE,
          res$EEtrue$local$KECE$SEE,
          res$EEtrue$local$IRTKEPSE$SEE,
          res$EEtrue$local$KEPSE$SEE
        )
        see_df <- data.frame(score, method, value = sees)
        bias_df$method <- factor(
          bias_df$method,
          levels = c("IRTKE EG", "LLKE EG", "IRTKE CE", "LLKE CE", "IRTKE PSE", "LLKE PSE")
        )
        see_df$method <- factor(
          see_df$method,
          levels = c("IRTKE EG", "LLKE EG", "IRTKE CE", "LLKE CE", "IRTKE PSE", "LLKE PSE")
        )

        y_max <- max(bias_df[!grepl("EG", bias_df$method), ]$value) +
          0.05 * max(bias_df[!grepl("EG", bias_df$method), ]$value)
        y_min <- min(bias_df[!grepl("EG", bias_df$method), ]$value) +
          0.05 * min(bias_df[!grepl("EG", bias_df$method), ]$value)
        y_max <- max(bias_df$value) + 0.05 * max(bias_df$value)
        y_min <- min(bias_df$value) + 0.05 * min(bias_df$value)

        pbiasirtke <- x_score_plot(bias_df[!grepl("LLKE", bias_df$method), ], y_min = y_min, y_max = y_max, "Bias")
        pbiasllke <- x_score_plot(bias_df[!grepl("IRTKE", bias_df$method), ], y_min = y_min, y_max = y_max, "Bias")

        y_max <- max(see_df$value) + 0.05 * max(see_df$value)
        pseeirtke <- x_score_plot(see_df[!grepl("LLKE", see_df$method), ], y_min = 0, y_max = y_max, "SEE")
        pseellke <- x_score_plot(see_df[!grepl("IRTKE", see_df$method), ], y_min = 0, y_max = y_max, "SEE")

        scenario <- paste(
          no_bin, "p",
          no_poly, " bA", no_bin_a, "pA",
          no_poly_a, " diffy ",
          test_scen, " diffpop ", diff_pop[pop_scen]
        )
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
  nrow = 4, ncol = 2, common.legend = TRUE, legend = "bottom"
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
  nrow = 4, ncol = 2, common.legend = TRUE, legend = "bottom"
)
