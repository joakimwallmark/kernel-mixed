library(mirt)
library(haven)
library(gridExtra)
library(ggplot2)
library(ggpubr)
library(gss)

# Note that the data files (NEATv13B2469.sav) cannot be shared for copyright reasons.
# Fil
# model scores by X score -------------------------------------------------
data1 <- read_sav("data/NEATv13B2469.sav")
data2 <- read_sav("data/NEATv14A2859.sav")
data1 <- as.data.frame(data1)
data2 <- as.data.frame(data2)
dim(data1)
dim(data2)
anchor_13 <- data1[, grep("UTP113BG", names(data1))] # verbal anchor items, 13B
anchor_14 <- data2[, grep("UTP114AG", names(data2))] # verbal anchor items, 13B
y_13_1 <- data1[, grep("VE113BG", names(data1))] # verbal 1 items, 13B
y_13_2 <- data1[, grep("VE213BG", names(data1))] # verbal 2 items, 13B
x_14_1 <- data2[, grep("VE114AG", names(data2))] # verbal 1 items, 14A
x_14_2 <- data2[, grep("VE214AG", names(data2))] # verbal 2 items, 14A
X_binform <- cbind(x_14_1, x_14_2)
Y_binform <- cbind(y_13_1, y_13_2)
XA_binform <- anchor_14

# Get 75% easiest items from Y and A13 forms
all_y_tests <- cbind(y_13_1, y_13_2, anchor_13)
easier_y <- sort(colMeans(all_y_tests))[31:120]
set.seed(28) # Easier Y
from_y_new <- sample(1:length(easier_y), size = 80) # 80 random of the easier ones from each form to new y
new_easier_y <- all_y_tests[names(easier_y[from_y_new])]
easy_y_tot <- rowSums(new_easier_y)
set.seed(231) # Easier A
from_a_new <- sample(1:length(easier_y), size = 40) # 40 random of the easier ones from each form to new y
new_easier_a <- all_y_tests[names(easier_y[from_a_new])]
easy_a_tot <- rowSums(new_easier_a)

x_tot <- rowSums(X_binform)
a_tot <- rowSums(XA_binform)
summary(x_tot)
summary(easy_y_tot)
summary(a_tot)
summary(easy_a_tot)
sd(x_tot)
sd(easy_y_tot)
sd(a_tot)
sd(easy_a_tot)

barplot(table(x_tot))
barplot(table(easy_y_tot))
barplot(table(a_tot))
barplot(table(easy_a_tot))

x_den <- ssden(~x_tot, domain = data.frame(x_tot = c(-0.5, 80.5)))
easy_y_den <- ssden(~easy_y_tot, domain = data.frame(easy_y_tot = c(-0.5, 80.5)))
a_den <- ssden(~a_tot, domain = data.frame(a_tot = c(-0.5, 40.5)))
easy_a_den <- ssden(~easy_a_tot, domain = data.frame(easy_a_tot = c(-0.5, 40.5)))

x_list <- list(x = x_den)
y_list <- list(y = x_den, easy_y = easy_y_den)
a_list <- list(a = a_den, easy_a = easy_a_den)
# save(x_list, y_list, a_list, file = "fitted_splines2.RData")
# load(file = "fitted_splines2.RData")

# Weights for data generation ---------------------------------------------
set.seed(24)
quantiles <- runif(100000) # generate responses from each density
totals_x <- round(qssden(x_list$x, p = quantiles))
totals_ey <- round(qssden(y_list$easy_y, p = quantiles))
totals_a <- round(qssden(a_list$a, p = quantiles))
totals_ea <- round(qssden(a_list$easy_a, p = quantiles))

hist(totals_x, breaks = 50)
hist(totals_ey, breaks = 50)
hist(totals_a, breaks = 50)
hist(totals_ea, breaks = 50)

# true probabilities
x_probs <- unname(colMeans(X_binform))
ey_probs <- unname(colMeans(new_easier_y))
a_probs <- unname(colMeans(anchor_14))
ea_probs <- unname(colMeans(new_easier_a))
# inital weights
init_wxy <- rep(0.5, 80)
init_wa <- rep(0.5, 40)

# Search for good weights
source("functions/optimize_fun.R")
opt_w_x <- optimize_fun(init_wxy, x_probs, totals_x)
opt_w_y <- opt_w_x
opt_w_a <- optimize_fun(init_wa, a_probs, totals_a)
opt_w_ey <- optimize_fun(init_wxy, ey_probs, totals_ey)
opt_w_ea <- optimize_fun(init_wa, ea_probs, totals_ea)
# save(opt_w_x, opt_w_y, opt_w_ey, opt_w_a, opt_w_ea, file = "weights2.RData")

# compare with real to decide if weights are good
sort(colMeans(get_item_resp_from_total(totals_x, opt_w_x)))-sort(x_probs)
sort(colMeans(get_item_resp_from_total(totals_ey, opt_w_ey)))-sort(ey_probs)
sort(colMeans(get_item_resp_from_total(totals_a, opt_w_a)))-sort(a_probs)
sort(colMeans(get_item_resp_from_total(totals_ea, opt_w_ea)))-sort(ea_probs)

# For plotting only ------------------------------------------------------------
df_x <- data.frame(score = x_tot)
df_e_y <- data.frame(score = easy_y_tot)
df_a <- data.frame(score = a_tot)
df_e_a <- data.frame(score = easy_a_tot)
x <- seq(-0.5, 80.5, 0.01)
a <- seq(-0.5, 40.5, 0.01)
xy_max <- 0.031
a_max <- 0.063
df_2x <- data.frame(x = x, y = sapply(x, function(x) { dssden(x_den, x) } ))
df_2e_y <- data.frame(x = x, y = sapply(x, function(x) { dssden(easy_y_den, x) } ))
df_2a <- data.frame(x = a, y = sapply(a, function(x) { dssden(a_den, x) } ))
df_2e_a <- data.frame(x = a, y = sapply(a, function(x) { dssden(easy_a_den, x) } ))

p1 <- ggplot(df_x, aes(x = score)) +
  geom_histogram(aes(y = ..density..), colour = "black", fill = "white", bins = 14) +
  geom_line(data = df_2x, mapping = aes(x, y)) +
  theme_bw() +
  xlab("score") +
  labs(title = "X/Y") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_continuous(limits = c(min(x), max(x)),
                     expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, xy_max), expand = c(0, 0, 0.05, 0))
p2 <- ggplot(df_e_y, aes(x = score)) +
  geom_histogram(aes(y = ..density..), colour = "black", fill = "white", bins = 14) +
  geom_line(data = df_2e_y, mapping = aes(x, y)) +
  theme_bw() +
  xlab("score") +
  labs(title = "Less difficult Y") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_continuous(limits = c(min(x), max(x)),
                     expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, xy_max), expand = c(0, 0, 0.05, 0))
p3 <- ggplot(df_a, aes(x = score)) +
  geom_histogram(aes(y = ..density..), colour = "black", fill = "white", bins = 14) +
  geom_line(data = df_2a, mapping = aes(x, y)) +
  theme_bw() +
  xlab("score") +
  labs(title = "A") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_continuous(limits = c(min(a), max(a)),
                     expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, a_max), expand = c(0, 0, 0.05, 0))
p4 <- ggplot(df_e_a, aes(x = score)) +
  geom_histogram(aes(y = ..density..), colour = "black", fill = "white", bins = 14) +
  geom_line(data = df_2e_a, mapping = aes(x, y)) +
  theme_bw() +
  xlab("score") +
  labs(title = "More able group A") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_continuous(limits = c(min(a), max(a)),
                     expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, a_max), expand = c(0, 0, 0.05, 0))

ggarrange(p1, p2, p3, p4, ncol = 2, nrow = 2, 
          labels = c("a)", "b)", "c)", "d)"), hjust = -1.2)

