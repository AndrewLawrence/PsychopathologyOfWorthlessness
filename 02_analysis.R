library(tidyverse)
library(cowplot)
library(bootnet)
library(NetworkComparisonTest)


# Load objects ------------------------------------------------------------

items <- read.csv("input/item_info.csv")

load("model_objects/s1s2_input.RData")
# [1] "study1" "study2"

load("model_objects/s1_bootstrap_D14_Ising_k10000_0.25.RData")
# [1] "bootD14_s1"

load("model_objects/s2_bootstrap_D14_Ising_k10000_0.25.RData")
# [1] "bootD14_s2"

load("model_objects/nctD14_s1s2_bootstrap_SCL90_Ising_k10000_0.25.RData")
# [1] "nctD14_s1s2"


# Supplementary Figure 3 --------------------------------------------------

# obtain phi coefficients:
cormats <- list(study1 = cor(study1),
                study2 = cor(study2))

# Order matrix elements by dendrogram:
matorder <- bootD14_s1$sample$graph %>%
  dist() %>%
  hclust() %>%
  as.dendrogram() %>%
  order.dendrogram()

matproc <- function(x) {
  # consistent format/names for coefficient matrices:
  diag(x) <- NA
  dimnames(x) <- list(items$tlname[matorder],
                      items$tlname[matorder])
  return(x)
}

# gather all matrix data:
plotmats <- list(
  # study 1 correlation:
  s1_cor = cormats$study1[matorder, matorder] %>% matproc(),
  # study 1 eLasso coefs:
  s1_bnet = bootD14_s1$sample$graph[matorder, matorder] %>% matproc(),
  # study 2 correlation:
  s2_cor = cormats$study2[matorder, matorder] %>% matproc(),
  # study 2 eLasso coefs:
  s2_bnet = bootD14_s2$sample$graph[matorder, matorder] %>% matproc()
)

bplotter <- function(mat, smin = -1, smax = 1) {
  df <- reshape2::melt(mat)
  df %>%
    ggplot(aes(y = Var1, x = Var2, fill = value)) +
    geom_tile() +
    scale_fill_gradient2(midpoint = 0, limits = round(c(smin, smax), 2),
                         low = "blue", mid = "white", high = "red",
                         na.value = "white", name = "",
                         oob = scales::squish,
                         breaks = seq(from = smin,
                                      to = smax,
                                      length.out = 7) %>% round(2),
    ) +
    theme_minimal() +
    theme(axis.title = element_blank())
}

# calculate a single plotting range over both studies:
cor_range <- c(plotmats$s1_cor, plotmats$s2_cor) %>% abs() %>% max(na.rm = TRUE)
# for the networks use the 99th percentile of data:
bnet_range <- 1.8

# the plots:
bplot1 <- bplotter(plotmats[[1]], smin = -cor_range, smax = cor_range)
bplot2 <- bplotter(plotmats[[2]], smin = -bnet_range, smax = bnet_range)
bplot3 <- bplotter(plotmats[[3]], smin = -cor_range, smax = cor_range)
bplot4 <- bplotter(plotmats[[4]], smin = -bnet_range, smax = bnet_range)

# save figure:
cowplot::save_plot(
  "output/SupplementaryFigure3.png",
  cowplot::plot_grid(
    bplot1,
    bplot2,
    bplot3,
    bplot4,
    labels = "AUTO",
    label_size = 14,
    scale = 0.98,
    hjust = 0.0
  ),
  base_height = 10 / 1.5,
  base_width = 14 / 1.5,
  dpi = 600
)

# Figure 1 ----------------------------------------------------------------

# create node positions averaged over both studies:
alayout <- qgraph::averageLayout(bootD14_s1$boots, bootD14_s2$boots)

# Plotting function:
plotnet <- function(net, title = "") {
  bootnet:::plot.bootnetResult(net,
                               title = title, title.cex = 2,
                               layout = alayout,
                               color = items$nodecol,
                               labels = items$tlname,
                               label.norm = "OOOO", label.cex = 2,
                               cut = 0, minimum = 0, maximum = 2.53,
                               details = FALSE, mar = c(2,2,2,2))
}


# A combined plot:
tiff(filename = "output/Figure1.tiff",
    width = 17/2, height = 6/2, units = "in", res = 600)
graphics::layout(mat = t(1:3), widths = t(c(2,2,1.5)))
# Study 1:
plotnet(bootD14_s1$sample, "a")
# Study 2:
plotnet(bootD14_s2$sample, "b")

# A key linking two-letter codes to questions:
par(mar = c(1,1,1,1))

leg.sc <- items$tlname
leg.lc <- paste0(items$idx, ": ",
                 gsub("distressed by ", "", items$fulltext))

plot(NULL, xaxt = 'n', yaxt = 'n', bty = 'n',
     ylab = '', xlab = '', xlim = 0:1, ylim = 0:1)
legend(x = 0, y = 0.5, bty = "n", cex = 0.9,
       xjust = 0, yjust = 0.5,
       legend = leg.sc[order(leg.sc)],
       text.font = 2, y.intersp = 1.2)
legend(x = 0.12, y = 0.5, bty = "n", cex = 0.9,
       xjust = 0, yjust = 0.5,
       legend = leg.lc[order(leg.sc)],
       text.font = 1, y.intersp = 1.2)

dev.off()

# revert graphics options to defaults:
par(mar = c(5.1, 4.1, 4.1, 2.1))
graphics::layout(mat = as.matrix(1))


# Hypothesis testing ------------------------------------------------------

# This function extracts and formats the worthlessness hypothesis:
extract_worthlessness_hypothesis <- function(x) {
  # x is a bootstrapped D14 network.
  RLH <- x$graph["worthlessn", c("selfblame", "hopeless", "guilt")]
  RPA <- x$graph["worthlessn", c("losssexint", "lowenergy", "nointerest")]
  return(as.data.frame(t(unlist(list(RLH = RLH,
                                     RPA = RPA,
                                     sumRLH = sum(RLH),
                                     sumRPA = sum(RPA),
                                     diff =  sum(RLH) - sum(RPA))))))
}

hypo_summary_table <- function(obs,
                               boots,
                               alpha = 0.05) {
  # obs is a n*(1+1) of quantities to summarise columns are Measure, Oserved
  # boots is a nboots * n quantities
  # alpha is the complement of the confidence interval width.

  loq <- alpha / 2
  upq <- 1 - loq

  # check input:
  stopifnot(identical(names(boots), obs$Measure))

  # functions for quantiles:
  qupper <- function(x) quantile(x, upq, type = 6)
  qlower <- function(x) quantile(x, loq, type = 6)

  bresult <- boots %>%
    summarise_all(.funs = list(upper = qupper,
                               lower = qlower)) %>%
    gather() %>%
    separate(key, into = c("Measure", "Metric"), sep = "_") %>%
    mutate(Metric = factor(Metric, levels = c("lower", "upper"))) %>%
    spread(key = Metric, value = value) %>%
    left_join(obs, by = "Measure") %>%
    mutate(excludeszero = ifelse(sign(upper) == sign(lower), "*", "")) %>%
    select(Measure, Observed, excludeszero, lower, upper)
  return(bresult)
}


# ~ Table 1(a) ------------------------------------------------------------
# Calculates the difference in the sum of RLH and RPA edges

# Extract the observed hypothesis results:
# study 1
hypo_obs_s1 <- extract_worthlessness_hypothesis(bootD14_s1$sample) %>%
  t() %>%  as.data.frame() %>% setNames("Observed") %>%
  rownames_to_column("Measure")
# study 2
hypo_obs_s2 <- extract_worthlessness_hypothesis(bootD14_s2$sample) %>%
  t() %>%  as.data.frame() %>% setNames("Observed") %>%
  rownames_to_column("Measure")

# Extract the bootstrapped hypothesis results:
hypo_boot_s1 <- map_dfr(bootD14_s1$boots, extract_worthlessness_hypothesis)
hypo_boot_s2 <- map_dfr(bootD14_s2$boots, extract_worthlessness_hypothesis)

# Summary tables for each study, corrected at alpha=0.05/2:
hypo_tab_s1 <- hypo_summary_table(hypo_obs_s1, hypo_boot_s1, alpha = 0.05/2)
hypo_tab_s2 <- hypo_summary_table(hypo_obs_s2, hypo_boot_s2, alpha = 0.05/2)

# Combine for Table1 part a:
tab1_a.rows <- c(1,8,9,2:7)
tab1_a <- cbind(hypo_tab_s1[tab1_a.rows,],
              hypo_tab_s2[tab1_a.rows,-1])


write.csv(tab1_a, file = "output/Table1_parta.csv", row.names = FALSE)

# ~ Table 1(b) ------------------------------------------------------------
# Contains all pair-wise difference tests
#   between the 6 edges involved in the hypothesis

# function to format input with brackets:
fmt_bsci <- function(x, lower, upper, digits = 3) {
  # input are vectors to format "x [lower, upper]"
  fmtstr <- paste0("%.",digits,"f")
  paste0(x %>% sprintf(fmt = fmtstr, .),
         " [",
         lower %>% sprintf(fmt = fmtstr, .),
         ", ",
         upper %>% sprintf(fmt = fmtstr, .),
         "]")
}

# obtain all unique pairwise combinations of x in a data.frame format:
make_allpairs <- function(x,
                              lower = TRUE,
                              sep = "_",
                              names = c("A", "B")) {
  # given vectors (x,y) generate all pairwise combinations and return a
  #   data.frame of either the lower or upper triangular orderings
  mat <- outer(x,x,function(a,b) { paste(a, b, sep = sep) } )
  tri <- ifelse(lower, "lower.tri", "upper.tri")
  vec <- mat[do.call(tri, list(mat))]
  return(
    data.frame(label = vec, stringsAsFactors = FALSE) %>%
      separate(col = label,
               into = names,
               sep = sep,
               remove = FALSE)
  )
}

hypo_pairwise.edges <- make_allpairs(
  c(
    "selfblame",
    "hopeless",
    "guilt",
    "losssexint",
    "lowenergy",
    "nointerest"
  ),
  names = c("x2", "y2")
) %>%
  mutate(y1 = "worthlessn",
         x1 = "worthlessn") %>%
  mutate(
    x2_set = ifelse(x2 %in% c("selfblame", "hopeless", "guilt"),
                    "RLH",
                    "RPA"),
    y2_set = ifelse(y2 %in% c("selfblame", "hopeless", "guilt"),
                    "RLH",
                    "RPA"),
    set = paste(x2_set, y2_set, sep = "_"),
    Type = ifelse(x2_set == y2_set,
                  yes = "Within",
                  no = "Between")
  ) %>%
  arrange(Type, set) %>%
  select(label, x1, x2, y1, y2, everything()) %>%
  # on reflection remove some of the info:
  filter(Type == "Between") %>%
  select(-Type, -x2_set, -y2_set, -set)

# The result:
#                  label         x1         x2         y1        y2
# 1 losssexint_selfblame worthlessn losssexint worthlessn selfblame
# 2  lowenergy_selfblame worthlessn  lowenergy worthlessn selfblame
# 3 nointerest_selfblame worthlessn nointerest worthlessn selfblame
# 4  losssexint_hopeless worthlessn losssexint worthlessn  hopeless
# 5   lowenergy_hopeless worthlessn  lowenergy worthlessn  hopeless
# 6  nointerest_hopeless worthlessn nointerest worthlessn  hopeless
# 7     losssexint_guilt worthlessn losssexint worthlessn     guilt
# 8      lowenergy_guilt worthlessn  lowenergy worthlessn     guilt
# 9     nointerest_guilt worthlessn nointerest worthlessn     guilt


# manipulate some tab1_a data to merge into tab1_b:
tab1_a_merge <- data.frame(
  Measure = tab1_a$Measure,
  raw_s1 = tab1_a[,2],
  raw_s2 = tab1_a[,6],
  study1 = fmt_bsci(tab1_a[,2], tab1_a[,4], tab1_a[,5], 3),
  study2 = fmt_bsci(tab1_a[,6], tab1_a[,8], tab1_a[,9], 3),
  stringsAsFactors = FALSE)[4:9, ]
# add merging label:
tab1_a_merge$lab <- c("guilt", "hopeless", "selfblame",
                      "losssexint", "lowenergy", "nointerest")

# do the merge:
hypo_pairwise.edges <- hypo_pairwise.edges %>%
  left_join(., tab1_a_merge %>%
              select(-Measure) %>%
              setNames(c("s1rawx", "s2rawx", "s1_x", "s2_x", "lab")),
            by = c("x2" = "lab")) %>%
  left_join(., tab1_a_merge %>%
              select(-Measure) %>%
              setNames(c("s1rawy", "s2rawy", "s1_y", "s2_y", "lab")),
            by = c("y2" = "lab")) %>%
  mutate(s1_yxdiff = s1rawy - s1rawx,
         s2_yxdiff = s2rawy - s2rawx) %>%
  select(1:5, s1_x, s1_y, s1_yxdiff, s2_x, s2_y, s2_yxdiff)


# function takes a bootnet and df containing edges and adds the bootnet
#   difference results to the df
get_pairwise_bootnetdiff <- function(edgedf,
                                     bootnetobj,
                                     alpha = 0.025) {
  # edgedf has a 'label' and
  #   4 node names making 2 cxns:
  #     connection x comprising nodes: 'x1', 'x2'
  #     connection y comprising ndoes: 'y1', 'y2'
  result <- map_dfr(setNames(1:nrow(edgedf),
                             edgedf$label),
                    function(i) {
                      bootnet::differenceTest(bootnetobj,
                                              x = edgedf$x1[[i]],
                                              x2 = edgedf$x2[[i]],
                                              y = edgedf$y1[[i]],
                                              y2 = edgedf$y2[[i]],
                                              measure = "edge",
                                              alpha = alpha,
                                              verbose = TRUE) } )
  return(result)
}

merge_pairwise_bootnetdiff <- function(edgedf,
                                       pwbnd_s1,
                                       pwbnd_s2) {
  pwbnd_s1 <- pwbnd_s1 %>%
    select(lower, upper, significant) %>%
    setNames(c("s1_yxdiff_lower", "s1_yxdiff_upper", "s1_yxdiff_sig"))

  pwbnd_s2 <- pwbnd_s2 %>%
    select(lower, upper, significant) %>%
    setNames(c("s2_yxdiff_lower", "s2_yxdiff_upper", "s2_yxdiff_sig"))

  result <- cbind(edgedf, pwbnd_s1, pwbnd_s2)
  result %>% select(label, x1, x2, y1, y2,
                    s1_x, s1_y,
                    s1_yxdiff, s1_yxdiff_lower, s1_yxdiff_upper, s1_yxdiff_sig,
                    s2_x, s2_y,
                    s2_yxdiff, s2_yxdiff_lower, s2_yxdiff_upper, s2_yxdiff_sig)
}

gu_pairwise_bootnetdiff <- function(edgedf,
                                    bootnetobj_s1,
                                    bootnetobj_s2,
                                    alpha = 0.025) {
  s1 <- get_pairwise_bootnetdiff(edgedf, bootnetobj_s1, alpha = alpha)
  s2 <- get_pairwise_bootnetdiff(edgedf, bootnetobj_s2, alpha = alpha)
  result <- merge_pairwise_bootnetdiff(edgedf, s1, s2)
  # additional formatting:
  result$s1_yxdiff_fmt <- fmt_bsci(x = result$s1_yxdiff,
                                   lower = result$s1_yxdiff_lower,
                                   upper = result$s1_yxdiff_upper,
                                   digits = 3)

  result$s2_yxdiff_fmt <- fmt_bsci(x = result$s2_yxdiff,
                                   lower = result$s2_yxdiff_lower,
                                   upper = result$s2_yxdiff_upper,
                                   digits = 3)
  result %>% select(-ends_with("upper"), -ends_with("lower"))
}

# Run the pairwise analysis:
hypo_pairwise <- gu_pairwise_bootnetdiff(hypo_pairwise.edges,
                                            bootD14_s1,
                                            bootD14_s2,
                                            alpha = 0.05/18)

# Extract for part b of the table:
tab1_b <- hypo_pairwise %>%
  transmute(label = paste0(y2, "-", x2),
            study1 = paste0(s1_yxdiff_fmt, ifelse(s1_yxdiff_sig, "*", " ")),
            study2 = paste0(s2_yxdiff_fmt, ifelse(s2_yxdiff_sig, "*", " ")))

write.csv(tab1_b, file = "output/Table1_partb.csv", row.names = FALSE)


# Network Comparison ------------------------------------------------------
#   permutation tests of differences between study1 and study2 networks

summary(nctD14_s1s2)

plot(nctD14_s1s2)

## corrected edge-wise invariance p-values:
nctD14_s1s2$einv.pvals %>%
  mutate(p.adj = p.adjust(`p-value`, method = "holm")) %>%
  slice_min(p.adj, n = 1)

