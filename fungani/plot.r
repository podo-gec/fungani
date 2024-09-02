library(ggplot2)
library(patchwork)

args <- commandArgs(trailingOnly = TRUE)
# TODO: check that on Windows this points to the same directory as Python
homedir <- path.expand("~")
outfile <- file.path(homedir, "fungani.pdf")
file_fwd <- args[1]
file_rev <- args[2]
# FIXME: use Fasta filename?
strain_ref <- "A"
strain_test <- "B"

d1 <- read.csv(file_fwd, header = FALSE)
nsample <- nrow(d1)
dd <- subset(d1, d1$V2 > 0.8)
NI <- cut(dd$V2,
  breaks = c(
    0.8, 0.85, 0.9, 0.95, 0.96,
    0.97, 0.98, 0.99, 0.995, 0.999, 1
  ),
  include.lowest = TRUE, labels = c(
    "0.80 ≤ % ≤ 0.85",
    "0.85 < % ≤ 0.90",
    "0.90 < % ≤ 0.95",
    "0.95 < % ≤ 0.96",
    "0.96 < % ≤ 0.97",
    "0.97 < % ≤ 0.98",
    "0.98 < % ≤ 0.99",
    "0.99 < % ≤ 0.995",
    "0.995 < % ≤ 0.999",
    "0.999 < % ≤ 1"
  )
)
dx <- as.data.frame(table(NI))
dx$Pct <- dx$Freq / sum(dx$Freq) * 100

d0 <- sum(d1$V2 == 0)
d0pct <- d0 / length(d1$V2) * 100
ani1 <- mean(dd$V2)

leg <- paste0(
  strain_test, ": ANI = ", round(ani1, 3),
  " (", round(nrow(d1) * 1000 / 2 / 1000 / 1000, 1),
  "Mb)", "\n",
  "Zero hits: ", round(d0pct, 2), "% or ",
  round(d0 / 2 / 1000, 2), "Mb", "\n", nsample,
  " Blast samples in total"
)

p1 <- ggplot(dx, aes(x = NI, y = Freq, fill = NI)) +
  geom_col() +
  geom_label(
    colour = "black", aes(label = Freq),
    position = position_stack(vjust = 0.5),
    size = 3, show.legend = FALSE
  ) +
  scale_x_discrete("", breaks = NULL) +
  scale_y_reverse() +
  labs(
    y = "", subtitle = leg,
    caption = paste0(
      "Reference: ", strain_ref,
      "\n", "Test: ", strain_test
    )
  ) +
  guides(fill = guide_legend(
    position = "top",
    title = paste("")
  )) +
  theme_minimal() +
  theme(plot.caption = element_text(hjust = 0)) +
  coord_flip()

d2 <- read.csv(file_rev, header = FALSE)
nsample <- nrow(d2)
dd <- subset(d2, d2$V2 > 0.8)
NI <- cut(dd$V2,
  breaks = c(
    0.8, 0.85, 0.9, 0.95, 0.96, 0.97,
    0.98, 0.99, 0.995, 0.999, 1
  ),
  include.lowest = TRUE, labels = c(
    "0.80 ≤ % ≤ 0.85",
    "0.85 < % ≤ 0.90",
    "0.90 < % ≤ 0.95",
    "0.95 < % ≤ 0.96",
    "0.96 < % ≤ 0.97",
    "0.97 < % ≤ 0.98",
    "0.98 < % ≤ 0.99",
    "0.99 < % ≤ 0.995",
    "0.995 < % ≤ 0.999",
    "0.999 < % ≤ 1"
  )
)
dx <- as.data.frame(table(NI))
dx$Pct <- signif(dx$Freq / sum(dx$Freq) * 100, 3)

d0 <- sum(d2$V2 == 0)
d0pct <- d0 / length(d2$V2) * 100
ani2 <- mean(dd$V2)

leg <- paste0(
  strain_ref, ": ANI = ", round(ani2, 3),
  " (", round(nrow(d2) * 1000 / 2 / 1000 / 1000, 1), "Mb)",
  "\n",
  "Zero hits: ", round(d0pct, 2), "% or ",
  round(d0 / 2 / 1000, 2), "Mb", "\n",
  nsample, " Blast samples in total"
)

p2 <- ggplot(dx, aes(x = NI, y = Freq, fill = NI)) +
  geom_col() +
  geom_label(
    colour = "black", aes(label = Freq),
    position = position_stack(vjust = 0.5),
    size = 3, show.legend = FALSE
  ) +
  scale_x_discrete("", breaks = NULL) +
  labs(
    y = "", subtitle = leg,
    caption = paste0(
      "Reference: ", strain_test,
      "\n", "Test: ", strain_ref
    )
  ) +
  guides(fill = "none") +
  theme_minimal() +
  theme(
    plot.caption = element_text(hjust = 0),
  ) +
  coord_flip()


d1$V3 <- as.numeric(d1$V2 == 0)
n <- 100
nr <- nrow(d1)
d1$groups <- rep(1:ceiling(nr / n), each = n, length.out = nr)
dd0_agg <- aggregate(V3 ~ groups, d1, sum)

p3 <- ggplot(data = dd0_agg, aes(x = groups, y = V3)) +
  geom_line(aes(group = 1)) +
  geom_hline(yintercept = 0, color = gray(.7)) +
  annotate(
    geom = "text", label = strain_test,
    x = -Inf, y = Inf, size = 3, hjust = 0, vjust = 2
  ) +
  labs(x = "", y = "") +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.title.y = element_text(angle = 0),
    panel.grid = element_blank(),
    plot.margin = unit(c(0, 0, 0, 0), "mm")
  )

d2$V3 <- as.numeric(d2$V2 == 0)
n <- 100
nr <- nrow(d2)
d2$groups <- rep(1:ceiling(nr / n), each = n, length.out = nr)
dd0_agg <- aggregate(V3 ~ groups, d2, sum)

p4 <- ggplot(data = dd0_agg, aes(x = groups, y = V3)) +
  geom_line(aes(group = 1)) +
  geom_hline(yintercept = 0, color = gray(.7)) +
  annotate(
    geom = "text", label = strain_ref,
    x = -Inf, y = Inf, size = 3, hjust = 0, vjust = 2
  ) +
  scale_x_discrete(
    breaks = seq(0, nrow(dd0_agg), by = 100),
    labels = paste(seq(0, nrow(dd0_agg), by = 100) / 20, "Mb")
  ) +
  labs(x = "", y = "") +
  theme_minimal() +
  theme(
    axis.text.y = element_blank(),
    axis.title.y = element_text(angle = 0),
    panel.grid = element_blank(),
    plot.margin = unit(c(0, 0, 0, 0), "mm")
  )

pp <- (p1 + p2) / plot_spacer() / p3 / p4 +
  plot_layout(heights = c(8, 0.5, 0.8, 0.8)) +
  plot_annotation(
    paste0(
      "AVERAGE NUCLEOTIDE IDENTITY = ",
      format(round(mean(ani1, ani2) * 100, 2), nsmall = 2), "%"
    ),
    theme = theme(plot.title = element_text(
      size = 16,
      face = "bold", hjust = 0.5
    ))
  )

ggsave(outfile, pp, width = 12, height = 14)
