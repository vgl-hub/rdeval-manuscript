library(tidyverse)
library(ggExtra)
library(bit64)
source('rdeval_interface.R')

# Set color palette (jco from ggsci package)
palette <- c("#0073C2", "#EFC000","#868686", "#CD534C", "#7AA6DC", 
             "#003C67", "#8F7700", "#3B3B3B", "#A73030", "#4A6990")

# Inverse Cumulative Distribution Function
icdf <- function(read_lengths) {
  sum <- sum(read_lengths)
  read_lengths <- sort(read_lengths)
  n <- length(read_lengths)
  vals <- unique(read_lengths)
  inv_cs <- sum - cumsum(as.numeric(tabulate(match(read_lengths, vals)) * vals))
  inv_cs_gb <- inv_cs/10^9
  return(tibble(read_length = vals, inv_cs = inv_cs, inv_cs_gb = inv_cs_gb))
}

reads <- tibble(read_length = numeric(), read_quality = numeric(), run = character())
inv_cs <- tibble(read_length = numeric(), inv_cs = numeric(), run = character())

files <- list.files(pattern = 'SRR[[:digit:]]*.rd')
n_runs <- length(files)
index <- 1
run_names <- c('D1 (Vipera latastei)', 'D2 (Amazona ochrocephala)', 'D3 (Ascaphus truei)', 'D4 (Gallus gallus)')
for (file in files) {
  rdFile <- generateRdFile(file)
  df <- tibble(read_length = rdFile$lengths, 
                   read_quality = rdFile$qualities)
  run <- strsplit(basename(file), '[.]')[[1]][1]
  df$run <- run_names[index]
  df_icdf <- icdf(df$read_length)
  df_icdf$run <- run_names[index]
  reads <- rbind(reads, df)
  inv_cs <- rbind(inv_cs, df_icdf)
  cat(printRdSummary(rdFile))
  index = index + 1
}


### Violin Plot
violin <- ggplot(reads, aes(x = run, y = read_length/1000, fill = run)) +
  geom_violin() +
  scale_x_discrete(name = ("Run"),
                   labels = c('D1', 'D2', 'D3', 'D4')) +
  scale_y_continuous(name = ("Read Length (kbp)"),
                     labels = scales::label_comma()) +
  scale_fill_manual(guide = "none", values = palette) + 
  theme(
    axis.title = element_text(family = "Arial", size = 10),
    axis.text = element_text(family = "Arial", size = 9)
  ) +
  theme_bw()

### Read lengths
read_lengths <- ggplot(reads, aes(x = read_length/1000, fill = run)) +
  geom_density(aes(y = after_stat(count)), alpha = 0.5) +
  scale_x_continuous(name = ("Read length (kbp)"),
                     breaks = scales::breaks_width(10),
                     labels = scales::label_comma(),
                     expand = expansion(mult = c(0, 0))) +
  scale_y_continuous(name = "Count", labels = scales::label_comma(),
                     expand = expansion(mult = c(0, 0.05))) +
  scale_fill_manual(guide = "none", values = palette) +
  theme(
    axis.title = element_text(family = "Arial", size = 10),
    axis.text = element_text(family = "Arial", size = 9)
  ) +
  theme_bw() 

### Inverse cumulative distribution plot
icd <- ggplot(inv_cs, aes(read_length/1000, inv_cs_gb, color = run)) +
  geom_point() +
  scale_x_continuous(name = "Read length (kbp)",
                     breaks = scales::breaks_width(10),
                     labels = scales::label_comma(),
                     expand = expansion(mult = c(0.05, 0.05))) +
  scale_y_continuous(name = "Cumulative yield (Gbp)",
                     expand = expansion(mult = c(0.05, 0.05))) +
  scale_color_manual(name = "Run", values = palette) +
  theme(
    axis.title = element_text(family = "Arial", size = 10),
    axis.text = element_text(family = "Arial", size = 9)
  ) +
  theme_bw() 

### Read lengths vs average quality
reads_by_avg_quality <- reads %>%
  group_by(run, read_length) %>%
  summarize_at("read_quality", mean)

g <- reads_by_avg_quality %>%
  ggplot(aes(x = read_length/1000, 
             y = read_quality, 
             )) +
  geom_density_2d(aes(color = ..level..)) +
  scale_color_viridis_c(guide = "none") +
  scale_x_continuous(name = ("Read length (kbp)"),
                     breaks = scales::breaks_width(10),
                     labels = scales::label_comma()) +
  scale_y_continuous(name = ("Average Read Quality")) +
  theme(
    axis.title = element_text(family = "Arial", size = 10),
    axis.text = element_text(family = "Arial", size = 9)
  ) +
  theme_bw() 

# Stitch 4 plots together
library(patchwork)
violin + read_lengths + icd + g +
  plot_layout(guides = 'collect') &
  theme(legend.position = 'top')
ggsave('Figure_1.png', width = 10, height = 6.5)

