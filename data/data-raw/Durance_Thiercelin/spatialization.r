rm(list = ls())
setwd("/home/famendezrios/Documents/These/VSCODE-R/HydroBayes/HydroBayes_git/data/processed_data/Durance_Thiercelin")
library(ggplot2)
library(dplyr)

# ── 1. Read the file ──────────────────────────────────────────────────────────
# Skip the comment line (starts with "*"), no header row
df <- read.table(
    "Thiercelin_K_distribution.RUG",
    comment.char = "*",
    header = FALSE,
    col.names = c("type", "reach", "x_start", "x_end", "K_start", "K_end")
)

# ── 2. Sort by x_start so segments appear in river order ─────────────────────
df <- df[order(df$x_start), ]

# ── 3. Plot ───────────────────────────────────────────────────────────────────
ggplot(df) +
    geom_segment(
        aes(
            x = x_start, xend = x_end,
            y = K_start, yend = K_start
        ), # horizontal segment at K value
        linewidth = 0.7,
        colour = "#2166ac"
    ) +
    # optional: vertical connectors between consecutive segments
    geom_segment(
        aes(
            x = x_end, xend = x_end,
            y = K_start, yend = lead(K_start)
        ),
        linewidth = 0.3,
        colour = "#2166ac",
        linetype = "dashed",
        na.rm = TRUE
    ) +
    labs(
        title    = "Spatially distributed friction",
        subtitle = "Piecewise-constant function along the channel",
        x        = "Streamwise position (m)",
        y        = "Strickler coefficient (m¹/³/s)"
    ) +
    theme_bw(base_size = 13) +
    theme(
        plot.title = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        panel.grid.minor = element_blank()
    )

# ── 4. Save ───────────────────────────────────────────────────────────────────
ggsave("K_thiercelin.png", width = 9, height = 5, dpi = 300)
