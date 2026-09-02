# =============================================================================
# theme_paper.R
# Unified greyscale ggplot2 style for all Section 5 figures, plus shared
# aesthetic scales so evaluation mode / metric / aggregation are encoded the
# same way in every figure.
# =============================================================================

library(ggplot2)

theme_paper <- function(base_size = 12) {
  theme_minimal(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      legend.position  = "top",
      legend.box       = "horizontal",
      strip.text       = element_text(face = "bold"),
      plot.title       = element_text(face = "bold", size = base_size + 1),
      plot.subtitle    = element_text(colour = "grey30")
    )
}

# Fixed encodings used across figures ------------------------------------------

EVAL_LEVELS  <- c("In-sample", "Held-out (reconstruction)", "Held-out (completion)")
EVAL_LINES   <- c("In-sample" = "solid",
                  "Held-out (reconstruction)" = "longdash",
                  "Held-out (completion)" = "dotted")
METRIC_LEVELS <- c("Deviance", "Pearson", "Squared-Error")
METRIC_SHAPES <- c("Deviance" = 16, "Pearson" = 17, "Squared-Error" = 15)
AGG_ALPHA     <- c("Micro" = 1.0, "Macro" = 0.45)

EVAL_SHAPES  <- c("In-sample" = 16,
                  "Held-out (reconstruction)" = 17,
                  "Held-out (completion)" = 15)

scale_eval_linetype <- function() {
  scale_linetype_manual(values = EVAL_LINES, name = "Evaluation", drop = TRUE)
}
#' Point-shape scale keyed to the same evaluation levels as the linetype scale;
#' same `name` so ggplot merges the two into a single legend.
scale_eval_shape <- function() {
  scale_shape_manual(values = EVAL_SHAPES, name = "Evaluation", drop = TRUE)
}
scale_metric_shape <- function() {
  scale_shape_manual(values = METRIC_SHAPES, name = "Metric", drop = TRUE)
}
scale_agg_alpha <- function() {
  scale_alpha_manual(values = AGG_ALPHA, name = "Aggregation", drop = TRUE)
}

metric_label <- function(x) {
  c(dev = "Deviance", chisq = "Pearson", se = "Squared-Error")[x]
}

eval_label <- function(x) {
  c(insample          = "In-sample",
    ho_reconstruction = "Held-out (reconstruction)",
    ho_completion     = "Held-out (completion)")[x]
}

#' Save a figure under Results/Figures/<label>/<exp>/ as PDF + PNG.
#' `label` defaults to the RUN_LABEL global set by make_figures.R.
save_fig <- function(p, name, exp, width = 8.2, height = 4.8,
                     label = get0("RUN_LABEL", ifnotfound = "adhoc")) {
  ggsave(p_figures(label, exp, paste0(name, ".pdf")), p,
         width = width, height = height)
  ggsave(p_figures(label, exp, paste0(name, ".png")), p,
         width = width, height = height, dpi = 200, bg = "white")
  invisible(name)
}
