## ==============================Intro=======================================
## Finalized analysis for the N-gradient manuscript.
## Experimental design:
##     * 6 N-treatments (0, 28, 56, 84, 140, 196) accross
##     * 4 randomized blocks (treatment block) and 
##     * 2 ingrowth core types (hyphae Enriched; hyphae suppressed)
##     * * 3 analytical reps of each ingrowth core per plot.
## 
## Response variables include:
##    * Diazotroph hill numbers (q=1) 
##    * AMF hill numbers (q=1)
##    * Diazotroph community composition (Bray-Curtis dissimilarity)
##    * AMF community composition (Brayh-Curtis dissimilarity)
##    * N-fixation potential (15-N incubations)
##    * Extra-radical hyphal (ERH) length 
##    * Soil chemistry 
##    * * Soil moisture
##    * * Total nitrogen
##    * * Soil pH
##    * * Soil Organic Matter (SOM)
##    * * Macro-nutrients (P, K)
##    * * Micro-nutrients (Mg, Ca, Mn, Zn, S, B, Cu, Fe)
##  
##  
##  Amplicon sequencing analyses include 
##    * Hill numbers (Alpha Diversity)
##    * Bray-Curtis dissimilarity (Beta Diversity)
##    * Mantel Correlations
##    * Differential Abundance Analysis (DESEQ2)
##    * Co-occurrence Network Analysis (SpiecEasi)


## For questions, please contact: kristybr@msu.edu 


## ==============================Data Objects=======================================
## All response variables can be found and extracted within the phyloseq
## object's metadata. Two pre-processed amplicon datasets are available in this GitHub.
## The raw sequencing data can be found on the NCBI Sequencing Archive: 
## Note: The sequencing data also includes roots and bulk soil samples, but these samples are not included in this
## manuscript analysis.

# Load phyloseq objects /w relevant metadata (Change to the relevant directory)
# AMF dataset (unrarefied)
lsu <- readRDS("amf_phylo.rds")

# NifH dataset (unrarefied)
nif <- readRDS("nifH_phylo.rds")

# Output folders: final figures/tables, and working figures not in the manuscript
for (d in c("Figures", "Figures/other", "Tables")) dir.create(d, showWarnings = FALSE)

# Create N Treatment labels 
N_treatment_labels <- c(expression("SF1 (0 kg N "* ha^{-1}*")"),
                        expression("SF2 (28 kg N "* ha^{-1}*")"),
                        expression("SF3 (56 kg N "* ha^{-1}*")"),
                        expression("SF4 (84 kg N "* ha^{-1}*")"),
                        expression("SF6 (140 kg N "* ha^{-1}*")"),
                        expression("SF8 (196 kg N "* ha^{-1}*")"))

## ==============================Rarefaction=======================================
## Both nifH and LSU datasets will be rarefied based on a sequencing cut-off that
## maintains sufficient depth while also maximizing the number of samples retained.
## 
## Multi-rarefaction will randomly subsample from the OTU over ~100 iterations 

## Necessary packages
library(phyloseq)
library(vegan)
library(tidyverse)
library(cowplot)
library(BRCore)


# Tidy rarefaction curves for a phyloseq object, joined to sample metadata.
# Samples are filtered *before* curves are computed (faster than filtering after).
rarecurve_data <- function(physeq,
                           step            = 100,
                           drop_habitats   = c("azotobacter", "root"),
                           drop_core_types = "Reference") {
  meta <- as(sample_data(physeq), "data.frame") %>%
    rownames_to_column("Site")
  
  keep <- meta %>%
    filter(!habitat %in% drop_habitats,
           !core_type %in% drop_core_types) %>%
    pull(Site)
  physeq <- prune_samples(keep, physeq)
  
  # vegan wants samples as rows
  otu <- as(otu_table(physeq), "matrix")
  if (taxa_are_rows(physeq)) otu <- t(otu)
  
  rarecurve(otu, step = step, tidy = TRUE) %>%
    left_join(meta, by = "Site") %>%
    mutate(core_type = factor(core_type, levels = c("Bulk", "Reference", "Fungi")),
           core_type = fct_recode(core_type,
                                  "Hyphae Suppressed" = "Bulk",
                                  "Hyphae Enriched"   = "Fungi"),
           core_type = fct_drop(core_type))
}

# Rarefaction curve plot, faceted by core type, with a dashed line at the
# chosen rarefaction depth. `xmax` zooms the x-axis without dropping data.
plot_rarecurve <- function(df,
                           depth,
                           xmax         = NULL,
                           legend_title = "N Fertilization",
                           labels       = N_treatment_labels) {
  p <- ggplot(df, aes(x = Sample, y = Species, colour = N_treatment)) +
    geom_point(size = 1) +
    geom_vline(xintercept = depth, linetype = "dashed") +
    facet_grid(~ core_type) +
    scale_color_viridis_d(name = legend_title, labels = labels) +
    labs(x = "Sequencing Depth", y = "Species Richness") +
    guides(colour = guide_legend(override.aes = list(size = 5))) +
    theme_bw() +
    theme(
      axis.title.y = element_text(size = 14),
      strip.text.x = element_text(size = 14),
      strip.text.y = element_text(size = 14),
      axis.text.y  = element_text(size = 12, color = "black"),
      legend.text  = element_text(size = 12),
      legend.title = element_text(size = 14),
      plot.margin  = unit(c(.5, 1, .5, .5), "cm")
    )
  
  if (!is.null(xmax)) p <- p + coord_cartesian(xlim = c(0, xmax))
  p
}

# Repeated rarefaction with BRCore; returns the updated phyloseq object.
rarefy_physeq <- function(physeq,
                          depth,
                          num_iter  = 100,
                          seed      = 37920,
                          iteration = num_iter) {
  rarefied <- multi_rarefy(
    physeq_obj  = physeq,
    depth_level = depth,
    num_iter    = num_iter,
    set_seed    = seed
  )
  update_otu_table(
    physeq_obj    = physeq,
    rarefied_otus = rarefied,
    iteration     = iteration
  )
}

# Specified depth(s)
nif_depth <- 3000
lsu_depth <- 500


## Rarefaction curves 
nif_rarecurve <- rarecurve_data(nif) %>%
  plot_rarecurve(depth = nif_depth)

lsu_rarecurve <- rarecurve_data(lsu) %>%
  plot_rarecurve(depth = lsu_depth, xmax = 3000)

# Optional BRCore diagnostic (works for either object)
nif_rarefaction_metrics <- add_rarefaction_metrics(data = nif)
plot_rarefaction_metrics(nif_rarefaction_metrics)


## Rarefy 
nif_rfy <- rarefy_physeq(nif, depth = nif_depth)
lsu_rfy <- rarefy_physeq(lsu, depth = lsu_depth)


## Combine rarefaction plots 
leg <- get_legend(nif_rarecurve)

rarecurves <- plot_grid(
  nif_rarecurve + theme(legend.position = "none"),
  lsu_rarecurve + theme(legend.position = "none"),
  ncol   = 1,
  labels = c("A", "B")
)

rarecurves_wleg <- plot_grid(
  rarecurves, leg,
  ncol       = 2,
  rel_widths = c(1, 0.25)
)
rarecurves_wleg

ggsave("Figures/Figure_S1.png",
       plot = rarecurves_wleg, dpi = 600, width = 10, height = 5)

## ==============================Core x N  split-plot models=======================================
## Alpha diversity, FLNF, ERH and soil chemistry.
##
## Design: block > plot (N treatment) > core pair (replicate_num) > core (core type)
## Inference: lmer + Kenward-Roger F tests (optionally on log(y));
##           glmmTMB fit kept only for DHARMa residual checks.

library(phyloseq)
library(tidyverse)
library(hillR)
library(glmmTMB)
library(lmerTest)   # lmer() with Kenward-Roger anova
library(DHARMa)
library(emmeans)
library(multcomp)   # cld()
library(ggsci)      # scale_*_npg()
library(cowplot)
library(ggrepel)


## Shared labels
# new label = old label
core_labels <- c("Hyphae Suppressed" = "Bulk",
                 "Reference"         = "Reference",
                 "Hyphae Enriched"   = "Fungi")
N_labels    <- c("0" = "SF1", "28" = "SF2", "56" = "SF3",
                 "84" = "SF4", "140" = "SF6", "196" = "SF8")

relabel_treatments <- function(df) {
  df %>%
    mutate(
      across(any_of("core_type"),
             ~ fct_drop(fct_recode(factor(.x, levels = core_labels), !!!core_labels))),
      across(any_of("N_treatment"),
             ~ fct_recode(factor(.x, levels = N_labels), !!!N_labels))
    )
}

## Data preparation 

# Sample metadata from a phyloseq object, with sample IDs in `Site`.
sample_df <- function(physeq) {
  as(sample_data(physeq), "data.frame") %>%
    rownames_to_column("Site")
}

# Filter samples, convert response/covariate columns to numeric, and build
# the design factors. Pass drop_core_types = NULL to keep the Reference cores.
prep_design <- function(df,
                        numeric_vars    = character(0),
                        drop_habitats   = c("azotobacter", "root"),
                        drop_core_types = "Reference") {
  df %>%
    filter(!habitat %in% drop_habitats,
           !core_type %in% drop_core_types) %>%
    mutate(
      across(any_of(numeric_vars), ~ as.numeric(as.character(.x))),
      plot            = factor(paste(N_treatment, treatment_block)),  # whole plot
      pair            = factor(paste(plot, replicate_num)),           # bulk/fungal core pair
      treatment_block = factor(treatment_block),
      across(any_of("data_harvested"), as.factor)
    ) %>%
    relabel_treatments()
}

# Hill numbers per sample, joined to metadata and run through prep_design().
hill_data <- function(physeq, q = 1, ...) {
  otu <- as(otu_table(physeq), "matrix")
  if (taxa_are_rows(physeq)) otu <- t(otu)          # hillR wants samples as rows
  hill <- hill_taxa(as.data.frame(otu), q = q)
  
  sample_df(physeq) %>%
    left_join(tibble(Site = names(hill), hill = unname(hill)), by = "Site") %>%
    prep_design(...)
}

# Expect n = 2 for every pair (one bulk + one fungal core).
check_pairs <- function(df) {
  dplyr::count(df, pair) %>% dplyr::count(n, name = "n_pairs")
}


## Models 

## Split-plot mixed models for one response.
##  response     column name of the response
##  covariates   extra fixed-effect columns (character vector)
##  log_response TRUE fits lmer on log(y) (emmeans back-transforms) and a
##               Gamma(log) glmmTMB for diagnostics
## Rows with NA in the response or covariates are dropped, and that same
## data is used for the post-hoc tests and the raw points in the figure.
fit_models <- function(df,
                       response,
                       covariates   = character(0),
                       log_response = FALSE,
                       family       = if (log_response) Gamma(link = "log") else gaussian(),
                       fit_glmm     = TRUE) {
  dat <- df %>%
    filter(if_all(all_of(c(response, covariates)), ~ !is.na(.x))) %>%
    droplevels()
  
  if (log_response && any(dat[[response]] <= 0)) {
    stop(response, " has values <= 0; set log_response = FALSE or handle zeros first.")
  }
  
  fixed <- paste(c("core_type * N_treatment", covariates, "treatment_block"),
                 collapse = " + ")
  rhs   <- paste(fixed, "+ (1 | plot) + (1 | pair)")
  lhs   <- if (log_response) sprintf("log(%s)", response) else response
  
  lmm  <- lmerTest::lmer(as.formula(paste(lhs, "~", rhs)), data = dat)
  glmm <- if (fit_glmm) {
    glmmTMB(as.formula(paste(response, "~", rhs)), data = dat, family = family)
  }
  
  list(
    response     = response,
    log_response = log_response,
    data         = dat,
    lmer         = lmm,
    glmm         = glmm,
    anova_t3     = anova(lmm, ddf = "Kenward-Roger"),            # lmerTest default
    anova_t2     = anova(lmm, type = 2, ddf = "Kenward-Roger")   # reported / label
  )
}

# Variance components and DHARMa residual plots for one fit.
check_fit <- function(fit) {
  print(VarCorr(fit$lmer))
  plot(simulateResiduals(fit$lmer))
  if (!is.null(fit$glmm)) plot(simulateResiduals(fit$glmm))
  invisible(fit)
}


## Post-hoc testing 

## Model estimates (core x N cells and core-type marginal means), N-treatment
## letters, the raw data for plotting, and a KR ANOVA text label.
posthoc <- function(fit) {
  dat       <- fit$data
  y_letters <- max(dat[[fit$response]], na.rm = TRUE) * 1.08  # outside mutate(): emmeans has a `df` column
  
  tidy_emm <- function(emm, ...) {
    emm %>%
      cld(Letters = letters, ...) %>%
      as.data.frame() %>%
      dplyr::rename(estimate = any_of(c("emmean", "response")))       # name differs if logged
  }
  
  cell <- emmeans(fit$lmer, ~ core_type * N_treatment, type = "response") %>%
    tidy_emm(alpha = 0.05)
  core <- emmeans(fit$lmer, ~ core_type, type = "response") %>%
    tidy_emm(alpha = 0.05)
  
  letters_N <- emmeans(fit$lmer, ~ N_treatment, lmer.df = "kenward-roger") %>%
    cld(Letters = letters, adjust = "tukey", reversed = TRUE) %>%
    as.data.frame() %>%
    mutate(.group = trimws(.group),
           y      = y_letters)
  
  at <- as.data.frame(fit$anova_t2)[c("N_treatment", "core_type",
                                      "core_type:N_treatment"), ]
  fmt_p <- function(p) ifelse(p < 0.001, "p < 0.001", sprintf("p = %.3f", p))
  stat_lab <- paste(sprintf("%s: F(%d, %.1f) = %.2f, %s",
                            c("N", "Core type", "N × Core type"),
                            as.integer(at$NumDF), at$DenDF, at$`F value`,
                            fmt_p(at$`Pr(>F)`)),
                    collapse = "\n")
  
  list(cell      = cell,
       core      = core,
       letters_N = letters_N,
       raw       = mutate(dat, value = .data[[fit$response]]),
       stat_lab  = stat_lab)
}


## Figure

# Raw points + model estimates with 95% CIs + N-treatment letters.
# facet = TRUE expects a `variable` column (see combine_posthoc()).
plot_response <- function(post, ylab, dodge = 0.35, facet = FALSE, ncol = 4) {
  p <- ggplot(post$cell, aes(x = N_treatment, y = estimate, group = core_type)) +
    geom_line(aes(col = core_type), linewidth = 1.15) +
    geom_point(data = post$raw, aes(y = value, col = core_type),
               alpha = 0.4, size = 1.8,
               position = position_jitterdodge(jitter.width = 0.1,
                                               dodge.width  = dodge,
                                               seed         = 1)) +
    geom_point(aes(fill = core_type), shape = 21, col = "black", size = 5.5,
               position = position_dodge(width = dodge)) +
    geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL, col = core_type),
                  position = position_dodge(width = dodge),
                  width = 0.2, linewidth = 0.75) +
    geom_text(data = post$letters_N, aes(N_treatment, y = y, label = .group),
              inherit.aes = FALSE, size = 5) +
    labs(x = "N Treatment", y = ylab) +
    scale_fill_npg(name = "Core Type") +
    scale_color_npg(name = "Core Type") +
    theme_bw() +
    theme(
      axis.text.x  = element_text(size = 16, color = "black"),
      axis.title.x = element_text(size = 16),
      axis.title.y = element_text(size = 16),
      strip.text.x = element_text(size = 16, color = "black"),
      strip.text.y = element_text(size = 16, color = "black"),
      axis.text.y  = element_text(size = 14, color = "black"),
      legend.text  = element_text(size = 14),
      legend.title = element_text(size = 16)
    )
  
  if (facet) p <- p + facet_wrap(~ variable, scales = "free_y", ncol = ncol)
  p
}

# Panels with shared legend (legend taken from the first plot).
panel_figure <- function(plots, ncol = 2, legend_width = 0.35,
                         labels = LETTERS[seq_along(plots)]) {
  leg    <- get_legend(plots[[1]])
  panels <- plot_grid(plotlist = purrr::map(plots, ~ .x + theme(legend.position = "none")),
                      ncol = ncol, labels = labels)
  plot_grid(panels, leg, ncol = 2, rel_widths = c(1, legend_width))
}


## One-call wrapper and multi-variable summaries 

# Fit, post-hoc and plot one response.
analyze_response <- function(df, response, ylab = response, ...) {
  fit  <- fit_models(df, response, ...)
  post <- posthoc(fit)
  list(fit = fit, post = post, plot = plot_response(post, ylab))
}

# Stack post-hoc output from several analyze_response() results for a
# faceted figure. `labels` is a named vector: names = response columns.
combine_posthoc <- function(results, labels = set_names(names(results))) {
  pick <- function(el) {
    purrr::map(results, ~ .x$post[[el]]) %>%
      bind_rows(.id = "variable") %>%
      mutate(variable = factor(variable, levels = names(labels), labels = labels))
  }
  list(cell = pick("cell"), raw = pick("raw"), letters_N = pick("letters_N"))
}

# KR (Type II) tests for the design terms across several responses, with
# Benjamini-Hochberg adjustment within each term across responses.
anova_table <- function(results,
                        terms = c("N_treatment", "core_type", "core_type:N_treatment")) {
  purrr::imap_dfr(results, ~ as.data.frame(.x$fit$anova_t2) %>%
                    rownames_to_column("term") %>%
                    filter(term %in% terms) %>%
                    mutate(variable = .y, .before = 1)) %>%
    group_by(term) %>%
    mutate(p_BH = p.adjust(`Pr(>F)`, method = "BH")) %>%
    ungroup()
}

# DHARMa uniformity / dispersion tests and singular-fit flags (lmer fits).
residual_checks <- function(results) {
  purrr::imap_dfr(results, function(res, nm) {
    sim <- simulateResiduals(res$fit$lmer)
    tibble(variable     = nm,
           uniformity_p = testUniformity(sim, plot = FALSE)$p.value,
           dispersion_p = testDispersion(sim, plot = FALSE)$p.value,
           singular     = lme4::isSingular(res$fit$lmer))
  })
}


## Alpha diversity analysis 


## nifH 
nif_df <- hill_data(nif_rfy)
check_pairs(nif_df)

nif_div_res <- analyze_response(nif_df, "hill",
                                ylab = "Diazotroph Hill Num. (q=1)")
check_fit(nif_div_res$fit)
nif_div_res$fit$anova_t3
nif_div_res$post$core          # core-type marginal means
nif_div_res$post$cell
nif_div_res$post$letters_N
nif_div <- nif_div_res$plot
nif_div


## LSU 
lsu_df <- hill_data(lsu_rfy)
check_pairs(lsu_df)

lsu_div_res <- analyze_response(lsu_df, "hill",
                                ylab         = "AMF Hill Num. (q=1)",
                                log_response = TRUE)
check_fit(lsu_div_res$fit)
lsu_div_res$fit$anova_t3
lsu_div_res$post$core
lsu_div_res$post$letters_N
lsu_div <- lsu_div_res$plot
lsu_div


## Combine nifH-LSU diversity plots 
div_plots_wleg <- panel_figure(list(nif_div, lsu_div))
div_plots_wleg

ggsave("Figures/other/hill_num_plots.png",
       plot = div_plots_wleg, dpi = 600, width = 10, height = 5)


# ERH, FLNF, and soil chemistry 


## FLNF and ERH 
# Full, unrarefied metadata (no diversity covariate needed)
fun_df <- sample_df(nif) %>%
  prep_design(numeric_vars = c("n_fixation", "erh", "soil_moisture"))
check_pairs(fun_df)

# Planned N contrasts (N levels are ordered 0, 28, 56, 84, 140, 196).
# Liming began at 84 kg N, so 28-56 = unlimed and 84-196 = limed fertilized plots.
fert_contrast <- list("Unfertilized / Fertilized" = c(1, rep(-1/5, 5)))
N_contrasts   <- list(
  "Unfertilized / Unlimed fertilized (28-56)" = c(1, -1/2, -1/2, 0, 0, 0),
  "Unfertilized / Limed fertilized (84-196)"  = c(1, 0, 0, -1/3, -1/3, -1/3)
)

# FLNF
flnf_res <- analyze_response(
  fun_df, "n_fixation",
  ylab         = expression(FLNF ~ ("µg" ~ N ~ g^-1 ~ dry ~ soil ~ day^-1)),
  log_response = TRUE
)
check_fit(flnf_res$fit)
flnf_res$fit$anova_t3
flnf_res$post$core

flnf <- flnf_res$plot
flnf


# Averaged over core types
emmeans(flnf_res$fit$lmer, ~ N_treatment, lmer.df = "kenward-roger") %>%
  contrast(fert_contrast) %>%
  summary(infer = TRUE, type = "response")

# Within each core type
flnf_fert_core <- emmeans(flnf_res$fit$lmer, ~ N_treatment | core_type,
                          lmer.df = "kenward-roger") %>%
  contrast(fert_contrast)
summary(flnf_fert_core, infer = TRUE, type = "response")


emmeans(flnf_res$fit$lmer, ~ N_treatment, lmer.df = "kenward-roger") %>%
  contrast(N_contrasts) %>%
  summary(infer = TRUE, type = "response")
# Does the fertilization effect differ between core types?
pairs(flnf_fert_core, by = NULL) %>%
  summary(infer = TRUE, type = "response")


# ERH
erh_res <- analyze_response(
  fun_df, "erh",
  ylab         = expression("ERH" * " (mm " * g^{-1} * " soil)"),
  covariates   = c("soil_moisture"),
  log_response = TRUE
)
check_fit(erh_res$fit)
erh_res$fit$anova_t3
erh_res$post$core

erh <- erh_res$plot
erh

# Averaged over core types
emmeans(erh_res$fit$lmer, ~ N_treatment, lmer.df = "kenward-roger") %>%
  contrast(fert_contrast) %>%
  summary(infer = TRUE, type = "response")

# Within each core type
erh_fert_core <- emmeans(erh_res$fit$lmer, ~ N_treatment | core_type,
                         lmer.df = "kenward-roger") %>%
  contrast(fert_contrast)
summary(erh_fert_core, infer = TRUE, type = "response")


emmeans(erh_res$fit$lmer, ~ N_treatment, lmer.df = "kenward-roger") %>%
  contrast(N_contrasts) %>%
  summary(infer = TRUE, type = "response")
# Does the fertilization effect differ between core types?
pairs(erh_fert_core, by = NULL) %>%
  summary(infer = TRUE, type = "response")


## Combine FLNF and ERH plots 
fun_plots_wleg <- panel_figure(list(flnf, erh))
fun_plots_wleg

ggsave("Figures/other/flnf_erh_plots.png",
       plot = fun_plots_wleg, dpi = 600, width = 10, height = 5)

# Diversity + function (Figure 1)
div_fun_wleg <- panel_figure(list(nif_div, lsu_div, flnf, erh))
div_fun_wleg

ggsave("Figures/Figure_1.png",
       plot = div_fun_wleg, dpi = 600, width = 10, height = 8)


## Soil chemistry: per-variable models
# log = TRUE: model log(y). Check residual_checks() / DHARMa and adjust.
chem_spec <- tribble(
  ~var,            ~label,          ~log_y,
  "soil_moisture", "Soil moisture", FALSE,
  "pH",            "pH",            FALSE,
  "BpH",           "Buffer pH",     FALSE,
  "OM",            "OM",            TRUE,
  "P",             "P",             TRUE,
  "K",             "K",             TRUE,
  "Mg",            "Mg",            TRUE,
  "Ca",            "Ca",            TRUE,
  "CEC",           "CEC",           TRUE,
  "K_Saturation",  "K saturation",  FALSE,
  "Mg_Saturation", "Mg saturation", FALSE,
  "Ca_Saturation", "Ca saturation", FALSE,
  "S",             "S",             TRUE,
  "B",             "B",             TRUE,
  "Cu",            "Cu",            TRUE,
  "Fe",            "Fe",            TRUE,
  "Mn",            "Mn",            TRUE,
  "Zn",            "Zn",            TRUE,
  "total_n",       "Total N",       TRUE
)

chem_df <- sample_df(nif) %>%
  prep_design(numeric_vars = chem_spec$var)

# One variable failing (e.g. zeros with log_y = TRUE) gives a warning, not a stop
chem_res <- purrr::pmap(chem_spec, function(var, label, log_y) {
  tryCatch(
    analyze_response(chem_df, var, ylab = label,
                     log_response = log_y, fit_glmm = FALSE),
    error = function(e) {
      warning(var, ": ", conditionMessage(e), call. = FALSE)
      NULL
    }
  )
}) %>%
  set_names(chem_spec$var) %>%
  purrr::compact()

# KR tests for all variables, BH-adjusted within each term
chem_stats <- anova_table(chem_res)
chem_stats

# Residual checks across variables (lmer fits)
residual_checks(chem_res)
# For a closer look at one variable: check_fit(chem_res$P$fit)

# Faceted figure of model estimates
chem_post   <- combine_posthoc(chem_res,
                               labels = set_names(chem_spec$label, chem_spec$var)[names(chem_res)])
chem_trends <- plot_response(chem_post, ylab = "Model estimate (95% CI)",
                             facet = TRUE, ncol = 4)
chem_trends

ggsave("Figures/Figure_S2.png",
       plot = chem_trends, dpi = 600, width = 15, height = 14)


## Soil chemistry: PCA
chem_pca_df <- sample_df(nif) %>%
  prep_design(numeric_vars = chem_spec$var) %>%   # Reference cores dropped by default
  dplyr::select(Site, N_treatment, core_type, all_of(chem_spec$var)) %>%
  drop_na()


pca_res <- prcomp(dplyr::select(chem_pca_df, all_of(chem_spec$var)),
                  center = TRUE, scale. = TRUE)
percent_var <- summary(pca_res)$importance[2, 1:2] * 100

# Mean +/- SE of site scores per treatment
chem_ord_sum <- chem_pca_df %>%
  dplyr::select(N_treatment, core_type) %>%
  bind_cols(as.data.frame(pca_res$x[, 1:2])) %>%
  group_by(N_treatment, core_type) %>%
  summarise(N        = n(),
            Mean_PC1 = mean(PC1), SE_PC1 = sd(PC1) / sqrt(N),
            Mean_PC2 = mean(PC2), SE_PC2 = sd(PC2) / sqrt(N),
            .groups  = "drop")

# Loadings (strongest vectors only), scaled for plotting
arrow_mult <- 1.8
loadings_filt <- as.data.frame(pca_res$rotation) %>%
  rownames_to_column("Variable") %>%
  mutate(arrow_length = sqrt(PC1^2 + PC2^2)) %>%
  filter(arrow_length > 0.25) %>%
  mutate(PC1_plot = PC1 * arrow_mult,
         PC2_plot = PC2 * arrow_mult)

chem_pca <- ggplot(chem_ord_sum, aes(x = Mean_PC1, y = Mean_PC2)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey70", linewidth = 0.4) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey70", linewidth = 0.4) +
  geom_point(aes(fill = N_treatment, shape = core_type),
             size = 5.5, stroke = 1, color = "black") +
  geom_segment(data = loadings_filt,
               aes(x = 0, y = 0, xend = PC1_plot, yend = PC2_plot),
               inherit.aes = FALSE,
               arrow = arrow(length = unit(0.16, "cm"), type = "closed"),
               linewidth = 0.6, color = "grey25", alpha = 0.85) +
  geom_text_repel(data = loadings_filt,
                  aes(x = PC1_plot, y = PC2_plot, label = Variable),
                  inherit.aes = FALSE, size = 4, color = "grey15",
                  box.padding = 0.25, point.padding = 0.15,
                  segment.color = "grey50", segment.size = 0.3,
                  max.overlaps = Inf) +
  scale_shape_manual(name = "Core Type", values = c(21, 24)) +
  scale_fill_viridis_d(name = "N Treatment") +
  labs(x = paste0("PC1 (", round(percent_var[1], 1), "%)"),
       y = paste0("PC2 (", round(percent_var[2], 1), "%)")) +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  theme_bw() +
  theme(
    axis.text.x  = element_text(size = 16, color = "black"),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text.y  = element_text(size = 14, color = "black"),
    legend.text  = element_text(size = 14),
    legend.title = element_text(size = 16)
  )
chem_pca

ggsave("Figures/other/soil_chem_PCA.png",
       plot = chem_pca, dpi = 600, width = 8, height = 5)


## 
## ============================== Beta Diversity Analysis =======================================
## Community composition: nifH (diazotrophs) and LSU (AMF)
##   * split-plot PERMANOVA (both datasets), accounting for core pairs
##   * NMDS ordinations + envfit vectors (both datasets)
##
## Sample sets
##   nif_ps          all nifH soil cores            -> nifH PERMANOVA
##   lsu_ps          all LSU soil cores - outliers  -> LSU PERMANOVA, LSU NMDS (panel B)
##   nif_ps_matched  nifH cores also in lsu_ps      -> nifH NMDS (panel A), so the
##                                                     AMF NMDS axes can be fitted
## Design (same hierarchy as the mixed models above)
##   treatment_block = 4 blocks (harvest date confounded with block)
##   plot            = N_treatment x block (whole plot)
##   pair            = bulk/fungal core pair within plot (replicate_num)
##   N_treatment     = whole-plot factor, tested against plot (block:N) by
##                     permuting N labels among plots within blocks
##   core_type       = applied within pairs; core_type and N x core_type are
##                     tested against the within-pair residual by permuting
##                     core labels within pairs

library(vegan)
library(permute)
library(patchwork)

SEED  <- 37920
NPERM <- 9999

## Helper functions

# Keep soil cores only (no roots, Azotobacter control, or Reference soils),
# optionally dropping named samples. Counts are left as they are.
soil_cores <- function(ps, drop = character(0)) {
  sd   <- data.frame(sample_data(ps))
  keep <- !(sd$habitat   %in% c("root", "azotobacter")) &
    !(sd$core_type %in% "Reference") &
    !(rownames(sd) %in% drop)
  prune_samples(rownames(sd)[keep], ps)
}

# Soil cores, empty taxa removed, converted to relative abundance.
prep_rel <- function(ps, drop = character(0)) {
  ps <- soil_cores(ps, drop)
  ps <- prune_taxa(taxa_sums(ps) > 0, ps)
  transform_sample_counts(ps, function(x) x / sum(x))
}

# Metadata with design factors (plot, pair) and the shared treatment labels.
# Row names stay as sample names so rows can be matched to distance labels.
get_meta <- function(ps) {
  sam <- data.frame(sample_data(ps))
  ids <- rownames(sam)
  sam <- sam %>%
    mutate(treatment_block = factor(treatment_block),
           plot            = factor(paste(N_treatment, treatment_block)),
           pair            = factor(paste(plot, replicate_num))) %>%
    relabel_treatments()
  rownames(sam) <- ids
  sam
}

# Cores sharing no taxa with the rest of the dataset (disconnected), or whose
# nearest neighbour is >= threshold Bray-Curtis. Documents the outlier rule.
flag_outliers <- function(dis, threshold = 0.99) {
  D <- as.matrix(dis); diag(D) <- NA
  grp  <- distconnected(dis, toolong = 1, trace = FALSE)
  main <- as.numeric(names(which.max(table(grp))))
  out  <- data.frame(sample    = rownames(D),
                     min_BC    = apply(D, 1, min, na.rm = TRUE),
                     connected = grp == main)
  out[!out$connected | out$min_BC >= threshold, ]
}

# Split-plot PERMANOVA with core pairs.
# Sequential terms: block -> N -> plot (= block:N) -> pair (within plot)
#                   -> core_type -> N:core_type
# `plot` is used instead of block:N so the plot stratum is fitted before
# `pair` (R would otherwise move the interaction after the main effects).
splitplot_permanova <- function(dis, sam, nperm = NPERM, seed = SEED) {
  sam <- droplevels(sam[labels(dis), ])
  stopifnot(identical(rownames(sam), labels(dis)),
            all(rowSums(table(sam$plot, sam$N_treatment)     > 0) == 1),
            all(rowSums(table(sam$plot, sam$treatment_block) > 0) == 1),
            all(rowSums(table(sam$pair, sam$plot)            > 0) == 1),  # pairs nested in plots
            all(table(sam$pair, sam$core_type) <= 1),                     # <= 1 core of each type per pair
            nlevels(sam$plot) == nlevels(sam$treatment_block) * nlevels(sam$N_treatment))
  
  # Within-pair tests: core_type and N x core_type, permuting cores within pairs
  set.seed(seed)
  fit <- adonis2(dis ~ treatment_block + N_treatment + plot + pair +
                   core_type + N_treatment:core_type,
                 data = sam, by = "terms",
                 permutations = how(nperm = nperm, blocks = sam$pair))
  tab <- as.data.frame(fit)
  int_term <- grep("^(core_type:N_treatment|N_treatment:core_type)$",
                   rownames(tab), value = TRUE)
  
  # Whole-plot test: N_treatment against plot (block:N), permuting N among plots
  D2 <- as.matrix(dis)^2
  n  <- nrow(D2)
  J  <- diag(n) - matrix(1 / n, n, n)
  G  <- -0.5 * J %*% D2 %*% J                       # Gower-centred matrix
  basis <- function(f, d) {
    q <- qr(model.matrix(f, data = d))
    qr.Q(q)[, seq_len(q$rank), drop = FALSE]
  }
  trHG <- function(Q) sum(Q * (G %*% Q))            # model sum of squares
  
  Q1 <- basis(~ treatment_block, sam)               # unchanged by N relabelling
  Q3 <- basis(~ plot, sam)                          # spans all plots: also unchanged
  s1 <- trHG(Q1); s3 <- trHG(Q3)
  
  F_N <- function(N_lab) {
    d  <- sam; d$N_treatment <- N_lab
    Q2 <- basis(~ treatment_block + N_treatment, d)
    s2 <- trHG(Q2)
    ((s2 - s1) / (ncol(Q2) - ncol(Q1))) / ((s3 - s2) / (ncol(Q3) - ncol(Q2)))
  }
  F_obs <- F_N(sam$N_treatment)
  stopifnot(isTRUE(all.equal(                       # hand calc matches adonis2
    F_obs,
    (tab["N_treatment", "SumOfSqs"] / tab["N_treatment", "Df"]) /
      (tab["plot",        "SumOfSqs"] / tab["plot",        "Df"]))))
  
  plots <- unique(sam[, c("plot", "treatment_block", "N_treatment")])
  set.seed(seed)
  perms <- shuffleSet(nrow(plots),
                      control = how(nperm = nperm, blocks = plots$treatment_block))
  idx    <- match(sam$plot, plots$plot)
  F_perm <- apply(perms, 1, function(p) F_N(plots$N_treatment[p][idx]))
  p_N    <- (sum(F_perm >= F_obs) + 1) / (length(F_perm) + 1)
  
  SS_tot <- tab["Total", "SumOfSqs"]
  row_of <- function(term, F, P, err) data.frame(
    Term = term, Df = tab[term, "Df"], SumOfSqs = tab[term, "SumOfSqs"],
    R2 = tab[term, "SumOfSqs"] / SS_tot, F = F, P = P, Error_term = err)
  
  list(
    adonis  = fit,   # block, plot and pair R2 from here; their p-values are not valid
    results = rbind(
      row_of("N_treatment", F_obs, p_N, "plot (whole-plot)"),
      row_of("core_type", tab["core_type", "F"], tab["core_type", "Pr(>F)"],
             "Residual (within-pair)"),
      row_of(int_term, tab[int_term, "F"], tab[int_term, "Pr(>F)"],
             "Residual (within-pair)")),
    sam = sam)
}

# Core-type effect within each N level (follow-up to a significant interaction),
# permuting cores within pairs
core_within_N <- function(dis, sam, nperm = NPERM, seed = SEED) {
  sam <- droplevels(sam[labels(dis), ])
  out <- do.call(rbind, lapply(levels(sam$N_treatment), function(lv) {
    keep <- sam$N_treatment == lv
    s    <- droplevels(sam[keep, ])
    d    <- as.dist(as.matrix(dis)[keep, keep])
    set.seed(seed)
    res  <- adonis2(d ~ pair + core_type, data = s, by = "terms",
                    permutations = how(nperm = nperm, blocks = s$pair))
    data.frame(N_treatment = lv, R2 = res["core_type", "R2"],
               F = res["core_type", "F"], P = res["core_type", "Pr(>F)"])
  }))
  out$P_holm <- p.adjust(out$P, method = "holm")
  out
}

# Seeded 2-D NMDS
run_nmds <- function(dis, k = 2, seed = SEED) {
  set.seed(seed)
  metaMDS(dis, k = k, trymax = 100, trace = FALSE)
}

# envfit with metadata rows matched to ordination sites BY NAME
# Plot-level envfit: site scores and env averaged per plot, permuted within blocks
fit_vectors_plot <- function(ord, meta, vars, alpha = 0.05,
                             nperm = NPERM, seed = SEED) {
  sc  <- as.data.frame(scores(ord, display = "sites"))
  dat <- cbind(sc, meta[rownames(sc), c("plot", "treatment_block", vars)])
  dat[vars] <- lapply(dat[vars], function(x) as.numeric(as.character(x)))
  
  agg <- dat %>%
    group_by(plot, treatment_block) %>%
    summarise(across(c(NMDS1, NMDS2, all_of(vars)), ~ mean(.x, na.rm = TRUE)),
              .groups = "drop")
  
  set.seed(seed)
  fit <- envfit(as.matrix(agg[, c("NMDS1", "NMDS2")]),
                as.data.frame(agg[vars]),
                permutations = how(nperm = nperm, blocks = agg$treatment_block),
                na.rm = TRUE)
  vec <- as.data.frame(scores(fit, display = "vectors"))
  vec$r2   <- fit$vectors$r
  vec$p    <- fit$vectors$pvals
  vec$name <- rownames(vec)
  list(fit = fit, table = vec, sig = vec[vec$p < alpha, ], plot_means = agg)
}

# Centroid (+/- SE) ordination plot with significant vectors
plot_ord <- function(ord, meta, vec_sig) {
  sc <- as.data.frame(scores(ord, display = "sites"))
  sc <- cbind(sc, meta[rownames(sc), c("N_treatment", "core_type")])
  summ <- sc %>%
    group_by(N_treatment, core_type) %>%
    summarise(across(c(NMDS1, NMDS2),
                     list(mean = ~ mean(.x), se = ~ sd(.x) / sqrt(length(.x)))),
              n = n(), .groups = "drop")
  
  ggplot(summ, aes(NMDS1_mean, NMDS2_mean)) +
    geom_hline(yintercept = 0, colour = "grey30", linetype = "dashed") +
    geom_vline(xintercept = 0, colour = "grey30", linetype = "dashed") +
    geom_segment(data = vec_sig, inherit.aes = FALSE,
                 aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
                 linewidth = 1.25, alpha = 0.7, colour = "darkgrey") +
    geom_errorbar(aes(ymin = NMDS2_mean - NMDS2_se, ymax = NMDS2_mean + NMDS2_se), width = 0) +
    geom_errorbar(aes(xmin = NMDS1_mean - NMDS1_se, xmax = NMDS1_mean + NMDS1_se), width = 0) +
    geom_point(aes(fill = N_treatment, shape = core_type), size = 5.5, stroke = 1) +
    geom_label_repel(data = vec_sig, inherit.aes = FALSE,
                     aes(x = NMDS1, y = NMDS2, label = name),
                     fontface = "bold", size = 3.5, min.segment.length = Inf) +
    scale_shape_manual(name = "Core Type", values = c(21, 22)) +
    scale_fill_viridis_d(name = "N Treatment", labels = N_treatment_labels) +
    guides(fill = guide_legend(override.aes = list(shape = 21))) +
    labs(x = "NMDS1", y = "NMDS2") +
    theme_bw() +
    theme(axis.text    = element_text(size = 14, colour = "black"),
          axis.title   = element_text(size = 16),
          legend.text  = element_text(size = 14),
          legend.title = element_text(size = 16))
}


## Sample sets
# LSU outliers: share no ASVs with any other core (see flag_outliers below)
lsu_outliers <- c("SF8_B4_Bul_R3", "SF3_B4_Bul_R3")

nif_ps <- prep_rel(nif_rfy)
lsu_ps <- prep_rel(lsu_rfy, drop = lsu_outliers)

common         <- base::intersect(sample_names(nif_ps), sample_names(lsu_ps))
nif_ps_matched <- prune_samples(common, nif_ps)

nif_meta   <- get_meta(nif_ps)
lsu_meta   <- get_meta(lsu_ps)
nif_meta_m <- get_meta(nif_ps_matched)

check_pairs(nif_meta)   # LSU: the two outliers leave their partner cores unpaired
check_pairs(lsu_meta)

nif_bray   <- phyloseq::distance(nif_ps,         method = "bray")
lsu_bray   <- phyloseq::distance(lsu_ps,         method = "bray")
nif_bray_m <- phyloseq::distance(nif_ps_matched, method = "bray")

c(nifH = nsamples(nif_ps), LSU = nsamples(lsu_ps), nifH_matched = nsamples(nif_ps_matched))

# Documents the outlier rule on the full LSU core set
flag_outliers(phyloseq::distance(prep_rel(lsu_rfy), method = "bray"))


## PERMANOVA
nif_sp <- splitplot_permanova(nif_bray, nif_meta)
nif_sp$adonis
nif_sp$results
nif_core <- core_within_N(nif_bray, nif_sp$sam)
nif_core

lsu_sp <- splitplot_permanova(lsu_bray, lsu_meta)
lsu_sp$adonis
lsu_sp$results
# Run if the LSU interaction is significant:
# lsu_core <- core_within_N(lsu_bray, lsu_sp$sam); lsu_core

# Sensitivity: LSU with outliers retained
lsu_ps_all <- prep_rel(lsu_rfy)
lsu_sp_all <- splitplot_permanova(phyloseq::distance(lsu_ps_all, method = "bray"),
                                  get_meta(lsu_ps_all))
lsu_sp_all$results


## Ordinations + vectors
soil_vars <- c("pH", "OM", "P", "K", "S", "B", "Cu", "Fe", "Mn", "Zn",
               "total_n", "n_fixation", "erh", "soil_moisture")

# LSU: one ordination, used for panel B AND for the AMF axes in panel A
lsu_nmds <- run_nmds(lsu_bray)
lsu_nmds$stress # 0.14
lsu_vec <- fit_vectors_plot(lsu_nmds, lsu_meta, soil_vars)
lsu_vec$table

# nifH (matched set), with AMF scores taken from lsu_nmds by sample name
lsu_sc <- scores(lsu_nmds, display = "sites")
nif_meta_m$AMF_NMDS1 <- lsu_sc[rownames(nif_meta_m), "NMDS1"]
nif_meta_m$AMF_NMDS2 <- lsu_sc[rownames(nif_meta_m), "NMDS2"]
stopifnot(!anyNA(nif_meta_m$AMF_NMDS2))

nif_nmds <- run_nmds(nif_bray_m)
nif_nmds$stress # 0.13
nif_vec <- fit_vectors_plot(nif_nmds, nif_meta_m, c(soil_vars, "AMF_NMDS1", "AMF_NMDS2"))
nif_vec$table

vector_labels <- c(soil_moisture = "Soil moisture", total_n = "Total N", erh = "ERH",
                   n_fixation = "FLNF", AMF_NMDS1 = "AMF NMDS1", AMF_NMDS2 = "AMF NMDS2")
relabel_vectors <- function(v) {
  i <- v$name %in% names(vector_labels)
  v$name[i] <- vector_labels[v$name[i]]
  v
}

p_nif <- plot_ord(nif_nmds, nif_meta_m, relabel_vectors(nif_vec$sig))
p_lsu <- plot_ord(lsu_nmds, lsu_meta,   relabel_vectors(lsu_vec$sig))

fig_ord <- (p_nif | p_lsu) + plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A")
fig_ord

ggsave("Figures/Figure_2.png",
       plot = fig_ord, dpi = 600, width = 12, height = 5)



## ============================== Differential Abundance (DESeq2) =======================================
## Raw (unrarefied) counts. DESeq2 has no random effects, and N is a
## whole-plot factor, so replicate cores are pseudobulked (counts summed)
## within each plot x core type: each plot is then one replicate of N.
##
## Core pairs: hyphae-enriched and -suppressed cores were installed as pairs
## within plots. Models that compare core types WITHIN plots (core-type term,
## Model B) use only cores in complete pairs, so a plot's two pseudobulk sums
## are built from the same pair locations and pair-to-pair variation cancels
## out of the core contrast. Model A compares plots within each core type, so
## it uses all cores.
##
## Models
##   Core-type term:  ~ plot + core                        (complete pairs)
##   Model A:         ~ treatment_block + core + core:N_20  (all cores)
##       N slope within each core type
##   Model B:         ~ plot + core + fungi_N               (complete pairs)
##       difference in N slopes between core types, tested within plots
##       (fungi_N = N_20 in hyphae-enriched cores, 0 in suppressed cores)
##   Sensitivity:     ~ pair + core on individual cores in complete pairs
##       (core type tested within pairs, as in the mixed models and PERMANOVA)
##
## N_20 = N rate / 20, so slopes are log2 fold change per 20 kg N ha-1.
## `core` (Suppressed = reference, Enriched) is a syntactic copy of core_type
## for DESeq2 coefficient names; core_type keeps the display labels.
## Block is fixed (it also absorbs harvest date, which is confounded with block).

library(DESeq2)
library(apeglm)     # lfcShrink(type = "apeglm")

# DESeq2's Bioconductor dependencies mask several dplyr verbs; restore them.
# Re-run these lines if an earlier section is re-run after DESeq2 is loaded.
select <- dplyr::select
rename <- dplyr::rename
desc   <- dplyr::desc
slice  <- dplyr::slice
first  <- dplyr::first
count  <- dplyr::count    # matrixStats::count (via MatrixGenerics)

## Helper functions

### Design and sample units

# get_meta() (shared labels, plot, pair) + N_20 and the DESeq2 core factor.
# Row names stay as sample names.
deseq_meta <- function(ps) {
  sam <- get_meta(ps)
  sam$N_20 <- as.numeric(as.character(sam$N_treatment)) / 20
  sam$core <- factor(ifelse(sam$core_type == "Hyphae Enriched", "Enriched", "Suppressed"),
                     levels = c("Suppressed", "Enriched"))
  sam
}

# TRUE for cores whose pair partner (the other core type) is also present
in_complete_pair <- function(sam) {
  ave(as.integer(sam$core), sam$pair, FUN = function(x) length(unique(x))) == 2
}

# taxa x samples count matrix
otu_counts <- function(ps) {
  otu <- as(otu_table(ps), "matrix")
  if (!taxa_are_rows(ps)) otu <- t(otu)
  otu
}

# Prevalence filter, then drop units left with no counts (their size factor
# would be NA). paired = TRUE re-applies the complete-pair rule afterwards,
# so a core whose partner was dropped is dropped too.
finish_units <- function(counts, meta, ps, min_prev, paired = FALSE) {
  counts <- counts[rowMeans(counts > 0) >= min_prev, , drop = FALSE]
  keep   <- colSums(counts) > 0
  if (paired) keep[keep] <- in_complete_pair(meta[keep, ])
  if (any(!keep)) message("Dropped ", sum(!keep), " unit(s) with no counts after filtering",
                          if (paired) " (or whose pair partner was dropped)", ": ",
                          paste(colnames(counts)[!keep], collapse = ", "))
  counts <- counts[, keep, drop = FALSE]
  meta   <- droplevels(meta[keep, ])
  storage.mode(counts) <- "integer"
  tax <- as.data.frame(tax_table(ps)) %>% rownames_to_column("OTU")
  list(counts = counts, meta = meta, tax = tax)
}

# Sum raw counts of replicate cores within each plot x core type.
#   paired_only: keep only cores in complete pairs (within-plot models)
#   min_prev:    drop taxa present in fewer than this fraction of units
pseudobulk <- function(ps, paired_only = FALSE, min_prev = 0.10) {
  otu <- otu_counts(ps)
  sd  <- deseq_meta(ps)[colnames(otu), ]
  if (paired_only) sd <- sd[in_complete_pair(sd), ]
  otu <- otu[, rownames(sd), drop = FALSE]
  
  sd$unit <- paste(sd$plot, sd$core, sep = "_")              # plot x core type
  counts  <- t(rowsum(t(otu), group = sd$unit))              # taxa x units
  
  meta <- sd[!duplicated(sd$unit),
             c("unit", "plot", "N_treatment", "N_20", "treatment_block",
               "core", "core_type")]
  meta$n_cores   <- as.vector(table(sd$unit)[meta$unit])     # cores summed per unit
  rownames(meta) <- meta$unit
  meta <- droplevels(meta[colnames(counts), ])
  meta$fungi_N <- as.numeric(meta$core == "Enriched") * meta$N_20   # Model B term
  
  finish_units(counts, meta, ps, min_prev)
}

# Individual cores in complete pairs (within-pair sensitivity model)
paired_cores <- function(ps, min_prev = 0.10) {
  otu <- otu_counts(ps)
  sd  <- deseq_meta(ps)[colnames(otu), ]
  sd  <- droplevels(sd[in_complete_pair(sd),
                       c("plot", "pair", "N_treatment", "N_20",
                         "treatment_block", "core", "core_type")])
  finish_units(otu[, rownames(sd), drop = FALSE], sd, ps, min_prev, paired = TRUE)
}

### Model fitting and results

run_deseq <- function(units, design) {
  meta <- units$meta
  for (v in base::intersect(c("plot", "pair", "treatment_block"), names(meta))) {
    levels(meta[[v]]) <- make.names(levels(meta[[v]]), unique = TRUE)  # syntactic levels for DESeq2
  }
  dds <- DESeqDataSetFromMatrix(countData = units$counts, colData = meta,
                                design = design)
  dds <- estimateSizeFactors(dds, type = "poscounts")   # robust to zero-heavy amplicon data
  DESeq(dds, sfType = "poscounts", quiet = TRUE)
}

# Shrunken log2 fold changes (apeglm) with the ordinary Wald p-values.
# Shrinkage tames extreme estimates (e.g. when a taxon disappears from part
# of the gradient). It changes effect sizes only; padj is unchanged.
shrink <- function(dds, coef) {
  lfcShrink(dds, coef = coef, res = results(dds, name = coef),
            type = "apeglm", quiet = TRUE)
}

# DESeq2 results -> data frame with taxonomy, a readable taxon label
# (lowest classified rank in `ranks`) and -log10(padj)
tidy_results <- function(res, tax, ranks) {
  df  <- as.data.frame(res) %>% rownames_to_column("OTU") %>% left_join(tax, by = "OTU")
  lab <- rep(NA_character_, nrow(df))
  for (r in ranks) {
    v  <- df[[r]]
    ok <- is.na(lab) & !is.na(v) & !str_detect(v, "Unclassified")
    lab[ok] <- v[ok]
  }
  last <- df[[ranks[length(ranks)]]]
  lab[is.na(lab)] <- last[is.na(lab)]                 # fall back to highest rank
  df$Taxon_label <- lab
  df$log10padj   <- -log10(df$padj + 1e-300)
  df
}

# Core-type term: within plots on paired pseudobulk units
core_type_da <- function(pb_paired, ranks) {
  dds <- run_deseq(pb_paired, ~ plot + core)
  list(dds = dds,
       res = tidy_results(shrink(dds, "core_Enriched_vs_Suppressed"), pb_paired$tax, ranks))
}

# Sensitivity: core type within pairs on individual cores
core_type_da_pairs <- function(cores, ranks) {
  dds <- run_deseq(cores, ~ pair + core)
  list(dds = dds,
       res = tidy_results(shrink(dds, "core_Enriched_vs_Suppressed"), cores$tax, ranks))
}

# Agreement between the pseudobulk and within-pair core-type results
compare_core_models <- function(res_pb, res_pair, alpha = 0.05) {
  j <- inner_join(dplyr::select(res_pb,   OTU, lfc_pb   = log2FoldChange, padj_pb   = padj),
                  dplyr::select(res_pair, OTU, lfc_pair = log2FoldChange, padj_pair = padj),
                  by = "OTU")
  list(significant = table(pseudobulk  = !is.na(j$padj_pb)   & j$padj_pb   < alpha,
                           within_pair = !is.na(j$padj_pair) & j$padj_pair < alpha),
       lfc_cor     = cor(j$lfc_pb, j$lfc_pair, use = "complete.obs"))
}

# Model A (N slope within each core type, all cores) +
# Model B (slope difference, tested within plots on complete pairs)
within_core_slopes <- function(pb_all, pb_paired, ranks) {
  ddsA <- run_deseq(pb_all, ~ treatment_block + core + core:N_20)
  stopifnot(all(c("coreSuppressed.N_20", "coreEnriched.N_20") %in% resultsNames(ddsA)))
  
  res_within_core <- bind_rows(
    tidy_results(shrink(ddsA, "coreSuppressed.N_20"), pb_all$tax, ranks) %>%
      mutate(core_type = "Hyphae Suppressed"),
    tidy_results(shrink(ddsA, "coreEnriched.N_20"),   pb_all$tax, ranks) %>%
      mutate(core_type = "Hyphae Enriched")
  )
  
  ddsB <- run_deseq(pb_paired, ~ plot + core + fungi_N)
  slope_diff <- tidy_results(shrink(ddsB, "fungi_N"), pb_paired$tax, ranks)
  
  list(res_within_core = res_within_core, slope_diff = slope_diff,
       ddsA = ddsA, ddsB = ddsB)
}

### Linking taxa to FLNF

# Summed relative abundance of a set of OTUs in each sample
group_rel_abund <- function(ps, otus) {
  rel <- transform_sample_counts(ps, function(x) x / sum(x))
  m   <- as(otu_table(rel), "matrix")
  if (!taxa_are_rows(rel)) m <- t(m)
  colSums(m[rownames(m) %in% otus, , drop = FALSE])
}

# Core-level data for FLNF models: design factors from prep_design(), plus the
# summed relative abundance (rarefied counts) of each named group of OTUs.
# Groups with no variation (e.g. no OTUs) are dropped.
otu_group_data <- function(ps_rfy, groups) {
  ps <- soil_cores(ps_rfy)
  df <- sample_df(ps) %>%
    prep_design(numeric_vars = "n_fixation") %>%
    filter(is.na(n_fixation) | n_fixation > 0)          # log-scale model needs FLNF > 0
  for (g in names(groups)) df[[g]] <- unname(group_rel_abund(ps, groups[[g]])[df$Site])
  keep <- vapply(names(groups), function(g) isTRUE(sd(df[[g]], na.rm = TRUE) > 0), logical(1))
  attr(df, "groups") <- names(groups)[keep]
  df
}

### Plotting helpers

da_theme <- theme_bw(base_size = 16) +
  theme(axis.text.x  = element_text(size = 16, colour = "black"),
        axis.text.y  = element_text(size = 14, colour = "black"),
        axis.title   = element_text(size = 16, colour = "black"),
        strip.text   = element_text(size = 16, colour = "black"),
        legend.text  = element_text(size = 14, colour = "black"),
        legend.title = element_text(size = 16, colour = "black"),
        panel.grid.major = element_line(colour = "grey85"),
        panel.grid.minor = element_blank())

# Classify each taxon row as consistently positive, divergent (OTUs with
# opposite signs) or consistently negative, and order rows top -> bottom.
# OTUs sharing a Taxon_label share one row. `group` (e.g. "core_type") orders
# rows separately within each facet. Adds `response` and the y variable `Taxon_y`.
order_by_response <- function(df, group = NULL,
                              labels = c("Positive", "Divergent", "Negative")) {
  grp <- c(group, "Taxon_label")
  df <- df %>%
    group_by(across(all_of(grp))) %>%
    mutate(response = case_when(all(log2FoldChange > 0) ~ labels[1],
                                all(log2FoldChange < 0) ~ labels[3],
                                TRUE                    ~ labels[2]),
           mean_lfc = mean(log2FoldChange)) %>%
    ungroup() %>%
    mutate(response = factor(response, levels = labels))
  
  df$.key <- if (is.null(group)) df$Taxon_label else
    paste(df$Taxon_label, df[[group]], sep = "___")
  
  lev <- df %>%                                        # first level drawn at the bottom
    distinct(.key, response, mean_lfc) %>%
    arrange(dplyr::desc(response), mean_lfc) %>%
    pull(.key)
  df$Taxon_y <- factor(df$.key, levels = unique(lev))
  dplyr::select(df, -.key)
}

# top_n strongest positive and negative significant N responders per core type
top_responders <- function(res, top_n = 20) {
  sig <- res %>% filter(!is.na(padj), padj < 0.05) %>% group_by(core_type)
  bind_rows(slice_max(sig, log2FoldChange, n = top_n, with_ties = FALSE),
            slice_min(sig, log2FoldChange, n = top_n, with_ties = FALSE)) %>%
    ungroup() %>%
    distinct(core_type, OTU, .keep_all = TRUE) %>%
    mutate(core_type = factor(core_type, levels = c("Hyphae Suppressed", "Hyphae Enriched")))
}

# Lollipop plot of significant taxa, ordered by response type.
#   facet:  column to facet by (e.g. "core_type"); NULL for a single panel
#   labels: names of the three response groups (top to bottom)
#   strips: single panel only - label response groups with side strips;
#           otherwise groups are separated by dashed lines
#   size_by_padj: scale point size by -log10(padj)
plot_lollipop <- function(df, fill_var, fill_scale, y_lab, fill_lab, x_lab = "",
                          facet = NULL,
                          labels = c("Positive", "Divergent", "Negative"),
                          strips = is.null(facet), size_by_padj = FALSE) {
  df <- order_by_response(df, group = facet, labels = labels)
  facet_vars <- if (is.null(facet)) character(0) else facet
  
  seps <- df %>%                                       # dashed separators between groups
    distinct(across(all_of(c(facet_vars, "Taxon_y", "response")))) %>%
    dplyr::count(across(all_of(c(facet_vars, "response"))), name = "n_rows") %>%
    group_by(across(all_of(facet_vars))) %>%
    arrange(dplyr::desc(response), .by_group = TRUE) %>%
    mutate(yintercept = cumsum(n_rows) + 0.5) %>%
    slice_head(n = -1) %>%
    ungroup()
  
  points <- if (size_by_padj) {
    geom_point(aes(fill = .data[[fill_var]], size = log10padj), shape = 21)
  } else {
    geom_point(aes(fill = .data[[fill_var]]), shape = 21, size = 4)
  }
  
  p <- ggplot(df, aes(x = log2FoldChange, y = Taxon_y)) +
    geom_vline(xintercept = 0, colour = "grey40") +
    geom_segment(aes(x = 0, xend = log2FoldChange, yend = Taxon_y),
                 linewidth = 1, alpha = 0.8, colour = "black") +
    points +
    fill_scale +
    guides(size = "none") +
    scale_y_discrete(labels = function(x) sub("___.*$", "", x)) +
    scale_x_continuous(expand = expansion(mult = c(0.1, 0.1))) +
    labs(x = x_lab, y = y_lab, fill = fill_lab) +
    da_theme
  
  if (!is.null(facet)) {
    p <- p + facet_wrap(vars(.data[[facet]]), scales = "free")
  } else if (strips) {
    p <- p + facet_grid(rows = vars(response), scales = "free_y", space = "free_y") +
      theme(strip.text.y = element_text(angle = 0))
  }
  if (!strips || !is.null(facet)) {
    p <- p + geom_hline(data = seps, aes(yintercept = yintercept),
                        linetype = 2, colour = "grey50")
  }
  p
}

core_da_labels <- c("Enriched", "Divergent", "Depleted")   # core-type lollipop groups
N_slope_lab    <- expression("Slope across nitrogen gradient (log"[2]*
                               " fold change per 20 kg N ha"^-1*")")


## nifH ######################################

### Sample units
nif_subset    <- soil_cores(nif)                          # raw counts, soil cores only
nif_pb        <- pseudobulk(nif_subset)                   # all cores (Model A)
nif_pb_paired <- pseudobulk(nif_subset, paired_only = TRUE)  # complete pairs (core term, Model B)
nif_tax_ranks <- c("Genus", "Family", "Order", "Class", "Phylum")

check_pairs(deseq_meta(nif_subset))                       # unpaired cores are left out of nif_pb_paired
table(nif_pb$meta$core_type);        table(nif_pb$meta$n_cores)
table(nif_pb_paired$meta$core_type); table(nif_pb_paired$meta$n_cores)
dim(nif_pb$counts); dim(nif_pb_paired$counts)              # taxa x pseudobulk units


### Core-type term
# Diazotroph differential abundance between hyphae-enriched and -suppressed
# cores, compared within plots on complete pairs, averaged across N
nif_core_da <- core_type_da(nif_pb_paired, nif_tax_ranks)
resultsNames(nif_core_da$dds)
res_coretype <- nif_core_da$res

nif_top_core <- res_coretype %>% filter(!is.na(padj), padj < 0.05)
nrow(nif_top_core)                                        # significant OTUs

# Sensitivity: within-pair model on individual cores
nif_cores_paired   <- paired_cores(nif_subset)
nif_core_da_pairs  <- core_type_da_pairs(nif_cores_paired, nif_tax_ranks)
compare_core_models(res_coretype, nif_core_da_pairs$res)

# Core-type responders (Figure 3A, assembled at the end of the script)
deseq_core_nif <- plot_lollipop(nif_top_core,
                                fill_var     = "Phylum",
                                fill_scale   = scale_fill_observable(),
                                y_lab        = "Diazotroph Taxon",
                                fill_lab     = "Diazotroph Phylum",
                                x_lab        = "log2 fold change (Hyphae Enriched vs Suppressed)",
                                labels       = core_da_labels,
                                size_by_padj = TRUE)
deseq_core_nif

# Do hyphae-enriched OTUs predict FLNF? (core level: plot and pair random effects)
nif_enriched_otus <- nif_top_core %>% filter(log2FoldChange > 0) %>% pull(OTU)
nif_enriched_df   <- otu_group_data(nif_rfy, list(Enriched_OTUs = nif_enriched_otus))

ggplot(nif_enriched_df, aes(x = Enriched_OTUs, y = n_fixation)) +
  geom_point() + geom_smooth(method = "lm")
ggplot(nif_enriched_df, aes(x = N_treatment, y = Enriched_OTUs, fill = core_type)) +
  geom_boxplot()

if (length(attr(nif_enriched_df, "groups")) > 0) {
  nif_enriched_flnf <- fit_models(nif_enriched_df, "n_fixation",
                                  covariates   = attr(nif_enriched_df, "groups"),
                                  log_response = TRUE, fit_glmm = FALSE)
  nif_enriched_flnf$anova_t2
}
# Previous result: the hyphae-enriched OTUs' relative abundance did not predict FLNF


### Core-type x N interaction

# Models A and B
nif_slopes      <- within_core_slopes(nif_pb, nif_pb_paired, nif_tax_ranks)
res_within_core <- nif_slopes$res_within_core   # N slope within each core type
nif_slope_diff  <- nif_slopes$slope_diff        # Enriched slope - Suppressed slope
resultsNames(nif_slopes$ddsA)
resultsNames(nif_slopes$ddsB)

# Significant responders by direction
nif_sig_N <- res_within_core %>%
  filter(!is.na(padj), padj < 0.05) %>%
  mutate(direction = if_else(log2FoldChange > 0, "Positive", "Negative"))
dplyr::count(nif_sig_N, core_type, direction)

# Top N responders per core type (Figure 4A)
deseq_nif <- plot_lollipop(top_responders(res_within_core),
                           fill_var   = "Phylum",
                           fill_scale = scale_fill_observable(),
                           y_lab      = "Diazotroph Taxon",
                           fill_lab   = "Diazotroph Phylum",
                           x_lab      = N_slope_lab,
                           facet      = "core_type")
deseq_nif


### Compare N responses between core types
# One row per OTU: N slope and padj in each core type (Model A), plus the
# within-plot slope-difference test (Model B)
nif_N_compare <- res_within_core %>%
  mutate(core = recode(core_type,
                       "Hyphae Enriched"   = "HyphaeEnriched",
                       "Hyphae Suppressed" = "HyphaeSuppressed")) %>%
  dplyr::select(OTU, core, log2FoldChange, padj, Phylum, Taxon_label) %>%
  pivot_wider(names_from = core, values_from = c(log2FoldChange, padj),
              names_sep = "__") %>%
  left_join(dplyr::select(nif_slope_diff, OTU, diff_lfc = log2FoldChange,
                          diff_pvalue = pvalue, diff_padj = padj), by = "OTU") %>%
  mutate(
    sig_hyphae   = !is.na(padj__HyphaeEnriched)   & padj__HyphaeEnriched   < 0.05,
    sig_nohyphae = !is.na(padj__HyphaeSuppressed) & padj__HyphaeSuppressed < 0.05,
    sig_group = case_when(
      sig_hyphae &  sig_nohyphae ~ "Significant in both",
      sig_hyphae & !sig_nohyphae ~ "Hyphae Enriched only",
      !sig_hyphae & sig_nohyphae ~ "Hyphae Suppressed only",
      TRUE                       ~ "Not significant"),
    min_padj       = pmin(padj__HyphaeEnriched, padj__HyphaeSuppressed, na.rm = TRUE),
    log10_min_padj = -log10(min_padj + 1e-300),
    Phylum = case_when(                                   # merge split phyla
      Phylum %in% c("Bacillota_A", "Bacillota_B", "Bacillota_G") ~ "Bacillota",
      Phylum %in% c("Desulfobacterota", "Desulfobacterota_C")     ~ "Desulfobacterota",
      TRUE ~ Phylum)
  )

nif_N_compare_sig <- filter(nif_N_compare, sig_group != "Not significant")  # significant in >= 1 core

# Exploratory: all significant responders
ggplot(nif_N_compare_sig, aes(x = log2FoldChange__HyphaeEnriched,
                              y = log2FoldChange__HyphaeSuppressed)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey60") +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey60") +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "black") +
  geom_point(aes(colour = sig_group), alpha = 0.7, size = 2.5) +
  labs(x = "N response in Hyphae Enriched\n(log2 fold change per 20 kg N ha⁻¹)",
       y = "N response in Hyphae Suppressed\n(log2 fold change per 20 kg N ha⁻¹)",
       colour = NULL) +
  theme_bw(base_size = 16)

# Negative responders in both core types
nif_N_neg <- nif_N_compare_sig %>%
  filter(log2FoldChange__HyphaeEnriched < 0, log2FoldChange__HyphaeSuppressed < 0)

ggplot(nif_N_neg, aes(x = log2FoldChange__HyphaeEnriched,
                      y = log2FoldChange__HyphaeSuppressed, fill = Phylum)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "black") +
  geom_point(shape = 21, colour = "black", stroke = 0.4, alpha = 0.9, size = 3.5) +
  scale_fill_manual(values = pal_observable()(11)) +
  labs(x = "N response in Hyphae Enriched",
       y = "N response in Hyphae Suppressed", fill = "Phylum") +
  theme_bw(base_size = 16) +
  theme(panel.grid.minor = element_blank(), legend.key = element_blank(),
        axis.text = element_text(colour = "black"),
        axis.title = element_text(colour = "black"))

### Classify differential N responders
# "test":        slope-difference test (Model B), BH across all taxa
# "test_subset": same test, BH re-applied within the plotted OTUs   <- used
# "nominal":     unadjusted p < 0.05 (exploratory)
# "z":           |z| > 1 for distance from the 1:1 line (descriptive)
outlier_method <- "test_subset"
fdr_cut <- 0.10   # FDR threshold for the test-based methods
lfc_cut <- 0      # optional minimum |slope difference| (log2 per 20 kg N)

grp_enr <- "More N-sensitive in Hyphae Enriched"    # above the 1:1 line
grp_sup <- "More N-sensitive in Hyphae Suppressed"  # below the 1:1 line
grp_std <- "Similar response"

# line_residual (and diff_lfc) = Enriched slope - Suppressed slope. Both slopes
# are negative here, so a NEGATIVE difference = steeper decline in
# hyphae-enriched cores; POSITIVE = steeper in hyphae-suppressed cores.
nif_N_neg <- nif_N_neg %>%
  mutate(
    line_residual    = log2FoldChange__HyphaeEnriched - log2FoldChange__HyphaeSuppressed,
    z_residual       = as.numeric(scale(line_residual)),
    diff_padj_subset = p.adjust(diff_pvalue, method = "BH"),   # BH within plotted OTUs
    flagged = case_when(
      outlier_method == "test"        ~ !is.na(diff_padj)        & diff_padj        < fdr_cut,
      outlier_method == "test_subset" ~ !is.na(diff_padj_subset) & diff_padj_subset < fdr_cut,
      outlier_method == "nominal"     ~ !is.na(diff_pvalue)      & diff_pvalue      < 0.05,
      outlier_method == "z"           ~ abs(z_residual) > 1
    ) & abs(line_residual) >= lfc_cut,
    direction     = if (outlier_method == "z") line_residual else diff_lfc,
    outlier_group = case_when(flagged & direction < 0 ~ grp_enr,
                              flagged & direction > 0 ~ grp_sup,
                              TRUE                    ~ grp_std),
    outlier_group = factor(outlier_group, levels = c(grp_std, grp_enr, grp_sup))
  )

table(nif_N_neg$outlier_group)

outliers_above <- nif_N_neg %>% filter(outlier_group == grp_enr) %>% arrange(line_residual)
outliers_below <- nif_N_neg %>% filter(outlier_group == grp_sup) %>% arrange(dplyr::desc(line_residual))
outliers_above
outliers_below

# N responses with labelled outliers (Figure S4A)
N_curve <- ggplot(nif_N_neg, aes(x = log2FoldChange__HyphaeEnriched,
                                 y = log2FoldChange__HyphaeSuppressed)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, linewidth = 1) +
  geom_point(aes(fill = outlier_group), shape = 21, size = 4, stroke = 1, alpha = 0.9) +
  geom_text_repel(data = filter(nif_N_neg, outlier_group != grp_std),
                  box.padding = 0.6, min.segment.length = 0,
                  aes(label = Taxon_label), size = 4, max.overlaps = Inf) +
  scale_fill_manual(values = setNames(c("white", "red", "blue"),
                                      c(grp_std, grp_enr, grp_sup)),
                    drop = TRUE) +
  labs(x = "N response in Hyphae Enriched (log2 fold change)",
       y = "N response in Hyphae Suppressed (log2 fold change)",
       fill = "Nitrogen Response Outliers") +
  theme_bw(base_size = 18)
N_curve

### Outlier-group abundance vs FLNF
# Core-level model with the same structure as the other FLNF models:
# core_type * N + block, random plot and pair, log(FLNF), Kenward-Roger tests.
nif_outlier_df <- otu_group_data(nif_rfy, list(
  Hyphae_Enriched_Sensitive   = outliers_above$OTU,
  Hyphae_Suppressed_Sensitive = outliers_below$OTU
))

nif_outlier_flnf <- fit_models(nif_outlier_df, "n_fixation",
                               covariates   = attr(nif_outlier_df, "groups"),
                               log_response = TRUE, fit_glmm = FALSE)
check_fit(nif_outlier_flnf)
nif_outlier_flnf$anova_t2

# Outlier abundance vs FLNF (Figure S4B)
nif_outlier_long <- nif_outlier_df %>%
  filter(!is.na(n_fixation)) %>%
  dplyr::select(Site, n_fixation, Hyphae_Enriched_Sensitive, Hyphae_Suppressed_Sensitive) %>%
  pivot_longer(c(Hyphae_Enriched_Sensitive, Hyphae_Suppressed_Sensitive),
               names_to = "Group", values_to = "Relative_Abundance") %>%
  mutate(Group = recode(Group,
                        "Hyphae_Enriched_Sensitive"   = "Hyphae Enriched Outliers",
                        "Hyphae_Suppressed_Sensitive" = "Hyphae Suppressed Outliers")) %>%
  filter(Relative_Abundance > 0)     # cores where the group was detected

# Spearman correlations for the panel labels. Descriptive only (cores are not
# independent); inference comes from nif_outlier_flnf above.
cor_stats <- nif_outlier_long %>%
  group_by(Group) %>%
  summarise(rho  = cor(Relative_Abundance, n_fixation, method = "spearman"),
            pval = cor.test(Relative_Abundance, n_fixation,
                            method = "spearman", exact = FALSE)$p.value,
            .groups = "drop") %>%
  mutate(label = paste0("rho = ", round(rho, 2), "\nP = ", signif(pval, 2)),
         x = 0.12, y = c(15, 17)[seq_len(n())])      # label positions: adjust as needed

outlier_cols <- c("Hyphae Enriched Outliers" = "red", "Hyphae Suppressed Outliers" = "blue")

flnf_outlier_plot <- ggplot(nif_outlier_long, aes(x = Relative_Abundance, y = n_fixation,
                                                  colour = Group)) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 1.5) +
  geom_point(aes(fill = Group), size = 3, alpha = 0.8, shape = 21, colour = "black") +
  geom_text(data = cor_stats, aes(x = x, y = y, label = label, colour = Group),
            hjust = 0, size = 5, show.legend = FALSE) +
  scale_fill_manual(values = outlier_cols, drop = TRUE) +
  scale_colour_manual(values = outlier_cols) +
  ylab(expression(FLNF ~ ("µg" ~ N ~ g^-1 ~ dry ~ soil ~ day^-1))) +
  labs(x = "Relative abundance of outlier OTUs", colour = NULL, fill = NULL) +
  theme_bw(base_size = 18) +
  theme(panel.grid.minor = element_blank())
flnf_outlier_plot

# Figure S4: N-response outliers and FLNF (legend below)
common_theme <- theme(plot.margin  = margin(10, 10, 10, 10),
                      legend.position = "none",
                      axis.title.y = element_text(margin = margin(r = 12)),
                      axis.title.x = element_text(margin = margin(t = 12)))

fig_N_outliers <- plot_grid(
  plot_grid(N_curve + common_theme, flnf_outlier_plot + common_theme,
            ncol = 2, labels = c("A", "B")),
  get_legend(N_curve),
  ncol = 1, rel_heights = c(1, 0.12)
)
fig_N_outliers

ggsave("Figures/Figure_S4.png", plot = fig_N_outliers, dpi = 600, width = 16, height = 7.5)


## LSU ######################################

### Sample units (same outlier cores dropped as in the beta diversity analysis)
lsu_subset    <- soil_cores(lsu, drop = lsu_outliers)     # raw counts, soil cores only
lsu_subset    <- prune_samples(sample_sums(lsu_subset) > 1, lsu_subset)
lsu_pb        <- pseudobulk(lsu_subset)
lsu_pb_paired <- pseudobulk(lsu_subset, paired_only = TRUE)
lsu_tax_ranks <- c("genus", "family", "phylum")           # LSU taxonomy uses lower case

check_pairs(deseq_meta(lsu_subset))
table(lsu_pb$meta$core_type);        table(lsu_pb$meta$n_cores)
table(lsu_pb_paired$meta$core_type); table(lsu_pb_paired$meta$n_cores)
dim(lsu_pb$counts); dim(lsu_pb_paired$counts)


### Core-type term
lsu_core_da <- core_type_da(lsu_pb_paired, lsu_tax_ranks)
resultsNames(lsu_core_da$dds)
res_coretype_lsu <- lsu_core_da$res

lsu_top_core <- res_coretype_lsu %>% filter(!is.na(padj), padj < 0.05)
lsu_top_core %>% dplyr::count(direction = if_else(log2FoldChange > 0, "Enriched", "Depleted"), genus)

# Sensitivity: within-pair model on individual cores
lsu_cores_paired  <- paired_cores(lsu_subset)
lsu_core_da_pairs <- core_type_da_pairs(lsu_cores_paired, lsu_tax_ranks)
compare_core_models(res_coretype_lsu, lsu_core_da_pairs$res)

# Core-type responders (Figure S3)
deseq_core_lsu <- plot_lollipop(lsu_top_core,
                                fill_var     = "family",
                                fill_scale   = scale_fill_npg(),
                                y_lab        = "AMF Taxon",
                                fill_lab     = "AMF Family",
                                x_lab        = "log2 fold change (Hyphae Enriched vs Suppressed)",
                                labels       = core_da_labels,
                                size_by_padj = TRUE)
deseq_core_lsu

ggsave("Figures/Figure_S3.png", plot = deseq_core_lsu, dpi = 600, width = 10, height = 7)


### Core-type x N interaction
lsu_slopes          <- within_core_slopes(lsu_pb, lsu_pb_paired, lsu_tax_ranks)
res_within_core_lsu <- lsu_slopes$res_within_core
lsu_slope_diff      <- lsu_slopes$slope_diff
resultsNames(lsu_slopes$ddsA)
resultsNames(lsu_slopes$ddsB)

lsu_sig_N <- res_within_core_lsu %>%
  filter(!is.na(padj), padj < 0.05) %>%
  mutate(direction = if_else(log2FoldChange > 0, "Positive", "Negative"))
dplyr::count(lsu_sig_N, core_type, direction)

# AMF taxa whose N response differs between core types
filter(lsu_slope_diff, !is.na(padj), padj < 0.05)

# Top N responders per core type (Figure 4B)
deseq_lsu <- plot_lollipop(top_responders(res_within_core_lsu),
                           fill_var   = "family",
                           fill_scale = scale_fill_npg(),
                           y_lab      = "AMF Taxon",
                           fill_lab   = "AMF Family",
                           x_lab      = N_slope_lab,
                           facet      = "core_type")
deseq_lsu


## Figure 4: N responders, nifH and LSU ######
# Remove legends, align panels vertically, then add legends back
aligned <- cowplot::align_plots(deseq_nif + theme(legend.position = "none"),
                                deseq_lsu + theme(legend.position = "none"),
                                align = "v", axis = "l")

nif_row <- plot_grid(aligned[[1]], cowplot::get_legend(deseq_nif),
                     ncol = 2, rel_widths = c(1, 0.25),
                     labels = c("A", ""), label_fontface = "bold", label_size = 16)
lsu_row <- plot_grid(aligned[[2]], cowplot::get_legend(deseq_lsu),
                     ncol = 2, rel_widths = c(1, 0.25),
                     labels = c("B", ""), label_fontface = "bold", label_size = 16)

fig_deseq_N <- plot_grid(nif_row, lsu_row, ncol = 1)
fig_deseq_N

ggsave("Figures/Figure_4.png", plot = fig_deseq_N, dpi = 600, width = 17, height = 10)



## ============================== Co-occurrence Networks (SpiecEasi) =======================================
## AMF-diazotroph networks. Diazotrophs: genus level; AMF: ASV level by
## default (NET_LSU_RANK = "genus" for genus level). Each domain has its own
## pool filter, so the finer AMF resolution does not inflate the taxon count.
##
##   A. Core type (all N pooled)        2 networks  (Hyphae Suppressed, Hyphae Enriched)
##   B. N group x core type             6 networks  (Low/Intermediate/High x core type)
##   C. N group (core types pooled)     3 networks  (Low/Intermediate/High)
##
## Common steps for every network
##   - soil cores in both datasets, above the rarefaction depths
##     (nif_depth, lsu_depth) and without the LSU outlier cores
##   - ONE shared taxon pool for every network and analysis, chosen from all
##     network cores (prevalence and mean relative abundance, set per
##     domain: NET_MIN_PREV / NET_MIN_ABUND). Each network is fitted on the pool taxa
##     present in >= NET_MIN_CORES of its cores; absent taxa have no edges.
##     Metrics are reported on the whole pool (turnover + rewiring) and on
##     the present taxa only (rewiring among taxa that co-occur).
##   - clr within each domain; group means removed before inference
##       A, B: plot means           C: plot means + core-type effect
##   - MB neighbourhood selection, lambda chosen by StARS (pulsar)
##   - edge stability = share of StARS subsamples that select the edge
##
## Permutation tests (networks refit at the observed lambdas)
##   core type: core-type labels swapped within core PAIRS        (A, B)
##              (same restriction as the PERMANOVA core-type test)
##   N group:   N treatments shuffled among plots within blocks   (B, C)
##
## N groups (NET_N_BREAKS): Low = 0-28, Intermediate = 56-84, High = 140-196.
## Liming began at 84 kg N, so "Intermediate" mixes unlimed (56) and limed (84)
## plots; change NET_N_BREAKS / NET_N_GROUPS to group by liming instead.
## Core-count check (N contrasts only): N groups differ in their number of
## cores, and the number of taxa present depends on it. Within each core type
## (and pooled), every N group is repeatedly subsampled to the smallest group's
## core count and the present taxa are recounted. Core-type contrasts need no
## check: swapping labels within pairs keeps each network's core count, so the
## permutation null already has the observed imbalance. Densities are not
## checked this way (refitting at a fixed lambda with fewer cores inflates them).
## Runtime: with NET_RUN_TESTS = TRUE and NET_N_PERM = 999, analysis B refits
## 6 networks x 999 x 2 test types. Set NET_RUN_TESTS <- FALSE while iterating.

library(SpiecEasi)
library(pulsar)
library(ggnewscale)

### Parameters
NET_LSU_RANK   <- "asv"  # AMF resolution: "asv" or "genus"
# Taxon pool filters, per domain: present in >= NET_MIN_PREV of all network cores
# and mean relative abundance (within domain) >= NET_MIN_ABUND. Keep the pool
# at <= ~100 taxa in total; ASVs need a stricter filter than genera.
NET_MIN_PREV   <- c(nif = 0.025, lsu = 0.10)
NET_MIN_ABUND  <- c(nif = 0.00025, lsu = 0.001)
NET_MIN_CORES  <- 3      # a pool taxon enters a network if present in >= this many of its cores
NET_LMIN_RATIO <- 1e-1   # end of lambda path (lmax x NET_LMIN_RATIO)
NET_NLAMBDA    <- 100
NET_REP_NUM    <- 100    # StARS subsamples
NET_EDGE_FREQ  <- 0.50   # stable edge: selected in >= 50% of subsamples
NET_RUN_TESTS  <- TRUE   # permutation tests
NET_RUN_SUBSAMPLE <- TRUE  # core-count check
NET_SUB_REPS   <- 999    # subsamples (counts only, so fast)
NET_N_PERM     <- 999
NET_N_BREAKS   <- c(-Inf, 28, 84, Inf)               # kg N ha-1
NET_N_GROUPS   <- c("Low", "Intermediate", "High")

net_core_types <- c("Hyphae Suppressed", "Hyphae Enriched")   # shared core labels
amf_group_rank <- if (NET_LSU_RANK == "genus") "family" else "genus"   # AMF node colour
metric_cols    <- c("n_present", "nif_present", "amf_present",   # turnover (total, by domain)
                    "density", "cross_density",                   # pool: turnover + rewiring
                    "density_present", "cross_density_present",   # present taxa: rewiring
                    "cross_frac", "amf_amf_frac", "diaz_diaz_frac", "positive_frac")

## Helper functions

### Data

# Aggregate to a taxonomic rank (taxonomy bookkeeping columns dropped first)
glom_rank <- function(ps, rank) {
  tax <- as.data.frame(as(tax_table(ps), "matrix"), stringsAsFactors = FALSE)
  tax <- tax[, !(colnames(tax) %in% c("Taxon", "Confidence")), drop = FALSE]
  tax_table(ps) <- tax_table(as.matrix(tax))
  tax_glom(ps, taxrank = rank)
}

# Taxa with a real name at `rank` (tax_glom lumps "Unclassified ..." strings
# into single artificial nodes, which become spurious hubs)
classified <- function(ps, rank) {
  g <- as.character(tax_table(ps)[, rank])
  taxa_names(ps)[!is.na(g) & g != "" & !grepl("unclassified|uncultured", g, ignore.case = TRUE)]
}

# N treatment (kg N ha-1 labels) -> N group
N_to_group <- function(N) {
  factor(cut(as.numeric(as.character(N)), breaks = NET_N_BREAKS, labels = NET_N_GROUPS),
         levels = NET_N_GROUPS)
}

fill_na <- function(x) ifelse(is.na(x) | x == "", "Unclassified", x)

# Shared taxon pool from all network cores: prevalence and mean relative
# abundance (relative to the domain's total reads in each core)
pool_taxa <- function(m, min_prev, min_abund) {
  rel <- m / pmax(rowSums(m), 1)
  colnames(m)[colMeans(m > 0) >= min_prev & colMeans(rel) >= min_abund]
}
clr <- function(m) { l <- log(m + 1); l - rowMeans(l) }

# Taxon pool, clr matrix, presence matrix and labels (shared by all analyses)
make_data <- function() {
  nk <- pool_taxa(net_nif_counts, NET_MIN_PREV[["nif"]], NET_MIN_ABUND[["nif"]])
  lk <- pool_taxa(net_lsu_counts, NET_MIN_PREV[["lsu"]], NET_MIN_ABUND[["lsu"]])
  nc <- clr(net_nif_counts[, nk, drop = FALSE]); colnames(nc) <- paste0("nif|", nk)
  lc <- clr(net_lsu_counts[, lk, drop = FALSE]); colnames(lc) <- paste0("lsu|", lk)
  X  <- cbind(nc, lc)
  P  <- cbind(net_nif_counts[, nk, drop = FALSE], net_lsu_counts[, lk, drop = FALSE]) > 0
  colnames(P) <- colnames(X)
  list(
    X      = X,
    P      = P,                                        # presence (count > 0)
    source = setNames(c(rep("nif", length(nk)), rep("lsu", length(lk))), colnames(X)),
    label  = c(setNames(fill_na(net_nif_tax[nk, "Genus"]), colnames(nc)),
               setNames(if (NET_LSU_RANK == "genus") fill_na(net_lsu_tax[lk, "genus"]) else
                 paste0(fill_na(net_lsu_tax[lk, "genus"]), " (", lk, ")"), colnames(lc))),
    group  = c(setNames(fill_na(net_nif_tax[nk, "Class"]), colnames(nc)),     # diazotroph class
               setNames(fill_na(net_lsu_tax[lk, amf_group_rank]), colnames(lc)))  # AMF family / genus
  )
}

### Network inference

# Remove group means (factors in df), within one network's samples
residualize <- function(X, df) {
  df <- droplevels(as.data.frame(df))
  df <- df[, vapply(df, function(x) nlevels(factor(x)) > 1, logical(1)), drop = FALSE]
  mm <- if (ncol(df) == 0) matrix(1, nrow(X), 1) else
    model.matrix(~ ., data = as.data.frame(lapply(df, factor)))
  qr.resid(qr(mm), X)
}

# MB network; lambda = NULL -> StARS selection (returns edge stability), else refit
fit_mb <- function(X, lambda = NULL) {
  tryCatch({
    X <- X[, apply(X, 2, sd) > 1e-10, drop = FALSE]       # guard: constant columns
    stability <- NA_real_; freq <- NULL; path_end <- NA
    if (is.null(lambda)) {
      C    <- cor(X)
      lmax <- max(abs(C[upper.tri(C)]), na.rm = TRUE)
      lams <- exp(seq(log(lmax), log(lmax * NET_LMIN_RATIO), length.out = NET_NLAMBDA))
      ps   <- pulsar::pulsar(X, fun = SpiecEasi::sparseiCov,
                             fargs = list(method = "mb", lambda = lams),
                             criterion = "stars", thresh = 0.05,
                             rep.num = NET_REP_NUM, seed = SEED, refit = FALSE)
      opt       <- ps$stars$opt.index
      lambda    <- lams[opt]
      stability <- ps$stars$summary[opt]
      path_end  <- opt == length(lams)
      freq      <- as.matrix(ps$stars$merge[[opt]])
      dimnames(freq) <- list(colnames(X), colnames(X))
    }
    est  <- SpiecEasi::sparseiCov(X, method = "mb", lambda = lambda)
    beta <- as.matrix(SpiecEasi::symBeta(est$beta[[1]], mode = "maxabs"))
    dimnames(beta) <- list(colnames(X), colnames(X))
    list(adj = beta, freq = freq, lambda = lambda, stability = stability, path_end = path_end)
  }, error = function(e) { message("Fit failed: ", conditionMessage(e)); NULL })
}

# One network per level of `net` (a factor with fixed levels), each fitted on
# the pool taxa present in >= NET_MIN_CORES of that network's cores
build_networks <- function(D, net, resid_vars, lambdas = NULL) {
  setNames(lapply(levels(net), function(k) {
    idx     <- which(net == k)
    present <- colnames(D$X)[colSums(D$P[idx, , drop = FALSE]) >= NET_MIN_CORES]
    fit     <- fit_mb(residualize(D$X[idx, present, drop = FALSE],
                                  net_meta[idx, resid_vars, drop = FALSE]),
                      lambda = lambdas[[k]])
    if (!is.null(fit)) fit$present <- colnames(fit$adj)   # taxa actually fitted
    fit
  }), levels(net))
}

# Network labels
net_core  <- function(core)    factor(as.character(core), levels = net_core_types)
net_N     <- function(N)       factor(as.character(N), levels = NET_N_GROUPS)
net_Ncore <- function(N, core) factor(paste(N, core, sep = " | "),
                                      levels = as.vector(outer(NET_N_GROUPS, net_core_types,
                                                               paste, sep = " | ")))

# Split network keys ("Low | Hyphae Suppressed") into columns
separate_keys <- function(m, key_cols) {
  parts <- strsplit(as.character(m$network), " | ", fixed = TRUE)
  for (j in seq_along(key_cols)) m[[key_cols[j]]] <- vapply(parts, `[`, character(1), j)
  if ("N_group"   %in% key_cols) m$N_group   <- factor(m$N_group,   levels = NET_N_GROUPS)
  if ("core_type" %in% key_cols) m$core_type <- factor(m$core_type, levels = net_core_types)
  m
}

### Edges and metrics

edge_table <- function(fit, key, D) {
  if (is.null(fit)) return(NULL)
  adj <- fit$adj
  idx <- which(adj != 0 & upper.tri(adj), arr.ind = TRUE)
  if (nrow(idx) == 0) return(NULL)
  from <- rownames(adj)[idx[, 1]]; to <- colnames(adj)[idx[, 2]]
  tibble(network = key, from = from, to = to, weight = adj[idx],
         freq = if (is.null(fit$freq)) NA_real_ else fit$freq[idx]) %>%
    mutate(type = case_when(D$source[from] != D$source[to] ~ "AMF-diazotroph",
                            D$source[from] == "lsu"        ~ "AMF-AMF",
                            TRUE                           ~ "Diazotroph-diazotroph"),
           sign       = ifelse(weight > 0, "Positive", "Negative"),
           stable     = freq >= NET_EDGE_FREQ,
           from_taxon = D$label[from], to_taxon = D$label[to],
           from_group = D$group[from], to_group = D$group[to])
}

all_edge_tables <- function(fits, D, key_cols) {
  purrr::imap_dfr(fits, ~ edge_table(.x, .y, D)) %>% separate_keys(key_cols)
}

# Network metrics. Pool-based densities use every pool taxon as a possible node
# (turnover + rewiring); *_present densities use only the taxa present in that
# network (rewiring among co-occurring taxa).
metrics_of <- function(fits, D, key_cols) {
  n_nif <- sum(D$source == "nif"); n_lsu <- sum(D$source == "lsu"); p <- n_nif + n_lsu
  purrr::imap_dfr(fits, function(f, k) {
    if (is.null(f)) return(tibble(network = k))
    e  <- edge_table(f, k, D)
    ne <- if (is.null(e)) 0L else nrow(e)
    tp <- if (is.null(e)) character(0) else e$type
    st <- if (is.null(e)) logical(0)   else e$stable
    pp <- length(f$present)
    np <- sum(D$source[f$present] == "nif"); lp <- sum(D$source[f$present] == "lsu")
    tibble(network        = k,
           n_present      = pp,
           nif_present    = np,
           amf_present    = lp,
           n_edges        = ne,
           cross_edges    = sum(tp == "AMF-diazotroph"),
           stable_edges   = sum(st, na.rm = TRUE),
           stable_cross   = sum(st & tp == "AMF-diazotroph", na.rm = TRUE),
           density        = ne / (p * (p - 1) / 2),
           cross_density  = sum(tp == "AMF-diazotroph") / (n_nif * n_lsu),
           density_present       = ne / (pp * (pp - 1) / 2),
           cross_density_present = if (np * lp) sum(tp == "AMF-diazotroph") / (np * lp) else NA_real_,
           cross_frac     = if (ne) mean(tp == "AMF-diazotroph") else NA_real_,
           amf_amf_frac   = if (ne) mean(tp == "AMF-AMF") else NA_real_,
           diaz_diaz_frac = if (ne) mean(tp == "Diazotroph-diazotroph") else NA_real_,
           positive_frac  = if (ne) mean(e$weight > 0) else NA_real_,
           stability      = f$stability,
           path_end       = f$path_end)
  }) %>% separate_keys(key_cols)
}

# Present taxa, edges and mean neighbours per present taxon for each network
network_size <- function(met) {
  met %>%
    dplyr::select(network, n_present, nif_present, amf_present,
                  n_edges, stable_edges, stable_cross, path_end) %>%
    mutate(mean_neighbours = 2 * n_edges / n_present)
}

### Test statistics and permutation helpers

# Hyphae Enriched - Hyphae Suppressed (optionally within N groups)
core_diff <- function(m, by = character(0)) {
  m %>%
    dplyr::select(all_of(by), core_type, all_of(metric_cols)) %>%
    pivot_longer(all_of(metric_cols), names_to = "metric") %>%
    pivot_wider(names_from = core_type, values_from = value) %>%
    transmute(across(all_of(by)), metric, contrast = "Enriched - Suppressed",
              value = `Hyphae Enriched` - `Hyphae Suppressed`)
}

# Each N group minus the first (lowest) N group (optionally within core types)
n_contrast <- function(m, by = character(0)) {
  ref <- NET_N_GROUPS[1]
  m %>%
    dplyr::select(all_of(by), N_group, all_of(metric_cols)) %>%
    pivot_longer(all_of(metric_cols), names_to = "metric") %>%
    group_by(across(all_of(c(by, "metric")))) %>%
    mutate(value = value - value[N_group == ref][1]) %>%
    ungroup() %>%
    filter(N_group != ref) %>%
    transmute(across(all_of(by)), metric,
              contrast = paste(N_group, "-", ref), value)
}

# N x core type: does the Enriched - Suppressed difference change between N groups?
nx_core <- function(m) {
  n_contrast(m, "core_type") %>%
    pivot_wider(names_from = core_type, values_from = value) %>%
    transmute(metric, contrast = paste0("N x core (Enriched - Suppressed): ", contrast),
              value = `Hyphae Enriched` - `Hyphae Suppressed`)
}

perm_p <- function(obs, null) (sum(abs(null) >= abs(obs), na.rm = TRUE) + 1) /
  (sum(!is.na(null)) + 1)

run_null <- function(perms, make_net, D, resid_vars, lambdas, key_cols, stat) {
  purrr::map_dfr(seq_len(nrow(perms)), function(i) {
    stat(metrics_of(build_networks(D, make_net(perms[i, ]), resid_vars, lambdas), D, key_cols))
  }, .progress = TRUE)
}

summarise_test <- function(obs, null, keys) {
  obs %>%
    dplyr::rename(observed = value) %>%
    left_join(null, by = keys, relationship = "one-to-many") %>%
    group_by(across(all_of(keys))) %>%
    summarise(observed = dplyr::first(observed),
              p_value  = perm_p(dplyr::first(observed), value),
              n_perm   = sum(!is.na(value)), .groups = "drop")
}

### Core-count check (N contrasts)

count_cols <- c("n_present", "nif_present", "amf_present")

# Pool taxa present (in >= NET_MIN_CORES cores) among a set of cores, by domain
present_counts <- function(D, idx) {
  pres <- colSums(D$P[idx, , drop = FALSE]) >= NET_MIN_CORES
  tibble(n_present   = sum(pres),
         nif_present = sum(pres & D$source == "nif"),
         amf_present = sum(pres & D$source == "lsu"))
}

# Strata for the check: core types (B) or one pooled stratum (C)
count_strata <- function(by, n) {
  if (is.null(by)) factor(rep("Pooled", n)) else droplevels(factor(by))
}

# Present-taxa counts with every N group cut to the smallest group's core
# count within its stratum, `reps` times
subsample_counts <- function(D, groups, by = NULL, reps = NET_SUB_REPS, seed = SEED) {
  strata <- count_strata(by, length(groups))
  set.seed(seed)
  purrr::map_dfr(seq_len(reps), function(r) {
    purrr::map_dfr(levels(strata), function(st) {
      in_st <- strata == st
      n_min <- min(table(droplevels(groups[in_st])))
      purrr::map_dfr(levels(groups), function(g) {
        i <- which(in_st & groups == g)
        if (length(i) > n_min) i <- sample(i, n_min)
        bind_cols(tibble(rep = r, stratum = st, N_group = g, n_cores = n_min),
                  present_counts(D, i))
      })
    })
  }) %>%
    mutate(N_group = factor(N_group, levels = levels(groups)))
}

# Each N group minus the first, per rep (or for the observed data)
count_contrasts <- function(df) {
  ref <- levels(df$N_group)[1]
  df %>%
    pivot_longer(all_of(count_cols), names_to = "metric") %>%
    group_by(across(any_of(c("rep", "stratum", "metric")))) %>%
    mutate(value = value - value[N_group == ref][1]) %>%
    ungroup() %>%
    filter(N_group != ref) %>%
    mutate(contrast = paste(N_group, "-", ref))
}

# Observed (all cores) vs subsampled (equal cores) N contrasts: subsample
# mean, 95% interval, and share of subsamples with the observed sign
summarise_counts <- function(sub, D, groups, by = NULL) {
  strata <- count_strata(by, length(groups))
  obs <- purrr::map_dfr(levels(strata), function(st) {
    purrr::map_dfr(levels(groups), function(g) {
      bind_cols(tibble(stratum = st, N_group = g),
                present_counts(D, which(strata == st & groups == g)))
    })
  }) %>%
    mutate(N_group = factor(N_group, levels = levels(groups)))
  
  obs_c <- count_contrasts(obs) %>%
    dplyr::select(stratum, metric, contrast, observed = value)
  
  count_contrasts(sub) %>%
    left_join(obs_c, by = c("stratum", "metric", "contrast")) %>%
    group_by(stratum, metric, contrast) %>%
    summarise(n_cores   = dplyr::first(n_cores),
              observed  = dplyr::first(observed),
              sub_mean  = mean(value),
              sub_lo    = quantile(value, 0.025),
              sub_hi    = quantile(value, 0.975),
              same_sign = mean(sign(value) == sign(observed)),
              .groups   = "drop")
}

# Mean and 95% interval of present-taxa counts per N group at equal cores
count_table <- function(sub) {
  sub %>%
    group_by(stratum, N_group, n_cores) %>%
    summarise(across(all_of(count_cols),
                     list(mean = ~ mean(.x),
                          lo   = ~ quantile(.x, 0.025),
                          hi   = ~ quantile(.x, 0.975))),
              .groups = "drop")
}

### Figures

make_pal <- function(levels, base) {
  levels <- sort(unique(levels))
  cols <- if (length(levels) <= length(base)) base[seq_along(levels)] else
    colorRampPalette(base)(length(levels))
  setNames(cols, levels)
}

# Bipartite AMF-diazotroph network (all selected edges): AMF left, diazotrophs
# right, edge width = stability (share of StARS subsamples selecting the edge),
# edge colour = sign, node colour = AMF genus/family and diazotroph class
# (amf_pal, diaz_pal). Nodes are unlabelled; the colour legends identify them.
# legend_all = TRUE lists every colour group in the palettes, not just those
# drawn, so one legend can be shared by several network figures.
plot_cross_network <- function(edges, D, facet, title = NULL, node_size = 4.5,
                               legend_all = FALSE) {
  lim <- function(pal) if (legend_all) names(pal) else NULL       # limits: NULL = from data
  brk <- function(pal) if (legend_all) names(pal) else waiver()    # breaks: waiver = default
  ce <- edges %>%
    filter(type == "AMF-diazotroph") %>%
    mutate(amf  = ifelse(D$source[from] == "lsu", from, to),
           diaz = ifelse(D$source[from] == "lsu", to, from))
  if (nrow(ce) == 0) { message("No AMF-diazotroph edges to plot"); return(NULL) }
  
  facet_info <- ce %>% distinct(network, across(any_of(c("N_group", "core_type"))))
  
  # Node positions, sorted by colour group
  nodes <- bind_rows(ce %>% transmute(network, node = amf,  domain = "AMF"),
                     ce %>% transmute(network, node = diaz, domain = "Diazotroph")) %>%
    distinct() %>%
    mutate(group = D$group[node], label = D$label[node]) %>%
    group_by(network, domain) %>%
    arrange(group, label, .by_group = TRUE) %>%
    mutate(y = if (n() == 1) 0 else seq(1, -1, length.out = n()),
           x = ifelse(domain == "AMF", 0, 1)) %>%
    ungroup() %>%
    left_join(facet_info, by = "network")
  
  seg <- ce %>%
    left_join(nodes %>% dplyr::select(network, amf = node, x, y),                by = c("network", "amf")) %>%
    left_join(nodes %>% dplyr::select(network, diaz = node, xend = x, yend = y), by = c("network", "diaz"))
  
  ggplot() +
    geom_segment(data = seg,
                 aes(x = x, y = y, xend = xend, yend = yend,
                     linewidth = freq, colour = sign),
                 alpha = 0.8, lineend = "round") +
    scale_linewidth_continuous(name = "Edge stability", range = c(0.2, 2.5), limits = c(0, 1)) +
    scale_colour_manual(name = "Association",
                        values = c(Positive = "grey25", Negative = "firebrick3")) +
    geom_point(data = filter(nodes, domain == "AMF"),
               aes(x = x, y = y, fill = group), shape = 24, size = node_size, colour = "black") +
    scale_fill_manual(name = paste("AMF", amf_group_rank), values = amf_pal,
                      limits = lim(amf_pal), breaks = brk(amf_pal)) +
    ggnewscale::new_scale_fill() +
    geom_point(data = filter(nodes, domain == "Diazotroph"),
               aes(x = x, y = y, fill = group), shape = 21, size = node_size, colour = "black") +
    scale_fill_manual(name = "Diazotroph class", values = diaz_pal,
                      limits = lim(diaz_pal), breaks = brk(diaz_pal)) +
    geom_text(data = data.frame(x = c(0, 1), y = 1.12, lab = c("AMF", "Diazotrophs")),
              aes(x = x, y = y, label = lab), inherit.aes = FALSE,
              fontface = "bold", size = 3.5) +
    facet +
    scale_x_continuous(limits = c(-0.3, 1.3)) +
    scale_y_continuous(limits = c(-1.05, 1.18)) +
    guides(colour = guide_legend(override.aes = list(linewidth = 1.2))) +
    labs(title = title) +
    theme_void(base_size = 12) +
    theme(strip.text      = element_text(size = 12, margin = margin(4, 4, 4, 4)),
          panel.border    = element_rect(colour = "grey70", fill = NA, linewidth = 0.4),
          panel.spacing   = unit(1, "lines"),
          legend.position = "right")
}


## Data
# Raw counts, soil cores only: diazotroph genera and AMF genera (or ASVs).
# The LSU outlier cores are dropped, as in the beta diversity and DESeq2 analyses.
net_nif_ps <- glom_rank(soil_cores(nif), "Genus")
net_lsu_ps <- soil_cores(lsu, drop = lsu_outliers)
if (NET_LSU_RANK == "genus") net_lsu_ps <- glom_rank(net_lsu_ps, "genus")

net_ids  <- base::intersect(sample_names(net_nif_ps), sample_names(net_lsu_ps))
nif_mat  <- t(otu_counts(net_nif_ps))[net_ids, , drop = FALSE]   # samples x taxa
lsu_mat  <- t(otu_counts(net_lsu_ps))[net_ids, , drop = FALSE]
depth_ok <- rowSums(lsu_mat) >= lsu_depth & rowSums(nif_mat) >= nif_depth
cat("Cores dropped for low depth:", sum(!depth_ok), "of", length(net_ids), "\n")
net_ids  <- net_ids[depth_ok]

# Lumped "Unclassified" genera are removed; unclassified ASVs are kept (distinct variants)
net_nif_counts <- nif_mat[net_ids, classified(net_nif_ps, "Genus"), drop = FALSE]
net_lsu_counts <- if (NET_LSU_RANK == "genus") {
  lsu_mat[net_ids, classified(net_lsu_ps, "genus"), drop = FALSE]
} else {
  lsu_mat[net_ids, , drop = FALSE]
}

net_meta <- get_meta(net_nif_ps)[net_ids, ]          # shared labels, plot and pair
net_meta$N_group <- N_to_group(net_meta$N_treatment)
stopifnot(!anyNA(net_meta$N_group), identical(rownames(net_meta), net_ids))

table(net_meta$N_group, net_meta$core_type)          # cores per N group x core type
check_pairs(net_meta)                                # unpaired cores stay put in core-type permutations

net_nif_tax <- as.data.frame(as(tax_table(net_nif_ps), "matrix"))
net_lsu_tax <- as.data.frame(as(tax_table(net_lsu_ps), "matrix"))

# One taxon pool for every network and analysis
net_D <- make_data()
c(diazotroph_genera = sum(net_D$source == "nif"), AMF_taxa = sum(net_D$source == "lsu"))


## Permutation sets (shared across analyses)
# Core type: labels swapped within pairs
set.seed(SEED)
perms_core <- shuffleSet(nrow(net_meta),
                         control = how(nperm = NET_N_PERM, blocks = net_meta$pair))

# N group: N treatments shuffled among plots within blocks
net_plots <- net_meta %>% distinct(plot, treatment_block, N_treatment)
set.seed(SEED)
perms_N   <- shuffleSet(nrow(net_plots),
                        control = how(nperm = NET_N_PERM, blocks = net_plots$treatment_block))
idx_plot  <- match(net_meta$plot, net_plots$plot)
permute_N <- function(p) N_to_group(net_plots$N_treatment[p])[idx_plot]


## A. Core type (all N)
set.seed(SEED)
fits_A  <- build_networks(net_D, net_core(net_meta$core_type), "plot")
lam_A   <- purrr::map(fits_A, "lambda")
met_A   <- metrics_of(fits_A, net_D, "core_type")
edges_A <- all_edge_tables(fits_A, net_D, "core_type")
met_A
network_size(met_A)   # mean neighbours should stay well below ~ n / log(taxa)

if (NET_RUN_TESTS) {
  null_A <- run_null(perms_core, function(p) net_core(net_meta$core_type[p]),
                     net_D, "plot", lam_A, "core_type", core_diff)
  test_A <- summarise_test(core_diff(met_A), null_A, c("metric", "contrast"))
  test_A
}


## B. N group x core type
set.seed(SEED)
fits_B  <- build_networks(net_D, net_Ncore(net_meta$N_group, net_meta$core_type), "plot")
lam_B   <- purrr::map(fits_B, "lambda")
met_B   <- metrics_of(fits_B, net_D, c("N_group", "core_type"))
edges_B <- all_edge_tables(fits_B, net_D, c("N_group", "core_type"))
met_B
network_size(met_B)   # mean neighbours should stay well below ~ n / log(taxa)

if (NET_RUN_TESTS) {
  # Core type within each N group
  stat_B_core <- function(m) core_diff(m, by = "N_group")
  null_B_core <- run_null(perms_core,
                          function(p) net_Ncore(net_meta$N_group, net_meta$core_type[p]),
                          net_D, "plot", lam_B, c("N_group", "core_type"), stat_B_core)
  test_B_core <- summarise_test(stat_B_core(met_B), null_B_core,
                                c("N_group", "metric", "contrast"))
  
  # N effect within each core type, and N x core type
  stat_B_N <- function(m) bind_rows(n_contrast(m, by = "core_type"), nx_core(m))
  null_B_N <- run_null(perms_N,
                       function(p) net_Ncore(permute_N(p), net_meta$core_type),
                       net_D, "plot", lam_B, c("N_group", "core_type"), stat_B_N)
  test_B_N <- summarise_test(stat_B_N(met_B), null_B_N, c("core_type", "metric", "contrast"))
  
  test_B_core
  test_B_N
}


## C. N group (core types pooled)
set.seed(SEED)
fits_C  <- build_networks(net_D, net_N(net_meta$N_group), c("plot", "core_type"))
lam_C   <- purrr::map(fits_C, "lambda")
met_C   <- metrics_of(fits_C, net_D, "N_group")
edges_C <- all_edge_tables(fits_C, net_D, "N_group")
met_C
network_size(met_C)   # mean neighbours should stay well below ~ n / log(taxa)

if (NET_RUN_TESTS) {
  null_C <- run_null(perms_N, function(p) net_N(permute_N(p)),
                     net_D, c("plot", "core_type"), lam_C, "N_group", n_contrast)
  test_C <- summarise_test(n_contrast(met_C), null_C, c("metric", "contrast"))
  test_C
}


## Core-count check: N contrasts in present taxa
# Within each core type (B) and pooled (C), every N group is cut to the
# smallest group's core count. Contrasts with same_sign near 1 and an interval
# excluding 0 hold at equal sampling effort.
if (NET_RUN_SUBSAMPLE) {
  sub_B_N <- subsample_counts(net_D, net_meta$N_group, by = net_meta$core_type)
  sub_C_N <- subsample_counts(net_D, net_meta$N_group)
  
  net_subsample <- bind_rows(
    `B: N within core type` = summarise_counts(sub_B_N, net_D, net_meta$N_group,
                                               by = net_meta$core_type),
    `C: N group, pooled`    = summarise_counts(sub_C_N, net_D, net_meta$N_group),
    .id = "analysis")
  
  net_subsample_groups <- bind_rows(
    `B: N within core type` = count_table(sub_B_N),
    `C: N group, pooled`    = count_table(sub_C_N),
    .id = "analysis")
  
  net_subsample
  net_subsample_groups
}


## Combined tables
dir.create("Tables", showWarnings = FALSE)

net_metrics <- bind_rows(`A: core type`       = met_A,
                         `B: N group x core`  = met_B,
                         `C: N group, pooled` = met_C, .id = "analysis")
net_edges   <- bind_rows(`A: core type`       = edges_A,
                         `B: N group x core`  = edges_B,
                         `C: N group, pooled` = edges_C, .id = "analysis")
write.csv(net_metrics, "Tables/network_metrics.csv", row.names = FALSE)
write.csv(net_edges,   "Tables/network_edges.csv",   row.names = FALSE)

if (NET_RUN_TESTS) {
  net_tests <- bind_rows(`A: core type`                = test_A,
                         `B: core type within N group` = test_B_core,
                         `B: N group and N x core`     = test_B_N,
                         `C: N group, pooled`          = test_C, .id = "analysis")
  write.csv(net_tests, "Tables/network_tests.csv", row.names = FALSE)
}
if (NET_RUN_SUBSAMPLE) {
  write.csv(net_subsample,        "Tables/network_subsample_contrasts.csv", row.names = FALSE)
  write.csv(net_subsample_groups, "Tables/network_subsample_groups.csv",    row.names = FALSE)
}

# Stable AMF-diazotroph edges, all analyses
net_edges %>%
  filter(type == "AMF-diazotroph", stable %in% TRUE) %>%
  arrange(analysis, network, dplyr::desc(freq)) %>%
  dplyr::select(analysis, network, from_taxon, to_taxon, sign, weight, freq)


## Bipartite AMF-diazotroph network figures
amf_pal  <- make_pal(net_D$group[net_D$source == "lsu"], ggsci::pal_npg()(10))
diaz_pal <- make_pal(net_D$group[net_D$source == "nif"], ggsci::pal_observable()(10))

net_fig_A <- plot_cross_network(edges_A, net_D, facet_wrap(~ core_type, nrow = 1),
                                "A. Core type (all N levels)")
net_fig_B <- plot_cross_network(edges_B, net_D, facet_grid(core_type ~ N_group),
                                "B. N group x core type")
net_fig_C <- plot_cross_network(edges_C, net_D, facet_wrap(~ N_group, nrow = 1),
                                "C. N group (core types pooled)")
net_fig_A; net_fig_B; net_fig_C

# Saved as Figure 5A (analysis B) and Figure S5 (analyses A and C) at the end of the script


## ============================== Network-linked diazotrophs and FLNF =======================================
## Do diazotroph genera in stable, positive AMF-diazotroph associations predict
## free-living N fixation? One genus set per N group x core type network.
##   1. For each network in analysis B, the diazotroph genera with >= 1 positive,
##      stable edge (selected in >= NET_EDGE_FREQ of StARS subsamples)
##   2. Each set's summed relative abundance in every core: all nifH OTUs in
##      those genera, as a share of the core's rarefied nifH reads (genus level)
##   3. One model per set: FLNF ~ abundance + core_type * N_treatment + block,
##      random plot and pair, log(FLNF). Kenward-Roger F tests (lmer), with a
##      Gamma(log) GLMM (glmmTMB) as a check; BH correction across the sets.
##      N and core type are in every model, so each test asks whether the set
##      predicts FLNF beyond the treatment effects.
##   4. Within-network test: each set only in its own network's cores, overall
##      and within plots (abundance centred on the plot mean), matching how the
##      networks were inferred.
## Sets share genera, so their abundances are correlated: they are tested in
## separate models, not together. Genera are chosen by co-occurrence with AMF,
## not by FLNF, so the selection does not use the response.

LINK_EDGE_TYPES <- "AMF-diazotroph"   # add "Diazotroph-diazotroph" to include those edges

## Network-linked genera, per network
link_genera <- edges_B %>%
  filter(type %in% LINK_EDGE_TYPES, sign == "Positive", stable %in% TRUE) %>%
  dplyr::select(network, N_group, core_type, from, to) %>%
  pivot_longer(c(from, to), values_to = "node") %>%
  filter(startsWith(node, "nif|")) %>%
  mutate(Genus = net_nif_tax[sub("^nif\\|", "", node), "Genus"],
         set   = paste(N_group, sub("Hyphae ", "", core_type), sep = "_")) %>%   # e.g. Low_Suppressed
  distinct(set, network, N_group, core_type, Genus)

link_sets <- split(link_genera$Genus, link_genera$set)
lengths(link_sets)                                   # genera per network

# Which genera are shared between networks
link_membership <- link_genera %>%
  mutate(present = 1L) %>%
  dplyr::select(Genus, set, present) %>%
  pivot_wider(names_from = set, values_from = present, values_fill = 0L)
link_membership

## Relative abundance of each set in every core (genus level, rarefied nifH)
nif_rfy_tax <- as.data.frame(as(tax_table(nif_rfy), "matrix"))
link_otus   <- lapply(link_sets, function(g) rownames(nif_rfy_tax)[nif_rfy_tax$Genus %in% g])

link_df   <- otu_group_data(nif_rfy, link_otus)      # core level, FLNF > 0
link_vars <- attr(link_df, "groups")                 # sets with variation across cores
setdiff(names(link_sets), link_vars)                 # sets dropped (no variation)

# Correlation among sets and with N rate
round(cor(link_df[link_vars], use = "pairwise.complete.obs", method = "spearman"), 2)
sapply(link_df[link_vars], function(x)
  cor(x, as.numeric(as.character(link_df$N_treatment)), method = "spearman", use = "complete.obs"))

## One model per set
link_fits <- set_names(link_vars) %>%
  purrr::map(~ fit_models(link_df, "n_fixation", covariates = .x, log_response = TRUE))

link_tests <- purrr::imap_dfr(link_fits, function(fit, v) {
  kr    <- as.data.frame(fit$anova_t2)[v, ]
  glmm  <- as.data.frame(car::Anova(fit$glmm, type = 2))[v, ]
  slope <- emtrends(fit$lmer, ~ 1, var = v, lmer.df = "kenward-roger") %>%
    summary(infer = TRUE) %>%
    as.data.frame()
  tibble(set              = v,
         n_genera         = length(link_sets[[v]]),
         F                = kr$`F value`,
         NumDF            = kr$NumDF,
         DenDF            = kr$DenDF,
         p_KR             = kr$`Pr(>F)`,
         p_GLMM           = glmm$`Pr(>Chisq)`,
         fold_per_10pp    = exp(0.1 * slope[[paste0(v, ".trend")]]),   # FLNF fold change per +10 pp
         fold_lo          = exp(0.1 * slope$lower.CL),
         fold_hi          = exp(0.1 * slope$upper.CL))
}) %>%
  mutate(p_BH = p.adjust(p_KR, method = "BH"), .after = p_KR)
link_tests

# Residual checks for each set's model
purrr::walk(link_fits, check_fit)

# Sensitivity: does a set's relationship differ between core types?
# fit_models(link_df, "n_fixation", covariates = c("Low_Enriched", "Low_Enriched:core_type"),
#            log_response = TRUE)$anova_t2

## Figure: cores + model-adjusted trend for each set (averaged over core type, N, block)
link_pred <- purrr::imap_dfr(link_fits, function(fit, v) {
  grid <- seq(min(fit$data[[v]]), max(fit$data[[v]]), length.out = 50)
  emmeans(fit$lmer, as.formula(paste("~", v)), at = setNames(list(grid), v),
          type = "response", lmer.df = "kenward-roger") %>%
    as.data.frame() %>%
    dplyr::rename(estimate = any_of(c("response", "emmean")), abund = all_of(v)) %>%
    mutate(set = v)
})

link_points <- link_df %>%
  filter(!is.na(n_fixation)) %>%
  dplyr::select(n_fixation, N_treatment, core_type, all_of(link_vars)) %>%
  pivot_longer(all_of(link_vars), names_to = "set", values_to = "abund")

# Facet order (Low -> High N; suppressed row first) and readable labels,
# e.g. "Low_Suppressed" -> "Low-N: Hyphae Suppressed"
set_levels <- as.vector(t(outer(c("Suppressed", "Enriched"), NET_N_GROUPS,
                                function(core, N) paste(N, core, sep = "_"))))
set_names_pretty <- setNames(sub("^(.*)_(.*)$", "\\1-N: Hyphae \\2", set_levels), set_levels)
set_labels <- function(df) {
  lev <- base::intersect(set_levels, unique(df$set))
  df %>% mutate(set = factor(set, levels = lev, labels = set_names_pretty[lev]))
}

link_plot <- ggplot() +
  geom_ribbon(data = set_labels(link_pred), aes(x = abund, ymin = lower.CL, ymax = upper.CL),
              fill = "grey80", alpha = 0.6) +
  geom_line(data = set_labels(link_pred), aes(x = abund, y = estimate), linewidth = 1.1) +
  geom_point(data = set_labels(link_points),
             aes(x = abund, y = n_fixation, fill = N_treatment, shape = core_type),
             size = 2.4, alpha = 0.85, colour = "black") +
  facet_wrap(~ set, scales = "free_x", nrow = 2) +
  scale_y_log10() +
  scale_fill_viridis_d(name = "N Treatment", labels = N_treatment_labels, drop = FALSE) +
  scale_shape_manual(name = "Core Type", values = c(21, 24)) +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  labs(x = "Relative abundance of the network's linked diazotroph genera",
       y = expression(FLNF ~ ("µg" ~ N ~ g^-1 ~ dry ~ soil ~ day^-1))) +
  theme_bw() +
  theme(axis.text    = element_text(size = 12, colour = "black"),
        axis.title   = element_text(size = 16),
        strip.text   = element_text(size = 14),
        legend.text  = element_text(size = 14),
        legend.title = element_text(size = 16))
link_plot


## Within-network test: each set only in the cores of its own network
# Matches how the networks were inferred: the same cores, one core type and two
# N levels. Model: log(FLNF) ~ abundance + N_treatment + block + (1 | plot).
# The within-plot version uses abundance centred on its plot mean (plot mean as
# a separate term), mirroring the networks' removal of plot means: do cores with
# more of these genera than their plot-mates fix more N? ~14-20 cores and 8
# plots per network, so power is limited.
fit_within_network <- function(df, v) {
  dat <- df %>%
    filter(!is.na(n_fixation), !is.na(.data[[v]])) %>%
    group_by(plot) %>%
    mutate(abund_plot   = mean(.data[[v]]),
           abund_within = .data[[v]] - abund_plot) %>%
    ungroup() %>%
    droplevels()
  fixed <- "N_treatment + treatment_block + (1 | plot)"
  list(data   = dat,
       lmer   = lmerTest::lmer(as.formula(paste("log(n_fixation) ~", v, "+", fixed)), data = dat),
       within = lmerTest::lmer(as.formula(paste("log(n_fixation) ~ abund_within + abund_plot +", fixed)),
                               data = dat))
}

net_cores <- function(set) {          # cores of the network a set came from
  key <- link_genera %>% filter(set == !!set) %>% distinct(N_group, core_type)
  rownames(net_meta)[as.character(net_meta$N_group)   == as.character(key$N_group) &
                       as.character(net_meta$core_type) == as.character(key$core_type)]
}

link_within_fits <- set_names(link_vars) %>%
  purrr::map(~ fit_within_network(filter(link_df, Site %in% net_cores(.x)), .x))

link_within_tests <- purrr::imap_dfr(link_within_fits, function(fit, v) {
  kr  <- as.data.frame(anova(fit$lmer,   type = 2, ddf = "Kenward-Roger"))[v, ]
  krw <- as.data.frame(anova(fit$within, type = 2, ddf = "Kenward-Roger"))["abund_within", ]
  b   <- summary(fit$lmer)$coefficients[v, "Estimate"]
  tibble(set            = v,
         n_cores        = nrow(fit$data),
         n_plots        = dplyr::n_distinct(fit$data$plot),
         abund_sd       = sd(fit$data[[v]]),
         F              = kr$`F value`,
         DenDF          = kr$DenDF,
         p_KR           = kr$`Pr(>F)`,
         fold_per_sd    = exp(b * sd(fit$data[[v]])),   # FLNF fold change per 1 SD of abundance
         F_within       = krw$`F value`,
         p_within       = krw$`Pr(>F)`)
}) %>%
  mutate(p_BH        = p.adjust(p_KR,     method = "BH"), .after = p_KR) %>%
  mutate(p_within_BH = p.adjust(p_within, method = "BH"))
link_within_tests

# Figure: each network's own cores with the fitted trend (N and block averaged)
link_within_pred <- purrr::imap_dfr(link_within_fits, function(fit, v) {
  grid <- seq(min(fit$data[[v]]), max(fit$data[[v]]), length.out = 50)
  emmeans(fit$lmer, as.formula(paste("~", v)), at = setNames(list(grid), v),
          type = "response", lmer.df = "kenward-roger") %>%
    as.data.frame() %>%
    dplyr::rename(estimate = any_of(c("response", "emmean")), abund = all_of(v)) %>%
    mutate(set = v)
})
link_within_points <- purrr::imap_dfr(link_within_fits, function(fit, v) {
  fit$data %>% transmute(set = v, abund = .data[[v]], n_fixation,
                         N_treatment = as.character(N_treatment))
}) %>%
  mutate(N_treatment = factor(N_treatment, levels = names(N_labels)))   # keep 0 -> 196 order and colours

link_within_plot <- ggplot() +
  geom_ribbon(data = set_labels(link_within_pred), aes(x = abund, ymin = lower.CL, ymax = upper.CL),
              fill = "grey80", alpha = 0.6) +
  geom_line(data = set_labels(link_within_pred), aes(x = abund, y = estimate), linewidth = 1.1) +
  geom_point(data = set_labels(link_within_points),
             aes(x = abund, y = n_fixation, fill = N_treatment),
             shape = 21, size = 2.8, alpha = 0.9, colour = "black") +
  facet_wrap(~ set, scales = "free_x", nrow = 2) +
  scale_y_log10() +
  scale_fill_viridis_d(name = "N Treatment", labels = N_treatment_labels, drop = FALSE) +
  labs(x = "Relative abundance of the network's linked diazotroph genera (network cores only)",
       y = expression(FLNF ~ ("µg" ~ N ~ g^-1 ~ dry ~ soil ~ day^-1))) +
  theme_bw() +
  theme(axis.text    = element_text(size = 12, colour = "black"),
        axis.title   = element_text(size = 16),
        strip.text   = element_text(size = 14),
        legend.text  = element_text(size = 14),
        legend.title = element_text(size = 16))
link_within_plot

write.csv(link_genera, "Tables/network_linked_genera.csv",       row.names = FALSE)
write.csv(link_tests,  "Tables/network_linked_genera_FLNF.csv",  row.names = FALSE)
write.csv(link_within_tests, "Tables/network_linked_genera_FLNF_within_network.csv", row.names = FALSE)
# Saved as Figure S6 (A: all cores, B: network cores) at the end of the script


## ============================== Manuscript figures and tables =======================================
## Assembled here because they combine results from several sections.
## Saved earlier in the script:
##   Figure 1   Figures/Figure_1.png    alpha diversity, FLNF, ERH
##   Figure 2   Figures/Figure_2.png    beta diversity (NMDS + plot-level vectors)
##   Figure 4   Figures/Figure_4.png    DESeq2 N responders, nifH and LSU
##   Figure S1  Figures/Figure_S1.png   rarefaction curves
##   Figure S2  Figures/Figure_S2.png   soil chemistry model estimates
##   Figure S3  Figures/Figure_S3.png   AMF core-type DESeq2
##   Figure S4  Figures/Figure_S4.png   DESeq2 N-response outliers vs FLNF
## Built here: Figures 3 and 5, Figures S5 and S6, Tables 1 and 2,
## Tables S2 (soil chemistry models) and S3 (envfit vectors).
## Table S1 (primers / PCR conditions) is in the manuscript, not generated here.
## Figure 3B and Table 1 need the permutation tests (NET_RUN_TESTS = TRUE);
## Figure 5B and the Table 1 subsample columns need NET_RUN_SUBSAMPLE = TRUE.

ms_theme <- theme_bw() +
  theme(axis.text    = element_text(size = 13, colour = "black"),
        axis.title   = element_text(size = 15),
        strip.text   = element_text(size = 13),
        legend.text  = element_text(size = 13),
        legend.title = element_text(size = 15))

## Figure 3: diazotroph core-type DESeq2 (A), core-type networks (B) and the
## cross-domain share vs its permutation null (C)
if (NET_RUN_TESTS) {
  cross_null <- null_A %>% filter(metric == "cross_frac")
  cross_obs  <- test_A %>% filter(metric == "cross_frac")
  obs_frac   <- met_A %>% dplyr::select(core_type, cross_frac)
  
  fig3_null <- ggplot(cross_null, aes(x = value)) +
    geom_histogram(bins = 40, fill = "grey75", colour = "grey40", linewidth = 0.2) +
    geom_vline(xintercept = cross_obs$observed, colour = "firebrick3", linewidth = 1.2) +
    annotate("label", x = cross_obs$observed, y = Inf, vjust = 1.3, hjust = 1.05, size = 4.2,
             label = sprintf("Observed = %.3f\np = %.3f", cross_obs$observed, cross_obs$p_value)) +
    labs(x = "Difference in AMF-diazotroph share of edges\n(Hyphae Enriched - Hyphae Suppressed)",
         y = "Permutations",
         caption = sprintf("Observed share: Suppressed %.2f, Enriched %.2f",
                           obs_frac$cross_frac[obs_frac$core_type == "Hyphae Suppressed"],
                           obs_frac$cross_frac[obs_frac$core_type == "Hyphae Enriched"])) +
    ms_theme
  
  # Core-type networks (analysis A), no legends: strips give the number and share
  # of AMF-diazotroph edges. Edge width = stability; grey = positive, red = negative.
  fig3_strips <- met_A %>%
    mutate(lab = sprintf("%s\n%d AMF-diazotroph edges (%.0f%%)",
                         core_type, cross_edges, 100 * cross_frac))
  fig3_net <- plot_cross_network(
    edges_A, net_D,
    facet_wrap(~ core_type, nrow = 1,
               labeller = as_labeller(setNames(fig3_strips$lab,
                                               as.character(fig3_strips$core_type))))
  ) + theme(legend.position = "none",
            plot.margin = margin(t = 5, r = 5, b = 15, l = 5))
  
  fig3_null <- fig3_null + theme(plot.margin = margin(t = 25, r = 10, b = 5, l = 20))
  
  fig3 <- plot_grid(
    deseq_core_nif,
    plot_grid(fig3_net, fig3_null, ncol = 1, rel_heights = c(1.3, 1),
              labels = c("B", "C"), label_fontface = "bold", label_size = 16),
    ncol = 2, rel_widths = c(1.4, 1),
    labels = c("A", ""), label_fontface = "bold", label_size = 16
  )
  fig3
  ggsave("Figures/Figure_3.png", plot = fig3, dpi = 600, width = 20, height = 12, bg = "white")
}

## Figure 5: N effects on networks
# A: analysis B networks; B: present taxa by domain at equal core counts;
# C: AMF-diazotroph edge counts and connectance among present taxa across N
fig5_A <- net_fig_B + labs(title = NULL)

if (NET_RUN_SUBSAMPLE) {
  fig5_counts <- net_subsample_groups %>%
    filter(analysis == "B: N within core type") %>%
    transmute(core_type = factor(stratum, levels = net_core_types),
              N_group   = factor(N_group, levels = NET_N_GROUPS),
              n_cores,
              AMF_mean = amf_present_mean, AMF_lo = amf_present_lo, AMF_hi = amf_present_hi,
              Diaz_mean = nif_present_mean, Diaz_lo = nif_present_lo, Diaz_hi = nif_present_hi) %>%
    pivot_longer(-c(core_type, N_group, n_cores),
                 names_to = c("domain", ".value"), names_sep = "_") %>%
    mutate(domain = factor(ifelse(domain == "AMF", "AMF ASVs", "Diazotroph genera"),
                           levels = c("AMF ASVs", "Diazotroph genera")))
  
  fig5_B <- ggplot(fig5_counts, aes(x = N_group, y = mean, colour = core_type, group = core_type)) +
    geom_line(linewidth = 1.1, position = position_dodge(width = 0.3)) +
    geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.15, linewidth = 0.8,
                  position = position_dodge(width = 0.3)) +
    geom_point(aes(fill = core_type), shape = 21, colour = "black", size = 4,
               position = position_dodge(width = 0.3)) +
    facet_wrap(~ domain, scales = "free_y") +
    scale_colour_npg(name = "Core Type") +
    scale_fill_npg(name = "Core Type") +
    labs(x = "N group", y = "Taxa present\n(equal cores per N group; mean, 95% interval)") +
    ms_theme
} else {
  fig5_B <- NULL
}

# C: AMF-diazotroph associations across N within core types, two facets:
#    top    - number of AMF-diazotroph edges (turnover + rewiring; tested as the
#             pool-based cross_density, which is the count / a fixed constant)
#    bottom - AMF-diazotroph connectance among the taxa present
#             (edges / (AMF present x diazotrophs present); rewiring only)
# Labels give the permutation p-value range for the N contrasts in each core type.
fig5_C_measures <- c(cross_edges           = "AMF-diazotroph edges",
                     cross_density_present = "Connectance among taxa present")
fig5_C_tested   <- c(cross_edges = "cross_density", cross_density_present = "cross_density_present")

fig5_C_df <- met_B %>%
  dplyr::select(N_group, core_type, all_of(names(fig5_C_measures))) %>%
  pivot_longer(all_of(names(fig5_C_measures)), names_to = "measure", values_to = "value") %>%
  mutate(measure = factor(fig5_C_measures[measure], levels = fig5_C_measures))

fig5_C_lab <- if (NET_RUN_TESTS) {
  net_tests %>%
    filter(analysis == "B: N group and N x core", !is.na(core_type),
           metric %in% fig5_C_tested) %>%
    mutate(measure = fig5_C_measures[names(fig5_C_tested)[match(metric, fig5_C_tested)]],
           core    = sub("Hyphae ", "", core_type)) %>%
    group_by(measure, core) %>%
    summarise(p_lo = min(p_value), p_hi = max(p_value), .groups = "drop") %>%
    mutate(txt = ifelse(sprintf("%.3f", p_lo) == sprintf("%.3f", p_hi),
                        sprintf("%s: p = %.3f", core, p_lo),
                        sprintf("%s: p = %.3f-%.3f", core, p_lo, p_hi))) %>%
    group_by(measure) %>%
    summarise(label = paste(c("N effect", rev(sort(txt))), collapse = "\n"), .groups = "drop") %>%
    mutate(measure = factor(measure, levels = fig5_C_measures))
} else NULL

fig5_C <- ggplot(fig5_C_df, aes(x = N_group, y = value, colour = core_type, group = core_type)) +
  geom_line(linewidth = 1.1) +
  geom_point(aes(fill = core_type), shape = 21, colour = "black", size = 4) +
  facet_wrap(~ measure, ncol = 1, scales = "free_y") +
  scale_colour_npg(name = "Core Type") +
  scale_fill_npg(name = "Core Type") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.35))) +   # room for the labels
  expand_limits(y = 0) +
  labs(x = "N group", y = NULL) +
  ms_theme
if (!is.null(fig5_C_lab)) {
  fig5_C <- fig5_C + geom_text(data = fig5_C_lab, aes(x = 3.4, y = Inf, label = label),
                               inherit.aes = FALSE, hjust = 1, vjust = 1.2, size = 3.6)
}

fig5_BC <- if (!is.null(fig5_B)) {
  plot_grid(fig5_B + theme(legend.position = "none"),
            fig5_C + theme(legend.position = "none"),
            cowplot::get_legend(fig5_C),
            ncol = 3, rel_widths = c(1.6, 1, 0.5),
            labels = c("B", "C", ""), label_fontface = "bold", label_size = 16)
} else {
  plot_grid(fig5_C, labels = "C", label_fontface = "bold", label_size = 16)
}
fig5 <- plot_grid(fig5_A, fig5_BC, ncol = 1, rel_heights = c(1.5, 1),
                  labels = c("A", ""), label_fontface = "bold", label_size = 16)
fig5
ggsave("Figures/Figure_5.png", plot = fig5, dpi = 600, width = 16, height = 16, bg = "white")

## Figure S5: networks from analysis A (core type) and C (N group, pooled)
# One shared legend listing every AMF and diazotroph group in either panel
figS5_A <- plot_cross_network(edges_A, net_D, facet_wrap(~ core_type, nrow = 1), legend_all = TRUE)
figS5_C <- plot_cross_network(edges_C, net_D, facet_wrap(~ N_group,   nrow = 1), legend_all = TRUE)

figS5 <- plot_grid(
  plot_grid(figS5_A + theme(legend.position = "none"),
            figS5_C + theme(legend.position = "none"),
            ncol = 1, labels = c("A", "B"), label_fontface = "bold", label_size = 16),
  cowplot::get_legend(figS5_A),
  ncol = 2, rel_widths = c(1, 0.22)
)
figS5
ggsave("Figures/Figure_S5.png", plot = figS5, dpi = 600, width = 14, height = 16, bg = "white")

## Figure S6: network-linked diazotroph genera vs FLNF (A: all cores, B: network cores)
figS6 <- plot_grid(link_plot, link_within_plot, ncol = 1,
                   labels = c("A", "B"), label_fontface = "bold", label_size = 16)
figS6
ggsave("Figures/Figure_S6.png", plot = figS6, dpi = 600, width = 16, height = 16)

## Table 1: network permutation tests, primary metrics
# BH adjustment within each analysis block. Count metrics carry the
# equal-core subsample interval for N contrasts (core-type contrasts need none:
# within-pair permutations keep each network's core count).
table1_metrics <- c("n_present", "amf_present", "nif_present",
                    "density_present", "cross_frac", "cross_density_present")

if (NET_RUN_TESTS) {
  sub_cols <- if (NET_RUN_SUBSAMPLE) {
    net_subsample %>%
      transmute(sub_key = ifelse(stratum == "Pooled", "C Pooled", paste("B", stratum)),
                metric, contrast, sub_cores = n_cores,
                sub_mean, sub_lo, sub_hi, sub_same_sign = same_sign)
  } else {
    tibble(sub_key = character(), metric = character(), contrast = character())
  }
  
  table1 <- net_tests %>%
    filter(metric %in% table1_metrics) %>%
    mutate(sub_key = case_when(analysis == "B: N group and N x core" & !is.na(core_type) ~
                                 paste("B", core_type),
                               analysis == "C: N group, pooled" ~ "C Pooled")) %>%
    left_join(sub_cols, by = c("sub_key", "metric", "contrast")) %>%
    group_by(analysis) %>%
    mutate(p_BH = p.adjust(p_value, method = "BH")) %>%
    ungroup() %>%
    mutate(metric = factor(metric, levels = table1_metrics)) %>%
    arrange(analysis, N_group, core_type, contrast, metric) %>%
    dplyr::select(analysis, N_group, core_type, contrast, metric, observed,
                  p_value, p_BH, n_perm, any_of(c("sub_cores", "sub_mean", "sub_lo",
                                                  "sub_hi", "sub_same_sign"))) %>%
    mutate(across(where(is.numeric) & !any_of(c("n_perm", "sub_cores")), ~ signif(.x, 3)))
  table1
  write.csv(table1, "Tables/Table_1_network_permutation_tests.csv", row.names = FALSE)
}

## Table 2: strongest AMF-diazotroph links
# Stable edges in the six N group x core type networks (analysis B), which use
# separate sets of cores, so the network count is a count of independent
# detections. Detection in the core-type networks (analysis A) is reported
# separately: A reuses B's cores, so it is not an extra replicate.
amf_diaz <- function(edges) {
  edges %>%
    filter(type == "AMF-diazotroph", stable %in% TRUE) %>%
    mutate(amf  = ifelse(startsWith(from, "lsu|"), from, to),
           diaz = ifelse(startsWith(from, "lsu|"), to, from))
}

links_A <- amf_diaz(edges_A) %>%
  group_by(amf, diaz) %>%
  summarise(A_core_types = paste(sort(as.character(core_type)), collapse = "; "),
            A_max_freq   = max(freq), .groups = "drop")

table2 <- amf_diaz(edges_B) %>%
  group_by(amf, diaz) %>%
  summarise(n_networks_B = dplyr::n_distinct(network),
            networks_B   = paste(network, collapse = "; "),
            sign         = if (dplyr::n_distinct(sign) == 1) dplyr::first(sign) else "Mixed",
            mean_freq    = mean(freq),
            max_freq     = max(freq),
            mean_weight  = mean(weight),
            .groups      = "drop") %>%
  left_join(links_A, by = c("amf", "diaz")) %>%
  mutate(AMF_taxon        = net_D$label[amf],
         AMF_group        = net_D$group[amf],
         Diazotroph_genus = net_D$label[diaz],
         Diazotroph_class = net_D$group[diaz],
         A_core_types     = tidyr::replace_na(A_core_types, "Not detected")) %>%
  arrange(dplyr::desc(n_networks_B), dplyr::desc(mean_freq)) %>%
  dplyr::select(AMF_taxon, AMF_group, Diazotroph_genus, Diazotroph_class, sign,
                n_networks_B, networks_B, mean_freq, max_freq, mean_weight,
                A_core_types, A_max_freq) %>%
  mutate(across(c(mean_freq, max_freq, mean_weight, A_max_freq), ~ round(.x, 2)))

table2 %>% slice_head(n = 20)                   # main-text table: top links
write.csv(table2, "Tables/Table_2_AMF_diazotroph_links.csv", row.names = FALSE)


## Table S2: soil chemistry models (Kenward-Roger Type II F tests)
# BH adjustment across chemistry variables within each term (from anova_table())
tableS2 <- chem_stats %>%
  mutate(variable = set_names(chem_spec$label, chem_spec$var)[variable],
         term     = recode(term, N_treatment = "N", core_type = "Core type",
                           `core_type:N_treatment` = "N x Core type"),
         term     = factor(term, levels = c("N", "Core type", "N x Core type"))) %>%
  arrange(factor(variable, levels = chem_spec$label), term) %>%
  dplyr::select(variable, term, NumDF, DenDF, F = `F value`, p = `Pr(>F)`, p_BH) %>%
  mutate(across(c(DenDF, F), ~ round(.x, 2)), across(c(p, p_BH), ~ signif(.x, 3)))
tableS2
write.csv(tableS2, "Tables/Table_S2_chemistry_models.csv", row.names = FALSE)

## Table S3: envfit vectors on the NMDS ordinations (plot means, n = 24 plots,
# permutations within blocks)
tableS3 <- bind_rows(`Diazotrophs (nifH)` = nif_vec$table,
                     `AMF (LSU)`          = lsu_vec$table, .id = "community") %>%
  mutate(variable    = dplyr::coalesce(vector_labels[name], name),
         significant = p < 0.05) %>%
  arrange(community, p) %>%
  dplyr::select(community, variable, NMDS1, NMDS2, r2, p, significant) %>%
  mutate(across(c(NMDS1, NMDS2, r2), ~ round(.x, 3)), p = signif(p, 3))
tableS3
write.csv(tableS3, "Tables/Table_S3_envfit_vectors.csv", row.names = FALSE)
