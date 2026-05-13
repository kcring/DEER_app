# Deer Density Lab – USCR vs REM vs TTE
# Models run via run_USCR(), run_REM(), run_TTE() in R/sim_and_models.R

# -------------------------------------------------------------------
# Packages
# -------------------------------------------------------------------

needs <- c(
  "shiny", "bslib", "DT", "ggplot2", "dplyr", "tidyr",
  "readr", "purrr", "stringr", "secr", "data.table",
  "leaflet", "ggrepel", "shinyjs",
  "nimble", "parallel", "MCMCvis", "lubridate", "sf", "tibble",
  "future", "promises"
)

# Redwood-inspired palette (hex only — no extra color package)
redwood_colors <- c(
  "#303018", "#604830", "#609048", "#90A860", "#786048", "#B8C4A8"
)

missing <- needs[!needs %in% rownames(installed.packages())]
if (length(missing)) {
  message("DEER App: installing missing packages: ", paste(missing, collapse = ", "))
  utils::install.packages(missing, dependencies = TRUE)
}
failed <- character()
for (pkg in needs) {
  if (!requireNamespace(pkg, quietly = TRUE)) failed <- c(failed, pkg)
}
if (length(failed)) {
  stop(
    "DEER App could not load: ", paste(failed, collapse = ", "),
    "\nInstall manually, e.g. install.packages(c(",
    paste0("\"", failed, "\"", collapse = ", "), "))",
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(shiny)
  library(bslib)
  library(DT)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(purrr)
  library(stringr)
  library(secr)
  library(data.table)
  library(leaflet)
  library(ggrepel)
  library(tibble)
  library(sf)
  library(nimble)
  library(shinyjs)
  library(future)
  library(promises)
})

# Conservative async worker pool for server-side background REM jobs.
# Keep this small because the model code can also use parallel chains internally.
available_cores <- suppressWarnings(parallel::detectCores(logical = TRUE))
if (is.na(available_cores) || available_cores < 2) {
  async_workers <- 1L
} else {
  async_workers <- min(2L, available_cores - 1L)
}
future::plan(future::multisession, workers = async_workers)

# -------------------------------------------------------------------
# Helper files
# -------------------------------------------------------------------

source("R/sim_and_models.R")   # run_USCR(), run_REM(), run_TTE()
source("R/data_checks.R")      # QC + summary helpers
source("R/camera_array_helpers.R")

# Pin the active USCR implementation to the current file contents so stale
# objects from prior sessions cannot route the app through an older wrapper.
uscr_model_env <- new.env(parent = globalenv())
sys.source("R/sim_and_models.R", envir = uscr_model_env)
run_USCR_app <- uscr_model_env$run_USCR

# -------------------------------------------------------------------
# Simulation helpers
# -------------------------------------------------------------------
# Simulation and model-input helpers now come directly from R/sim_and_models.R:
#   simulate_camera_counts(), get_counts_matrix(), capthist_to_events(),
#   sim_model_inputs(), build_nps_model_inputs()

# -------------------------------------------------------------------
# Generic helpers for summaries & WAIC-combo (like Rmd)
# -------------------------------------------------------------------

get_waic_value <- function(fit) {
  if (is.null(fit$waic)) return(NA_real_)
  if (is.list(fit$waic) && !is.null(fit$waic$WAIC)) {
    as.numeric(fit$waic$WAIC)
  } else {
    as.numeric(fit$waic)
  }
}

summarize_method <- function(fit) {
  D_mi2 <- fit$samples_all[, "D_mi2"]
  D_km2 <- D_mi2 / 2.59
  
  list(
    mean_km2  = mean(D_km2),
    sd_km2    = sd(D_km2),
    q2.5_km2  = quantile(D_km2, 0.025),
    q97.5_km2 = quantile(D_km2, 0.975),
    mean_mi2  = mean(D_mi2),
    q2.5_mi2  = quantile(D_mi2, 0.025),
    q97.5_mi2 = quantile(D_mi2, 0.975),
    waic      = get_waic_value(fit)
  )
}

build_combo_table_from_fits <- function(fits) {
  model_order <- c("REM", "TTE", "USCR")
  fits <- fits[model_order[model_order %in% names(fits)]]
  fits <- fits[!vapply(fits, is.null, logical(1))]
  if (!length(fits)) return(NULL)

  model_names <- names(fits)
  waic_tbl <- tibble::tibble(
    model = model_names,
    waic = vapply(fits, get_waic_value, numeric(1))
  )

  if (length(model_names) == 1L) {
    waic_tbl <- waic_tbl %>%
      dplyr::mutate(deltaWAIC = 0, rel_lik = 1, w = 1)
  } else if (any(!is.finite(waic_tbl$waic))) {
    waic_tbl <- waic_tbl %>%
      dplyr::mutate(
        deltaWAIC = NA_real_,
        rel_lik = NA_real_,
        w = 1 / length(model_names)
      )
  } else {
    waic_tbl <- waic_tbl %>%
      dplyr::mutate(
        deltaWAIC = waic - min(waic),
        rel_lik = exp(-0.5 * deltaWAIC),
        w = rel_lik / sum(rel_lik)
      )
  }

  density_draws <- lapply(fits, function(fit) {
    x <- as.numeric(fit$samples_all[, "D_mi2"])
    sort(x[is.finite(x)])
  })
  n_draws <- min(lengths(density_draws))
  if (!is.finite(n_draws) || n_draws < 1L) return(NULL)

  mat <- do.call(
    cbind,
    lapply(density_draws, function(x) x[seq_len(n_draws)])
  )
  colnames(mat) <- model_names

  density_est <- mat
  summary_rows <- character()
  if (ncol(mat) > 1L) {
    density_est <- cbind(
      mat,
      `unweighted mean` = rowMeans(mat),
      `weighted mean` = as.numeric(mat %*% waic_tbl$w)
    )
    summary_rows <- c("unweighted mean", "weighted mean")
  }

  means <- apply(density_est, 2, mean)
  lower <- apply(density_est, 2, stats::quantile, probs = 0.025)
  upper <- apply(density_est, 2, stats::quantile, probs = 0.975)
  prob20 <- apply(density_est, 2, function(x) mean(x > 20))

  table_out <- tibble::tibble(
    Method = c(model_names, summary_rows),
    `Mean density (animals/mi²)` = as.numeric(means),
    `Lower 2.5%` = as.numeric(lower),
    `Upper 97.5%` = as.numeric(upper),
    `Prob > 20 DPSM` = as.numeric(prob20),
    `WAIC weight` = c(waic_tbl$w, rep(NA_real_, length(summary_rows)))
  )

  list(
    table = table_out,
    waic = waic_tbl,
    models_available = model_names
  )
}

#' Simulated data: only USCR is fit — single-model summary table (no WAIC averaging)
build_sim_combo_table_uscr_only <- function(uscr_fit) {
  if (is.null(uscr_fit) || is.null(uscr_fit$samples_all)) return(NULL)
  D <- uscr_fit$samples_all[, "D_mi2"]
  waic_val <- get_waic_value(uscr_fit)
  waic_tbl <- tibble::tibble(
    model     = "USCR",
    waic      = waic_val,
    deltaWAIC = 0,
    rel_lik   = 1,
    w         = 1
  )
  table_out <- tibble::tibble(
    Method                    = "USCR",
    `Mean density (animals/mi²)` = mean(D),
    `Lower 2.5%`              = stats::quantile(D, 0.025),
    `Upper 97.5%`             = stats::quantile(D, 0.975),
    `Prob > 20 DPSM`          = mean(D > 20),
    `WAIC weight`             = 1
  )
  list(table = table_out, waic = waic_tbl)
}

#' Posterior summary for CSV export (parameter names, mean, 95% CI)
posterior_summary_df <- function(fit) {
  if (is.null(fit) || is.null(fit$samples_all)) return(NULL)
  m <- as.data.frame(fit$samples_all)
  tibble::tibble(
    parameter = names(m),
    mean      = sapply(m, mean),
    q025      = sapply(m, function(x) stats::quantile(x, 0.025)),
    q975      = sapply(m, function(x) stats::quantile(x, 0.975))
  )
}

# -------------------------------------------------------------------
# UI
# -------------------------------------------------------------------

ui <- page_fillable(
  shinyjs::useShinyjs(),
  theme = bs_theme(
    version = 5, 
    bootswatch = "flatly",
    primary = redwood_colors[3],  # #609048 - green
    secondary = redwood_colors[2], # #604830 - medium brown
    success = redwood_colors[4],   # #90A860 - light green
    bg = "#fafafa",
    fg = redwood_colors[1]         # #303018 - dark brown/green
  ),
  tags$head(
    tags$script(src = "https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-mml-chtml.js"),
    tags$script(HTML("
      window.MathJax = { 
        tex: { 
          inlineMath: [['$', '$'], ['\\\\(', '\\\\)']],
          displayMath: [['$$', '$$'], ['\\\\[', '\\\\]']],
          processEscapes: true,
          processEnvironments: true
        } 
      };
    ")),
    tags$style(HTML("
      :root {
        --rw1: #A45A52;
        --rw2: #52372E;
        --rw3: #8C9B7A;
        --rw4: #EDE8E2;
        --rw5: #3C4B3E;
        --ink: #213026;
        --muted: #5c675d;
        --radius: 14px;
        --shadow: 0 10px 26px rgba(33,48,38,.10);
      }
      h1, h2, h3, h4 {
        color: var(--ink) !important;
      }
      .hero h1, .hero h2, .hero h3, .hero h4 {
        color: #609048 !important;
      }
      h1 {
        font-size: clamp(1.9rem, 2.4vw + 1rem, 2.6rem);
        line-height: 1.15;
        margin: 0.25rem 0;
      }
      h2 {
        font-size: clamp(1.2rem, 1.2vw + .9rem, 1.6rem);
        margin: 1.2rem 0 0.5rem;
      }
      h3 {
        font-size: 1.05rem;
        margin: 0.8rem 0 0.35rem;
      }
      .nav-link {
        color: var(--rw2) !important;
      }
      .nav-link.active {
        color: var(--rw1) !important;
        font-weight: 600;
      }
      .btn-primary {
        background-color: #609048 !important;
        border-color: #609048 !important;
        color: white !important;
      }
      .btn-primary:hover {
        background-color: #90A860 !important;
        border-color: #90A860 !important;
        color: white !important;
      }
      .btn-danger {
        background-color: #786048 !important;
        border-color: #786048 !important;
        color: white !important;
      }
      .btn-danger:hover {
        background-color: #604830 !important;
        border-color: #604830 !important;
        color: white !important;
      }
      code {
        background-color: var(--rw4) !important;
        color: var(--rw2) !important;
        padding: 2px 6px !important;
        border-radius: 3px !important;
      }
      .card {
        border-left: 3px solid var(--rw1) !important;
      }
      details summary {
        background: var(--rw4) !important;
        border-left: 3px solid var(--rw1) !important;
      }
      details summary:hover {
        background: #e8e4df !important;
      }
      hr {
        border-color: var(--rw3) !important;
        opacity: 0.3;
      }
      .hero {
        background: transparent;
        color: var(--ink);
        padding: 44px 0 24px;
        margin: -20px -15px 24px -15px;
      }
      .badge {
        display: inline-block;
        padding: 0.18rem 0.55rem;
        border-radius: 999px;
        background: rgba(255,255,255,.16);
        border: 1px solid rgba(255,255,255,.28);
        font-size: 0.8rem;
        letter-spacing: 0.02em;
        margin-bottom: 0.5rem;
      }
      .divider {
        height: 2px;
        background: linear-gradient(90deg, transparent, var(--rw3), transparent);
        margin: 16px 0;
        border: none;
      }
      .about-card {
        background: #ffffff;
        border: 1px solid #e7e6e2;
        border-radius: var(--radius);
        padding: 14px 16px;
        box-shadow: var(--shadow);
        margin-bottom: 14px;
      }
      .tag {
        display: inline-block;
        background: var(--rw1);
        color: #fff;
        padding: 0.22rem 0.55rem;
        border-radius: 10px;
        font-size: 0.8rem;
        margin: 0.05rem 0.4rem 0.35rem 0;
      }
      .tag-uscr {
        background: #609048 !important;
      }
      .tag-rem {
        background: #786048 !important;
      }
      .tag-tte {
        background: #604830 !important;
      }
      .pcbox {
        background: #fff;
        border: 1px solid #e5e3df;
        border-radius: 12px;
        padding: 10px 12px;
        box-shadow: 0 8px 20px rgba(33,48,38,.06);
        margin-top: 0.5rem;
      }
      .pcbox .pcsec + .pcsec {
        margin-top: 0.5rem;
        padding-top: 0.5rem;
        border-top: 1px dashed #ddd6cd;
      }
      .pcbox h4 {
        margin: 0.1rem 0 0.25rem;
        font-size: 0.95rem;
        color: #3C4B3E;
      }
      .pcbox ul {
        margin: 0.1rem 0 0.1rem 1rem;
      }
      .kicker {
        font-weight: 600;
        color: var(--rw5);
        letter-spacing: 0.02em;
      }
      .ref {
        display: block;
        font-size: 0.86rem;
        color: var(--muted);
        margin-top: 0.25rem;
      }
      .lead {
        font-size: 1.05rem;
      }
      .small {
        font-size: 0.92rem;
        color: var(--muted);
      }
      ol.simple {
        padding-left: 1rem;
        margin: 0.6rem 0;
      }
      ol.simple > li {
        margin: 0.35rem 0;
        padding-left: 0.3rem;
      }
      ol.simple strong {
        color: var(--rw5);
      }
      a {
        color: var(--rw1);
        text-decoration: none;
        border-bottom: 1px dotted rgba(0,0,0,.18);
      }
      .nav-tabs {
        justify-content: center !important;
      }
      .banner-logos {
        text-align: center;
        padding: 10px 0 18px;
        display: flex;
        align-items: center;
        justify-content: center;
        gap: clamp(12px, 2vw, 28px);
        flex-wrap: wrap;
      }
      .banner-logo-side {
        height: clamp(60px, 7vw, 90px);
        width: auto;
        flex: 0 0 auto;
      }
      .banner-logo-nps {
        height: clamp(88px, 10vw, 132px);
      }
      .banner-logo-main {
        height: clamp(320px, 34vw, 430px);
        width: auto;
        flex: 0 1 auto;
      }
    "))
  ),
  tags$div(
    class = "banner-logos",
    tagList(
      tags$img(src = "wvu_logo.png", alt = "WVU", class = "banner-logo-side"),
      tags$img(src = "deer_app_logo.png", alt = "DEER App", class = "banner-logo-main"),
      tags$img(src = "nps_logo.png", alt = "NPS", class = "banner-logo-side banner-logo-nps"),
      if (file.exists(file.path("www", "usgs_logo.png"))) {
        tags$img(src = "usgs_logo.png", alt = "USGS", class = "banner-logo-side")
      } else {
        tags$span(style = "display:none;")
      }
    )
  ),
  
  navset_tab(
    id = "deer_tabs",
        
        # ------------------------- ABOUT ------------------------------
        nav_panel(
          "About",
          tags$div(
            class = "hero",
            tags$div(
              style = "max-width: 980px; margin: 0 auto; padding: 0 16px;",
              tags$h1(style = "color: var(--ink); margin-top: 0.5rem; text-align: center;",
                "DEER App"
              ),
              tags$h2(style = "color: var(--ink); margin-top: 0.5rem; font-size: 1.2rem; font-weight: 400; text-align: center;",
                tags$span(style = "color: var(--rw3);", "D"), tags$span(style = "color: #666;", "ensity "),
                tags$span(style = "color: var(--rw3);", "E"), tags$span(style = "color: #666;", "stimation "),
                tags$span(style = "color: var(--rw3);", "f"), tags$span(style = "color: #666;", "rom "),
                tags$span(style = "color: var(--rw3);", "E"), tags$span(style = "color: #666;", "ncounter "),
                tags$span(style = "color: var(--rw3);", "R"), tags$span(style = "color: #666;", "ates")
              ),
              tags$p(class = "small", style = "opacity: 0.9; color: var(--muted); text-align: center;",
                "uSCR · REM · TTE — three unmarked camera methods for estimating animal density, reported as animals/mi² by default."
              )
            )
          ),
          tags$div(
            style = "max-width: 980px; margin: 0 auto; padding: 0 16px 42px;",
            
            # Welcome section
            tags$section(
              tags$h1("Welcome to the DEER app!"),
              tags$h2(style = "font-size: 1.1rem; margin-top: 0.5rem; margin-bottom: 1rem; color: var(--muted); font-weight: normal;", "What it does"),
              tags$p(class = "lead",
                "The DEER app helps estimate ", tags$strong("animal density at sites"), " from ",
                tags$strong("unmarked camera detections"), ". It fits ", tags$strong("three complementary models"),
                ", reports ", tags$strong("credible intervals and variability"), ", and lets users ",
                tags$strong("compare, diagnose, and average results"), ", rather than relying on a point estimate alone. The current upload workflow expects ",
                "data on ", tags$strong("camera deployment"), " (e.g. where and when cameras were set and began recording) and ",
                tags$strong("images CSVs"), ". The broader goal is a ", tags$strong("consistent, reproducible workflow"),
                " that estimates ", tags$strong("animal densities"), " from camera trap data without ",
                tags$strong("tedious individual marking strategies"), "."
              ),
              tags$div(class = "divider"),
              tags$div(
                style = "display: grid; grid-template-columns: repeat(auto-fit, minmax(280px, 1fr)); gap: 14px;",
                tags$div(
                  class = "about-card",
                  tags$h3("Efficient!"),
                  tags$ul(
                    tags$li("Uses unmarked detections rather than requiring marked animals or individual IDs."),
                    tags$li("Pairs field data with a streamlined upload, QC, and analysis workflow that can benefit from AI-assisted image tagging.")
                  )
                ),
                tags$div(
                  class = "about-card",
                  tags$h3("Reproducible!"),
                  tags$ul(
                    tags$li("Standardized CSV checks and shared defaults make analyses easier to repeat."),
                    tags$li("Downloadable per-model summaries support site-to-site or project-to-project comparisons.")
                  )
                ),
                tags$div(
                  class = "about-card",
                  tags$h3("Rigorous"),
                  tags$ul(
                    tags$li("Reports animal density with uncertainty, including credible intervals."),
                    tags$li("Uses three model types and a WAIC-weighted model-averaged estimate for uploaded field data."),
                    tags$li("Brings together spatial and encounter-rate approaches in one transparent workflow.")
                  )
                )
              )
            ),
            
            tags$div(class = "divider"),
            
            # How to use
            tags$section(
              tags$h2("How to use the DEER App"),
              tags$div(
                class = "about-card",
                tags$h2(style = "font-size: 1.5rem; font-weight: 600; margin-top: 1rem;", "Step 1: Simulate or Upload data"),
                tags$p("You can explore the workflow with simulated data or analyze field data."),
                tags$div(
                  style = "margin: 1rem 0;",
                  tags$h4(style = "font-size: 1.1rem; font-weight: 500;", "Option 1: Simulate data"),
                  tags$p(
                    "Go to the ", tags$strong("'Simulate data'"), " tab. The app uses ",
                    tags$strong("one shared spatial simulator"),
                    " that can feed all three model tabs."
                  ),
                  tags$ul(
                    tags$li(tags$strong("Shared spatial simulator"), " — spatial toy data generated under the uSCR assumptions, then reformatted so USCR, REM, and TTE can all run on the same simulated dataset.")
                  ),
                  tags$p("For the shared simulator you can adjust:"),
                  tags$ul(
                    tags$li("Grid dimension (n × n cameras)"),
                    tags$li("Camera spacing (meters)"),
                    tags$li("Number of days"),
                    tags$li("True density (animals/km²)"),
                    tags$li("Detection parameters (sigma, lambda0)"),
                    tags$li("Random seed for reproducibility")
                  ),
                  tags$p(
                    "Press ", tags$strong("'Simulate data'"), " to generate the shared simulated dataset, then run ",
                    tags$strong("USCR"), ", ", tags$strong("REM"), ", and ", tags$strong("TTE"), " from their own tabs."
                  ),
                  tags$p(
                    "The ", tags$strong("Compare & combine"), " tab summarizes whichever simulated model fits have finished from that shared dataset."
                  )
                ),
                tags$div(
                  style = "margin: 1rem 0;",
                  tags$h4(style = "font-size: 1.1rem; font-weight: 500;", "Option 2: Upload your field data"),
                  tags$p("You will need two data files:"),
                  tags$ul(
                    tags$li(tags$strong("Deployment CSV"), " — camera deployment information such as locations, dates, and detection distances."),
                    tags$li(tags$strong("Images CSV"), " — detection records such as timestamps, species, and ", tags$strong("Cluster ID"), " values. In this workflow, ", tags$strong("Cluster ID"), " means the unique identifier for an independent encounter event.")
                  ),
                  tags$p(
                    "The current upload pipeline expects the ", tags$strong("exact column names"), " described in the Add your data tab. Recommended camera spacing and camera counts come from the protocol, but the app analyzes the data you provide rather than requiring one exact array design or minimum camera count to function."
                  ),
                  tags$p(
                    "Click the ", tags$strong("'Add your data'"), " tab for column requirements and examples. The app will automatically:"
                  ),
                  tags$ul(
                    tags$li("Standardize column names, whitespace, and common logical values such as Yes/No or TRUE/FALSE"),
                    tags$li("Run quality checks on the deployment and images files"),
                    tags$li("Optionally trim each camera to the first 56 deployed days to meet the closed-population assumption used by all three models")
                  )
                ),
                tags$h2(style = "font-size: 1.5rem; font-weight: 600; margin-top: 1rem;", "Step 2: Adjust settings (optional)"),
                tags$p(
                  "The default settings work for many use cases. To change MCMC settings, priors, or detection geometry, open the ",
                  tags$strong("'Model settings'"), " tab, choose ", tags$strong("Advanced"), ", and edit the fields."
                ),
                tags$ul(
                  tags$li("Adjust MCMC settings such as iterations, burn-in, thinning, and number of chains."),
                  tags$li("Modify priors for movement speed, viewshed or detection parameters, and camera heterogeneity."),
                  tags$li("Change camera detection angle (default: 55°)."),
                  tags$li("Remember that priors matter most for REM and TTE because speed and viewshed assumptions directly affect the density calculation.")
                ),
                tags$h2(style = "font-size: 1.5rem; font-weight: 600; margin-top: 1rem;", "Step 3: Run the models"),
                tags$p(
                  tags$strong("Simulated data:"), " run ", tags$strong("USCR"), ", ", tags$strong("REM"), ", and ", tags$strong("TTE"), " from their own tabs after generating the shared simulated dataset."
                ),
                tags$p(
                  tags$strong("Uploaded field data:"), " use the ", tags$strong("USCR"), ", ", tags$strong("REM"), ", and ", tags$strong("TTE"),
                  " tabs. Each tab has its own run button and troubleshooting panel."
                ),
                tags$ul(
                  tags$li("Click the ", tags$strong("'Run'"), " button for your data type."),
                  tags$li("Watch progress bars, stage labels, and troubleshooting text where available."),
                  tags$li("Use the red ", tags$strong("'Stop'"), " button when a model supports stopping."),
                  tags$li("Results appear below the buttons once each model completes.")
                ),
                tags$p(
                  tags$em("Note:"), " USCR can take longer than 60 minutes for high-density sites, larger arrays, or larger ", tags$code("M"),
                  " settings, so plan to be patient and, if needed, leave the computer running. REM and TTE are usually faster."
                ),
                tags$h2(style = "font-size: 1.5rem; font-weight: 600; margin-top: 0.5rem;", "Step 4: Compare results"),
                tags$p("Use the ", tags$strong("'Compare & combine'"), " tab to:"),
                tags$ul(
                  tags$li("For ", tags$strong("uploaded field data"), ": compare whichever model fits have finished so far, and compute WAIC-weighted summaries once multiple models are available."),
                  tags$li("For ", tags$strong("simulated data"), ": compare whichever shared-simulation model fits have finished so far, just like the uploaded-data workflow."),
                  tags$li("Compare model performance using WAIC values for uploaded field data."),
                  tags$li("Download posterior summaries with parameter names, means, and 95% credible intervals.")
                ),
                tags$p(
                  "The app currently reports ", tags$strong("animals per square mile (mi²)"), " by default, with credible intervals. Some inputs and simulations still use ",
                  tags$strong("animals/km²"), " where noted."
                ),
                tags$p(
                  tags$strong("Unit convention:"), " 1 animal/km² = 2.59 animals/mi², and 1 animal/mi² = 0.386 animals/km². The current app does not yet have a global unit toggle, so values are labeled explicitly where they appear."
                )
              )
            ),
            
            tags$div(class = "divider"),
            
            # The three models
            tags$section(
              tags$h2("Quick model descriptions"),
              tags$div(
                style = "display: grid; grid-template-columns: repeat(auto-fit, minmax(280px, 1fr)); gap: 14px;",
                tags$div(
                  class = "about-card",
                  tags$span(class = "tag tag-uscr", "uSCR"),
                  tags$p(
                    tags$strong("Unmarked Spatial Capture–Recapture (uSCR)"), " — ", tags$em("Chandler & Royle, 2013"), ".",
                    tags$br(),
                    tags$strong("Spatial pattern → density."), " Where detections fall across the array informs density estimates inside the defined state-space."
                  ),
                  tags$div(
                    class = "pcbox",
                    tags$div(
                      class = "pcsec",
                      tags$h4("Pros"),
                      tags$ul(
                        tags$li("Uses spatial information across the array and estimates density inside an explicit state-space."),
                        tags$li("Produces density with uncertainty without requiring individual IDs.")
                      )
                    ),
                    tags$div(
                      class = "pcsec",
                      tags$h4("Cons"),
                      tags$ul(
                        tags$li("Most computationally intensive option in the app."),
                        tags$li("Sensitive to camera spacing, movement scale, and prior settings.")
                      )
                    )
                  ),
                  tags$span(
                    class = "ref",
                    "Reference: ",
                    tags$a(
                      href = "https://doi.org/10.1214/12-AOAS610",
                      target = "_blank",
                      style = "color: var(--rw1); text-decoration: underline;",
                      "Chandler, R.B. & Royle, J.A. (2013)"
                    )
                  )
                ),
                tags$div(
                  class = "about-card",
                  tags$span(class = "tag tag-rem", "REM"),
                  tags$p(
                    tags$strong("Random Encounter Model (REM)"), " — ", tags$em("Rowcliffe et al., 2008"), ".",
                    tags$br(),
                    tags$strong("Events per unit time → density."), " Encounter rates are converted to density using animal speed and the camera viewshed."
                  ),
                  tags$div(
                    class = "pcbox",
                    tags$div(
                      class = "pcsec",
                      tags$h4("Pros"),
                      tags$ul(
                        tags$li("Works with unmarked detections and is faster than uSCR."),
                        tags$li("Provides a practical encounter-rate approach when cameras are randomly placed and unbaited.")
                      )
                    ),
                    tags$div(
                      class = "pcsec",
                      tags$h4("Cons"),
                      tags$ul(
                        tags$li("Priors for movement speed and viewshed matter."),
                        tags$li("Assumes random, unbaited camera placement and provides less spatial detail than uSCR.")
                      )
                    )
                  ),
                  tags$span(
                    class = "ref",
                    "Reference: ",
                    tags$a(
                      href = "https://doi.org/10.1111/j.1365-2664.2008.01473.x",
                      target = "_blank",
                      style = "color: var(--rw1); text-decoration: underline;",
                      "Rowcliffe, J.M. et al. (2008)"
                    )
                  )
                ),
                tags$div(
                  class = "about-card",
                  tags$span(class = "tag tag-tte", "TTE"),
                  tags$p(
                    tags$strong("Time‑to‑Event (TTE)"), " — ", tags$em("Moeller et al., 2018"), ".",
                    tags$br(),
                    tags$strong("Waiting time → density."), " Shorter effective time between encounters implies higher density when speed and viewshed are held fixed."
                  ),
                  tags$div(
                    class = "pcbox",
                    tags$div(
                      class = "pcsec",
                      tags$h4("Pros"),
                      tags$ul(
                        tags$li("Works with unmarked detections and is faster than uSCR."),
                        tags$li("Provides a second encounter-rate perspective that complements REM.")
                      )
                    ),
                    tags$div(
                      class = "pcsec",
                      tags$h4("Cons"),
                      tags$ul(
                        tags$li("Priors for movement speed and viewshed matter here too."),
                        tags$li("Shares many of the same placement and calibration assumptions as REM.")
                      )
                    )
                  ),
                  tags$span(
                    class = "ref",
                    "Reference: ",
                    tags$a(
                      href = "https://doi.org/10.1002/ecs2.2331",
                      target = "_blank",
                      style = "color: var(--rw1); text-decoration: underline;",
                      "Moeller, A.K. et al. (2018, Ecosphere)"
                    )
                  )
                )
              )
            ),
            
            tags$div(class = "divider"),
            
            # What you get
            tags$section(
              tags$h2("What you get"),
              tags$div(
                class = "about-card",
                tags$p(
                  tags$span(class = "kicker", "Current default output:"),
                  " ",
                  tags$strong("Animals per square mile (mi²)"),
                  " with credible intervals, plus per-model estimates and downloadable posterior summaries for transparency."
                )
              )
            ),
            
            tags$div(class = "divider"),
            
            # Acknowledgments
            tags$section(
              tags$h2("Acknowledgments"),
              tags$div(
                class = "about-card",
                tags$h3(style = "font-size: 1.2rem; font-weight: 600; margin-top: 0.5rem;", "Model Development"),
                tags$p(
                  "The underlying models and code were created by ",
                  tags$strong("Dr. Amanda Van Buskirk"),
                  " with guidance from ",
                  tags$a(
                    href = "https://www.davis.wvu.edu/faculty-staff/directory/christopher-rota",
                    target = "_blank",
                    style = "color: var(--rw1); text-decoration: underline;",
                    tags$strong("Dr. Christopher Rota")
                  ),
                  " at ",
                  tags$strong("West Virginia University, Davis College of Agriculture and Natural Resources"),
                  "."
                ),
                tags$p(
                  style = "margin-top: 0.75rem;",
                  "Scientific collaboration and feedback on model integration and quality control also came from ",
                  tags$strong("Dr. Laura C. Gigliotti"),
                  ", ",
                  tags$strong("U.S. Geological Survey, West Virginia Cooperative Fish and Wildlife Research Unit, West Virginia University"),
                  "."
                )
              ),
              tags$div(
                class = "about-card",
                tags$h3(style = "font-size: 1.2rem; font-weight: 600; margin-top: 0.5rem;", "Shiny App Development"),
                tags$p(
                  "This Shiny application was developed as part of the ", 
                  tags$a(href = "https://esa.org/programs/scip/", target = "_blank", 
                         style = "color: var(--rw1); text-decoration: underline;",
                         tags$strong("Science in the Parks Communications Fellowship")), 
                  ", a collaborative effort between the ", 
                  tags$strong("Ecological Society of America (ESA)"), 
                  " and the ", 
                  tags$strong("National Park Service (NPS)"), 
                  "."
                ),
                tags$p(
                  style = "margin-top: 0.75rem;",
                  tags$strong("Fellowship Support:"), tags$br(),
                  "• ", tags$strong("Dr. Brian Mitchell"), " (NPS) — Fellowship Liaison", tags$br(),
                  "• ", tags$strong("Jasjeet Dhanota"), " (ESA) — Mentor", tags$br(),
                  "• ", tags$strong("Mary Joy Mulumba"), " (ESA) — Mentor"
                )
              )
            ),
            
            tags$div(class = "divider")
          )
        ),
        
        # ---------------------- SIMULATE DATA -------------------------
        nav_panel(
          "Simulate data",
          markdown(paste(
            "Simulate a toy camera grid under **SECR** (*spatially explicit capture–recapture*) and use that **shared spatial simulator** to run **USCR, REM, and TTE** from their model tabs for compare-ready simulated analyses.",
            "",
            "Adjust the simulation parameters below, then choose the model tab you want to test.",
            sep = "\n"
          )),
          
          h3("Simulation parameters"),
          fluidRow(
            column(4,
              sliderInput("n_side", "Grid dimension (n x n cameras)",
                          min = 3, max = 8, value = 5, step = 1),
              sliderInput("spacing", "Camera spacing (m)",
                          min = 150, max = 500, value = 300, step = 10),
              sliderInput("days", "Survey length (days)",
                          min = 7, max = 60, value = 21, step = 1)
            ),
            column(4,
              sliderInput("Dtrue", "True density for simulation (animals/km²)",
                          min = 5, max = 80, value = 25, step = 1),
              sliderInput("sigma", "Home-range scale σ (m)",
                          min = 80, max = 300, value = 150, step = 5),
              sliderInput("lambda0", "Baseline detection rate λ₀ (per occasion)",
                          min = 0.05, max = 0.6, value = 0.20, step = 0.01)
            ),
            column(4,
              numericInput("seed", "Random seed", value = 1, min = 1),
              h4("Camera geometry (simulation)"),
              sliderInput("r_m_sim", "Detection radius r (m) for simulated grid",
                          min = 8, max = 25, value = 12, step = 1),
              p(class = "small", style = "color: var(--muted);",
                "Detection angle θ for REM/TTE/USCR is set under ", tags$strong("Model settings"), "."),
              br(),
              actionButton("run_sim", "Simulate data", class = "btn-primary")
            )
          ),
          tags$details(
            class = "about-card",
            style = "margin-top: 1rem;",
            tags$summary(
              style = "cursor: pointer; font-weight: 600;",
              "Simulation parameter guide"
            ),
            tags$div(
              style = "margin-top: 0.85rem;",
              tags$p(
                tags$strong("Grid dimension"),
                " sets the number of cameras along one side of the toy grid. A 5 x 5 grid gives 25 cameras."
              ),
              tags$p(
                tags$strong("Camera spacing"),
                " is the distance between neighboring cameras. Wider spacing makes detections sparser."
              ),
              tags$p(
                tags$strong("Survey length"),
                " is the number of deployed days simulated for each camera."
              ),
              tags$p(
                tags$strong("True density"),
                " is the underlying animal density used to generate the toy data. This input is in animals/km^2, while app summaries are reported in animals/mi^2."
              ),
              tags$p(
                tags$strong("Home-range scale sigma"),
                " is the spatial scale parameter passed into the uSCR detection function. Larger values spread detections across more cameras and correspond to broader space use."
              ),
              tags$p(
                tags$strong("Baseline detection rate lambda0"),
                " is the expected number of detections per occasion if an animal's activity center were at the camera. Larger values produce more detections everywhere in the array."
              ),
              tags$p(
                tags$strong("Random seed"),
                " reproduces the same simulated dataset when you rerun the grid with identical settings."
              ),
              tags$p(
                tags$strong("Detection radius r"),
                " is not used to generate the SECR detections themselves. It is carried into the model-input table used by REM and TTE."
              )
            )
          ),
          
          hr(),
          h3("Simulated visualizations"),
          fluidRow(
            column(6,
              h4("Simulated field grid"),
              plotOutput("grid_plot", height = "420px")
            ),
            column(6,
              h4("Simulated animal distribution"),
              plotOutput("sim_deer_distribution_plot", height = "420px")
            )
          ),
          h3("Camera data table"),
          DTOutput("camera_table"),
        ),

        # ------------------- CAMERA ARRAY DESIGN --------------------
        nav_panel(
          "Camera array design",
          markdown(paste(
            "Design a camera layout for a site from uploaded spatial layers.",
            "Upload a boundary and any optional exclusion layers, choose spacing and buffer settings, and the app will suggest a camera array you can preview and export.",
            "",
            "This first in-app version uses the collaborator workflow as a base, but keeps everything inside the DEER app instead of launching a separate selector app.",
            sep = "\n"
          )),
          fluidRow(
            column(
              4,
              h3("Step 1: Upload site layers"),
              p(
                "Boundary is required. Optional layers can keep cameras away from roads, trails, buildings, parking lots, or other excluded areas."
              ),
              fileInput(
                "camera_boundary_files",
                "Boundary file(s)",
                multiple = TRUE,
                accept = c(".zip", ".shp", ".dbf", ".shx", ".prj", ".cpg", ".kml", ".kmz", ".gpkg")
              ),
              fileInput(
                "camera_roads_files",
                "Roads layer (optional)",
                multiple = TRUE,
                accept = c(".zip", ".shp", ".dbf", ".shx", ".prj", ".cpg", ".kml", ".kmz", ".gpkg")
              ),
              fileInput(
                "camera_trails_files",
                "Trails layer (optional)",
                multiple = TRUE,
                accept = c(".zip", ".shp", ".dbf", ".shx", ".prj", ".cpg", ".kml", ".kmz", ".gpkg")
              ),
              fileInput(
                "camera_buildings_files",
                "Buildings layer (optional)",
                multiple = TRUE,
                accept = c(".zip", ".shp", ".dbf", ".shx", ".prj", ".cpg", ".kml", ".kmz", ".gpkg")
              ),
              fileInput(
                "camera_parking_files",
                "Parking lots layer (optional)",
                multiple = TRUE,
                accept = c(".zip", ".shp", ".dbf", ".shx", ".prj", ".cpg", ".kml", ".kmz", ".gpkg")
              ),
              fileInput(
                "camera_exclusion_files",
                "Additional exclusion layer(s) (optional)",
                multiple = TRUE,
                accept = c(".zip", ".shp", ".dbf", ".shx", ".prj", ".cpg", ".kml", ".kmz", ".gpkg")
              )
            ),
            column(
              4,
              h3("Step 2: Choose layout settings"),
              textInput("camera_prefix", "Camera name prefix", value = "CAM"),
              sliderInput(
                "camera_spacing_m",
                "Maximum target spacing between nearby cameras (m)",
                min = 100, max = 600, value = 400, step = 10
              ),
              sliderInput(
                "camera_buffer_m",
                "Buffer around roads / trails / buildings / parking (m)",
                min = 0, max = 100, value = 10, step = 1
              ),
              numericInput(
                "camera_budget",
                "Camera budget (0 = keep all suggested cameras)",
                value = 0,
                min = 0,
                step = 1
              ),
              numericInput(
                "camera_alternates",
                "Alternate cameras to keep when budget is smaller than the suggested layout",
                value = 2,
                min = 0,
                step = 1
              ),
              actionButton("run_camera_array", "Generate camera array", class = "btn-primary")
            ),
            column(
              4,
              h3("What this produces"),
              tags$ul(
                tags$li("A projected design area with optional exclusion buffers removed."),
                tags$li("Suggested camera points spaced across the allowable area."),
                tags$li("A final camera set, plus alternates if you enter a smaller budget."),
                tags$li("Downloadable camera coordinates in latitude/longitude and UTM.")
              ),
              tags$p(
                class = "small",
                style = "color: var(--muted);",
                "Tip: if the suggested layout is denser than you can deploy, enter your budget and the app will keep a well-spread subset plus optional alternates."
              )
            )
          ),
          hr(),
          h3("Camera-array log"),
          verbatimTextOutput("camera_array_log"),
          fluidRow(
            column(
              5,
              h3("Layout summary"),
              verbatimTextOutput("camera_array_summary"),
              downloadButton("download_camera_array_csv", "Download camera CSV"),
              tags$span(" "),
              downloadButton("download_camera_array_gpkg", "Download camera GeoPackage"),
              tags$span(" "),
              downloadButton("download_camera_area_gpkg", "Download design-area GeoPackage")
            ),
            column(
              7,
              h3("Camera-array map"),
              leafletOutput("camera_array_map", height = "520px")
            )
          ),
          h3("Camera coordinates"),
          DTOutput("camera_array_table")
        ),
        
        # ---------------------- ADD YOUR DATA -------------------------
        nav_panel(
          "Add your data",
          markdown(paste(
            "Upload a **deployment CSV** (where and when cameras were set and recording) and an **images CSV** (timestamps, species, counts, and Cluster IDs).",
            "The current upload checker expects the exact column names listed below.",
            "",
            sep = "\n"
          )),
          h3("The app will:"),
          markdown(paste(
            "1. Standardize column names, whitespace, and common logical values;",
            "2. Check required columns, data types, matching site names, and deployment/image consistency;",
            "3. Flag image timestamps outside each camera's deployment window (±3 days);",
            "4. Optionally trim images at each camera to the first 56 deployed days to meet the closed-population assumption used by all three models.",
            "",
            "Download processed CSVs if needed, then configure MCMC and priors in the **Model settings** tab and run models in the USCR/REM/TTE tabs.",
            "",
            "---",
            "",
            "## Required data columns",
            "",
            "### Deployment file columns (**bold = required**, others optional for quality control):",
            "",
            "- **`Site Name`** — Camera site identifier",
            "- `Site` — Optional site identifier (the app also accepts `Park` for legacy files)",
            "- `Camera ID` — Optional camera identifier",
            "- `SD Card ID` — Optional SD card identifier",
            "- **`Start Date`** — Deployment start date (`MM/DD/YYYY`)",
            "- **`Start Time`** — Deployment start time (`HH:MM:SS` when available)",
            "- **`End Date`** — Deployment end date (`MM/DD/YYYY`)",
            "- **`End Time`** — Deployment end time (`HH:MM:SS` when available)",
            "- **`Latitude`** — Camera latitude (decimal degrees)",
            "- **`Longitude`** — Camera longitude (decimal degrees)",
            "- `Camera Height` — Optional camera height",
            "- `Camera Orientation` — Optional cardinal direction or 0-359 degrees",
            "- **`Camera Functioning`** — Camera status; common `Yes/No`, `TRUE/FALSE`, `T/F`, and `1/0` values are normalized during import",
            "- **`Camera Malfunction Date`** — Keep this column in the file; fill it when `Camera Functioning = No` for a site with images",
            "- **`Detection Distance`** — Effective detection radius in meters",
            "- `Notes` — Optional notes column",
            "",
            "### **Images file** required columns:",
            "",
            "- **`Site Name`** — Must match deployment file",
            "- **`Latitude`** — Image/site latitude in decimal degrees",
            "- **`Longitude`** — Image/site longitude in decimal degrees",
            "- **`Timestamp`** — Image timestamp (date and time)",
            "- **`Species`** — Species identifier (e.g., your target species name)",
            "- **`Cluster ID`** — Unique identifier for independent detection events",
            "- **`Sighting Count`** — Number of individuals in the image; pipe-delimited values such as `1|1` are allowed for multi-species rows",
            "- `Image URL` — Optional image reference/link column that can help with later validation",
            "",
            "**Notes:** Deployment QC expects date-only fields in `MM/DD/YYYY`. Image timestamps are parsed more flexibly, including 2-digit years (for example `2/3/25`).",
            "The app expects the exact column names listed above. TrapTagger exports already use the images-column names expected here, and files from other tagging software can usually be renamed to match.",
            "Cross-year winter surveys (for example `12/2025` to `01/2026`) are supported; the app uses the actual deployment dates/times and image timestamps, so no separate `Survey Year` field is required.",
            sep = "\n"
          )),
          
          h3("Step 1: Deployment file"),
          fileInput("deployment_csv", "Upload deployment CSV", accept = ".csv"),
          h4("Deployment check log"),
          verbatimTextOutput("deployment_check_log"),
          h4("Preview of cleaned deployment data"),
          DTOutput("deployment_preview"),
          downloadButton("download_deployment_checked", "Download cleaned deployment CSV"),
          
          hr(),
          
          h3("Step 2: Images file"),
          fileInput("images_csv", "Upload images CSV", accept = ".csv"),
          checkboxInput(
            "apply_56day_trim",
            "After validation, trim images at each camera to the first 56 deployed days (recommended for meeting model assumptions). If no data were collected after 56 days, the dataset is left unchanged.",
            value = TRUE
          ),
          h4("Images check log"),
          verbatimTextOutput("images_check_log"),
          h4("Preview of processed images data"),
          DTOutput("images_preview"),
          downloadButton("download_images_checked", "Download processed images CSV"),
          
          hr(),
          
          h3("Model settings"),
          p(
            "Configure MCMC iterations, priors, and detection-angle settings (",
            tags$strong("θ"), ") in the ",
            tags$strong("Model settings"), " tab."
          )
        ),
        
        nav_panel(
          "Model settings",
          markdown(paste(
            "### MCMC settings",
            "",
            "**Recommended starting values in the app**:",
            "",
            "- **REM/TTE**: 6000 iterations, 1000 burn-in, thin = 5",
            "- **USCR**: 6000 iterations, 1000 burn-in, thin = 5, M = 300",
            "- **Number of chains**: 3 by default; the app enforces a minimum of 3 chains for model runs",
            "- **Convergence criteria**: R̂ < 1.1 for all monitored parameters when multiple chains are used",
            "",
            "USCR now does a short adaptive tuning phase before the final run. If convergence is poor, doubling iterations is usually more helpful than increasing thinning alone. If the posterior `psi` stays high or posterior `N` approaches `M`, increase `M` and rerun.",
            "",
            "### Meta-analysis and pooled priors (not implemented in-app)",
            "",
            "Each model run uses its own informative priors. The current defaults are based on the winter home-range and daily movement literature review used to build this app.",
            "If you use the app for other species or different field conditions (for example season, biome, or sex-age class), adjust the `sigma` and `v` priors under **Advanced** to reflect your system. If you are unsure, you can widen the priors, but uncertainty in the density estimate will usually increase.",
            "The app does **not** currently pool information across sites or studies.",
            "If you combine results outside the app, use outputs that keep uncertainty, such as posterior samples or credible intervals, rather than point estimates alone.",
            "",
            "### Camera geometry",
            "",
            "Default detection angle is 55° (based on Browning-style camera specifications). Adjust below if you used a different camera or measured detection angle in the field.",
            "",
            sep = "\n"
          )),
          
          h4("Camera geometry"),
          fluidRow(
            column(12,
              sliderInput("theta", "Detection angle θ (degrees)",
                          min = 20, max = 80, value = 55, step = 1)
            )
          ),
          
          hr(),
          
          radioButtons(
            "mode", "Settings mode",
            choices  = c("Default", "Advanced"),
            selected = "Default",
            inline   = TRUE
          ),
          
          conditionalPanel(
            "input.mode == 'Default'",
            markdown(paste(
              "**Using speed-aware app defaults** for MCMC and priors.",
              "",
              "These settings are lighter than a full offline analysis and are meant to keep the app usable.",
              "Switch to 'Advanced' to edit priors, data augmentation M, and MCMC settings.",
              sep = "\n"
            ))
          ),
          
          conditionalPanel(
            "input.mode == 'Advanced'",
            tagList(
              h4("MCMC settings"),
              fluidRow(
                column(4,
                  numericInput("n_chains", "Number of chains",
                               value = 3, min = 3, max = 8, step = 1)
                ),
                column(4,
                  numericInput("iter_rem_tte", "REM/TTE iterations",
                               value = 6000, min = 1000, step = 1000),
                  numericInput("burnin_rem_tte", "REM/TTE burn-in",
                               value = 1000, min = 0, step = 500),
                  numericInput("thin_rem_tte", "REM/TTE thinning",
                               value = 5, min = 1, step = 1)
                ),
                column(4,
                  numericInput("iter_uscr", "USCR iterations",
                               value = 6000, min = 1000, step = 1000),
                  numericInput("burnin_uscr", "USCR burn-in",
                               value = 1000, min = 0, step = 500),
                  numericInput("thin_uscr", "USCR thinning",
                               value = 5, min = 1, step = 1),
                  numericInput("M_uscr", "USCR M (data augmentation)",
                               value = 300, min = 100, step = 100)
                )
              ),
              
              markdown(paste(
                "**Note on M (data augmentation)**: The app starts lower for speed. If posterior N approaches M,",
                "or if posterior psi concentrates near 1, increase M and rerun.",
                "",
                "**Note on USCR iterations**: If convergence warnings appear, try doubling iterations before changing thinning.",
                "The app now does a short adaptive tuning phase before the final run.",
                sep = "\n"
              )),
              
              hr(),
              
              h4("Priors: REM & TTE"),
              markdown(paste(
                "**Movement speed is the main informative prior here; the density and heterogeneity bounds are mainly safeguards:**",
                "",
                "- **D** ~ Uniform(0, D_max): Upper bound on density, not meant to be strongly informative. If the posterior presses against `D_max`, increase it.",
                "- **log(v)** ~ Normal(mean, SD): Informative prior on animal movement speed (log scale). The current default corresponds to an average winter daily movement rate of `3.09 km/day` with an approximate 95% interval of `1.60` to `6.00 km/day`.",
                "- **sd_eps** ~ Uniform(0, max): Upper bound on camera-to-camera heterogeneity. If the posterior presses against the maximum, increase the bound.",
                "",
                "**Example modifications**:",
                "- Wider uncertainty in movement speed: log_v ~ Normal(1.130, 0.45)",
                "- Allowing higher densities: D ~ Uniform(0, 400)",
                "- Allowing more camera heterogeneity: sd_eps ~ Uniform(0, 20)",
                "",
                sep = "\n"
              )),
              
              fluidRow(
                column(6,
                  numericInput("D_max", "Max density prior (D max, animals/km²)",
                               value = 200, min = 10, step = 10),
                  numericInput("log_v_mean", "log(v) mean (km/day)",
                               value = 1.130, step = 0.1),
                  numericInput("log_v_sd", "log(v) SD",
                               value = 0.3372, min = 0.01, step = 0.05)
                ),
                column(6,
                  numericInput("sd_eps_max", "sd_eps upper bound",
                               value = 10, min = 1, step = 1)
                )
              ),
              
              hr(),
              
              h4("Priors: USCR"),
              markdown(paste(
                "**USCR uses a biologically informative prior on space use plus weaker priors on other terms:**",
                "",
                "- **log(σ)** ~ Normal(mean, SD): Informative prior on the space-use scale. The current default is `log(σ) ~ Normal(-1.5269, 0.1535)`, which corresponds to about `σ = 0.22 km` and is based on an average winter 95% home-range size of `0.89 km²` with an approximate interval of `0.48` to `1.62 km²`.",
                "- **log(λ₀)** ~ Normal(0, SD): Baseline detection intensity. Because `exp(0) = 1`, the default prior is centered near one expected detection per deployed day if an animal's activity center were at the camera.",
                "- **psi** ~ Uniform(0, 1): Data-augmentation inclusion probability. This is usually left alone in the app.",
                "- **sd_eps** ~ Gamma(shape, rate): Camera-level random effect. The app exposes both shape and rate; change these only if you have a specific prior in mind.",
                "",
                "**Example modifications**:",
                "- Looser movement range prior: log_sigma ~ Normal(-1.5269, 0.30)",
                "- Higher uncertainty in baseline detection: log_lam_0 ~ Normal(0, 2)",
                "- Lighter penalty on large camera effects: sd_eps ~ Gamma(1, 0.5)",
                "",
                "⚠️ **Model performance warning**: If the posterior `psi` peaks near 1 or posterior `N` approaches `M`, increase `M` and/or run longer chains.",
                "",
                sep = "\n"
              )),
              
              fluidRow(
                column(6,
                  numericInput("log_sigma_mean", "log(σ) mean",
                               value = -1.5269, step = 0.1),
                  numericInput("log_sigma_sd", "log(σ) SD",
                               value = 0.1535, min = 0.01, step = 0.05),
                  numericInput("log_lam0_sd", "log(λ₀) SD",
                               value = 1, min = 0.1, step = 0.1)
                ),
                column(6,
                  numericInput("sd_eps_shape", "sd_eps gamma shape",
                               value = 1, min = 0.1, step = 0.1),
                  numericInput("sd_eps_rate", "sd_eps gamma rate",
                               value = 1, min = 0.1, step = 0.1)
                )
              )
            )
          ),
          
          # Hidden defaults for Default mode
          conditionalPanel(
            "input.mode == 'Default'",
            tags$div(
              style = "display:none;",
              numericInput("n_chains", NULL, value = 3),
              numericInput("iter_rem_tte", NULL, value = 6000),
              numericInput("burnin_rem_tte", NULL, value = 1000),
              numericInput("thin_rem_tte", NULL, value = 5),
              numericInput("iter_uscr", NULL, value = 6000),
              numericInput("burnin_uscr", NULL, value = 1000),
              numericInput("thin_uscr", NULL, value = 5),
              numericInput("M_uscr", NULL, value = 300),
              numericInput("D_max", NULL, value = 200),
              numericInput("log_v_mean", NULL, value = 1.130),
              numericInput("log_v_sd", NULL, value = 0.3372),
              numericInput("sd_eps_max", NULL, value = 10),
              numericInput("log_sigma_mean", NULL, value = -1.5269),
              numericInput("log_sigma_sd", NULL, value = 0.1535),
              numericInput("log_lam0_sd", NULL, value = 1),
              numericInput("sd_eps_shape", NULL, value = 1),
              numericInput("sd_eps_rate", NULL, value = 1)
            )
          )
        ),
        
        # ---------------------- DATA SUMMARY --------------------------
        nav_panel(
          "Data summary",
          h4("Deployment summary"),
          DTOutput("deploy_summary_table"),
          leafletOutput("camera_map", height = "300px"),
          hr(),
          h4("Image summary by species"),
          DTOutput("image_summary_table"),
          plotOutput("species_bar_plot", height = "300px"),
          hr(),
          h4("Detection details by species"),
          selectInput(
            "summary_species",
            "Species for site table and plots",
            choices = character(0)
          ),
          DTOutput("deer_summary_table"),
          plotOutput("deer_bubble_plot", height = "300px"),
          plotOutput("deer_daily_plot", height = "300px")
        ),
        
        # ---------------------- USCR MODEL TAB ------------------------
        nav_panel(
          "USCR model",
          tags$div(
            id = "uscr-content",
            HTML('
            <h2>Model 1 — uSCR (Unmarked Spatial Capture–Recapture)</h2>
            <p><strong>The gist:</strong> Estimates <strong>animal density in the site/state‑space</strong> by learning from <strong>where and when</strong> animals are detected <strong>across a camera array</strong>.</p>
            <details>
              <summary style="cursor: pointer; font-weight: 600; margin: 1rem 0; padding: 0.5rem; background: #f0f4e8; border-left: 3px solid #609048; border-radius: 6px;"><strong>Under the hood (equations)</strong></summary>
              <div style="margin: 1rem 0; padding-left: 1rem;">
                <p>Per-camera observation model:</p>
                <p>$$y_j \\sim \\mathrm{Poisson}(\\mu_j), \\qquad \\log \\mu_j = \\log(\\mathrm{days}_j) + \\log(\\Lambda_j) + \\epsilon_j$$</p>
                <p>Expected encounter rate from augmented individuals:</p>
                <p>$$\\Lambda_j = \\sum_{i=1}^{M} z_i\\,\\lambda_0\\,\\exp\\!\\Big(-\\frac{d_{ij}^2}{2\\sigma^2}\\Big), \\qquad z_i \\sim \\mathrm{Bernoulli}(\\psi)$$</p>
                <p>Activity centers \\(s_i\\) are uniform over the buffered state-space, with a constraint that each center remains within the allowed buffer distance of at least one camera. Camera effects follow \\(\\epsilon_j \\sim \\mathcal{N}(0, sd_\\epsilon)\\).</p>
                <p>Density is obtained from \\(N = \\sum_i z_i\\) divided by the study-area size in mi². With real coordinates, area comes from the buffered camera state-space; on the toy simulator, area comes from the rectangular simulated state-space. The app reports \\(D_{\\mathrm{mi}^2}\\).</p>
                <h3>Variables</h3>
                <ul>
                  <li>\\(y_j\\) — total independent animal detections at camera \\(j\\) across the deployment, built from unique <code>Cluster ID</code> values.</li>
                  <li>\\(\\mathrm{days}_j\\) — number of deployed days for camera \\(j\\).</li>
                  <li>\\(\\mu_j\\) — expected number of animal detections at camera \\(j\\) during the deployment.</li>
                  <li>\\(\\lambda_0\\) — baseline expected detections per deployed day if an animal\'s activity center were at the camera.</li>
                  <li>\\(\\sigma\\) — spatial scale parameter (km), related to home-range size and how quickly detection falls with distance.</li>
                  <li>\\(d_{ij}\\) — projected distance (km) from augmented individual \\(i\\)\'s activity center to camera \\(j\\).</li>
                  <li>\\(M\\) — augmented population size used to estimate density.</li>
                  <li>\\(z_i\\) — indicator that augmented individual \\(i\\) is part of the real population.</li>
                  <li>\\(\\epsilon_j\\) — camera-level random effect, with \\(\\epsilon_j \\sim \\mathcal{N}(0, sd_\\epsilon)\\).</li>
                  <li>\\(N = \\sum_i z_i\\) — estimated abundance inside the state-space.</li>
                  <li>\\(D_{\\mathrm{mi}^2}\\) — density, animals per square mile (reported).</li>
                </ul>
                <h3>Priors used in the app</h3>
                <ul>
                  <li>\\(\\log \\sigma \\sim \\mathcal{N}(-1.5269,\\,0.1535)\\), which anchors \\(\\sigma\\) near 0.22 km.</li>
                  <li>\\(\\log \\lambda_0 \\sim \\mathcal{N}(0,\\,1)\\), centered near one expected detection per day at distance 0.</li>
                  <li>\\(\\psi \\sim \\mathcal{U}(0,1)\\) for data augmentation.</li>
                  <li>\\(sd_\\epsilon \\sim \\mathrm{Gamma}(\\text{shape},\\text{rate})\\), with the default app values set to 1 and 1.</li>
                </ul>
              </div>
            </details>
            '),
            tags$script(HTML("
              function processMathJax(element) {
                if (window.MathJax) {
                  if (window.MathJax.typesetPromise) {
                    MathJax.typesetPromise([element]).catch(function(err) {
                      console.log('MathJax error:', err);
                    });
                  } else if (window.MathJax.typeset) {
                    MathJax.typeset([element]);
                  }
                }
              }
              setTimeout(function() {
                var content = document.getElementById('uscr-content');
                if (content) {
                  processMathJax(content);
                  var details = content.querySelector('details');
                  if (details) {
                    details.addEventListener('toggle', function(e) {
                      if (e.target.open) {
                        setTimeout(function() {
                          processMathJax(content);
                        }, 300);
                      }
                    });
                  }
                }
              }, 500);
            "))
          ),
          hr(),
          h4("Simulated data"),
          div(
            actionButton("run_uscr_sim", "Run USCR on simulated data", class = "btn-primary"),
            actionButton("stop_uscr_sim", "Stop", class = "btn-danger", style = "margin-left: 10px;"),
            style = "margin-bottom: 10px;"
          ),
          br(),
          verbatimTextOutput("uscr_sim_text"),
          h5("Run status & troubleshooting"),
          verbatimTextOutput("uscr_sim_debug"),
          
          hr(),
          h4("NPS data"),
          div(
            actionButton("run_uscr_nps", "Run USCR on NPS data", class = "btn-primary"),
            actionButton("stop_uscr_nps", "Stop", class = "btn-danger", style = "margin-left: 10px;"),
            style = "margin-bottom: 10px;"
          ),
          br(),
          verbatimTextOutput("uscr_nps_text"),
          h5("Run status & troubleshooting"),
          verbatimTextOutput("uscr_nps_debug"),
          tags$div(
            style = "text-align: center; margin-top: 40px; padding-top: 20px; border-top: 1px solid #e0e0e0;",
            tags$h3(style = "margin: 0; color: var(--rw3);", "DEER App"),
            tags$p(style = "margin: 5px 0 0 0; color: var(--muted);", "Density Estimation from Encounter Rates")
          )
        ),
        
        # ---------------------- REM MODEL TAB -------------------------
        nav_panel(
          "REM model",
          tags$div(
            id = "rem-content",
            HTML('
            <h2>Model 2 — REM (Random Encounter Model)</h2>
            <p><strong>The gist:</strong> Converts <strong>how often animals pass a camera</strong> into density, correcting for <strong>movement speed</strong> and <strong>camera view geometry</strong>.</p>
            <details>
              <summary style="cursor: pointer; font-weight: 600; margin: 1rem 0; padding: 0.5rem; background: #f0f4e8; border-left: 3px solid #609048; border-radius: 6px;"><strong>Under the hood (equations)</strong></summary>
              <div style="margin: 1rem 0; padding-left: 1rem;">
                <p>Per-camera Poisson model:</p>
                <p>$$y_j \\sim \\mathrm{Poisson}(\\lambda_j)$$</p>
                <p>$$\\log \\lambda_j = \\log D + \\log(\\mathrm{camera\\_days}_j) + \\log v + \\log r_j + \\log\\!\\Big(\\frac{2 + \\theta_{\\mathrm{rad}}}{\\pi}\\Big) + \\epsilon_j$$</p>
                <p>The app takes the user-entered camera angle in degrees and converts it internally to \\(\\theta_{\\mathrm{rad}} = \\theta_{\\mathrm{deg}}\\pi/180\\). Camera effects follow \\(\\epsilon_j \\sim \\mathcal{N}(0, sd_\\epsilon)\\), and the reported density is \\(D_{\\mathrm{mi}^2} = 2.59\\,D_{\\mathrm{km}^2}\\).</p>
                <h3>Variables</h3>
                <ul>
                  <li>\\(y_j\\) — number of independent animal detection events at camera \\(j\\) (unique <code>Cluster ID</code>).</li>
                  <li>\\(\\mathrm{camera\\_days}_j\\) — deployed days for camera \\(j\\).</li>
                  <li>\\(D\\) — animal density in animals/km² before converting to animals/mi² for display.</li>
                  <li>\\(v\\) — animal movement speed (km/day).</li>
                  <li>\\(r_j\\) — effective detection radius in km; the app converts meters from <code>Detection Distance</code> to km.</li>
                  <li>\\(\\theta_{\\mathrm{deg}}\\) — full detection angle supplied by the user (default \\(55^\\circ\\)); converted internally to radians.</li>
                  <li>\\(\\epsilon_j\\) — camera-level random effect, with \\(\\epsilon_j \\sim \\mathcal{N}(0, sd_\\epsilon)\\).</li>
                  <li>\\(D_{\\mathrm{km}^2}, D_{\\mathrm{mi}^2}\\) — density in animals/km² and animals/mi² (reported in mi²).</li>
                </ul>
                <h3>Priors/inputs used in the app</h3>
                <ul>
                  <li>\\(D \\sim \\mathcal{U}(0, D_{\\max})\\), which acts as an upper bound rather than a strongly informative prior.</li>
                  <li>\\(\\log v \\sim \\mathcal{N}(1.130,\\,0.3372)\\), an informative movement prior that currently uses the app winter default values.</li>
                  <li>\\(\\theta_{\\mathrm{deg}} = 55^\\circ\\) by default, and cameras are assumed to be randomly placed and unbaited.</li>
                  <li>\\(sd_\\epsilon \\sim \\mathcal{U}(0, 10)\\), an upper bound on camera heterogeneity.</li>
                </ul>
              </div>
            </details>
            '),
            tags$script(HTML("
              function processMathJax(element) {
                if (window.MathJax) {
                  if (window.MathJax.typesetPromise) {
                    MathJax.typesetPromise([element]).catch(function(err) {
                      console.log('MathJax error:', err);
                    });
                  } else if (window.MathJax.typeset) {
                    MathJax.typeset([element]);
                  }
                }
              }
              setTimeout(function() {
                var content = document.getElementById('rem-content');
                if (content) {
                  processMathJax(content);
                  var details = content.querySelector('details');
                  if (details) {
                    details.addEventListener('toggle', function(e) {
                      if (e.target.open) {
                        setTimeout(function() {
                          processMathJax(content);
                        }, 300);
                      }
                    });
                  }
                }
              }, 500);
            "))
          ),
          hr(),
          h4("Simulated data"),
          p(
            style = "max-width: 52rem;",
            "REM runs on the",
            tags$strong("shared spatial simulator"),
            " generated in the Simulate data tab, so it can be compared and combined with simulated USCR and TTE."
          ),
          div(
            actionButton("run_rem_sim", "Run REM on simulated data", class = "btn-primary"),
            style = "margin-bottom: 10px;"
          ),
          br(),
          verbatimTextOutput("rem_sim_text"),
          h5("Run status & troubleshooting"),
          verbatimTextOutput("rem_sim_debug"),
          
          hr(),
          h4("NPS data"),
          div(
            actionButton("run_rem_nps", "Run REM on NPS data", class = "btn-primary"),
            actionButton("stop_rem_nps", "Stop", class = "btn-danger", style = "margin-left: 10px;"),
            style = "margin-bottom: 10px;"
          ),
          p(
            style = "max-width: 52rem; color: var(--muted); margin-bottom: 0.75rem;",
            "Server note: uploaded-data REM runs still start in a background worker so other sessions can keep using the app. This path now honors the app minimum of 3 chains, so server runs may take longer than before. The Stop button cannot cancel an already-started background REM job yet."
          ),
          br(),
          verbatimTextOutput("rem_nps_text"),
          h5("Run status & troubleshooting"),
          verbatimTextOutput("rem_nps_debug"),
          tags$div(
            style = "text-align: center; margin-top: 40px; padding-top: 20px; border-top: 1px solid #e0e0e0;",
            tags$h3(style = "margin: 0; color: var(--rw3);", "DEER App"),
            tags$p(style = "margin: 5px 0 0 0; color: var(--muted);", "Density Estimation from Encounter Rates")
          )
        ),
        
        # ---------------------- TTE MODEL TAB -------------------------
        nav_panel(
          "TTE model",
          tags$div(
            id = "tte-content",
            HTML('
            <h2>Model 3 — TTE (Time‑to‑Event)</h2>
            <p><strong>The gist:</strong> Uses <strong>animal detection events per camera</strong>, scaled by <strong>movement-based time units</strong> and the <strong>viewshed area</strong>. Shorter effective time between encounters implies higher density.</p>
            <details>
              <summary style="cursor: pointer; font-weight: 600; margin: 1rem 0; padding: 0.5rem; background: #f0f4e8; border-left: 3px solid #609048; border-radius: 6px;"><strong>Under the hood (equations)</strong></summary>
              <div style="margin: 1rem 0; padding-left: 1rem;">
                <p>Per-camera Poisson model:</p>
                <p>$$y_j \\sim \\mathrm{Poisson}(\\lambda_j), \\qquad \\log \\lambda_j \\;=\\; \\log D \\;+\\; \\log(\\mathrm{tte\\_units}_j)\\;+\\;\\log(A_j)\\;+\\;\\epsilon_j$$</p>
                <p>Movement‑based time‑units and viewshed area:</p>
                <p>$$\\mathrm{tte\\_units}_j \\;=\\; \\frac{\\mathrm{camera\\_days}_j}{\\mathrm{time\\_unit}}, \\qquad \\mathrm{time\\_unit} \\;\\approx\\; \\frac{0.59\\, r_j}{v}$$</p>
                <p>$$A_j \\;=\\; \\pi r_j^2 \\frac{\\theta_{\\mathrm{deg}}}{360}, \\qquad \\epsilon_j \\sim \\mathcal{N}(0, sd_\\epsilon)$$</p>
                <p>Here \\(r_j\\) is in km and \\(\\theta_{\\mathrm{deg}}\\) is the user-entered detection angle in degrees. The app reports \\(D_{\\mathrm{mi}^2} = 2.59\\,D_{\\mathrm{km}^2}\\).</p>
                <h3>Variables</h3>
                <ul>
                  <li>\\(y_j\\) — number of animal detection events for camera \\(j\\) in the current app workflow.</li>
                  <li>\\(\\mathrm{camera\\_days}_j\\) — total deployed days for camera \\(j\\) from <code>Start/End</code>.</li>
                  <li>\\(\\mathrm{time\\_unit}\\) — expected time to traverse the viewshed once; \\(\\approx 0.59\\,r_j/v\\) (days).</li>
                  <li>\\(v\\) — movement speed (km/day); prior on \\(\\log v\\) as below.</li>
                  <li>\\(r_j\\) — effective detection radius (km), converted from meters in <code>Detection Distance</code>.</li>
                  <li>\\(\\theta_{\\mathrm{deg}}\\) — full detection angle supplied by the user (default \\(55^\\circ\\)).</li>
                  <li>\\(A_j\\) — viewshed area for camera \\(j\\) (km²).</li>
                  <li>\\(\\epsilon_j\\) — camera-level random effect, with \\(\\epsilon_j \\sim \\mathcal{N}(0, sd_\\epsilon)\\).</li>
                  <li>\\(D_{\\mathrm{mi}^2}\\) — density (animals/mi²), reported.</li>
                </ul>
                <h3>Priors/inputs used in the app</h3>
                <ul>
                  <li>\\(D \\sim \\mathcal{U}(0, D_{\\max})\\), which acts as an upper bound rather than a strongly informative prior.</li>
                  <li>\\(\\log v \\sim \\mathcal{N}(1.130,\\,0.3372)\\), which calibrates the movement-based time unit.</li>
                  <li>\\(\\theta_{\\mathrm{deg}} = 55^\\circ\\) by default, with \\(r_j\\) taken from field-measured <code>Detection Distance</code>.</li>
                  <li>\\(sd_\\epsilon \\sim \\mathcal{U}(0, 10)\\), an upper bound on camera heterogeneity.</li>
                </ul>
              </div>
            </details>
            '),
            tags$script(HTML("
              function processMathJax(element) {
                if (window.MathJax) {
                  if (window.MathJax.typesetPromise) {
                    MathJax.typesetPromise([element]).catch(function(err) {
                      console.log('MathJax error:', err);
                    });
                  } else if (window.MathJax.typeset) {
                    MathJax.typeset([element]);
                  }
                }
              }
              setTimeout(function() {
                var content = document.getElementById('tte-content');
                if (content) {
                  processMathJax(content);
                  var details = content.querySelector('details');
                  if (details) {
                    details.addEventListener('toggle', function(e) {
                      if (e.target.open) {
                        setTimeout(function() {
                          processMathJax(content);
                        }, 300);
                      }
                    });
                  }
                }
              }, 500);
            "))
          ),
          hr(),
          h4("Simulated data"),
          p(
            style = "max-width: 52rem;",
            "TTE runs on the",
            tags$strong("shared spatial simulator"),
            " generated in the Simulate data tab, so it can be compared and combined with simulated USCR and REM."
          ),
          div(
            actionButton("run_tte_sim", "Run TTE on simulated data", class = "btn-primary"),
            style = "margin-bottom: 10px;"
          ),
          br(),
          verbatimTextOutput("tte_sim_text"),
          h5("Run status & troubleshooting"),
          verbatimTextOutput("tte_sim_debug"),
          
          hr(),
          h4("NPS data"),
          div(
            actionButton("run_tte_nps", "Run TTE on NPS data", class = "btn-primary"),
            actionButton("stop_tte_nps", "Stop", class = "btn-danger", style = "margin-left: 10px;"),
            style = "margin-bottom: 10px;"
          ),
          br(),
          verbatimTextOutput("tte_nps_text"),
          h5("Run status & troubleshooting"),
          verbatimTextOutput("tte_nps_debug"),
          tags$div(
            style = "text-align: center; margin-top: 40px; padding-top: 20px; border-top: 1px solid #e0e0e0;",
            tags$h3(style = "margin: 0; color: var(--rw3);", "DEER App"),
            tags$p(style = "margin: 5px 0 0 0; color: var(--muted);", "Density Estimation from Encounter Rates")
          )
        ),
        
        # ---------------------- COMPARE & COMBINE --------------------
        nav_panel(
          "Compare & combine",
          markdown(paste(
            "**Uploaded field data:** the table updates as REM, TTE, and USCR finish. If only some models have completed, the table will still summarize the completed fits.",
            "",
            "1. Compute ΔWAIC and WAIC weights when multiple models are available;",
            "2. Combine posterior draws of `D_mi²` (animals per square mile) across the completed fits;",
            "3. Report model-specific densities, and when possible also report unweighted and WAIC-weighted summaries plus P(D > 20 DPSM).",
            "",
            "**Simulated data:** if you run USCR, REM, and TTE from the **shared spatial simulator**, the table below will compare and combine those completed simulated fits too.",
            "",
            "Run the models from their tabs first, then check the tables.",
            sep = "\n"
          )),
          
          h4("Simulated data – shared spatial simulator comparison (animals/mi²)"),
          DTOutput("sim_combo_table"),
          p(class = "small", style = "color: var(--muted);",
            "Only fits from the current shared spatial simulator are combined here."),
          downloadButton("dl_sim_uscr_csv", "Download available shared-simulation posterior summaries (CSV)"),
          
          hr(),
          h4("Uploaded field data – model comparison and combined results (animals/mi²)"),
          DTOutput("nps_combo_table"),
          p("Posterior summaries (all monitored parameters, mean, 95% CI) for each completed model:"),
          downloadButton("dl_nps_all_csv", "Download available uploaded-data posterior summaries (CSV)")
        )
  )
)

# -------------------------------------------------------------------
# SERVER
# -------------------------------------------------------------------

server <- function(input, output, session) {
  
  observeEvent(input$deer_tabs, {
    shinyjs::runjs("window.scrollTo(0, 0);")
  }, ignoreInit = TRUE)
  
  # ================================================================
  # SIMULATION: grid + inputs for models
  # ================================================================
  
  sim <- eventReactive(input$run_sim, {
    simulate_camera_counts(
      n_side    = input$n_side,
      spacing_m = input$spacing,
      days      = input$days,
      D_per_km2 = input$Dtrue,
      lambda0   = input$lambda0,
      sigma_m   = input$sigma,
      seed      = input$seed
    )
  }, ignoreInit = TRUE)
  
  sim_data <- reactive({
    req(sim())
    sim_model_inputs(
      sim = sim(),
      detection_radius_m = input$r_m_sim
    )
  })
  
  observeEvent(input$run_sim, {
    sim_id <- paste0("shared_sim_", format(Sys.time(), "%Y%m%d%H%M%S"))
    current_shared_sim_id(sim_id)
    uscr_sim_fit(NULL)
    rem_sim_fit(NULL)
    tte_sim_fit(NULL)
    uscr_sim_dataset_id(NULL)
    rem_sim_dataset_id(NULL)
    tte_sim_dataset_id(NULL)
  }, ignoreInit = TRUE)
  
  output$grid_plot <- renderPlot({
    req(sim())
    traps_df <- as.data.frame(secr::traps(sim()$ch)) %>%
      mutate(camera = row_number())
    mask_df <- as.data.frame(sim()$mask)
    
    ggplot() +
      geom_point(data = mask_df, aes(x = x, y = y), alpha = 0.1, color = redwood_colors[5]) +
      geom_point(data = traps_df, aes(x = x, y = y), size = 3, color = redwood_colors[3]) +
      geom_label(
        data = traps_df,
        aes(x = x, y = y, label = camera),
        nudge_y = 18,
        size = 3.8,
        color = redwood_colors[1],
        fill = "white",
        label.size = 0.2,
        label.padding = grid::unit(0.12, "lines")
      ) +
      coord_equal() +
      theme_minimal() +
      theme(axis.text = element_text(size = 11)) +
      labs(
        x = "m", y = "m",
        title = sprintf(
          "%d cameras, spacing %dm (simulated density: %.1f animals/km²)",
          sim()$truth$n_cams,
          sim()$truth$spacing_m,
          sim()$truth$D_per_km2
        )
      )
  })
  
  output$sim_deer_distribution_plot <- renderPlot({
    req(sim_data())
    d <- sim_data()
    
    # Combine camera locations with counts
    deer_dist <- d$out %>%
      mutate(
        x = utm_e * 1000,  # Convert km to m for consistency with grid plot
        y = utm_n * 1000,
        total_animals = d$camera_counts
      )
    
    td_max <- max(deer_dist$total_animals, na.rm = TRUE)
    size_breaks <- pretty(c(0, max(1, td_max)), n = 4)
    size_breaks <- unique(size_breaks[size_breaks >= 0 & size_breaks <= max(1, td_max)])
    ggplot(deer_dist, aes(x = x, y = y)) +
      geom_point(
        aes(size = total_animals),
        shape = 21,
        fill = redwood_colors[3],
        color = redwood_colors[1],
        stroke = 0.35,
        alpha = 0.75
      ) +
      ggrepel::geom_label_repel(
        aes(label = Site),
        size = 3,
        max.overlaps = Inf,
        fill = "white",
        label.size = 0.2,
        box.padding = 0.25,
        point.padding = 0.2
      ) +
      scale_size_area(
        name   = "Total detections",
        max_size = 14,
        breaks = size_breaks
      ) +
      coord_equal() +
      labs(
        title = "Total Animal Detections by Camera",
        x = "m",
        y = "m"
      ) +
      theme_minimal() +
      theme(legend.position = "right")
  })
  
  output$camera_table <- renderDT({
    req(sim())
    ch <- sim()$ch
    cams <- as.data.frame(secr::traps(ch)) %>%
      mutate(camera_id = paste0("C", row_number()))
    
    counts_mat <- get_counts_matrix(ch)
    totcounts  <- rowSums(counts_mat, na.rm = TRUE)
    
    cams %>%
      mutate(total_counts = totcounts) %>%
      select(camera_id, x, y, total_counts) %>%
      datatable(options = list(pageLength = 8))
  })
  
  # ================================================================
  # CAMERA ARRAY DESIGN
  # ================================================================

  camera_array_result <- reactiveVal(NULL)
  camera_array_log <- reactiveVal("No camera array has been generated yet.")

  observeEvent(input$run_camera_array, {
    req(input$camera_boundary_files)

    log_lines <- character()
    add_log <- function(...) {
      log_lines <<- c(log_lines, paste0(...))
      camera_array_log(paste(log_lines, collapse = "\n"))
    }

    camera_array_result(NULL)
    camera_array_log("Starting camera-array generation...")

    res <- tryCatch(
      {
        withProgress(
          message = "Generating camera array",
          detail = "Reading uploaded spatial layers...",
          value = 0,
          {
            add_log("Reading boundary file(s)...")
            boundary_sf <- read_uploaded_spatial_group(
              input$camera_boundary_files,
              label = "boundary",
              required = TRUE
            )

            setProgress(0.15, detail = "Reading optional layers...")
            roads_sf <- read_uploaded_spatial_group(input$camera_roads_files, "roads layer")
            trails_sf <- read_uploaded_spatial_group(input$camera_trails_files, "trails layer")
            buildings_sf <- read_uploaded_spatial_group(input$camera_buildings_files, "buildings layer")
            parking_sf <- read_uploaded_spatial_group(input$camera_parking_files, "parking layer")
            exclusion_sf <- read_uploaded_spatial_group(input$camera_exclusion_files, "exclusion layer")

            add_log("Boundary loaded successfully.")
            optional_counts <- c(
              roads = if (is.null(roads_sf)) 0 else nrow(roads_sf),
              trails = if (is.null(trails_sf)) 0 else nrow(trails_sf),
              buildings = if (is.null(buildings_sf)) 0 else nrow(buildings_sf),
              parking = if (is.null(parking_sf)) 0 else nrow(parking_sf),
              exclusions = if (is.null(exclusion_sf)) 0 else nrow(exclusion_sf)
            )
            add_log(
              "Optional layer features loaded: ",
              paste(names(optional_counts), optional_counts, sep = "=", collapse = ", ")
            )

            setProgress(0.4, detail = "Building allowable camera area...")
            design_area <- build_camera_design_area(
              boundary_sf = boundary_sf,
              roads = roads_sf,
              trails = trails_sf,
              buildings = buildings_sf,
              parking_lots = parking_sf,
              exclusion_areas = exclusion_sf,
              buffer_dist_m = input$camera_buffer_m
            )

            boundary_area_ha <- as.numeric(sum(sf::st_area(design_area$boundary))) / 10000
            allowed_area_ha <- as.numeric(sum(sf::st_area(design_area$allowed_area))) / 10000
            excluded_area_ha <- if (!is.null(design_area$excluded_area) && nrow(design_area$excluded_area)) {
              as.numeric(sum(sf::st_area(design_area$excluded_area))) / 10000
            } else {
              0
            }

            add_log("Projected design area created.")
            add_log("Boundary area: ", format_num(boundary_area_ha, 1), " ha")
            add_log("Allowable camera area: ", format_num(allowed_area_ha, 1), " ha")
            if (excluded_area_ha > 0) {
              add_log("Excluded / buffered area removed: ", format_num(excluded_area_ha, 1), " ha")
            }

            setProgress(0.7, detail = "Generating candidate camera locations...")
            prefix <- trimws(input$camera_prefix)
            if (!nzchar(prefix)) prefix <- "CAM"

            candidates <- generate_camera_candidates(
              allowed_area = design_area$allowed_area,
              spacing_m = input$camera_spacing_m,
              camera_prefix = prefix
            )
            selected <- select_camera_subset(
              candidates_sf = candidates,
              budget = input$camera_budget,
              n_alternates = input$camera_alternates
            )

            add_log("Suggested camera candidates: ", nrow(candidates))
            add_log(
              "Final cameras retained: ",
              sum(selected$status == "final"),
              if (any(selected$status == "alternate")) {
                paste0(" (plus ", sum(selected$status == "alternate"), " alternates)")
              } else {
                ""
              }
            )

            setProgress(0.9, detail = "Preparing exports and map outputs...")
            candidates_export <- camera_points_for_export(candidates)
            selected_export <- camera_points_for_export(selected)

            list(
              boundary = design_area$boundary,
              allowed_area = design_area$allowed_area,
              excluded_area = design_area$excluded_area,
              candidates = candidates_export,
              selected = selected_export,
              summary = list(
                spacing_m = input$camera_spacing_m,
                buffer_m = input$camera_buffer_m,
                budget = input$camera_budget,
                alternates = input$camera_alternates,
                boundary_area_ha = boundary_area_ha,
                allowed_area_ha = allowed_area_ha,
                excluded_area_ha = excluded_area_ha
              )
            )
          }
        )
      },
      error = function(e) {
        add_log("ERROR: ", e$message)
        showNotification(
          paste("Camera array design failed:", e$message),
          type = "error",
          duration = NULL
        )
        NULL
      }
    )

    if (!is.null(res)) {
      camera_array_result(res)
      add_log("Camera-array design completed successfully.")
      showNotification("Camera array generated successfully.", type = "message", duration = 6)
    }
  })

  output$camera_array_log <- renderText({
    camera_array_log()
  })

  output$camera_array_summary <- renderText({
    res <- camera_array_result()
    if (is.null(res)) {
      return("No camera array has been generated yet.")
    }

    summary <- res$summary
    selected <- res$selected

    paste(
      c(
        paste("Suggested spacing:", summary$spacing_m, "m"),
        paste("Buffer around optional layers:", summary$buffer_m, "m"),
        paste("Boundary area:", format_num(summary$boundary_area_ha, 1), "ha"),
        paste("Allowable area:", format_num(summary$allowed_area_ha, 1), "ha"),
        paste("Candidates generated:", nrow(res$candidates)),
        paste("Final cameras:", sum(selected$status == "final")),
        paste("Alternate cameras:", sum(selected$status == "alternate"))
      ),
      collapse = "\n"
    )
  })

  output$camera_array_map <- renderLeaflet({
    res <- camera_array_result()
    req(res)

    boundary_ll <- sf::st_transform(res$boundary, 4326)
    allowed_ll <- sf::st_transform(res$allowed_area, 4326)
    excluded_ll <- if (!is.null(res$excluded_area) && nrow(res$excluded_area)) {
      sf::st_transform(res$excluded_area, 4326)
    } else {
      NULL
    }
    candidates_ll <- sf::st_transform(res$candidates, 4326)
    selected_ll <- sf::st_transform(res$selected, 4326)

    m <- leaflet() |>
      addTiles() |>
      addPolygons(
        data = boundary_ll,
        color = redwood_colors[1],
        weight = 2,
        fill = FALSE,
        group = "Boundary"
      ) |>
      addPolygons(
        data = allowed_ll,
        color = redwood_colors[3],
        weight = 1,
        fillColor = redwood_colors[4],
        fillOpacity = 0.2,
        group = "Allowed area"
      ) |>
      addCircleMarkers(
        data = candidates_ll,
        radius = 3,
        color = "#666666",
        stroke = FALSE,
        fillOpacity = 0.7,
        group = "Candidates",
        popup = ~paste0("<strong>", camera_id, "</strong>")
      )

    if (!is.null(excluded_ll)) {
      m <- m |>
        addPolygons(
          data = excluded_ll,
          color = "#8b5a2b",
          weight = 1,
          fillColor = "#c59d74",
          fillOpacity = 0.3,
          group = "Excluded area"
        )
    }

    final_pts <- selected_ll[selected_ll$status == "final", , drop = FALSE]
    alt_pts <- selected_ll[selected_ll$status == "alternate", , drop = FALSE]

    if (nrow(final_pts)) {
      m <- m |>
        addCircleMarkers(
          data = final_pts,
          radius = 6,
          color = redwood_colors[3],
          fillColor = redwood_colors[3],
          fillOpacity = 1,
          weight = 1,
          group = "Final cameras",
          popup = ~paste0("<strong>", camera_id, "</strong><br>Status: final")
        )
    }

    if (nrow(alt_pts)) {
      m <- m |>
        addCircleMarkers(
          data = alt_pts,
          radius = 6,
          color = "#d97904",
          fillColor = "#f4b266",
          fillOpacity = 1,
          weight = 1,
          group = "Alternates",
          popup = ~paste0("<strong>", camera_id, "</strong><br>Status: alternate")
        )
    }

    m |>
      addLayersControl(
        overlayGroups = c("Boundary", "Allowed area", "Excluded area", "Candidates", "Final cameras", "Alternates"),
        options = layersControlOptions(collapsed = FALSE)
      )
  })

  output$camera_array_table <- DT::renderDT({
    res <- camera_array_result()
    req(res)

    export_df <- res$selected |>
      sf::st_drop_geometry() |>
      dplyr::select(camera_id, status, longitude, latitude, utm_easting_m, utm_northing_m)

    DT::datatable(export_df, options = list(pageLength = 10))
  })

  output$download_camera_array_csv <- downloadHandler(
    filename = function() {
      paste0("camera_array_", Sys.Date(), ".csv")
    },
    content = function(file) {
      res <- camera_array_result()
      req(res)
      readr::write_csv(sf::st_drop_geometry(res$selected), file)
    }
  )

  output$download_camera_array_gpkg <- downloadHandler(
    filename = function() {
      paste0("camera_array_", Sys.Date(), ".gpkg")
    },
    content = function(file) {
      res <- camera_array_result()
      req(res)
      sf::st_write(res$selected, file, layer = "camera_points", delete_dsn = TRUE, quiet = TRUE)
    }
  )

  output$download_camera_area_gpkg <- downloadHandler(
    filename = function() {
      paste0("camera_design_area_", Sys.Date(), ".gpkg")
    },
    content = function(file) {
      res <- camera_array_result()
      req(res)
      sf::st_write(res$boundary, file, layer = "boundary", delete_dsn = TRUE, quiet = TRUE)
      sf::st_write(res$allowed_area, file, layer = "allowed_area", append = TRUE, quiet = TRUE)
      if (!is.null(res$excluded_area) && nrow(res$excluded_area)) {
        sf::st_write(res$excluded_area, file, layer = "excluded_area", append = TRUE, quiet = TRUE)
      }
    }
  )
  
  # ================================================================
  # NPS QC: deployment + images
  # ================================================================
  
  deployment_checked <- reactiveVal(NULL)
  deployment_issues  <- reactiveVal(NULL)
  
  observeEvent(input$deployment_csv, {
    req(input$deployment_csv)
    
    deployment_raw <- readr::read_csv(
      input$deployment_csv$datapath,
      show_col_types = FALSE
    )
    
    deployment_raw <- clean_deployment_import(deployment_raw)
    
    msgs  <- character()
    warns <- character()
    
    res <- tryCatch(
      {
        withCallingHandlers(
          {
            check_deployments(deployment_raw)
          },
          message = function(m) {
            msgs  <<- c(msgs, m$message)
            invokeRestart("muffleMessage")
          },
          warning = function(w) {
            warns <<- c(warns, w$message)
            invokeRestart("muffleWarning")
          }
        )
      },
      error = function(e) {
        deployment_checked(NULL)
        deployment_issues(list(
          messages = msgs,
          warnings = c(warns, paste0("ERROR: ", e$message))
        ))
        return(NULL)
      }
    )
    
    if (!is.null(res)) {
      deployment_checked(res)
      deployment_issues(list(
        messages = msgs,
        warnings = warns
      ))
    }
  })
  
  output$deployment_check_log <- renderText({
    issues <- deployment_issues()
    req(issues)
    
    lines <- c()
    if (length(issues$messages)) {
      lines <- c(lines, "Messages:",
                 paste0("  - ", trimws(issues$messages)), "")
    }
    if (length(issues$warnings)) {
      lines <- c(lines, "Warnings / issues:",
                 paste0("  - ", trimws(issues$warnings)))
    }
    if (!length(lines)) lines <- "No issues reported."
    paste(lines, collapse = "\n")
  })
  
  output$deployment_preview <- DT::renderDT({
    req(deployment_checked())
    DT::datatable(deployment_checked(), options = list(pageLength = 10))
  })
  
  output$download_deployment_checked <- downloadHandler(
    filename = function() {
      base <- tools::file_path_sans_ext(input$deployment_csv$name)
      paste0(base, "_CHECKED.csv")
    },
    content = function(file) {
      req(deployment_checked())
      readr::write_csv(deployment_checked(), file)
    }
  )
  
  # ---- Images QC ----
  
  images_checked <- reactiveVal(NULL)
  images_issues  <- reactiveVal(NULL)
  
  observeEvent(input$images_csv, {
    req(input$images_csv)
    
    if (is.null(deployment_checked())) {
      showNotification(
        "Please upload and check a deployment CSV before checking images.",
        type = "error"
      )
      return(NULL)
    }
    
    images_raw <- readr::read_csv(
      input$images_csv$datapath,
      show_col_types = FALSE,
      col_types = readr::cols(
        `Sighting Count` = readr::col_character(),
        Species = readr::col_character()
      )
    )
    
    images_raw <- clean_images_import(images_raw)
    
    # Re-run deployment QC with images (malfunction date rules use site overlap)
    dep_msgs  <- character()
    dep_warns <- character()
    dep_res <- tryCatch(
      {
        withCallingHandlers(
          {
            check_deployments(deployment_checked(), images_raw)
          },
          message = function(m) {
            dep_msgs <<- c(dep_msgs, m$message)
            invokeRestart("muffleMessage")
          },
          warning = function(w) {
            dep_warns <<- c(dep_warns, w$message)
            invokeRestart("muffleWarning")
          }
        )
      },
      error = function(e) {
        showNotification(
          paste("Deployment re-check failed:", e$message),
          type = "error"
        )
        images_checked(NULL)
        images_issues(list(
          messages = dep_msgs,
          warnings = c(dep_warns, paste0("ERROR: ", e$message))
        ))
        return(NULL)
      }
    )
    
    if (is.null(dep_res)) return(NULL)
    
    deployment_checked(dep_res)
    deployment_issues(list(messages = dep_msgs, warnings = dep_warns))
    
    msgs  <- dep_msgs
    warns <- dep_warns
    
    checked <- NULL
    
    res <- tryCatch(
      {
        withCallingHandlers(
          {
            checked <- check_images(
              images      = images_raw,
              deployments = dep_res
            )
          },
          message = function(m) {
            msgs  <<- c(msgs, m$message)
            invokeRestart("muffleMessage")
          },
          warning = function(w) {
            warns <<- c(warns, w$message)
            invokeRestart("muffleWarning")
          }
        )
        checked
      },
      error = function(e) {
        images_checked(NULL)
        images_issues(list(
          messages = msgs,
          warnings = c(warns, paste0("ERROR: ", e$message))
        ))
        return(NULL)
      }
    )
    
    if (is.null(res)) return(NULL)
    
    if (isTRUE(input$apply_56day_trim)) {
      trim_msgs <- character()
      
      trimmed <- withCallingHandlers(
        {
          trim_images_56days(res)
        },
        message = function(m) {
          trim_msgs <<- c(trim_msgs, m$message)
          invokeRestart("muffleMessage")
        },
        warning = function(w) {
          warns <<- c(warns, w$message)
          invokeRestart("muffleWarning")
        }
      )
      
      images_checked(trimmed)
      images_issues(list(
        messages = c(msgs, trim_msgs),
        warnings = warns
      ))
    } else {
      images_checked(res)
      images_issues(list(
        messages = c(
          msgs,
          "Skipped the optional 56-day trim; using all validated images that fall within the deployment windows."
        ),
        warnings = warns
      ))
    }
  })
  
  output$images_check_log <- renderText({
    issues <- images_issues()
    req(issues)
    
    lines <- c()
    if (length(issues$messages)) {
      lines <- c(lines, "Messages:",
                 paste0("  - ", trimws(issues$messages)), "")
    }
    if (length(issues$warnings)) {
      lines <- c(lines, "Warnings / issues:",
                 paste0("  - ", trimws(issues$warnings)))
    }
    if (!length(lines)) lines <- "No issues reported."
    paste(lines, collapse = "\n")
  })
  
  output$images_preview <- DT::renderDT({
    req(images_checked())
    DT::datatable(images_checked(), options = list(pageLength = 10))
  })
  
  output$download_images_checked <- downloadHandler(
    filename = function() {
      base <- tools::file_path_sans_ext(input$images_csv$name)
      paste0(base, "_CHECKED_TRIM56.csv")
    },
    content = function(file) {
      req(images_checked())
      readr::write_csv(images_checked(), file)
    }
  )
  
  # ================================================================
  # DATA SUMMARY (NPS)
  # ================================================================
  
  deploy_summary <- reactive({
    req(deployment_checked())
    summarize_deployments(deployment_checked())
  })
  
  images_standardized <- reactive({
    req(images_checked())
    standardize_deer_species(images_checked())
  })
  
  image_summary <- reactive({
    summarize_images_by_species(images_standardized())
  })
  
  available_species <- reactive({
    req(images_standardized())
    species <- sort(unique(images_standardized()$Species))
    species <- species[!is.na(species) & nzchar(species)]
    species
  })
  
  observe({
    species <- available_species()
    current <- isolate(input$summary_species)
    selected <- if (!is.null(current) && current %in% species) {
      current
    } else if ("Deer" %in% species) {
      "Deer"
    } else if (length(species) > 0) {
      species[[1]]
    } else {
      character(0)
    }
    
    updateSelectInput(
      session,
      "summary_species",
      choices = species,
      selected = selected
    )
  })
  
  species_objects <- reactive({
    req(images_standardized(), input$summary_species)
    
    list(
      species_summary = species_summary_per_site(images_standardized(), input$summary_species),
      species_counts  = species_counts_per_camera(images_standardized(), input$summary_species),
      daily_species   = species_daily_detections(images_standardized(), input$summary_species)
    )
  })
  
  output$deploy_summary_table <- DT::renderDT({
    req(deploy_summary())
    DT::datatable(deploy_summary(), options = list(dom = "t", paging = FALSE))
  })
  
  output$camera_map <- leaflet::renderLeaflet({
    req(deployment_checked())
    dep <- deployment_checked()
    
    validate(
      need(all(c("Longitude", "Latitude", "Site Name") %in% names(dep)),
           "Deployment file must contain Longitude, Latitude, and Site Name.")
    )
    
    icon.fa <- leaflet::makeAwesomeIcon(
      icon = "camera", markerColor = "green",
      library = "fa", iconColor = "black"
    )
    
    leaflet::leaflet(dep) %>%
      leaflet::addTiles() %>%
      leaflet::addAwesomeMarkers(
        lng   = ~Longitude,
        lat   = ~Latitude,
        icon  = icon.fa,
        popup = ~as.character(`Site Name`),
        label = ~as.character(`Site Name`)
      )
  })
  
  output$image_summary_table <- DT::renderDT({
    req(image_summary())
    DT::datatable(image_summary(), options = list(pageLength = 10))
  })
  
  output$species_bar_plot <- renderPlot({
    df <- image_summary()
    validate(
      need(nrow(df) > 0, "No images available to summarize.")
    )
    
    ggplot(df, aes(x = reorder(Species, -total_detections),
                   y = total_detections)) +
      geom_bar(stat = "identity", fill = redwood_colors[3]) +
      labs(x = "Species", y = "Total detections") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  })
  
  output$deer_summary_table <- DT::renderDT({
    species_name <- req(input$summary_species)
    species_df <- species_objects()$species_summary
    validate(
      need(nrow(species_df) > 0, paste0("No '", species_name, "' images found in dataset."))
    )
    DT::datatable(species_df, options = list(pageLength = 10))
  })
  
  output$deer_bubble_plot <- renderPlot({
    species_name <- req(input$summary_species)
    dc <- species_objects()$species_counts
    validate(
      need(nrow(dc) > 0, paste0("No '", species_name, "' images found in dataset."))
    )
    
    td_max <- max(dc$total_detections, na.rm = TRUE)
    ggplot(dc, aes(x = Longitude, y = Latitude)) +
      geom_point(aes(size = total_detections),
                 color = redwood_colors[3], alpha = 0.7) +
      ggrepel::geom_text_repel(
        aes(label = `Site Name`),
        size = 3,
        max.overlaps = Inf
      ) +
      scale_size_continuous(
        name   = paste("Total", species_name),
        range  = c(2, 14),
        limits = c(0, NA),
        breaks = pretty(c(0, max(1, td_max)))
      ) +
      coord_fixed() +
      labs(
        title = paste("Total", species_name, "Detections by Camera"),
        x = "Longitude",
        y = "Latitude"
      ) +
      theme_minimal() +
      theme(legend.position = "right")
  })
  
  output$deer_daily_plot <- renderPlot({
    species_name <- req(input$summary_species)
    dd <- species_objects()$daily_species
    validate(
      need(nrow(dd) > 0, paste0("No '", species_name, "' images found in dataset."))
    )
    
    ggplot(dd, aes(x = Date, y = `Site Name`, size = detections)) +
      geom_point(
        color = redwood_colors[3], alpha = 0.7,
        position = position_jitter(width = 0, height = 0.1)
      ) +
      scale_size_continuous(name = paste(species_name, "count"), range = c(2, 8)) +
      labs(
        title = paste("Daily", species_name, "Detections by Camera"),
        x = "Date",
        y = "Site"
      ) +
      theme_minimal() +
      theme(
        axis.text.y     = element_text(size = 8),
        legend.position = "right"
      )
  })
  
  # ================================================================
  # MODEL FIT OBJECTS (reactiveVal) – SIM + NPS PER MODEL
  # ================================================================
  
  uscr_sim_fit <- reactiveVal(NULL)
  rem_sim_fit  <- reactiveVal(NULL)
  tte_sim_fit  <- reactiveVal(NULL)
  current_shared_sim_id <- reactiveVal(NULL)
  uscr_sim_dataset_id <- reactiveVal(NULL)
  rem_sim_dataset_id  <- reactiveVal(NULL)
  tte_sim_dataset_id  <- reactiveVal(NULL)

  uscr_nps_fit <- reactiveVal(NULL)
  rem_nps_fit  <- reactiveVal(NULL)
  tte_nps_fit  <- reactiveVal(NULL)
  
  make_model_debug <- function(model, data_source, supported_sources, preprocess_note) {
    list(
      model = model,
      data_source = data_source,
      supported_sources = supported_sources,
      preprocess_note = preprocess_note,
      status = "idle",
      stage = "Not started",
      started_at = NULL,
      finished_at = NULL,
      guidance = "No run has been started yet.",
      raw_error = NULL,
      context = character(),
      history = character()
    )
  }
  
  update_model_debug <- function(rv,
                                 status = NULL,
                                 stage = NULL,
                                 started_at = NULL,
                                 finished_at = NULL,
                                 guidance = NULL,
                                 raw_error = NULL,
                                 context = NULL,
                                 log_entry = NULL) {
    state <- rv()
    if (is.null(state)) state <- list()
    
    if (!is.null(status)) state$status <- status
    if (!is.null(stage)) state$stage <- stage
    if (!is.null(started_at)) state$started_at <- started_at
    if (!is.null(finished_at)) state$finished_at <- finished_at
    if (!is.null(guidance)) state$guidance <- guidance
    if (!is.null(raw_error)) state$raw_error <- raw_error
    if (!is.null(context)) state$context <- context
    if (!is.null(log_entry)) {
      state$history <- c(
        state$history,
        paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", log_entry)
      )
    }
    
    rv(state)
    invisible(state)
  }
  
  format_num <- function(x, digits = 2) {
    if (length(x) != 1 || is.na(x) || !is.finite(x)) return("NA")
    format(round(x, digits), nsmall = digits, trim = TRUE)
  }
  
  summarize_uscr_context <- function(d, source_label = "NPS", sim_truth = NULL) {
    det_dist <- if ("Detection Distance" %in% names(d$out)) {
      paste0(
        format_num(min(d$out$`Detection Distance`, na.rm = TRUE), 1),
        " to ",
        format_num(max(d$out$`Detection Distance`, na.rm = TRUE), 1),
        " m"
      )
    } else {
      "Not available"
    }
    
    lines <- c(
      paste("Cameras:", nrow(d$out)),
      paste("Total animal events:", sum(d$camera_counts, na.rm = TRUE)),
      paste("Total camera-days:", format_num(sum(d$camera_days, na.rm = TRUE), 1)),
      paste("Detection distance range:", det_dist),
      if (source_label == "NPS") {
        "Supported sources here: uploaded NPS data and simulated grid data."
      } else {
        "Supported sources here: simulated grid data and uploaded NPS data."
      },
      "Preprocessing: total animal events per camera, camera-days, and buffered spatial state space."
    )
    
    if (!is.null(sim_truth)) {
      lines <- c(lines, paste("True simulated density:", format_num(sim_truth$D_per_km2, 1), "animals/km^2"))
    }
    
    lines
  }
  
  summarize_rem_context <- function(d) {
    c(
      paste("Cameras:", nrow(d$out)),
      paste("Total animal events:", sum(d$camera_counts, na.rm = TRUE)),
      paste("Total camera-days:", format_num(sum(d$camera_days, na.rm = TRUE), 1)),
      paste(
        "Detection distance range:",
        paste0(
          format_num(min(d$out$`Detection Distance`, na.rm = TRUE), 1),
          " to ",
          format_num(max(d$out$`Detection Distance`, na.rm = TRUE), 1),
          " m"
        )
      ),
      "Supported sources here: uploaded NPS data only.",
      "Preprocessing: REM uses total animal events per camera and camera-days."
    )
  }
  
  summarize_tte_context <- function(d) {
    c(
      paste("Cameras:", nrow(d$out)),
      paste("Total animal events passed to TTE:", sum(d$camera_counts, na.rm = TRUE)),
      paste("Total camera-days:", format_num(sum(d$camera_days, na.rm = TRUE), 1)),
      paste(
        "Detection distance range:",
        paste0(
          format_num(min(d$out$`Detection Distance`, na.rm = TRUE), 1),
          " to ",
          format_num(max(d$out$`Detection Distance`, na.rm = TRUE), 1),
          " m"
        )
      ),
      "Supported sources here: uploaded NPS data only.",
      "Current app preprocessing: TTE is receiving total animal events per camera and camera-days."
    )
  }
  
  summarize_shared_sim_context <- function(d, model_label, sim_truth = NULL) {
    c(
      paste("Shared simulator model fit:", model_label),
      summarize_uscr_context(d, source_label = "Simulated", sim_truth = sim_truth),
      "This fit comes from the shared spatial simulator and is eligible for simulated Compare & combine."
    )
  }
  
  friendly_model_error <- function(model, data_source, error_text) {
    tips <- c(
      paste0(model, " on ", data_source, " did not finish."),
      "Read the raw error below, then try the model-specific guidance."
    )
    
    if (grepl("unused arguments", error_text, ignore.case = TRUE)) {
      tips <- c(tips, "This usually means an older function definition is still loaded. Restart R and launch this app from the current project folder.")
    }
    if (grepl("zero camera-days", error_text, ignore.case = TRUE)) {
      tips <- c(tips, "One or more cameras ended up with zero effort days. Check Start Date/End Date, malfunction handling, and the optional 56-day trimming choice.")
    }
    if (grepl("camera_counts and camera_days must have length nrow\\(out\\)", error_text)) {
      tips <- c(tips, "The camera table, counts, and camera-days vectors got out of sync. Re-upload the data and check the NPS model-input step.")
    }
    if (grepl("cannot allocate|vector memory|std::bad_alloc", error_text, ignore.case = TRUE)) {
      tips <- c(tips, "This looks like a memory issue. Try fewer chains, fewer iterations, or a smaller USCR M for a test run.")
    }
    if (grepl("Package '.*' is required|there is no package called", error_text, ignore.case = TRUE)) {
      tips <- c(tips, "A required R package is missing in this session. Install the package named in the raw error, restart R, and rerun.")
    }
    if (grepl("NA|NaN|Inf", error_text)) {
      tips <- c(tips, "Missing or invalid numeric values reached the model. Check detection distances, camera-days, timestamps, and coordinates in the upload summary tabs.")
    }
    
    if (identical(model, "USCR")) {
      tips <- c(tips, "USCR does an adaptive tuning run before the final run. If it is too slow for testing, reduce chains, iterations, or M.")
    }
    if (identical(model, "REM")) {
      tips <- c(tips, "REM expects per-camera animal events, camera-days, and detection distances.")
      tips <- c(tips, "REM now doubles iterations automatically until the maximum Rhat is acceptable, so longer runs can be expected when convergence is slow.")
    }
    if (identical(model, "TTE")) {
      tips <- c(tips, "TTE expects per-camera animal events, camera-days, and detection distances. The current build passes total animal events per camera to TTE.")
      tips <- c(tips, "TTE now doubles iterations automatically until the maximum Rhat is acceptable, so longer runs can be expected when convergence is slow.")
    }
    
    paste(tips, collapse = "\n")
  }

  make_encounter_status_callback <- function(rv, model_label, set_progress = NULL) {
    function(stage, detail = NULL, value = NULL) {
      stage_label <- switch(
        stage,
        setup = "Preparing model code",
        tuning = "Adaptive tuning / convergence checks",
        final_run = "Final MCMC run",
        stage
      )

      if (is.function(set_progress)) {
        progress_detail <- detail %||% stage_label
        if (!is.null(value)) {
          set_progress(value = value, detail = progress_detail)
        } else {
          set_progress(detail = progress_detail)
        }
      }

      update_model_debug(
        rv,
        stage = stage_label,
        log_entry = detail %||% paste(model_label, "stage:", stage_label)
      )
    }
  }
  
  format_model_debug <- function(state, fit = NULL) {
    lines <- c(
      paste("Model:", state$model),
      paste("Data source:", state$data_source),
      paste("Supported in app:", state$supported_sources),
      paste("Status:", state$status),
      paste("Stage:", state$stage)
    )
    
    if (!is.null(state$started_at)) {
      lines <- c(lines, paste("Started:", format(state$started_at, "%Y-%m-%d %H:%M:%S")))
    }
    if (!is.null(state$finished_at)) {
      lines <- c(lines, paste("Finished:", format(state$finished_at, "%Y-%m-%d %H:%M:%S")))
    }
    if (!is.null(state$started_at)) {
      end_time <- if (identical(state$status, "running")) Sys.time() else state$finished_at %||% Sys.time()
      elapsed_min <- as.numeric(difftime(end_time, state$started_at, units = "mins"))
      lines <- c(lines, paste("Elapsed:", format_num(elapsed_min, 1), "minutes"))
    }
    
    lines <- c(lines, "", "Guidance:", state$guidance)
    
    if (length(state$context)) {
      lines <- c(lines, "", "Input summary:", paste0("- ", state$context))
    }
    
    if (!is.null(fit) && !is.null(fit$settings)) {
      setting_lines <- c(
        paste0("- Final iterations: ", fit$settings$iter),
        paste0("- Burn-in: ", fit$settings$burnin),
        paste0("- Thin: ", fit$settings$thin),
        paste0("- Chains: ", fit$settings$n_chains)
      )
      if (!is.null(fit$settings$M)) {
        setting_lines <- c(paste0("- Final M: ", fit$settings$M), setting_lines)
      }
      lines <- c(lines, "", "Run settings used:", setting_lines)
      if (!is.null(fit$final_rhat_max) && is.finite(fit$final_rhat_max)) {
        lines <- c(lines, paste0("- Final max Rhat: ", format_num(fit$final_rhat_max, 3)))
      }
    }
    
    if (!is.null(fit) && !is.null(fit$tuning_history) && nrow(fit$tuning_history) > 0) {
      tuning_df <- fit$tuning_history
      tuning_lines <- apply(
        tuning_df,
        1,
        function(row) {
          parts <- c(
            paste0("round ", row[["round"]]),
            paste0("iter=", row[["niter"]]),
            paste0("thin=", row[["thin"]]),
            paste0("chains=", row[["n_chains"]]),
            paste0("rhat=", format_num(as.numeric(row[["rhat_max"]]), 3))
          )
          if ("M" %in% names(tuning_df)) {
            parts <- append(parts, paste0("M=", row[["M"]]), after = 1)
          }
          if ("M_too_small" %in% names(tuning_df)) {
            parts <- c(parts, paste0("M_too_small=", row[["M_too_small"]]))
          }
          paste(parts, collapse = ", ")
        }
      )
      lines <- c(lines, "", paste0(state$model, " tuning history:"), paste0("- ", tuning_lines))
    }
    
    if (!is.null(state$raw_error)) {
      lines <- c(lines, "", "Raw error:", state$raw_error)
    }
    
    if (length(state$history)) {
      lines <- c(lines, "", "Recent log:", paste0("- ", tail(state$history, 8)))
    }
    
    paste(lines, collapse = "\n")
  }
  
  `%||%` <- function(x, y) {
    if (is.null(x)) y else x
  }
  
  uscr_sim_debug <- reactiveVal(
    make_model_debug("USCR", "Simulated grid", "Simulated grid; uploaded NPS data", "Total animal events per camera, camera-days, and buffered state space.")
  )
  rem_sim_debug <- reactiveVal(
    make_model_debug("REM", "Shared spatial simulator", "Shared spatial simulator", "REM uses total animal events per camera and camera-days derived from the shared spatial simulator.")
  )
  tte_sim_debug <- reactiveVal(
    make_model_debug("TTE", "Shared spatial simulator", "Shared spatial simulator", "TTE uses total animal events per camera and camera-days derived from the shared spatial simulator.")
  )
  uscr_nps_debug <- reactiveVal(
    make_model_debug("USCR", "Uploaded NPS data", "Uploaded NPS data; simulated grid", "Total animal events per camera, camera-days, and buffered state space.")
  )
  rem_nps_debug <- reactiveVal(
    make_model_debug("REM", "Uploaded NPS data", "Uploaded NPS data only", "REM uses total animal events per camera and camera-days.")
  )
  tte_nps_debug <- reactiveVal(
    make_model_debug("TTE", "Uploaded NPS data", "Uploaded NPS data only", "Current build passes total animal events per camera and camera-days.")
  )
  
  # Status tracking for model runs
  uscr_sim_running <- reactiveVal(FALSE)
  rem_sim_running <- reactiveVal(FALSE)
  tte_sim_running <- reactiveVal(FALSE)
  uscr_nps_running <- reactiveVal(FALSE)
  rem_nps_running <- reactiveVal(FALSE)
  tte_nps_running <- reactiveVal(FALSE)
  rem_nps_background <- reactiveVal(FALSE)
  
  # Stop flags for interrupting model runs
  stop_uscr_sim <- reactiveVal(FALSE)
  stop_uscr_nps <- reactiveVal(FALSE)
  stop_rem_nps <- reactiveVal(FALSE)
  stop_tte_nps <- reactiveVal(FALSE)
  
  nps_model_inputs <- reactive({
    req(deployment_checked(), images_checked())
    build_nps_model_inputs(deployment_checked(), images_checked())
  })
  
  # --- Stop handlers ---
  
  # Use a global environment to store stop flags that can be checked during execution
  stop_flags_env <- new.env()
  
  observeEvent(input$stop_uscr_sim, {
    stop_uscr_sim(TRUE)
    stop_flags_env$stop_uscr_sim <- TRUE
    showNotification("Stopping USCR (sim) run...", 
                     type = "warning", duration = 3)
  })
  
  observeEvent(input$stop_uscr_nps, {
    stop_uscr_nps(TRUE)
    stop_flags_env$stop_uscr_nps <- TRUE
    showNotification("Stopping USCR (NPS) run...", 
                     type = "warning", duration = 3)
  })
  
  observeEvent(input$stop_rem_nps, {
    if (isTRUE(rem_nps_background())) {
      update_model_debug(
        rem_nps_debug,
        status = "running",
        stage = "Background run in progress",
        guidance = "This REM run is already running in a background worker. It will keep running until it finishes or fails; server-side cancellation is not available yet for this path.",
        log_entry = "Stop requested, but background REM cancellation is not available."
      )
      showNotification(
        "REM background jobs cannot be cancelled after launch yet. This run will keep going in the background.",
        type = "warning",
        duration = 8
      )
      return(NULL)
    }
    stop_rem_nps(TRUE)
    stop_flags_env$stop_rem_nps <- TRUE
    showNotification("Stopping REM (NPS) run...", 
                     type = "warning", duration = 3)
  })
  
  observeEvent(input$stop_tte_nps, {
    stop_tte_nps(TRUE)
    stop_flags_env$stop_tte_nps <- TRUE
    showNotification("Stopping TTE (NPS) run...", 
                     type = "warning", duration = 3)
  })
  
  # Helper function to check stop flag and throw error if set
  check_stop_flag <- function(flag_name) {
    if (exists(flag_name, envir = stop_flags_env) && 
        isTRUE(get(flag_name, envir = stop_flags_env))) {
      stop("Model execution stopped by user", call. = FALSE)
    }
  }

  app_n_chains <- function() {
    max(3L, as.integer(input$n_chains))
  }
  
  uscr_run_args <- function(waic = TRUE, status_callback = NULL) {
    args <- list(
      iter = input$iter_uscr,
      burnin = input$burnin_uscr,
      thin = input$thin_uscr,
      n_chains = app_n_chains(),
      M = input$M_uscr,
      log_sigma_mean = input$log_sigma_mean,
      log_sigma_sd = input$log_sigma_sd,
      log_lam0_sd = input$log_lam0_sd,
      sd_eps_shape = input$sd_eps_shape,
      sd_eps_rate = input$sd_eps_rate,
      adaptive = TRUE,
      compute_WAIC = waic,
      diagnostic_mode = FALSE,
      tuning_n_chains = min(2L, app_n_chains()),
      parallel_chains = isTRUE(app_n_chains() > 1),
      status_callback = status_callback,
      verbose = FALSE
    )
    
    # Keep the app compatible with whichever run_USCR() definition is
    # currently loaded, including older variants without these extras.
    valid_names <- names(formals(run_USCR_app))
    args[names(args) %in% valid_names]
  }

  # --- USCR: simulated ---
  
  observeEvent(input$run_uscr_sim, {
    req(sim_data())
    d <- sim_data()
    uscr_sim_debug(make_model_debug(
      "USCR",
      "Simulated grid",
      "Simulated grid; uploaded NPS data",
      "Total animal events per camera, camera-days, and buffered state space."
    ))
    update_model_debug(
      uscr_sim_debug,
      status = "running",
      stage = "Preflight checks",
      started_at = Sys.time(),
      guidance = "USCR supports simulated and uploaded NPS data. Watch this panel for tuning rounds and the final run stage.",
      raw_error = NULL,
      context = summarize_uscr_context(d, source_label = "Simulated", sim_truth = sim()$truth),
      log_entry = "Simulated USCR run requested."
    )
    
    # Reset stop flag
    stop_uscr_sim(FALSE)
    stop_flags_env$stop_uscr_sim <- FALSE
    uscr_sim_running(TRUE)
    uscr_sim_fit(NULL)  # Clear previous results
    uscr_sim_dataset_id(current_shared_sim_id())
    
    showNotification("Running USCR on simulated data.", type = "message", duration = 6)
    
    fit <- tryCatch(
      {
        # Check if stopped before starting
        check_stop_flag("stop_uscr_sim")
        
        withProgress(
          message = "Running USCR model",
          detail = "Checking inputs and preparing the USCR run...",
          value = 0,
          {
            setProgress(0.05, detail = "Checking inputs and preparing the USCR run...")
            update_model_debug(
              uscr_sim_debug,
              stage = "Preparing run",
              log_entry = "Input checks passed. Starting USCR setup."
            )
            
            # Check stop flag again before running
            check_stop_flag("stop_uscr_sim")
            
            status_callback <- function(stage, detail = NULL, value = NULL) {
              stage_label <- switch(
                stage,
                setup = "Preparing state space and model code",
                tuning = "Adaptive tuning",
                final_run = "Final MCMC run",
                stage
              )
              update_model_debug(
                uscr_sim_debug,
                status = "running",
                stage = stage_label,
                log_entry = detail %||% paste("USCR stage:", stage_label)
              )
              if (!is.null(value)) {
                setProgress(value = value, detail = detail)
              } else if (!is.null(detail)) {
                setProgress(detail = detail)
              }
            }
            
            fit_result <- do.call(
              run_USCR_app,
              c(
                list(
                  out = d$out,
                  camera_counts = d$camera_counts,
                  camera_days = d$camera_days
                ),
                uscr_run_args(waic = FALSE, status_callback = status_callback)
              )
            )
            setProgress(1.0, detail = "Complete!")
            fit_result
          }
        )
      },
      error = function(e) {
        if (stop_uscr_sim() || grepl("stopped by user", e$message, ignore.case = TRUE)) {
          update_model_debug(
            uscr_sim_debug,
            status = "stopped",
            stage = "Stopped by user",
            finished_at = Sys.time(),
            guidance = "The simulated USCR run was stopped manually before completion.",
            raw_error = e$message,
            log_entry = "Simulated USCR run stopped by user."
          )
          showNotification("USCR (sim) run stopped by user.", type = "warning")
        } else {
          update_model_debug(
            uscr_sim_debug,
            status = "error",
            stage = "Failed",
            finished_at = Sys.time(),
            guidance = friendly_model_error("USCR", "simulated data", e$message),
            raw_error = e$message,
            log_entry = paste("Simulated USCR failed:", e$message)
          )
          showNotification(
            paste("USCR (sim) failed:", e$message),
            type = "error", duration = NULL
          )
        }
        # If nimble compilation failed, errors can be inspected via nimble::printErrors() in the console
        return(NULL)
      },
      finally = {
        uscr_sim_running(FALSE)
        stop_uscr_sim(FALSE)  # Reset stop flag
        stop_flags_env$stop_uscr_sim <- FALSE
      }
    )
    uscr_sim_fit(fit)
    if (!is.null(fit) && !stop_uscr_sim()) {
      update_model_debug(
        uscr_sim_debug,
        status = "success",
        stage = "Complete",
        finished_at = Sys.time(),
        guidance = "Simulated USCR completed successfully. Review the summary above and the tuning history below.",
        log_entry = "Simulated USCR run completed."
      )
      showNotification("USCR (sim) complete!", type = "message")
    }
  })
  
  # --- REM: simulated ---
  
  observeEvent(input$run_rem_sim, {
    shared_truth <- sim()
    shared_d <- if (!is.null(shared_truth)) sim_data() else NULL

    if (is.null(shared_d)) {
      showNotification(
        "Run the shared spatial simulator in the Simulate data tab first.",
        type = "error",
        duration = 6
      )
      return(NULL)
    }
    
    rem_sim_dataset_id(current_shared_sim_id())
    rem_sim_debug(make_model_debug(
      "REM",
      "Shared spatial simulator",
      "Shared spatial simulator",
      "REM uses total animal events per camera and camera-days derived from the shared spatial simulator."
    ))
    update_model_debug(
      rem_sim_debug,
      status = "running",
      stage = "Preflight checks",
      started_at = Sys.time(),
      guidance = "REM is running on the shared spatial simulator so it can be compared and combined with simulated USCR and TTE.",
      raw_error = NULL,
      context = summarize_shared_sim_context(shared_d, "REM", sim_truth = shared_truth$truth),
      log_entry = "Simulated REM shared-spatial run requested."
    )
    
    rem_sim_running(TRUE)
    rem_sim_fit(NULL)
    
    showNotification("Running REM on shared simulated spatial data...", type = "message", duration = 6)
    
    fit <- tryCatch(
      {
        withProgress(
          message = "Running REM model",
          detail = "Preparing REM inputs from the shared spatial simulator...",
          value = 0,
          {
            status_callback <- make_encounter_status_callback(
              rem_sim_debug,
              "REM",
              set_progress = setProgress
            )
            status_callback(
              stage = "setup",
              detail = "Preparing REM inputs from the shared spatial simulator and initial model build...",
              value = 0.1
            )
            
            fit_result <- run_REM(
              y            = shared_d$camera_counts,
              r_km         = shared_d$out$`Detection Distance` / 1000,
              camera_days  = shared_d$camera_days,
              theta_deg    = input$theta,
              iter         = input$iter_rem_tte,
              burnin       = input$burnin_rem_tte,
              thin         = input$thin_rem_tte,
              n_chains     = app_n_chains(),
              D_max        = input$D_max,
              log_v_mean   = input$log_v_mean,
              log_v_sd     = input$log_v_sd,
              sd_eps_max   = input$sd_eps_max,
              status_callback = status_callback
            )
            setProgress(1.0, detail = "Complete!")
            fit_result
          }
        )
      },
      error = function(e) {
        update_model_debug(
          rem_sim_debug,
          status = "error",
          stage = "Failed",
          finished_at = Sys.time(),
          guidance = friendly_model_error("REM", "simulated shared spatial data", e$message),
          raw_error = e$message,
          log_entry = paste("Simulated REM run failed:", e$message)
        )
        showNotification(
          paste("REM (simulated) failed:", e$message),
          type = "error",
          duration = NULL
        )
        return(NULL)
      },
      finally = {
        rem_sim_running(FALSE)
      }
    )
    rem_sim_fit(fit)
    if (!is.null(fit)) {
      update_model_debug(
        rem_sim_debug,
        status = "success",
        stage = "Complete",
        finished_at = Sys.time(),
        guidance = "REM shared-simulator run completed successfully. Review the summary above and use Compare & combine after the other simulated models finish.",
        log_entry = "Simulated REM run completed."
      )
      showNotification("REM (simulated) complete!", type = "message")
    }
  })
  
  # --- TTE: simulated ---
  
  observeEvent(input$run_tte_sim, {
    shared_truth <- sim()
    shared_d <- if (!is.null(shared_truth)) sim_data() else NULL

    if (is.null(shared_d)) {
      showNotification(
        "Run the shared spatial simulator in the Simulate data tab first.",
        type = "error",
        duration = 6
      )
      return(NULL)
    }
    
    tte_sim_dataset_id(current_shared_sim_id())
    tte_sim_debug(make_model_debug(
      "TTE",
      "Shared spatial simulator",
      "Shared spatial simulator",
      "TTE uses total animal events per camera and camera-days derived from the shared spatial simulator."
    ))
    update_model_debug(
      tte_sim_debug,
      status = "running",
      stage = "Preflight checks",
      started_at = Sys.time(),
      guidance = "TTE is running on the shared spatial simulator so it can be compared and combined with simulated USCR and REM.",
      raw_error = NULL,
      context = summarize_shared_sim_context(shared_d, "TTE", sim_truth = shared_truth$truth),
      log_entry = "Simulated TTE shared-spatial run requested."
    )
    
    tte_sim_running(TRUE)
    tte_sim_fit(NULL)
    
    showNotification("Running TTE on shared simulated spatial data...", type = "message", duration = 6)
    
    fit <- tryCatch(
      {
        withProgress(
          message = "Running TTE model",
          detail = "Preparing TTE inputs from the shared spatial simulator...",
          value = 0,
          {
            status_callback <- make_encounter_status_callback(
              tte_sim_debug,
              "TTE",
              set_progress = setProgress
            )
            status_callback(
              stage = "setup",
              detail = "Preparing TTE inputs from the shared spatial simulator and initial model build...",
              value = 0.1
            )
            
            fit_result <- run_TTE(
              y            = shared_d$camera_counts,
              r_km         = shared_d$out$`Detection Distance` / 1000,
              camera_days  = shared_d$camera_days,
              theta_deg    = input$theta,
              iter         = input$iter_rem_tte,
              burnin       = input$burnin_rem_tte,
              thin         = input$thin_rem_tte,
              n_chains     = app_n_chains(),
              D_max        = input$D_max,
              log_v_mean   = input$log_v_mean,
              log_v_sd     = input$log_v_sd,
              sd_eps_max   = input$sd_eps_max,
              status_callback = status_callback
            )
            setProgress(1.0, detail = "Complete!")
            fit_result
          }
        )
      },
      error = function(e) {
        update_model_debug(
          tte_sim_debug,
          status = "error",
          stage = "Failed",
          finished_at = Sys.time(),
          guidance = friendly_model_error("TTE", "simulated shared spatial data", e$message),
          raw_error = e$message,
          log_entry = paste("Simulated TTE run failed:", e$message)
        )
        showNotification(
          paste("TTE (simulated) failed:", e$message),
          type = "error",
          duration = NULL
        )
        return(NULL)
      },
      finally = {
        tte_sim_running(FALSE)
      }
    )
    tte_sim_fit(fit)
    if (!is.null(fit)) {
      update_model_debug(
        tte_sim_debug,
        status = "success",
        stage = "Complete",
        finished_at = Sys.time(),
        guidance = "TTE shared-simulator run completed successfully. Review the summary above and use Compare & combine after the other simulated models finish.",
        log_entry = "Simulated TTE run completed."
      )
      showNotification("TTE (simulated) complete!", type = "message")
    }
  })
  
  # --- USCR: NPS ---
  
  observeEvent(input$run_uscr_nps, {
    req(nps_model_inputs())
    d <- nps_model_inputs()
    uscr_nps_debug(make_model_debug(
      "USCR",
      "Uploaded NPS data",
      "Uploaded NPS data; simulated grid",
      "Total animal events per camera, camera-days, and buffered state space."
    ))
    update_model_debug(
      uscr_nps_debug,
      status = "running",
      stage = "Preflight checks",
      started_at = Sys.time(),
      guidance = "USCR supports uploaded NPS data. This panel will show whether the run is still in setup, adaptive tuning, or the final run.",
      raw_error = NULL,
      context = summarize_uscr_context(d, source_label = "NPS"),
      log_entry = "NPS USCR run requested."
    )
    validate(
      need(all(d$camera_days > 0),
           "Some cameras have zero camera-days; check deployment dates.")
    )
    
    # Reset stop flag
    stop_uscr_nps(FALSE)
    uscr_nps_running(TRUE)
    uscr_nps_fit(NULL)  # Clear previous results
    
    showNotification("Running USCR on NPS data.", type = "message", duration = 6)
    
    fit <- tryCatch(
      {
        # Check if stopped before starting
        if (stop_uscr_nps()) {
          showNotification("USCR (NPS) run cancelled.", type = "warning")
          return(NULL)
        }
        
        withProgress(
          message = "Running USCR model",
          detail = "Checking inputs and preparing the USCR run...",
          value = 0,
          {
            setProgress(0.05, detail = "Checking inputs and preparing the USCR run...")
            update_model_debug(
              uscr_nps_debug,
              stage = "Preparing run",
              log_entry = "Input checks passed. Starting USCR setup."
            )
            
            # Check stop flag again before running
            if (stop_uscr_nps()) {
              showNotification("USCR (NPS) run cancelled.", type = "warning")
              return(NULL)
            }
            
            status_callback <- function(stage, detail = NULL, value = NULL) {
              stage_label <- switch(
                stage,
                setup = "Preparing state space and model code",
                tuning = "Adaptive tuning",
                final_run = "Final MCMC run",
                stage
              )
              update_model_debug(
                uscr_nps_debug,
                status = "running",
                stage = stage_label,
                log_entry = detail %||% paste("USCR stage:", stage_label)
              )
              if (!is.null(value)) {
                setProgress(value = value, detail = detail)
              } else if (!is.null(detail)) {
                setProgress(detail = detail)
              }
            }
            
            fit_result <- do.call(
              run_USCR_app,
              c(
                list(
                  out = d$out,
                  camera_counts = d$camera_counts,
                  camera_days = d$camera_days
                ),
                uscr_run_args(waic = TRUE, status_callback = status_callback)
              )
            )
            setProgress(1.0, detail = "Complete!")
            fit_result
          }
        )
      },
      error = function(e) {
        if (stop_uscr_nps()) {
          update_model_debug(
            uscr_nps_debug,
            status = "stopped",
            stage = "Stopped by user",
            finished_at = Sys.time(),
            guidance = "The NPS USCR run was stopped manually before completion.",
            raw_error = e$message,
            log_entry = "NPS USCR run stopped by user."
          )
          showNotification("USCR (NPS) run stopped by user.", type = "warning")
        } else {
          update_model_debug(
            uscr_nps_debug,
            status = "error",
            stage = "Failed",
            finished_at = Sys.time(),
            guidance = friendly_model_error("USCR", "uploaded NPS data", e$message),
            raw_error = e$message,
            log_entry = paste("NPS USCR failed:", e$message)
          )
          showNotification(
            paste("USCR (NPS) failed:", e$message),
            type = "error", duration = NULL
          )
        }
        return(NULL)
      },
      finally = {
        uscr_nps_running(FALSE)
        stop_uscr_nps(FALSE)  # Reset stop flag
      }
    )
    uscr_nps_fit(fit)
    if (!is.null(fit) && !stop_uscr_nps()) {
      update_model_debug(
        uscr_nps_debug,
        status = "success",
        stage = "Complete",
        finished_at = Sys.time(),
        guidance = "NPS USCR completed successfully. Review the summary above and the tuning history below.",
        log_entry = "NPS USCR run completed."
      )
      showNotification("USCR (NPS) complete!", type = "message")
    }
  })
  
  # --- REM: NPS ---
  
  observeEvent(input$run_rem_nps, {
    if (rem_nps_running()) {
      showNotification(
        "REM (NPS) is already running in this session. Wait for it to finish before starting another REM job.",
        type = "warning",
        duration = 6
      )
      return(NULL)
    }
    req(nps_model_inputs())
    d <- nps_model_inputs()
    requested_chains <- app_n_chains()
    rem_context <- c(
      summarize_rem_context(d),
      paste("Chains requested in UI:", requested_chains),
      "Background worker run: yes; the app enforces a minimum of 3 chains."
    )
    rem_args <- list(
      y            = d$camera_counts,
      r_km         = d$out$`Detection Distance` / 1000,
      camera_days  = d$camera_days,
      theta_deg    = input$theta,
      iter         = input$iter_rem_tte,
      burnin       = input$burnin_rem_tte,
      thin         = input$thin_rem_tte,
      n_chains     = requested_chains,
      D_max        = input$D_max,
      log_v_mean   = input$log_v_mean,
      log_v_sd     = input$log_v_sd,
      sd_eps_max   = input$sd_eps_max
    )
    rem_nps_debug(make_model_debug(
      "REM",
      "Uploaded NPS data",
      "Uploaded NPS data only",
      "REM uses total animal events per camera and camera-days."
    ))
    update_model_debug(
      rem_nps_debug,
      status = "running",
      stage = "Queued background run",
      started_at = Sys.time(),
      guidance = "REM on uploaded NPS data now runs in a background worker so other sessions stay responsive while still using the app minimum of 3 chains.",
      raw_error = NULL,
      context = rem_context,
      log_entry = "NPS REM run requested."
    )
    validate(
      need(all(d$camera_days > 0),
           "Some cameras have zero camera-days; check deployment dates.")
    )
    
    # Reset stop flag
    stop_rem_nps(FALSE)
    stop_flags_env$stop_rem_nps <- FALSE
    rem_nps_running(TRUE)
    rem_nps_background(TRUE)
    rem_nps_fit(NULL)  # Clear previous results
    
    update_model_debug(
      rem_nps_debug,
      stage = "Background run in progress",
      log_entry = paste(
        "Submitted REM background job.",
        "Chains used =", requested_chains,
        "(minimum enforced by the app).",
        "REM may internally increase iterations until max Rhat is acceptable."
      )
    )
    
    showNotification(
      "REM on NPS data is running in the background. This session should stay responsive while it runs.",
      type = "message",
      duration = NULL,
      id = "rem_nps_progress"
    )
    
    rem_future <- future::future({
      do.call(run_REM, rem_args)
    })
    
    promises::as.promise(rem_future) %...>% (function(fit) {
      removeNotification("rem_nps_progress")
      rem_nps_fit(fit)
      rem_nps_running(FALSE)
      rem_nps_background(FALSE)
      stop_rem_nps(FALSE)
      stop_flags_env$stop_rem_nps <- FALSE
      update_model_debug(
        rem_nps_debug,
        status = "success",
        stage = "Complete",
        finished_at = Sys.time(),
        guidance = "REM completed successfully in a background worker. Review the summary above and the tuning history below.",
        log_entry = "NPS REM background run completed."
      )
      showNotification("REM (NPS) complete!", type = "message", duration = 5)
      invisible(NULL)
    }) %...!% (function(e) {
      removeNotification("rem_nps_progress")
      rem_nps_fit(NULL)
      rem_nps_running(FALSE)
      rem_nps_background(FALSE)
      stop_rem_nps(FALSE)
      stop_flags_env$stop_rem_nps <- FALSE
      update_model_debug(
        rem_nps_debug,
        status = "error",
        stage = "Failed",
        finished_at = Sys.time(),
        guidance = friendly_model_error("REM", "uploaded NPS data", conditionMessage(e)),
        raw_error = conditionMessage(e),
        log_entry = paste("NPS REM background run failed:", conditionMessage(e))
      )
      showNotification(
        paste("REM (NPS) failed:", conditionMessage(e)),
        type = "error",
        duration = NULL
      )
      invisible(NULL)
    })
  })
  
  # --- TTE: NPS ---
  
  observeEvent(input$run_tte_nps, {
    req(nps_model_inputs())
    d <- nps_model_inputs()
    tte_nps_debug(make_model_debug(
      "TTE",
      "Uploaded NPS data",
      "Uploaded NPS data only",
      "Current build passes total animal events per camera and camera-days."
    ))
    update_model_debug(
      tte_nps_debug,
      status = "running",
      stage = "Preflight checks",
      started_at = Sys.time(),
      guidance = "TTE only runs on uploaded NPS data in this app.",
      raw_error = NULL,
      context = summarize_tte_context(d),
      log_entry = "NPS TTE run requested."
    )
    validate(
      need(all(d$camera_days > 0),
           "Some cameras have zero camera-days; check deployment dates.")
    )
    
    # Reset stop flag
    stop_tte_nps(FALSE)
    tte_nps_running(TRUE)
    tte_nps_fit(NULL)  # Clear previous results
    
    showNotification("Running TTE on NPS data... This may take several minutes.", 
                     type = "message", duration = 10)
    
    fit <- tryCatch(
      {
        # Check if stopped before starting
        if (stop_tte_nps()) {
          showNotification("TTE (NPS) run cancelled.", type = "warning")
          return(NULL)
        }
        
        withProgress(
          message = "Running TTE model",
          detail = "Preparing TTE inputs from uploaded NPS data...",
          value = 0,
          {
            status_callback <- make_encounter_status_callback(
              tte_nps_debug,
              "TTE",
              set_progress = setProgress
            )
            status_callback(
              stage = "setup",
              detail = "Preparing TTE inputs from uploaded data and initial model build...",
              value = 0.1
            )
            
            # Check stop flag again before running
            if (stop_tte_nps()) {
              showNotification("TTE (NPS) run cancelled.", type = "warning")
              return(NULL)
            }
            
            fit_result <- run_TTE(
              y            = d$camera_counts,
              r_km         = d$out$`Detection Distance` / 1000,
              camera_days  = d$camera_days,
              theta_deg    = input$theta,
              iter         = input$iter_rem_tte,
              burnin       = input$burnin_rem_tte,
              thin         = input$thin_rem_tte,
              n_chains     = app_n_chains(),
              D_max        = input$D_max,
              log_v_mean   = input$log_v_mean,
              log_v_sd     = input$log_v_sd,
              sd_eps_max   = input$sd_eps_max,
              status_callback = status_callback
            )
            setProgress(1.0, detail = "Complete!")
            fit_result
          }
        )
      },
      error = function(e) {
        if (stop_tte_nps()) {
          update_model_debug(
            tte_nps_debug,
            status = "stopped",
            stage = "Stopped by user",
            finished_at = Sys.time(),
            guidance = "The TTE run was stopped manually before completion.",
            raw_error = e$message,
            log_entry = "NPS TTE run stopped by user."
          )
          showNotification("TTE (NPS) run stopped by user.", type = "warning")
        } else {
          update_model_debug(
            tte_nps_debug,
            status = "error",
            stage = "Failed",
            finished_at = Sys.time(),
            guidance = friendly_model_error("TTE", "uploaded NPS data", e$message),
            raw_error = e$message,
            log_entry = paste("NPS TTE failed:", e$message)
          )
          showNotification(
            paste("TTE (NPS) failed:", e$message),
            type = "error", duration = NULL
          )
        }
        return(NULL)
      },
      finally = {
        tte_nps_running(FALSE)
        stop_tte_nps(FALSE)  # Reset stop flag
      }
    )
    tte_nps_fit(fit)
    if (!is.null(fit) && !stop_tte_nps()) {
      update_model_debug(
        tte_nps_debug,
        status = "success",
        stage = "Complete",
        finished_at = Sys.time(),
        guidance = "TTE completed successfully. Review the summary above.",
        log_entry = "NPS TTE run completed."
      )
      showNotification("TTE (NPS) complete!", type = "message")
    }
  })
  
  # ================================================================
  # MODEL TAB OUTPUTS (USCR / REM / TTE)
  # ================================================================
  
  # USCR summaries
  output$uscr_sim_text <- renderPrint({
    dbg <- uscr_sim_debug()
    if (uscr_sim_running()) {
      cat("⏳ USCR model is running...\n")
      cat("Current stage:", dbg$stage, "\n")
      cat("Open 'Run status & troubleshooting' below for more detail.\n")
      return(invisible(NULL))
    }
    if (identical(dbg$status, "error")) {
      cat("USCR (simulated): the last run failed.\n")
      cat("See 'Run status & troubleshooting' below for the raw error and guidance.\n")
      return(invisible(NULL))
    }
    if (identical(dbg$status, "stopped")) {
      cat("USCR (simulated): the last run was stopped before completion.\n")
      return(invisible(NULL))
    }
    fit <- uscr_sim_fit()
    if (is.null(fit)) {
      cat("USCR (simulated): not run yet. Click 'Run USCR on simulated data'.")
      return(invisible(NULL))
    }
    s <- summarize_method(fit)
    list(
      dataset                   = "Shared spatial simulator",
      mean_density_animals_per_km2 = round(s$mean_km2, 2),
      CI95_km2                  = c(round(s$q2.5_km2, 2), round(s$q97.5_km2, 2)),
      mean_density_animals_per_mi2 = round(s$mean_mi2, 2),
      CI95_mi2                  = c(round(s$q2.5_mi2, 2), round(s$q97.5_mi2, 2)),
      WAIC                      = round(s$waic, 2),
      note                      = if (!is.null(sim()))
        sprintf("True simulated density = %.1f animals/km²", sim()$truth$D_per_km2)
    )
  })
  
  output$uscr_nps_text <- renderPrint({
    dbg <- uscr_nps_debug()
    if (uscr_nps_running()) {
      cat("⏳ USCR model is running...\n")
      cat("Current stage:", dbg$stage, "\n")
      cat("Open 'Run status & troubleshooting' below for more detail.\n")
      return(invisible(NULL))
    }
    if (identical(dbg$status, "error")) {
      cat("USCR (NPS): the last run failed.\n")
      cat("See 'Run status & troubleshooting' below for the raw error and guidance.\n")
      return(invisible(NULL))
    }
    if (identical(dbg$status, "stopped")) {
      cat("USCR (NPS): the last run was stopped before completion.\n")
      return(invisible(NULL))
    }
    fit <- uscr_nps_fit()
    if (is.null(fit)) {
      cat("USCR (NPS): not run yet. Click 'Run USCR on NPS data'.")
      return(invisible(NULL))
    }
    s <- summarize_method(fit)
    list(
      dataset                   = "NPS data",
      mean_density_animals_per_km2 = round(s$mean_km2, 2),
      CI95_km2                  = c(round(s$q2.5_km2, 2), round(s$q97.5_km2, 2)),
      mean_density_animals_per_mi2 = round(s$mean_mi2, 2),
      CI95_mi2                  = c(round(s$q2.5_mi2, 2), round(s$q97.5_mi2, 2)),
      WAIC                      = round(s$waic, 2)
    )
  })
  
  # REM summaries
  output$rem_sim_text <- renderPrint({
    dbg <- rem_sim_debug()
    if (rem_sim_running()) {
      cat("⏳ REM model is running on the shared spatial simulator...\n")
      cat("Current stage:", dbg$stage, "\n")
      cat("Open 'Run status & troubleshooting' below for more detail.\n")
      return(invisible(NULL))
    }
    if (identical(dbg$status, "error")) {
      cat("REM (simulated): the last run failed.\n")
      cat("See 'Run status & troubleshooting' below for the raw error and guidance.\n")
      return(invisible(NULL))
    }
    fit <- rem_sim_fit()
    if (is.null(fit)) {
      cat("REM (simulated): not run yet. Generate the shared simulated dataset first, then click 'Run REM on simulated data'.")
      return(invisible(NULL))
    }
    s <- summarize_method(fit)
    truth <- sim()
    list(
      dataset                   = "Shared spatial simulator",
      mean_density_animals_per_km2 = round(s$mean_km2, 2),
      CI95_km2                  = c(round(s$q2.5_km2, 2), round(s$q97.5_km2, 2)),
      mean_density_animals_per_mi2 = round(s$mean_mi2, 2),
      CI95_mi2                  = c(round(s$q2.5_mi2, 2), round(s$q97.5_mi2, 2)),
      WAIC                      = round(s$waic, 2),
      note                      = if (!is.null(truth))
        sprintf(
          "%s true simulated density = %.1f animals/km²",
          "Shared spatial simulator",
          truth$truth$D_per_km2
        )
    )
  })
  
  output$rem_nps_text <- renderPrint({
    dbg <- rem_nps_debug()
    if (rem_nps_running()) {
      cat("⏳ REM model is running in the background...\n")
      cat("Current stage:", dbg$stage, "\n")
      cat("This session should stay responsive while the REM job runs.\n")
      cat("Open 'Run status & troubleshooting' below for more detail.\n")
      return(invisible(NULL))
    }
    if (identical(dbg$status, "error")) {
      cat("REM (NPS): the last run failed.\n")
      cat("See 'Run status & troubleshooting' below for the raw error and guidance.\n")
      return(invisible(NULL))
    }
    if (identical(dbg$status, "stopped")) {
      cat("REM (NPS): the last run was stopped before completion.\n")
      return(invisible(NULL))
    }
    fit <- rem_nps_fit()
    if (is.null(fit)) {
      cat("REM (NPS): not run yet. Click 'Run REM on NPS data'.")
      return(invisible(NULL))
    }
    s <- summarize_method(fit)
    list(
      dataset                   = "NPS data",
      mean_density_animals_per_km2 = round(s$mean_km2, 2),
      CI95_km2                  = c(round(s$q2.5_km2, 2), round(s$q97.5_km2, 2)),
      mean_density_animals_per_mi2 = round(s$mean_mi2, 2),
      CI95_mi2                  = c(round(s$q2.5_mi2, 2), round(s$q97.5_mi2, 2)),
      WAIC                      = round(s$waic, 2)
    )
  })
  
  # TTE summaries
  output$tte_sim_text <- renderPrint({
    dbg <- tte_sim_debug()
    if (tte_sim_running()) {
      cat("⏳ TTE model is running on the shared spatial simulator...\n")
      cat("Current stage:", dbg$stage, "\n")
      cat("Open 'Run status & troubleshooting' below for more detail.\n")
      return(invisible(NULL))
    }
    if (identical(dbg$status, "error")) {
      cat("TTE (simulated): the last run failed.\n")
      cat("See 'Run status & troubleshooting' below for the raw error and guidance.\n")
      return(invisible(NULL))
    }
    fit <- tte_sim_fit()
    if (is.null(fit)) {
      cat("TTE (simulated): not run yet. Generate the shared simulated dataset first, then click 'Run TTE on simulated data'.")
      return(invisible(NULL))
    }
    s <- summarize_method(fit)
    truth <- sim()
    list(
      dataset                   = "Shared spatial simulator",
      mean_density_animals_per_km2 = round(s$mean_km2, 2),
      CI95_km2                  = c(round(s$q2.5_km2, 2), round(s$q97.5_km2, 2)),
      mean_density_animals_per_mi2 = round(s$mean_mi2, 2),
      CI95_mi2                  = c(round(s$q2.5_mi2, 2), round(s$q97.5_mi2, 2)),
      WAIC                      = round(s$waic, 2),
      note                      = if (!is.null(truth))
        sprintf(
          "%s true simulated density = %.1f animals/km²",
          "Shared spatial simulator",
          truth$truth$D_per_km2
        )
    )
  })
  
  output$tte_nps_text <- renderPrint({
    dbg <- tte_nps_debug()
    if (tte_nps_running()) {
      cat("⏳ TTE model is running...\n")
      cat("Current stage:", dbg$stage, "\n")
      cat("Open 'Run status & troubleshooting' below for more detail.\n")
      return(invisible(NULL))
    }
    if (identical(dbg$status, "error")) {
      cat("TTE (NPS): the last run failed.\n")
      cat("See 'Run status & troubleshooting' below for the raw error and guidance.\n")
      return(invisible(NULL))
    }
    if (identical(dbg$status, "stopped")) {
      cat("TTE (NPS): the last run was stopped before completion.\n")
      return(invisible(NULL))
    }
    fit <- tte_nps_fit()
    if (is.null(fit)) {
      cat("TTE (NPS): not run yet. Click 'Run TTE on NPS data'.")
      return(invisible(NULL))
    }
    s <- summarize_method(fit)
    list(
      dataset                   = "NPS data",
      mean_density_animals_per_km2 = round(s$mean_km2, 2),
      CI95_km2                  = c(round(s$q2.5_km2, 2), round(s$q97.5_km2, 2)),
      mean_density_animals_per_mi2 = round(s$mean_mi2, 2),
      CI95_mi2                  = c(round(s$q2.5_mi2, 2), round(s$q97.5_mi2, 2)),
      WAIC                      = round(s$waic, 2)
    )
  })
  
  output$uscr_sim_debug <- renderText({
    state <- uscr_sim_debug()
    if (identical(state$status, "running")) invalidateLater(1000, session)
    format_model_debug(state, uscr_sim_fit())
  })
  
  output$uscr_nps_debug <- renderText({
    state <- uscr_nps_debug()
    if (identical(state$status, "running")) invalidateLater(1000, session)
    format_model_debug(state, uscr_nps_fit())
  })
  
  output$rem_sim_debug <- renderText({
    state <- rem_sim_debug()
    if (identical(state$status, "running")) invalidateLater(1000, session)
    format_model_debug(state, rem_sim_fit())
  })
  
  output$rem_nps_debug <- renderText({
    state <- rem_nps_debug()
    if (identical(state$status, "running")) invalidateLater(1000, session)
    format_model_debug(state, rem_nps_fit())
  })
  
  output$tte_sim_debug <- renderText({
    state <- tte_sim_debug()
    if (identical(state$status, "running")) invalidateLater(1000, session)
    format_model_debug(state, tte_sim_fit())
  })
  
  output$tte_nps_debug <- renderText({
    state <- tte_nps_debug()
    if (identical(state$status, "running")) invalidateLater(1000, session)
    format_model_debug(state, tte_nps_fit())
  })
  
  # ================================================================
  # COMPARE & COMBINE (WAIC-based)
  # ================================================================
  
  shared_sim_fits <- reactive({
    current_id <- current_shared_sim_id()
    if (is.null(current_id)) return(list())

    fits <- list(
      REM = if (identical(rem_sim_dataset_id(), current_id)) rem_sim_fit() else NULL,
      TTE = if (identical(tte_sim_dataset_id(), current_id)) tte_sim_fit() else NULL,
      USCR = if (identical(uscr_sim_dataset_id(), current_id)) uscr_sim_fit() else NULL
    )
    fits[!vapply(fits, is.null, logical(1))]
  })

  sim_combo <- reactive({
    fits <- shared_sim_fits()
    if (!length(fits)) return(NULL)
    build_combo_table_from_fits(fits)
  })
  
  nps_combo <- reactive({
    build_combo_table_from_fits(
      list(
        REM = rem_nps_fit(),
        TTE = tte_nps_fit(),
        USCR = uscr_nps_fit()
      )
    )
  })
  
  output$sim_combo_table <- renderDT({
    combo <- sim_combo()
    if (is.null(combo)) {
      return(DT::datatable(
        data.frame(Note = "Run one or more models from the shared spatial simulator first. Teaching-simulator runs are shown only in their own model tabs."),
        options = list(dom = "t", paging = FALSE),
        rownames = FALSE
      ))
    }
    df <- combo$table %>%
      mutate(across(where(is.numeric), ~ round(.x, 3)))
    datatable(df, options = list(pageLength = 5, dom = "t"))
  })
  
  output$dl_sim_uscr_csv <- downloadHandler(
    filename = function() {
      paste0("DEER_shared_sim_posterior_summary_", Sys.Date(), ".csv")
    },
    content = function(file) {
      fits <- shared_sim_fits()
      req(length(fits) > 0)
      out <- purrr::imap_dfr(fits, function(fit, model_name) {
        df <- posterior_summary_df(fit)
        if (is.null(df)) return(NULL)
        dplyr::mutate(df, model = model_name, .before = 1)
      })
      req(nrow(out) > 0)
      readr::write_csv(out, file)
    }
  )
  
  output$dl_nps_all_csv <- downloadHandler(
    filename = function() {
      paste0("DEER_NPS_posterior_summary_", Sys.Date(), ".csv")
    },
    content = function(file) {
      fits <- list(
        REM = rem_nps_fit(),
        TTE = tte_nps_fit(),
        USCR = uscr_nps_fit()
      )
      fits <- fits[!vapply(fits, is.null, logical(1))]
      req(length(fits) > 0)
      df <- dplyr::bind_rows(
        lapply(names(fits), function(model_name) {
          dplyr::mutate(
            posterior_summary_df(fits[[model_name]]),
            model = model_name,
            .before = 1
          )
        })
      )
      readr::write_csv(df, file)
    }
  )
  
  output$nps_combo_table <- renderDT({
    combo <- nps_combo()
    if (is.null(combo)) {
      return(DT::datatable(
        data.frame(
          Note = "Run at least one uploaded-data model from its tab first. The table will update as REM, TTE, and USCR finish."
        ),
        options = list(dom = "t", paging = FALSE),
        rownames = FALSE
      ))
    }
    df <- combo$table %>%
      mutate(across(where(is.numeric), ~ round(.x, 3)))
    datatable(df, options = list(pageLength = 5, dom = "t"))
  })
}

# -------------------------------------------------------------------
# Run the app
# -------------------------------------------------------------------

shinyApp(ui, server)
