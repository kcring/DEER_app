# app.R
# Deer Density Lab – USCR vs REM vs TTE
# Models run via run_USCR(), run_REM(), run_TTE() in R/sim_and_models.R

# -------------------------------------------------------------------
# Packages
# -------------------------------------------------------------------

needs <- c(
  "shiny", "bslib", "DT", "ggplot2", "dplyr", "tidyr",
  "readr", "purrr", "stringr", "secr", "data.table",
  "leaflet", "ggrepel",
  "nimble", "nimbleHMC", "parallel", "MCMCvis", "lubridate",
  "calecopal"
)

missing <- needs[!needs %in% installed.packages()[, "Package"]]
if (length(missing)) install.packages(missing, dependencies = TRUE)

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
  library(nimbleHMC)
  library(calecopal)
})

# Get redwood1 palette colors
redwood_colors <- cal_palette("redwood1")

# -------------------------------------------------------------------
# Helper files
# -------------------------------------------------------------------

source("R/sim_and_models.R")   # run_USCR(), run_REM(), run_TTE()
source("R/data_checks.R")      # QC + summary helpers

# -------------------------------------------------------------------
# Simulation helpers (secr + model input builders)
# -------------------------------------------------------------------

simulate_camera_counts <- function(n_side = 5, spacing_m = 300, days = 21,
                                   D_per_km2 = 25, lambda0 = 0.20,
                                   sigma_m = 150, seed = 1) {
  set.seed(seed)
  
  traps <- secr::make.grid(
    nx       = n_side,
    ny       = n_side,
    spacing  = spacing_m,
    detector = "count"
  )
  
  mask <- secr::make.mask(
    traps   = traps,
    buffer  = 4 * sigma_m,
    spacing = spacing_m / 3
  )
  
  # secr uses animals/ha when coords are in metres
  D_ha <- D_per_km2 / 100
  
  ch <- secr::sim.capthist(
    traps      = traps,
    popn       = list(D = D_ha),
    detectfn   = "HHN",
    detectpar  = list(lambda0 = lambda0, sigma = sigma_m),
    noccasions = days,
    renumber   = FALSE
  )
  
  list(
    ch    = ch,
    traps = traps,
    mask  = mask,
    truth = list(
      D_per_km2 = D_per_km2,
      lambda0   = lambda0,
      sigma_m   = sigma_m,
      days      = days,
      spacing_m = spacing_m,
      n_cams    = nrow(traps)
    )
  )
}

get_counts_matrix <- function(ch) {
  dims <- dim(ch)
  if (length(dims) == 3) {
    apply(ch, c(3, 2), sum, na.rm = TRUE)  # [trap, occasion]
  } else {
    as.matrix(ch)
  }
}

build_sim_data_for_nimble <- function(ch, detection_radius_m) {
  tr <- as.data.frame(secr::traps(ch))
  n_traps <- nrow(tr)
  
  counts_mat <- get_counts_matrix(ch)
  n_occ <- ncol(counts_mat)
  
  camera_counts <- rowSums(counts_mat, na.rm = TRUE)
  camera_days   <- rep(n_occ, n_traps)  # all cameras active all days
  
  out <- data.frame(
    Site                 = paste0("SimCam_", seq_len(n_traps)),
    `Detection Distance` = rep(detection_radius_m, n_traps),  # one value per camera
    utm_e                = tr$x / 1000,   # m -> km
    utm_n                = tr$y / 1000
  )
  
  list(out = out, camera_counts = camera_counts, camera_days = camera_days)
}

# --- Note: build_nps_model_inputs() is defined in R/sim_and_models.R ---
# This function matches the Rmd workflow exactly and handles:
# - Date parsing (including 2-digit years like "2/3/25")
# - UTM conversion and centering
# - Daily detection matrices
# - Start/End Index creation
# - Camera counts and camera days calculation

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
  waic_tbl <- tibble::tibble(
    model = c("REM", "TTE", "USCR"),
    waic  = c(
      get_waic_value(fits$REM),
      get_waic_value(fits$TTE),
      get_waic_value(fits$USCR)
    )
  ) %>%
    dplyr::mutate(
      deltaWAIC = waic - min(waic),
      rel_lik   = exp(-0.5 * deltaWAIC),
      w         = rel_lik / sum(rel_lik)
    ) %>%
    dplyr::arrange(model)  # REM, TTE, USCR
  
  rem_D  <- fits$REM$samples_all[,  "D_mi2"]
  tte_D  <- fits$TTE$samples_all[,  "D_mi2"]
  uscr_D <- fits$USCR$samples_all[, "D_mi2"]
  
  mat <- cbind(sort(rem_D), sort(tte_D), sort(uscr_D))
  colnames(mat) <- c("REM", "TTE", "USCR")
  
  unweighted_mean <- rowMeans(mat)
  weighted_mean   <- as.numeric(mat %*% waic_tbl$w)
  
  density_est <- cbind(mat,
                       unweighted_mean = unweighted_mean,
                       weighted_mean   = weighted_mean)
  
  means <- apply(density_est, 2, mean)
  lower <- apply(density_est, 2, stats::quantile, probs = 0.025)
  upper <- apply(density_est, 2, stats::quantile, probs = 0.975)
  prob20 <- apply(density_est, 2, function(x) mean(x > 20))
  
  table_out <- tibble::tibble(
    Method                    = c("REM", "TTE", "USCR",
                                  "unweighted mean", "weighted mean"),
    `Mean density (deer/mi²)` = as.numeric(means),
    `Lower 2.5%`              = as.numeric(lower),
    `Upper 97.5%`             = as.numeric(upper),
    `Prob > 20 DPSM`          = as.numeric(prob20),
    `WAIC weight`             = c(waic_tbl$w, NA, NA)
  )
  
  list(table = table_out, waic = waic_tbl)
}

# -------------------------------------------------------------------
# UI
# -------------------------------------------------------------------

ui <- page_fillable(
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
    "))
  ),
  tags$div(
    style = "text-align: center; padding: 5px 0; display: flex; align-items: center; justify-content: center; gap: 5px;",
    tags$img(src = "wvu_logo.png", alt = "WVU", style = "height: 50px; width: auto;"),
    tags$img(src = "deer_app_logo.png", alt = "DEER App", style = "height: 300px; width: auto;"),
    tags$img(src = "nps_logo.png", alt = "NPS", style = "height: 50px; width: auto;")
  ),
  
  navset_tab(
        
        # ------------------------- ABOUT ------------------------------
        nav_panel(
          "About",
          tags$div(
            class = "hero",
            tags$div(
              style = "max-width: 980px; margin: 0 auto; padding: 0 16px;",
              tags$span(class = "badge", style = "background: var(--rw4); border: 1px solid var(--rw3); color: var(--ink);", "National Parks • Winter camera survey"),
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
                "uSCR · REM · TTE — unmarked camera methods, model‑averaged to deer/mi²."
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
                "We use", tags$strong("three complementary models"), "to count and estimate deer in national parks.",
                "We then", tags$strong("average the three estimates"), "to report", tags$strong("density (deer/mi²)"), "with uncertainty.",
                "No model is perfect—so we", tags$strong("balance the three results"), "to reduce bias from any one method.",
                "The purpose of this app is to let NPS staff", tags$strong("calculate deer densities in their own park"), "using consistent methods,",
                "so", tags$strong("populations are comparable across parks at large scale"), "."
              ),
              tags$div(class = "divider"),
              tags$div(
                style = "display: grid; grid-template-columns: repeat(auto-fit, minmax(280px, 1fr)); gap: 14px;",
                tags$div(
                  class = "about-card",
                  tags$h3("Efficient!"),
                  tags$ul(
                    tags$li("No marking or individual ID; count independent events via", tags$em("Cluster ID"), "."),
                    tags$li("One winter deployment (≈2–8 weeks) meets closure assumptions.")
                  )
                ),
                tags$div(
                  class = "about-card",
                  tags$h3("Reproducible!"),
                  tags$ul(
                    tags$li("Arrays of ~16–25 cameras, ≤400 m spacing, unbaited placement."),
                    tags$li("Standardized CSVs and workflow make park‑to‑park comparisons possible.")
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
                tags$p(
                  "You have two options for data input:"
                ),
                tags$div(
                  style = "margin: 1rem 0;",
                  tags$h4(style = "font-size: 1.1rem; font-weight: 500;", "Option 1: Simulate data"),
                  tags$p(
                    "Go to the", tags$strong("'Simulate data'"), "tab. Here you can adjust the simulation parameters:",
                    tags$ul(
                      tags$li("Grid dimension (n × n cameras)"),
                      tags$li("Camera spacing (meters)"),
                      tags$li("Number of days"),
                      tags$li("True density (deer/km²)"),
                      tags$li("Detection parameters (sigma, lambda0)"),
                      tags$li("Random seed for reproducibility")
                    ),
                    "Press the", tags$strong("'Simulate grid'"), "button and you will have toy data to run the models on. This is useful for testing the app and understanding how the models work."
                  )
                ),
                tags$div(
                  style = "margin: 1rem 0;",
                  tags$h4(style = "font-size: 1.1rem; font-weight: 500;", "Option 2: Upload your NPS data"),
                  tags$p(
                    "You will need two data files:",
                    tags$ul(
                      tags$li(tags$strong("Deployment CSV"), "— Contains camera deployment information (locations, dates, detection distances)"),
                      tags$li(tags$strong("Images CSV"), "— Contains detection data (timestamps, species, cluster IDs)")
                    ),
                    "It is", tags$strong("crucially important"), "that the structure of your dataframes matches the requirements. Click the", 
                    tags$strong("'Add your data'"), "tab to get more information on the required columns and data format. The app will automatically:",
                    tags$ul(
                      tags$li("Pre-clean column names and whitespace"),
                      tags$li("Run quality checks on your data"),
                      tags$li("Trim images to the first 56 days per camera")
                    )
                  )
                ),
                tags$h2(style = "font-size: 1.5rem; font-weight: 600; margin-top: 1rem;", "Step 2: Adjust settings (optional)"),
                tags$p(
                  "The models have default settings that work well for most cases. However, if you need to change something, you can:",
                  tags$ul(
                    tags$li("Go to the", tags$strong("'Add your data'"), "tab"),
                    tags$li("Scroll to the", tags$strong("'Settings mode'"), "section"),
                    tags$li("Click the", tags$strong("'Advanced'"), "button to reveal additional options"),
                    tags$li("Adjust MCMC settings (iterations, burn-in, thinning, number of chains)"),
                    tags$li("Modify priors for movement speed, detection parameters, or camera heterogeneity"),
                    tags$li("Change camera detection angle (default: 55°)")
                  )
                ),
                tags$h2(style = "font-size: 1.5rem; font-weight: 600; margin-top: 1rem;", "Step 3: Run the models"),
                tags$p(
                  "Use the", tags$strong("USCR"), ",", tags$strong("REM"), ", and", tags$strong("TTE"), "tabs to run each model:",
                  tags$ul(
                    tags$li("Each model tab has", tags$strong("separate buttons"), "for simulated vs NPS data"),
                    tags$li("Click the", tags$strong("'Run'"), "button for your chosen data type"),
                    tags$li("Models show", tags$strong("progress bars"), "and", tags$strong("status messages"), "while running"),
                    tags$li("You can", tags$strong("stop a model run"), "using the red 'Stop' button if needed"),
                    tags$li("Results appear below the buttons once each model completes")
                  ),
                  tags$em("Note:"), "USCR typically takes 30-60 minutes, while REM and TTE are faster (several minutes)."
                ),
                tags$h2(style = "font-size: 1.5rem; font-weight: 600; margin-top: 1rem;", "Step 4: Compare results"),
                tags$p(
                  "Use the", tags$strong("'Compare & combine'"), "tab to:",
                  tags$ul(
                    tags$li("View per‑model results and density estimates"),
                    tags$li("See the", tags$strong("WAIC‑weighted average"), ", matching the Rmarkdown workflow"),
                    tags$li("Compare model performance using WAIC values"),
                    tags$li("Export the summary table for reporting")
                  ),
                  "The final output is", tags$strong("deer per square mile (mi²)"), "with credible intervals, providing a robust estimate that balances the strengths of all three models."
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
                    tags$strong("Unmarked Spatial Capture–Recapture (uSCR)"), "—", tags$em("Chandler & Royle, 2013"), ".",
                    tags$br(),
                    tags$strong("Spatial pattern → density."), "Where detections fall across the array tells us density inside your park's defined state‑space."
                  ),
                  tags$div(
                    class = "pcbox",
                    tags$div(
                      class = "pcsec",
                      tags$h4("Pros"),
                      tags$ul(
                        tags$li("Density within a clearly defined area (state‑space)."),
                        tags$li("No individual IDs; uses spatial correlation across cameras.")
                      )
                    ),
                    tags$div(
                      class = "pcsec",
                      tags$h4("Cons"),
                      tags$ul(
                        tags$li("Sensitive to camera spacing vs. movement scale; priors matter."),
                        tags$li("Heavier compute/time than REM/TTE.")
                      )
                    )
                  ),
                  tags$span(class = "ref", "Reference: Chandler, R.B. & Royle, J.A. (2013).")
                ),
                tags$div(
                  class = "about-card",
                  tags$span(class = "tag tag-rem", "REM"),
                  tags$p(
                    tags$strong("Random Encounter Model (REM)"), "—", tags$em("Rowcliffe et al., 2008"), ".",
                    tags$br(),
                    tags$strong("Encounters → density."), "Events per time become density using animal speed and the camera viewshed (radius & angle)."
                  ),
                  tags$div(
                    class = "pcbox",
                    tags$div(
                      class = "pcsec",
                      tags$h4("Pros"),
                      tags$ul(
                        tags$li("Simple and fast once speed & viewshed are known."),
                        tags$li("Well‑defined area (collective viewsheds of the array).")
                      )
                    ),
                    tags$div(
                      class = "pcsec",
                      tags$h4("Cons"),
                      tags$ul(
                        tags$li("Requires accurate speed & viewshed; assumes random, unbaited placement and independence."),
                        tags$li("Less spatial detail than uSCR.")
                      )
                    )
                  ),
                  tags$span(class = "ref", "Reference: Rowcliffe, J.M. et al. (2008).")
                ),
                tags$div(
                  class = "about-card",
                  tags$span(class = "tag tag-tte", "TTE"),
                  tags$p(
                    tags$strong("Time‑to‑Event (TTE)"), "—", tags$em("Moeller et al., 2018"), ".",
                    tags$br(),
                    tags$strong("Waiting time → density."), "Shorter time (or detection days) to first encounter implies higher density for the same viewshed & speed."
                  ),
                  tags$div(
                    class = "pcbox",
                    tags$div(
                      class = "pcsec",
                      tags$h4("Pros"),
                      tags$ul(
                        tags$li("Handles many zero‑days; complements REM/uSCR."),
                        tags$li("Well‑matched to short winter closure periods.")
                      )
                    ),
                    tags$div(
                      class = "pcsec",
                      tags$h4("Cons"),
                      tags$ul(
                        tags$li("Shares REM assumptions; needs a speed prior to set time‑units."),
                        tags$li("Sensitive to viewshed and time‑unit calibration.")
                      )
                    )
                  ),
                  tags$span(class = "ref", "Reference: Moeller, A.K. et al. (2018).")
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
                  tags$span(class = "kicker", "Final output:"),
                  tags$strong("Deer per square mile (mi²)"), "with credible intervals, plus per‑model estimates for transparency."
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
                  "The underlying models and code were created by", 
                  tags$strong("Dr. Amanda Van Buskirk"), 
                  "under the advisement of", 
                  tags$strong("Dr. Christopher Rota"), 
                  "in the", 
                  tags$a(href = "https://sites.google.com/mix.wvu.edu/rotalab/home", target = "_blank",
                         style = "color: var(--rw1); text-decoration: underline;",
                         tags$strong("Rota Quantitative Ecology Lab")), 
                  "within the", 
                  tags$strong("Davis College of Agriculture and Natural Resources at West Virginia University"), 
                  "."
                )
              ),
              tags$div(
                class = "about-card",
                tags$h3(style = "font-size: 1.2rem; font-weight: 600; margin-top: 0.5rem;", "Shiny App Development"),
                tags$p(
                  "This Shiny application was developed as part of the", 
                  tags$a(href = "https://esa.org/programs/scip/", target = "_blank", 
                         style = "color: var(--rw1); text-decoration: underline;",
                         tags$strong("Science in the Parks Communications Fellowship")), 
                  ", a collaborative effort between the", 
                  tags$strong("Ecological Society of America (ESA)"), 
                  "and the", 
                  tags$strong("National Park Service (NPS)"), 
                  "."
                ),
                tags$p(
                  style = "margin-top: 0.75rem;",
                  tags$strong("Fellowship Support:"), tags$br(),
                  "•", tags$strong("Dr. Brian Mitchell"), "(NPS) — Fellowship Liaison", tags$br(),
                  "•", tags$strong("Jasjeet Dhanota"), "(ESA) — Mentor", tags$br(),
                  "•", tags$strong("Mary Joy Mulumba"), "(ESA) — Mentor"
                )
              )
            ),
            
            tags$div(class = "divider"),
            
            tags$div(
              class = "small",
              style = "margin-top: 22px; color: var(--muted);",
              tags$p("Palette inspired by", tags$code("calecopal::redwood1"), "with California redwood forest colors.")
            )
          )
        ),
        
        # ---------------------- SIMULATE DATA -------------------------
        nav_panel(
          "Simulate data",
          markdown(paste(
            "Simulate a camera grid with secr; models are run from the USCR, REM, and TTE tabs.",
            "",
            "Adjust the simulation parameters below, then click 'Simulate grid' to generate data.",
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
              sliderInput("Dtrue", "True density for simulation (deer/km²)",
                          min = 5, max = 80, value = 25, step = 1),
              sliderInput("sigma", "Home-range scale σ (m)",
                          min = 80, max = 300, value = 150, step = 5),
              sliderInput("lambda0", "Baseline detection rate λ₀ (per occasion)",
                          min = 0.05, max = 0.6, value = 0.20, step = 0.01)
            ),
            column(4,
              numericInput("seed", "Random seed", value = 1, min = 1),
              h4("Camera geometry"),
              sliderInput("r_m", "Detection radius r (m)",
                          min = 8, max = 25, value = 12, step = 1),
              sliderInput("theta", "Detection angle θ (degrees)",
                          min = 20, max = 80, value = 55, step = 1),
              br(),
              actionButton("run_sim", "Simulate grid", class = "btn-primary")
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
              h4("Simulated deer distribution"),
              plotOutput("sim_deer_distribution_plot", height = "420px")
            )
          ),
          h3("Camera data table"),
          DTOutput("camera_table")
        ),
        
        # ---------------------- ADD YOUR DATA -------------------------
        nav_panel(
          "Add your data",
          markdown(paste(
            "Upload your **deployment** and **images** CSVs in NPS SOP format.",
            "",
            sep = "\n"
          )),
          h3("The app will:"),
          markdown(paste(
            "1. Pre-clean column names and whitespace;",
            "2. Run `check_deployments()` and `check_images()`;",
            "3. Trim images to the first 56 days per camera using `trim_images_56days()`.",
            "",
            "Download cleaned/trimmed CSVs if needed, then run models in the USCR/REM/TTE tabs.",
            "",
            "---",
            "",
            "## Required data columns",
            "",
            "### **Deployment file** required columns:",
            "",
            "- **`Park`** — 4-character park code",
            "- **`Site Name`** — Camera site identifier",
            "- **`Camera ID`** — Camera identifier",
            "- **`Start Date`** — Deployment start date (MM/DD/YYYY or MM/DD/YY)",
            "- **`Start Time`** — Deployment start time",
            "- **`End Date`** — Deployment end date (MM/DD/YYYY or MM/DD/YY)",
            "- **`End Time`** — Deployment end time",
            "- **`Latitude`** — Camera latitude (decimal degrees)",
            "- **`Longitude`** — Camera longitude (decimal degrees)",
            "- **`Detection Distance`** — Effective detection radius in meters",
            "",
            "### **Images file** required columns:",
            "",
            "- **`Site Name`** — Must match deployment file",
            "- **`Timestamp`** — Image timestamp (date and time)",
            "- **`Species`** — Species identifier (e.g., 'Deer')",
            "- **`Cluster ID`** — Unique identifier for independent detection events",
            "- **`Sighting Count`** — (Optional) Number of individuals in the image",
            "- **`Image URL`** — (Optional) Link to image file",
            "",
            "**Note:** The app automatically handles date parsing, including 2-digit years (e.g., '2/3/25').",
            sep = "\n"
          )),
          
          h3("Step 1: Deployment file"),
          fileInput("deployment_csv", "Upload deployment CSV (NPS format)", accept = ".csv"),
          h4("Deployment check log"),
          verbatimTextOutput("deployment_check_log"),
          h4("Preview of cleaned deployment data"),
          DTOutput("deployment_preview"),
          downloadButton("download_deployment_checked", "Download cleaned deployment CSV"),
          
          hr(),
          
          h3("Step 2: Images file"),
          fileInput("images_csv", "Upload images CSV (NPS format)", accept = ".csv"),
          numericInput("survey_year", "Survey year (yyyy)",
                       value = 2025, min = 2000, max = 2100),
          h4("Images check log"),
          verbatimTextOutput("images_check_log"),
          h4("Preview of cleaned + trimmed images data"),
          DTOutput("images_preview"),
          downloadButton("download_images_checked", "Download cleaned images CSV"),
          
          hr(),
          
          h3("Step 3: Model settings and parameters"),
          markdown(paste(
            "### MCMC Settings",
            "",
            "**Recommended starting values** (from NPS guidance):",
            "",
            "- **REM/TTE**: Start with 6000 iterations, 1000 burn-in, thin = 5",
            "- **USCR**: Start with 11000 iterations, 1000 burn-in, thin = 10",
            "- **Number of chains**: 3 (for convergence diagnostics)",
            "- **Convergence criteria**: R̂ < 1.1 for all monitored parameters",
            "",
            "If convergence is poor (R̂ > 1.1), increase iterations by iter × 2 or increase thinning.",
            "",
            "### Camera Geometry",
            "",
            "Current camera detection angle is 55 degrees based on Browning camera model.",
            "Adjust if using different camera model or taking measurements in the field.",
            "",
            sep = "\n"
          )),
          
          h4("Camera geometry"),
          fluidRow(
            column(6,
              sliderInput("r_m", "Detection radius r (m)",
                          min = 8, max = 25, value = 12, step = 1)
            ),
            column(6,
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
              "**Using recommended defaults** for MCMC and priors (from NPS guidance).",
              "",
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
                               value = 3, min = 1, max = 8, step = 1)
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
                               value = 11000, min = 1000, step = 1000),
                  numericInput("burnin_uscr", "USCR burn-in",
                               value = 1000, min = 0, step = 500),
                  numericInput("thin_uscr", "USCR thinning",
                               value = 10, min = 1, step = 1),
                  numericInput("M_uscr", "USCR M (data augmentation)",
                               value = 1000, min = 200, step = 100)
                )
              ),
              
              markdown(paste(
                "**Note on M (data augmentation)**: Start at M = 1000. If posterior N approaches M,",
                "or if psi concentrates near 1, increase M to 2000-4000.",
                "",
                "**Note on USCR iterations**: If convergence warnings appear, increase iterations by iter × 2.",
                sep = "\n"
              )),
              
              hr(),
              
              h4("Priors: REM & TTE"),
              markdown(paste(
                "**Informative priors** for animal movement speed and camera-level noise:",
                "",
                "- **D** ~ Uniform(0, D_max): Population density prior. Increase D_max if extremely high densities expected.",
                "- **log(v)** ~ Normal(mean, SD): Animal movement speed (log scale). Mean ≈ 3.99 km/day.",
                "  Increase SD to weaken constraint or change mean to reflect regionally-specific deer movement estimates.",
                "- **sd_eps** ~ Uniform(0, max): Random overdispersion among cameras. Adjust upper bound based on variability.",
                "",
                "**Example modifications**:",
                "- Wider uncertainty in movement speed: log_v ~ Normal(1.339, 0.40)",
                "- Allowing higher densities: D ~ Uniform(0, 400)",
                "",
                sep = "\n"
              )),
              
              fluidRow(
                column(6,
                  numericInput("D_max", "Max density prior (D max, deer/km²)",
                               value = 200, min = 10, step = 10),
                  numericInput("log_v_mean", "log(v) mean (km/day)",
                               value = 1.339, step = 0.1),
                  numericInput("log_v_sd", "log(v) SD",
                               value = 0.2955, min = 0.01, step = 0.05)
                ),
                column(6,
                  numericInput("sd_eps_max", "sd_eps upper bound",
                               value = 10, min = 1, step = 1)
                )
              ),
              
              hr(),
              
              h4("Priors: USCR"),
              markdown(paste(
                "**Informative priors** for spatial scale and detection intensity:",
                "",
                "- **log(σ)** ~ Normal(mean, SD): Spatial scale parameter (home range) on log scale.",
                "  Based on average home range size of 1.1 km². Increase SD to weaken constraint if poor mixing,",
                "  or modify mean based on regionally-specific home range size estimates.",
                "- **log(λ₀)** ~ Normal(0, SD): Baseline detection intensity. Increase SD for weaker prior.",
                "- **psi** ~ Uniform(0, 1): Proportion of augmented individuals actually present. Do not modify.",
                "- **sd_eps** ~ Gamma(shape, rate): Camera-level random effect. Adjust shape to reduce prior influence.",
                "",
                "**Example modifications**:",
                "- Looser movement range prior: log_sigma ~ Normal(-1.4442, 0.30)",
                "- Higher uncertainty in baseline detection: log_lam_0 ~ Normal(0, 2)",
                "",
                "⚠️ **Model performance warning**: If psi peaks near 1 or N approaches M, increase M or increase iterations.",
                "",
                sep = "\n"
              )),
              
              fluidRow(
                column(6,
                  numericInput("log_sigma_mean", "log(σ) mean",
                               value = -1.4442, step = 0.1),
                  numericInput("log_sigma_sd", "log(σ) SD",
                               value = 0.1451, min = 0.01, step = 0.05),
                  numericInput("log_lam0_sd", "log(λ₀) SD",
                               value = 1, min = 0.1, step = 0.1)
                ),
                column(6,
                  numericInput("sd_eps_shape", "sd_eps gamma shape",
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
              numericInput("iter_uscr", NULL, value = 11000),
              numericInput("burnin_uscr", NULL, value = 1000),
              numericInput("thin_uscr", NULL, value = 10),
              numericInput("M_uscr", NULL, value = 1000),
              numericInput("D_max", NULL, value = 200),
              numericInput("log_v_mean", NULL, value = 1.339),
              numericInput("log_v_sd", NULL, value = 0.2955),
              numericInput("sd_eps_max", NULL, value = 10),
              numericInput("log_sigma_mean", NULL, value = -1.4442),
              numericInput("log_sigma_sd", NULL, value = 0.1451),
              numericInput("log_lam0_sd", NULL, value = 1),
              numericInput("sd_eps_shape", NULL, value = 1)
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
          h4("Deer detections"),
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
            <p><strong>The gist:</strong> Estimates <strong>deer density in the park/state‑space</strong> by learning from <strong>where and when</strong> deer are detected <strong>across a camera array</strong>.</p>
            <details>
              <summary style="cursor: pointer; font-weight: 600; margin: 1rem 0; padding: 0.5rem; background: #f0f4e8; border-left: 3px solid #609048; border-radius: 6px;"><strong>Under the hood (equations)</strong></summary>
              <div style="margin: 1rem 0; padding-left: 1rem;">
                <p>Observation model per camera \\(j\\) and occasion/day \\(t\\):</p>
                <p>$$Y_{jt} \\sim \\mathrm{Poisson}(\\lambda_{jt})$$</p>
                <p>Encounter rate with distance decay and camera heterogeneity:</p>
                <p>$$\\lambda_{jt} = \\mathrm{effort}_{jt}\\;\\lambda_0 \\sum_{i=1}^{M} z_i \\exp\\!\\Big(-\\frac{d_{ij}^2}{2\\sigma^2}\\Big)\\; \\exp(\\epsilon_j)$$</p>
                <p>Density is obtained by integrating over activity centers and dividing the inferred abundance by the state‑space area; the app reports \\(D_{\\mathrm{mi}^2}\\).</p>
                <h3>Variables</h3>
                <ul>
                  <li>\\(Y_{jt}\\) — count of independent deer events at camera \\(j\\) on day \\(t\\) (from <code>Cluster ID</code>).</li>
                  <li>\\(\\lambda_{jt}\\) — expected event rate for \\(j,t\\).</li>
                  <li>\\(\\mathrm{effort}_{jt}\\) — effort indicator (1 for deployed days from <code>Start/End</code>).</li>
                  <li>\\(\\lambda_0\\) — baseline encounter rate at ~0 m (internal).</li>
                  <li>\\(\\sigma\\) — space‑use scale (m), controls decay with distance (internal).</li>
                  <li>\\(d_{ij}\\) — distance (m; projected) from animal \\(i\\)\'s activity center to camera \\(j\\).</li>
                  <li>\\(z_i\\) — augmentation indicator (1 if individual \\(i\\) is real; internal).</li>
                  <li>\\(\\epsilon_j\\) — camera random effect (mean 0; internal).</li>
                  <li>\\(D_{\\mathrm{mi}^2}\\) — density, deer per square mile (reported).</li>
                </ul>
                <h3>Priors used in the app</h3>
                <ul>
                  <li>\\(\\log \\sigma \\sim \\mathcal{N}(-1.4442,\\,0.1451)\\) (anchors HR95 ≈ 1.1 km² ⇒ σ ≈ 236 m).</li>
                  <li>\\(\\log \\lambda_0 \\sim \\mathcal{N}(0,\\,1)\\) (weak baseline).</li>
                  <li>\\(\\psi \\sim \\mathcal{U}(0,1)\\) (augmentation inclusion).</li>
                  <li>\\(\\mathrm{sd}_\\epsilon \\sim \\mathrm{Gamma}(1,1)\\) (camera heterogeneity).</li>
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
          
          hr(),
          h4("NPS data"),
          div(
            actionButton("run_uscr_nps", "Run USCR on NPS data", class = "btn-primary"),
            actionButton("stop_uscr_nps", "Stop", class = "btn-danger", style = "margin-left: 10px;"),
            style = "margin-bottom: 10px;"
          ),
          br(),
          verbatimTextOutput("uscr_nps_text"),
          tags$div(
            style = "text-align: center; margin-top: 40px; padding-top: 20px; border-top: 1px solid #e0e0e0;",
            tags$h3(style = "margin: 0; color: var(--rw3);", "DEER App"),
            tags$p(style = "margin: 5px 0 0 0; color: var(--muted);", "Density Estimation from encounter rates")
          )
        ),
        
        # ---------------------- REM MODEL TAB -------------------------
        nav_panel(
          "REM model",
          tags$div(
            id = "rem-content",
            HTML('
            <h2>Model 2 — REM (Random Encounter Model)</h2>
            <p><strong>The gist:</strong> Converts <strong>how often deer pass a camera</strong> into density, correcting for <strong>movement speed</strong> and <strong>camera view geometry</strong>.</p>
            <details>
              <summary style="cursor: pointer; font-weight: 600; margin: 1rem 0; padding: 0.5rem; background: #f0f4e8; border-left: 3px solid #609048; border-radius: 6px;"><strong>Under the hood (equations)</strong></summary>
              <div style="margin: 1rem 0; padding-left: 1rem;">
                <p>Closed‑form estimator:</p>
                <p>$$D_{\\mathrm{km}^2} = \\frac{(y/t)}{\\,v \\; r_{\\mathrm{km}} \\; \\frac{(2+\\theta_{\\mathrm{rad}})}{\\pi}} \\qquad\\Rightarrow\\qquad D_{\\mathrm{mi}^2} = 2.59\\,D_{\\mathrm{km}^2}$$</p>
                <h3>Variables</h3>
                <ul>
                  <li>\\(y\\) — number of independent deer events per camera (unique <code>Cluster ID</code>).</li>
                  <li>\\(t\\) — camera‑days of effort from <code>Start/End</code>.</li>
                  <li>\\(v\\) — movement speed (km/day), prior on \\(\\log v\\) as below.</li>
                  <li>\\(r_{\\mathrm{km}}\\) — effective detection radius in km; \\(r_{\\mathrm{km}} = r_m/1000\\).</li>
                  <li>\\(r_m\\) — effective detection radius (meters), from <code>Detection Distance</code>.</li>
                  <li>\\(\\theta_{\\mathrm{rad}}\\) — camera half‑angle in radians (default \\(55^\\circ\\)).</li>
                  <li>\\(D_{\\mathrm{km}^2}, D_{\\mathrm{mi}^2}\\) — density in deer/km² and deer/mi² (reported in mi²).</li>
                </ul>
                <h3>Priors/inputs used in the app</h3>
                <ul>
                  <li>\\(\\log v \\sim \\mathcal{N}(1.339,\\,0.2955)\\) (median ≈ 3.8–4.0 km/day; literature).</li>
                  <li>\\(\\theta = 55^\\circ\\) by default (can be changed if measured).</li>
                  <li>\\(r_m\\) from field <code>Detection Distance</code> (m).</li>
                  <li>\\(\\mathrm{sd}_\\epsilon \\sim \\mathcal{U}(0,10)\\) (optional camera heterogeneity).</li>
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
          div(
            actionButton("run_rem_sim", "Run REM on simulated data", class = "btn-primary"),
            actionButton("stop_rem_sim", "Stop", class = "btn-danger", style = "margin-left: 10px;"),
            style = "margin-bottom: 10px;"
          ),
          br(),
          verbatimTextOutput("rem_sim_text"),
          
          hr(),
          h4("NPS data"),
          div(
            actionButton("run_rem_nps", "Run REM on NPS data", class = "btn-primary"),
            actionButton("stop_rem_nps", "Stop", class = "btn-danger", style = "margin-left: 10px;"),
            style = "margin-bottom: 10px;"
          ),
          br(),
          verbatimTextOutput("rem_nps_text"),
          tags$div(
            style = "text-align: center; margin-top: 40px; padding-top: 20px; border-top: 1px solid #e0e0e0;",
            tags$h3(style = "margin: 0; color: var(--rw3);", "DEER App"),
            tags$p(style = "margin: 5px 0 0 0; color: var(--muted);", "Density Estimation from encounter rates")
          )
        ),
        
        # ---------------------- TTE MODEL TAB -------------------------
        nav_panel(
          "TTE model",
          tags$div(
            id = "tte-content",
            HTML('
            <h2>Model 3 — TTE (Time‑to‑Event, periodized)</h2>
            <p><strong>The gist:</strong> Uses <strong>how many days with ≥1 deer detection</strong> each camera has, scaled by <strong>movement‑based time‑units</strong> and the <strong>viewshed area</strong>. More detection days per time on → higher density.</p>
            <details>
              <summary style="cursor: pointer; font-weight: 600; margin: 1rem 0; padding: 0.5rem; background: #f0f4e8; border-left: 3px solid #609048; border-radius: 6px;"><strong>Under the hood (equations)</strong></summary>
              <div style="margin: 1rem 0; padding-left: 1rem;">
                <p>Periodized Poisson per camera \\(j\\):</p>
                <p>$$y_j \\sim \\mathrm{Poisson}(\\lambda_j), \\qquad \\log \\lambda_j \\;=\\; \\log D \\;+\\; \\log(\\mathrm{tte\\_units}_j)\\;+\\;\\log(A_j)\\;+\\;\\epsilon_j$$</p>
                <p>Movement‑based time‑units and viewshed area:</p>
                <p>$$\\mathrm{tte\\_units}_j \\;=\\; \\frac{\\mathrm{camera\\_days}_j}{\\mathrm{time\\_unit}}, \\qquad \\mathrm{time\\_unit} \\;\\approx\\; \\frac{0.59\\, r_m}{v}$$</p>
                <p>$$A_j \\;=\\; \\frac{\\theta_{\\mathrm{rad}}}{2} \\; r_m^2 \\quad (\\mathrm{m}^2) \\qquad\\Rightarrow\\qquad A_j\\;(\\mathrm{km}^2) = \\frac{A_j\\;(\\mathrm{m}^2)}{10^6}$$</p>
                <p>Report \\(D_{\\mathrm{mi}^2} = 2.59\\,D_{\\mathrm{km}^2}\\).</p>
                <h3>Variables</h3>
                <ul>
                  <li>\\(y_j\\) — number of <strong>days with ≥1 deer detection</strong> for camera \\(j\\).</li>
                  <li>\\(\\mathrm{camera\\_days}_j\\) — total deployed days for camera \\(j\\) from <code>Start/End</code>.</li>
                  <li>\\(\\mathrm{time\\_unit}\\) — expected time to traverse the viewshed once; \\(\\approx 0.59\\,r_m/v\\) (days).</li>
                  <li>\\(v\\) — movement speed (km/day); prior on \\(\\log v\\) as below.</li>
                  <li>\\(r_m\\) — effective detection radius (m) from <code>Detection Distance</code>.</li>
                  <li>\\(\\theta_{\\mathrm{rad}}\\) — camera half‑angle in radians (default \\(55^\\circ\\)).</li>
                  <li>\\(A_j\\) — viewshed area for camera \\(j\\) (km²).</li>
                  <li>\\(\\epsilon_j\\) — camera random effect (mean 0; optional).</li>
                  <li>\\(D_{\\mathrm{mi}^2}\\) — density (deer/mi²), reported.</li>
                </ul>
                <h3>Priors/inputs used in the app</h3>
                <ul>
                  <li>\\(\\log v \\sim \\mathcal{N}(1.339,\\,0.2955)\\) (used to define the time‑unit; calibrates the clock).</li>
                  <li>\\(\\theta = 55^\\circ\\) by default; \\(r_m\\) from <code>Detection Distance</code>.</li>
                  <li>\\(\\mathrm{sd}_\\epsilon \\sim \\mathcal{U}(0,10)\\) (optional camera heterogeneity).</li>
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
          div(
            actionButton("run_tte_sim", "Run TTE on simulated data", class = "btn-primary"),
            actionButton("stop_tte_sim", "Stop", class = "btn-danger", style = "margin-left: 10px;"),
            style = "margin-bottom: 10px;"
          ),
          br(),
          verbatimTextOutput("tte_sim_text"),
          
          hr(),
          h4("NPS data"),
          div(
            actionButton("run_tte_nps", "Run TTE on NPS data", class = "btn-primary"),
            actionButton("stop_tte_nps", "Stop", class = "btn-danger", style = "margin-left: 10px;"),
            style = "margin-bottom: 10px;"
          ),
          br(),
          verbatimTextOutput("tte_nps_text"),
          tags$div(
            style = "text-align: center; margin-top: 40px; padding-top: 20px; border-top: 1px solid #e0e0e0;",
            tags$h3(style = "margin: 0; color: var(--rw3);", "DEER App"),
            tags$p(style = "margin: 5px 0 0 0; color: var(--muted);", "Density Estimation from encounter rates")
          )
        ),
        
        # ---------------------- COMPARE & COMBINE --------------------
        nav_panel(
          "Compare & combine",
          markdown(paste(
            "Uses **WAIC-based model averaging** from the three models.",
            "",
            "For each dataset (simulated or NPS):",
            "",
            "1. Compute ΔWAIC and WAIC weights for REM, TTE, USCR;",
            "2. Combine posterior draws of `D_mi²` (deer per square mile);",
            "3. Report unweighted and WAIC-weighted mean densities, 95% CIs, and P(D > 20 DPSM).",
            "",
            "Run USCR, REM, and TTE from their tabs first, then check tables below.",
            sep = "\n"
          )),
          
          h4("Simulated data – WAIC-weighted results (deer/mi²)"),
          DTOutput("sim_combo_table"),
          
          hr(),
          h4("NPS data – WAIC-weighted results (deer/mi²)"),
          DTOutput("nps_combo_table")
        )
  )
)

# -------------------------------------------------------------------
# SERVER
# -------------------------------------------------------------------

server <- function(input, output, session) {
  
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
    build_sim_data_for_nimble(
      ch                 = sim()$ch,
      detection_radius_m = input$r_m
    )
  })
  
  output$grid_plot <- renderPlot({
    req(sim())
    traps_df <- as.data.frame(secr::traps(sim()$ch)) %>%
      mutate(camera = row_number())
    mask_df <- as.data.frame(sim()$mask)
    
    ggplot() +
      geom_point(data = mask_df, aes(x = x, y = y), alpha = 0.1, color = redwood_colors[5]) +
      geom_point(data = traps_df, aes(x = x, y = y), size = 3, color = redwood_colors[3]) +
      geom_text(data = traps_df, aes(x = x, y = y, label = camera),
                nudge_y = 15, size = 3, color = redwood_colors[1]) +
      coord_equal() +
      theme_minimal() +
      labs(
        x = "m", y = "m",
        title = sprintf(
          "%d cameras, spacing %dm (simulated density: %.1f deer/km²)",
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
        total_deer = d$camera_counts
      )
    
    ggplot(deer_dist, aes(x = x, y = y)) +
      geom_point(aes(size = total_deer),
                 color = redwood_colors[3], alpha = 0.7) +
      ggrepel::geom_text_repel(
        aes(label = Site),
        size = 3,
        max.overlaps = Inf
      ) +
      scale_size_continuous(name = "Total deer") +
      coord_equal() +
      labs(
        title = "Total Deer Detections by Camera",
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
      show_col_types = FALSE
    )
    
    images_raw <- clean_images_import(images_raw)
    
    msgs  <- character()
    warns <- character()
    
    checked <- NULL
    
    res <- tryCatch(
      {
        withCallingHandlers(
          {
            checked <- check_images(
              images      = images_raw,
              deployments = deployment_checked(),
              survey_year = input$survey_year
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
    
    # Trim to 56 days
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
  
  image_summary <- reactive({
    req(images_checked())
    imgs <- images_checked()
    imgs_std <- standardize_deer_species(imgs)
    summarize_images_by_species(imgs_std)
  })
  
  deer_objects <- reactive({
    req(images_checked())
    imgs <- images_checked()
    imgs_std <- standardize_deer_species(imgs)
    
    list(
      deer_summary = deer_summary_per_site(imgs_std),
      deer_counts  = deer_counts_per_camera(imgs_std),
      daily_deer   = deer_daily_detections(imgs_std)
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
    deer <- deer_objects()$deer_summary
    validate(
      need(nrow(deer) > 0, "No 'Deer' images found in dataset.")
    )
    DT::datatable(deer, options = list(pageLength = 10))
  })
  
  output$deer_bubble_plot <- renderPlot({
    dc <- deer_objects()$deer_counts
    validate(
      need(nrow(dc) > 0, "No 'Deer' images found in dataset.")
    )
    
    ggplot(dc, aes(x = Longitude, y = Latitude)) +
      geom_point(aes(size = total_deer),
                 color = redwood_colors[3], alpha = 0.7) +
      ggrepel::geom_text_repel(
        aes(label = `Site Name`),
        size = 3,
        max.overlaps = Inf
      ) +
      scale_size_continuous(name = "Total deer") +
      coord_fixed() +
      labs(
        title = "Total Deer Detections by Camera",
        x = "Longitude",
        y = "Latitude"
      ) +
      theme_minimal() +
      theme(legend.position = "right")
  })
  
  output$deer_daily_plot <- renderPlot({
    dd <- deer_objects()$daily_deer
    validate(
      need(nrow(dd) > 0, "No 'Deer' images found in dataset.")
    )
    
    ggplot(dd, aes(x = Date, y = `Site Name`, size = detections)) +
      geom_point(
        color = redwood_colors[3], alpha = 0.7,
        position = position_jitter(width = 0, height = 0.1)
      ) +
      scale_size_continuous(name = "Deer count", range = c(2, 8)) +
      labs(
        title = "Daily Deer Detections by Camera",
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
  
  uscr_nps_fit <- reactiveVal(NULL)
  rem_nps_fit  <- reactiveVal(NULL)
  tte_nps_fit  <- reactiveVal(NULL)
  
  # Status tracking for model runs
  uscr_sim_running <- reactiveVal(FALSE)
  uscr_nps_running <- reactiveVal(FALSE)
  rem_nps_running <- reactiveVal(FALSE)
  tte_nps_running <- reactiveVal(FALSE)
  rem_sim_running <- reactiveVal(FALSE)
  tte_sim_running <- reactiveVal(FALSE)
  
  # Stop flags for interrupting model runs
  stop_uscr_sim <- reactiveVal(FALSE)
  stop_uscr_nps <- reactiveVal(FALSE)
  stop_rem_sim <- reactiveVal(FALSE)
  stop_rem_nps <- reactiveVal(FALSE)
  stop_tte_sim <- reactiveVal(FALSE)
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
  
  observeEvent(input$stop_rem_sim, {
    stop_rem_sim(TRUE)
    stop_flags_env$stop_rem_sim <- TRUE
    showNotification("Stopping REM (sim) run...", 
                     type = "warning", duration = 3)
  })
  
  observeEvent(input$stop_rem_nps, {
    stop_rem_nps(TRUE)
    stop_flags_env$stop_rem_nps <- TRUE
    showNotification("Stopping REM (NPS) run...", 
                     type = "warning", duration = 3)
  })
  
  observeEvent(input$stop_tte_sim, {
    stop_tte_sim(TRUE)
    stop_flags_env$stop_tte_sim <- TRUE
    showNotification("Stopping TTE (sim) run...", 
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
  
  # --- USCR: simulated ---
  
  observeEvent(input$run_uscr_sim, {
    req(sim_data())
    d <- sim_data()
    
    # Reset stop flag
    stop_uscr_sim(FALSE)
    stop_flags_env$stop_uscr_sim <- FALSE
    uscr_sim_running(TRUE)
    uscr_sim_fit(NULL)  # Clear previous results
    
    showNotification("Running USCR on simulated data... This may take 30-60 minutes.", 
                     type = "message", duration = 10)
    
    fit <- tryCatch(
      {
        # Check if stopped before starting
        check_stop_flag("stop_uscr_sim")
        
        withProgress(
          message = "Running USCR model",
          detail = "Compiling model and running MCMC chains...",
          value = 0,
          {
            setProgress(0.1, detail = "Compiling NIMBLE model...")
            
            # Check stop flag again before running
            check_stop_flag("stop_uscr_sim")
            
            fit_result <- run_USCR(
              out            = d$out,
              camera_counts  = d$camera_counts,
              camera_days    = d$camera_days,
              iter           = input$iter_uscr,
              burnin         = input$burnin_uscr,
              thin           = input$thin_uscr,
              n_chains       = input$n_chains,
              M              = input$M_uscr,
              log_sigma_mean = input$log_sigma_mean,
              log_sigma_sd   = input$log_sigma_sd,
              log_lam0_sd    = input$log_lam0_sd,
              sd_eps_shape   = input$sd_eps_shape
            )
            setProgress(1.0, detail = "Complete!")
            fit_result
          }
        )
      },
      error = function(e) {
        if (stop_uscr_sim() || grepl("stopped by user", e$message, ignore.case = TRUE)) {
          showNotification("USCR (sim) run stopped by user.", type = "warning")
        } else {
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
      showNotification("USCR (sim) complete!", type = "message")
    }
  })
  
  # --- USCR: NPS ---
  
  observeEvent(input$run_uscr_nps, {
    req(nps_model_inputs())
    d <- nps_model_inputs()
    validate(
      need(all(d$camera_days > 0),
           "Some cameras have zero camera-days; check deployment dates.")
    )
    
    # Reset stop flag
    stop_uscr_nps(FALSE)
    uscr_nps_running(TRUE)
    uscr_nps_fit(NULL)  # Clear previous results
    
    showNotification("Running USCR on NPS data... This may take 30-60 minutes.", 
                     type = "message", duration = 10)
    
    fit <- tryCatch(
      {
        # Check if stopped before starting
        if (stop_uscr_nps()) {
          showNotification("USCR (NPS) run cancelled.", type = "warning")
          return(NULL)
        }
        
        withProgress(
          message = "Running USCR model",
          detail = "Compiling model and running MCMC chains...",
          value = 0,
          {
            setProgress(0.1, detail = "Compiling NIMBLE model...")
            
            # Check stop flag again before running
            if (stop_uscr_nps()) {
              showNotification("USCR (NPS) run cancelled.", type = "warning")
              return(NULL)
            }
            
            fit_result <- run_USCR(
              out            = d$out,
              camera_counts  = d$camera_counts,
              camera_days    = d$camera_days,
              iter           = input$iter_uscr,
              burnin         = input$burnin_uscr,
              thin           = input$thin_uscr,
              n_chains       = input$n_chains,
              M              = input$M_uscr,
              log_sigma_mean = input$log_sigma_mean,
              log_sigma_sd   = input$log_sigma_sd,
              log_lam0_sd    = input$log_lam0_sd,
              sd_eps_shape   = input$sd_eps_shape
            )
            setProgress(1.0, detail = "Complete!")
            fit_result
          }
        )
      },
      error = function(e) {
        if (stop_uscr_nps()) {
          showNotification("USCR (NPS) run stopped by user.", type = "warning")
        } else {
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
      showNotification("USCR (NPS) complete!", type = "message")
    }
  })
  
  # --- REM: simulated ---
  
  observeEvent(input$run_rem_sim, {
    req(sim_data())
    d <- sim_data()
    
    # Reset stop flag
    stop_rem_sim(FALSE)
    rem_sim_running(TRUE)
    
    showNotification("Running REM on simulated data...", type = "message")
    fit <- tryCatch(
      {
        # Check if stopped before starting
        if (stop_rem_sim()) {
          showNotification("REM (sim) run cancelled.", type = "warning")
          return(NULL)
        }
        
        run_REM(
          y            = d$camera_counts,
          r_km         = d$out$`Detection Distance` / 1000,
          camera_days  = d$camera_days,
          theta_deg    = input$theta,
          iter         = input$iter_rem_tte,
          burnin       = input$burnin_rem_tte,
          thin         = input$thin_rem_tte,
          n_chains     = input$n_chains,
          D_max        = input$D_max,
          log_v_mean   = input$log_v_mean,
          log_v_sd     = input$log_v_sd,
          sd_eps_max   = input$sd_eps_max
        )
      },
      error = function(e) {
        if (stop_rem_sim()) {
          showNotification("REM (sim) run stopped by user.", type = "warning")
        } else {
          showNotification(
            paste("REM (sim) failed:", e$message),
            type = "error", duration = NULL
          )
        }
        return(NULL)
      },
      finally = {
        rem_sim_running(FALSE)
        stop_rem_sim(FALSE)  # Reset stop flag
      }
    )
    rem_sim_fit(fit)
    if (!is.null(fit) && !stop_rem_sim()) {
      showNotification("REM (sim) complete.", type = "message")
    }
  })
  
  # --- REM: NPS ---
  
  observeEvent(input$run_rem_nps, {
    req(nps_model_inputs())
    d <- nps_model_inputs()
    validate(
      need(all(d$camera_days > 0),
           "Some cameras have zero camera-days; check deployment dates.")
    )
    
    # Reset stop flag
    stop_rem_nps(FALSE)
    rem_nps_running(TRUE)
    rem_nps_fit(NULL)  # Clear previous results
    
    showNotification("Running REM on NPS data... This may take several minutes.", 
                     type = "message", duration = 10)
    
    fit <- tryCatch(
      {
        # Check if stopped before starting
        if (stop_rem_nps()) {
          showNotification("REM (NPS) run cancelled.", type = "warning")
          return(NULL)
        }
        
        withProgress(
          message = "Running REM model",
          detail = "Compiling model and running MCMC chains...",
          value = 0,
          {
            setProgress(0.1, detail = "Compiling NIMBLE model...")
            
            # Check stop flag again before running
            if (stop_rem_nps()) {
              showNotification("REM (NPS) run cancelled.", type = "warning")
              return(NULL)
            }
            
            fit_result <- run_REM(
              y            = d$camera_counts,
              r_km         = d$out$`Detection Distance` / 1000,
              camera_days  = d$camera_days,
              theta_deg    = input$theta,
              iter         = input$iter_rem_tte,
              burnin       = input$burnin_rem_tte,
              thin         = input$thin_rem_tte,
              n_chains     = input$n_chains,
              D_max        = input$D_max,
              log_v_mean   = input$log_v_mean,
              log_v_sd     = input$log_v_sd,
              sd_eps_max   = input$sd_eps_max
            )
            setProgress(1.0, detail = "Complete!")
            fit_result
          }
        )
      },
      error = function(e) {
        if (stop_rem_nps()) {
          showNotification("REM (NPS) run stopped by user.", type = "warning")
        } else {
          showNotification(
            paste("REM (NPS) failed:", e$message),
            type = "error", duration = NULL
          )
        }
        return(NULL)
      },
      finally = {
        rem_nps_running(FALSE)
        stop_rem_nps(FALSE)  # Reset stop flag
      }
    )
    rem_nps_fit(fit)
    if (!is.null(fit) && !stop_rem_nps()) {
      showNotification("REM (NPS) complete!", type = "message")
    }
  })
  
  # --- TTE: simulated ---
  
  observeEvent(input$run_tte_sim, {
    req(sim_data())
    d <- sim_data()
    
    # Reset stop flag
    stop_tte_sim(FALSE)
    tte_sim_running(TRUE)
    
    showNotification("Running TTE on simulated data...", type = "message")
    fit <- tryCatch(
      {
        # Check if stopped before starting
        if (stop_tte_sim()) {
          showNotification("TTE (sim) run cancelled.", type = "warning")
          return(NULL)
        }
        
        run_TTE(
          y            = d$camera_counts,
          r_km         = d$out$`Detection Distance` / 1000,
          camera_days  = d$camera_days,
          theta_deg    = input$theta,
          iter         = input$iter_rem_tte,
          burnin       = input$burnin_rem_tte,
          thin         = input$thin_rem_tte,
          n_chains     = input$n_chains,
          D_max        = input$D_max,
          log_v_mean   = input$log_v_mean,
          log_v_sd     = input$log_v_sd,
          sd_eps_max   = input$sd_eps_max
        )
      },
      error = function(e) {
        if (stop_tte_sim()) {
          showNotification("TTE (sim) run stopped by user.", type = "warning")
        } else {
          showNotification(
            paste("TTE (sim) failed:", e$message),
            type = "error", duration = NULL
          )
        }
        return(NULL)
      },
      finally = {
        tte_sim_running(FALSE)
        stop_tte_sim(FALSE)  # Reset stop flag
      }
    )
    tte_sim_fit(fit)
    if (!is.null(fit) && !stop_tte_sim()) {
      showNotification("TTE (sim) complete.", type = "message")
    }
  })
  
  # --- TTE: NPS ---
  
  observeEvent(input$run_tte_nps, {
    req(nps_model_inputs())
    d <- nps_model_inputs()
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
          detail = "Compiling model and running MCMC chains...",
          value = 0,
          {
            setProgress(0.1, detail = "Compiling NIMBLE model...")
            
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
              n_chains     = input$n_chains,
              D_max        = input$D_max,
              log_v_mean   = input$log_v_mean,
              log_v_sd     = input$log_v_sd,
              sd_eps_max   = input$sd_eps_max
            )
            setProgress(1.0, detail = "Complete!")
            fit_result
          }
        )
      },
      error = function(e) {
        if (stop_tte_nps()) {
          showNotification("TTE (NPS) run stopped by user.", type = "warning")
        } else {
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
      showNotification("TTE (NPS) complete!", type = "message")
    }
  })
  
  # ================================================================
  # MODEL TAB OUTPUTS (USCR / REM / TTE)
  # ================================================================
  
  # USCR summaries
  output$uscr_sim_text <- renderPrint({
    if (uscr_sim_running()) {
      cat("⏳ USCR model is running...\n")
      cat("This may take 30-60 minutes depending on your settings.\n")
      cat("Please wait for the progress bar to complete.\n")
      return(invisible(NULL))
    }
    fit <- uscr_sim_fit()
    if (is.null(fit)) {
      cat("USCR (simulated): not run yet. Click 'Run USCR on simulated data'.")
      return(invisible(NULL))
    }
    s <- summarize_method(fit)
    list(
      dataset                   = "Simulated",
      mean_density_deer_per_km2 = round(s$mean_km2, 2),
      CI95_km2                  = c(round(s$q2.5_km2, 2), round(s$q97.5_km2, 2)),
      mean_density_deer_per_mi2 = round(s$mean_mi2, 2),
      CI95_mi2                  = c(round(s$q2.5_mi2, 2), round(s$q97.5_mi2, 2)),
      WAIC                      = round(s$waic, 2),
      note                      = if (!is.null(sim()))
        sprintf("True simulated density = %.1f deer/km²", sim()$truth$D_per_km2)
    )
  })
  
  output$uscr_nps_text <- renderPrint({
    if (uscr_nps_running()) {
      cat("⏳ USCR model is running...\n")
      cat("This may take 30-60 minutes depending on your settings.\n")
      cat("Please wait for the progress bar to complete.\n")
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
      mean_density_deer_per_km2 = round(s$mean_km2, 2),
      CI95_km2                  = c(round(s$q2.5_km2, 2), round(s$q97.5_km2, 2)),
      mean_density_deer_per_mi2 = round(s$mean_mi2, 2),
      CI95_mi2                  = c(round(s$q2.5_mi2, 2), round(s$q97.5_mi2, 2)),
      WAIC                      = round(s$waic, 2)
    )
  })
  
  # REM summaries
  output$rem_sim_text <- renderPrint({
    fit <- rem_sim_fit()
    if (is.null(fit)) {
      cat("REM (simulated): not run yet. Click 'Run REM on simulated data'.")
      return(invisible(NULL))
    }
    s <- summarize_method(fit)
    list(
      dataset                   = "Simulated",
      mean_density_deer_per_km2 = round(s$mean_km2, 2),
      CI95_km2                  = c(round(s$q2.5_km2, 2), round(s$q97.5_km2, 2)),
      mean_density_deer_per_mi2 = round(s$mean_mi2, 2),
      CI95_mi2                  = c(round(s$q2.5_mi2, 2), round(s$q97.5_mi2, 2)),
      WAIC                      = round(s$waic, 2)
    )
  })
  
  output$rem_nps_text <- renderPrint({
    if (rem_nps_running()) {
      cat("⏳ REM model is running...\n")
      cat("This may take 5-15 minutes depending on your settings.\n")
      cat("Please wait for the progress bar to complete.\n")
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
      mean_density_deer_per_km2 = round(s$mean_km2, 2),
      CI95_km2                  = c(round(s$q2.5_km2, 2), round(s$q97.5_km2, 2)),
      mean_density_deer_per_mi2 = round(s$mean_mi2, 2),
      CI95_mi2                  = c(round(s$q2.5_mi2, 2), round(s$q97.5_mi2, 2)),
      WAIC                      = round(s$waic, 2)
    )
  })
  
  # TTE summaries
  output$tte_sim_text <- renderPrint({
    fit <- tte_sim_fit()
    if (is.null(fit)) {
      cat("TTE (simulated): not run yet. Click 'Run TTE on simulated data'.")
      return(invisible(NULL))
    }
    s <- summarize_method(fit)
    list(
      dataset                   = "Simulated",
      mean_density_deer_per_km2 = round(s$mean_km2, 2),
      CI95_km2                  = c(round(s$q2.5_km2, 2), round(s$q97.5_km2, 2)),
      mean_density_deer_per_mi2 = round(s$mean_mi2, 2),
      CI95_mi2                  = c(round(s$q2.5_mi2, 2), round(s$q97.5_mi2, 2)),
      WAIC                      = round(s$waic, 2)
    )
  })
  
  output$tte_nps_text <- renderPrint({
    if (tte_nps_running()) {
      cat("⏳ TTE model is running...\n")
      cat("This may take 5-15 minutes depending on your settings.\n")
      cat("Please wait for the progress bar to complete.\n")
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
      mean_density_deer_per_km2 = round(s$mean_km2, 2),
      CI95_km2                  = c(round(s$q2.5_km2, 2), round(s$q97.5_km2, 2)),
      mean_density_deer_per_mi2 = round(s$mean_mi2, 2),
      CI95_mi2                  = c(round(s$q2.5_mi2, 2), round(s$q97.5_mi2, 2)),
      WAIC                      = round(s$waic, 2)
    )
  })
  
  # ================================================================
  # COMPARE & COMBINE (WAIC-based)
  # ================================================================
  
  sim_combo <- reactive({
    req(rem_sim_fit(), tte_sim_fit(), uscr_sim_fit())
    build_combo_table_from_fits(
      list(
        REM  = rem_sim_fit(),
        TTE  = tte_sim_fit(),
        USCR = uscr_sim_fit()
      )
    )
  })
  
  nps_combo <- reactive({
    if (is.null(rem_nps_fit()) || is.null(tte_nps_fit()) || is.null(uscr_nps_fit())) {
      return(NULL)
    }
    build_combo_table_from_fits(
      list(
        REM  = rem_nps_fit(),
        TTE  = tte_nps_fit(),
        USCR = uscr_nps_fit()
      )
    )
  })
  
  output$sim_combo_table <- renderDT({
    combo <- sim_combo()
    df <- combo$table %>%
      mutate(across(where(is.numeric), ~ round(.x, 3)))
    datatable(df, options = list(pageLength = 5, dom = "t"))
  })
  
  output$nps_combo_table <- renderDT({
    combo <- nps_combo()
    if (is.null(combo)) {
      return(DT::datatable(
        data.frame(
          Note = "For NPS data: run REM, TTE, and USCR from their tabs first."
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
