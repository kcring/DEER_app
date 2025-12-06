# DEER App

**Density Estimation from Encounter Rates** — A Shiny application for estimating deer density using unmarked camera trap methods.

![DEER App](www/deer_app_logo.png)

## Overview

The DEER App provides a user-friendly interface for running three Bayesian models (uSCR, REM, and TTE) to estimate deer density from camera trap data. The app supports both simulated data and real-world National Park Service (NPS) data formats, with model-averaged results using WAIC weights.

## Features

- **Three Bayesian Models:**
  - **uSCR** (Unmarked Spatial Capture–Recapture): Estimates deer density by learning from spatial detection patterns across a camera array
  - **REM** (Random Encounter Model): Converts encounter rates into density, correcting for movement speed and camera view geometry
  - **TTE** (Time-to-Event): Models the time between detections to estimate density

- **Data Options:**
  - Simulate camera trap data with customizable parameters
  - Upload NPS-format deployment and images CSV files

- **Model Averaging:**
  - WAIC-based model averaging across all three models
  - Weighted and unweighted density estimates with 95% credible intervals

- **Interactive Visualizations:**
  - Camera grid plots
  - Deer distribution maps
  - Interactive leaflet maps
  - Species summary plots
  - Daily detection time series

## Installation

### Prerequisites

- R (version 4.0 or higher)
- RStudio (recommended)

### Required R Packages

The app will automatically install missing packages, but you can also install them manually:

```r
install.packages(c(
  "shiny", "bslib", "DT", "ggplot2", "dplyr", "tidyr",
  "readr", "purrr", "stringr", "secr", "data.table",
  "leaflet", "ggrepel",
  "nimble", "nimbleHMC", "parallel", "MCMCvis", "lubridate",
  "calecopal"
))
```

**Note:** Some packages may require additional system dependencies:
- `nimble` and `nimbleHMC` require a C++ compiler
- `leaflet` requires system libraries for spatial data

### Running the App

1. Clone or download this repository
2. Open `app.R` in RStudio
3. Click "Run App" or run:
   ```r
   shiny::runApp()
   ```

Alternatively, from the command line:
```bash
Rscript -e "shiny::runApp()"
```

## Usage

### Step 1: Simulate or Upload Data

**Option 1: Simulate Data**
- Navigate to the "Simulate data" tab
- Adjust simulation parameters (grid size, spacing, days, density, etc.)
- Click "Simulate grid" to generate toy data

**Option 2: Upload NPS Data**
- Prepare two CSV files:
  - **Deployment file**: Camera deployment information
  - **Images file**: Detection records
- See the "Add your data" tab for required column specifications
- Upload files using the file input controls

### Step 2: Review Data Summary

- Check the "Data summary" tab to verify your data
- Review deployment summaries, species detections, and spatial distributions

### Step 3: Run Models

Navigate to the model tabs (USCR, REM, TTE):
- Separate "Run" buttons for simulated vs. NPS data
- Progress bars show model compilation and MCMC sampling status
- "Stop" buttons allow you to terminate long-running models
- Results appear below the buttons when complete

**Advanced Settings:**
- Click "Advanced" under "Settings mode" to adjust:
  - MCMC iterations, burnin, thinning
  - Number of chains
  - Prior parameters
  - Model-specific settings

### Step 4: Compare & Combine Results

- Use the "Compare & combine" tab to view:
  - WAIC values and weights for each model
  - Model-averaged density estimates (deer/mi²)
  - 95% credible intervals
  - Probability of exceeding threshold densities

## File Structure

```
deer_app_v2/
├── app.R                    # Main Shiny application
├── R/
│   ├── sim_and_models.R     # Model functions (USCR, REM, TTE)
│   └── data_checks.R        # Data validation and summary functions
├── www/
│   ├── deer_app_logo.png    # Main app logo
│   ├── wvu_logo.png         # WVU logo
│   └── nps_logo.png         # NPS logo
├── README.md                # This file
└── deer_app_v2.Rproj        # RStudio project file
```

## Data Format Requirements

### Deployment File Required Columns:
- **Site Name**: Camera location identifier
- **Latitude**: Decimal degrees
- **Longitude**: Decimal degrees
- **Start Date**: Camera deployment start (MM/DD/YYYY or MM/DD/YY)
- **End Date**: Camera deployment end (MM/DD/YYYY or MM/DD/YY)
- **Detection Distance**: Detection radius in meters

### Images File Required Columns:
- **Site Name**: Must match deployment file
- **Date**: Detection date (MM/DD/YYYY or MM/DD/YY)
- **Species**: Species name (deer species will be standardized)
- **Time**: Detection time (optional, for temporal analysis)

See the "Add your data" tab in the app for detailed column specifications.

## Model Details

### uSCR (Unmarked Spatial Capture–Recapture)
- **Pros:** Accounts for spatial detection patterns, handles camera heterogeneity
- **Cons:** Computationally intensive, requires spatial array design
- Estimates density by modeling activity centers and detection probability as a function of distance

### REM (Random Encounter Model)
- **Pros:** Fast, simple, works with single cameras
- **Cons:** Assumes random movement, requires movement speed estimates
- Converts encounter rates to density using movement speed and detection area

### TTE (Time-to-Event)
- **Pros:** Uses temporal information, accounts for detection heterogeneity
- **Cons:** Requires detection distance measurements, assumes Poisson process
- Models time between detections to estimate encounter rates and density

## Performance Notes

- **USCR** is the most computationally intensive model. For faster demos, reduce:
  - `M_uscr` (state-space size)
  - `iter_uscr` (MCMC iterations)
  - `n_chains` (number of chains)

- **REM** and **TTE** are generally faster but still benefit from reduced iterations for quick tests

- All models use parallel processing when multiple chains are specified

## Troubleshooting

**Models not running:**
- Check that all required packages are installed
- Verify data format matches requirements
- Ensure camera deployment dates are valid
- Check that detection distances are provided

**LaTeX equations not rendering:**
- The app uses MathJax for equation rendering
- If equations don't appear, try refreshing the page
- Check browser console for JavaScript errors

**Memory issues:**
- Reduce `M_uscr` for USCR model
- Reduce number of chains
- Close other R sessions

## Citation

If you use this app in your research, please cite the underlying methods:
- uSCR: [Chandler & Royle 2013]
- REM: [Rowcliffe et al. 2008]
- TTE: [Moeller et al. 2018]

## License

License to be determined.

## Acknowledgments

### Model Development
The underlying models and code were created by **Dr. Amanda Van Buskirk** under the advisement of **Dr. Christopher Rota** in the [**Rota Quantitative Ecology Lab**](https://sites.google.com/mix.wvu.edu/rotalab/home) within the **Davis College of Agriculture and Natural Resources at West Virginia University**.

### Shiny App Development
This Shiny application was developed as part of the **Science in the Parks Communications Fellowship**, a collaborative effort between the **Ecological Society of America (ESA)** and the **National Park Service (NPS)**. Learn more: [https://esa.org/programs/scip/](https://esa.org/programs/scip/)

**Fellowship Support:**
- **Dr. Brian Mitchell** (NPS) - Fellowship Liaison
- **Jasjeet Dhanota** (ESA) - Mentor
- **Mary Joy Mulumba** (ESA) - Mentor

### Technical Acknowledgments
- Built with R Shiny, NIMBLE, and the `calecopal` color palette
- West Virginia University (WVU)
- National Park Service (NPS)

## Contact

**Kacie Ring**  
University of California, Santa Barbara  
Website: [kaciering.com](https://kaciering.com)  
GitHub: [@kcring](https://github.com/kcring)

---

**DEER App** — Density Estimation from Encounter Rates  
uSCR · REM · TTE — unmarked camera methods, model‑averaged to deer/mi².

