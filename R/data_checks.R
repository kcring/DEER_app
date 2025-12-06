# R/data_checks.R
# Cleaning, QC, trimming, and summary helpers for NPS camera data

# -------------------------------------------------------------------
# PRE-QC CLEANING HELPERS
# -------------------------------------------------------------------

clean_deployment_import <- function(deployments) {
  # Standardize column names: strip dots & weird whitespace, collapse spaces
  nm <- names(deployments)
  nm <- gsub("\\.+$", "", nm)
  nm <- gsub("[\\s\\p{Z}]+$", "", nm, perl = TRUE)
  nm <- gsub("^[\\s\\p{Z}]+", "", nm, perl = TRUE)
  nm <- gsub("[\\s\\p{Z}]+", " ", nm, perl = TRUE)
  nm <- trimws(nm)
  names(deployments) <- nm
  
  # Trim leading/trailing whitespace from all character columns
  deployments[] <- lapply(deployments, function(x) {
    if (is.character(x)) {
      gsub("^[\\s\\p{Z}]+|[\\s\\p{Z}]+$", "", x, perl = TRUE)
    } else x
  })
  
  deployments
}

clean_images_import <- function(images) {
  # Check for NULL or empty input
  if (is.null(images)) {
    stop("clean_images_import: images is NULL")
  }
  if (!is.data.frame(images)) {
    stop("clean_images_import: images must be a data frame")
  }
  if (nrow(images) == 0) {
    warning("clean_images_import: images data frame is empty")
    return(images)
  }
  
  # Standardize column names
  nm <- names(images)
  if (is.null(nm) || length(nm) == 0) {
    stop("clean_images_import: images has no column names")
  }
  nm <- gsub("\\.+$", "", nm)
  nm <- trimws(nm)
  names(images) <- nm
  
  # Strip explicit " UTC" suffix if present in Timestamp
  if ("Timestamp" %in% names(images)) {
    images$Timestamp <- sub(" UTC$", "", images$Timestamp)
  }
  
  # Split multi-species rows like "Deer|Raccoon" into separate rows
  if (all(c("Species", "Sighting Count") %in% names(images))) {
    images <- images %>%
      dplyr::mutate(`Sighting Count` = as.character(`Sighting Count`)) %>%
      tidyr::separate_rows(Species, `Sighting Count`, sep = "\\|") %>%
      dplyr::mutate(`Sighting Count` = suppressWarnings(as.numeric(`Sighting Count`)))
  }
  
  images
}

# -------------------------------------------------------------------
# QC: check_deployments (as before)
# -------------------------------------------------------------------

check_deployments <- function(deployment) {
  issues <- list()
  
  # ---- Column presence ----
  required_cols <- c(
    "Park", "Site Name", "Camera ID", "SD Card ID", 
    "Start Date", "Start Time", "End Date", "End Time", 
    "Latitude", "Longitude", "Camera Height", 
    "Camera Orientation", "Camera Functioning", 
    "Camera Malfunction Date", "Detection Distance", "Notes"
  )
  
  missing_cols <- setdiff(required_cols, names(deployment))
  if (length(missing_cols) > 0) {
    cat(
      "❌ Missing required column(s):", paste(missing_cols, collapse = ", "), "\n\n",
      "💡 Guidance:\n",
      "• If these columns are truly missing from your dataset (for example, no Start Time or End Time was ever recorded),\n",
      "  please **add the column names** to your CSV file and leave the cells blank.\n",
      "• If you believe these columns exist but are labeled slightly differently (e.g., 'Start_Date' instead of 'Start Date'),\n",
      "  please rename them to match exactly the expected column names above.\n\n",
      "After fixing the file, re-run `check_deployments()`.\n"
    )
    stop("Deployment file is missing one or more required columns.")
  }
  
  # ---- Auto-fix Site Names missing leading zero ----
  deployment$`Site Name` <- mapply(function(sn, park) {
    if (is.na(sn) || sn == "") return(sn)
    
    # Match pattern like HOFU_1 or MORR_9 (no leading zero)
    if (grepl(paste0("^", park, "_[0-9]$"), sn)) {
      num <- sub(".*_", "", sn)
      corrected <- paste0(park, "_0", num)
      message("🛠 Fixed Site Name: ", sn, " → ", corrected)
      return(corrected)
    } else {
      return(sn)
    }
  }, deployment$`Site Name`, deployment$Park, USE.NAMES = FALSE)
  
  # --- Site Name validation ---
  for (i in seq_len(nrow(deployment))) {
    park <- deployment$Park[i]
    sn   <- trimws(deployment$`Site Name`[i]) 
    
    if (!is.na(sn) && sn != "") {
      # Suffix: 01–09, 10–99, or 100–999
      re_suffix <- "(0[1-9]|[1-9][0-9]|[1-9][0-9]{2})$"
      
      # Build patterns
      re_park_only        <- paste0("^", park, "_", re_suffix)              # FRSP_01
      re_unit_only        <- paste0("^[A-Za-z]{2,}_", re_suffix)            # WILD_01
      re_park_unit_us     <- paste0("^", park, "_[A-Za-z]{2,}_", re_suffix) # FRSP_WILD_01
      re_park_unit_concat <- paste0("^", park, "[A-Za-z]{2,}_", re_suffix)  # FRSPWILD_01
      
      patterns <- c(re_park_only, re_unit_only, re_park_unit_us, re_park_unit_concat)
      
      is_valid <- any(vapply(
        patterns,
        function(p) grepl(p, sn, ignore.case = TRUE),
        logical(1)
      ))
      
      if (!is_valid) {
        issues <- c(issues, paste0(
          "Invalid Site_Name: ", sn,
          " → must match one of: 'PARK_##/PARK_###' (e.g., FRSP_09 or FRSP_101), ",
          "'UNIT_##/UNIT_###' (e.g., WILD_01), ",
          "'PARK_UNIT_##' (e.g., FRSP_WILD_01), or 'PARKUNIT_##' (e.g., FRSPWILD_01)."
        ))
      }
    }
  }
  
  # ---- Per-row checks ----
  deployment_cols <- setdiff(required_cols, c("Park", "Site Name", "Notes"))
  
  for (i in seq_len(nrow(deployment))) {
    row_values <- deployment[i, deployment_cols]
    
    # Skip row if all deployment columns except Park, Site Name, Notes are blank
    if (all(is.na(row_values) | row_values == "")) next
    
    cam <- deployment$`Camera ID`[i]
    
    # --- Camera ID ---
    if (is.na(cam) || cam == "") {
      issues <- c(issues, paste("❌ Missing Camera ID in row", i, " — must have a value"))
    }
    
    # --- SD Card ID ---
    sd_val <- deployment$`SD Card ID`[i]
    if (is.na(sd_val) || sd_val == "") {
      issues <- c(issues, paste("❌ Missing SD Card ID in row", i, " — must have a value"))
    }
    
    # --- Dates ---
    date_cols <- c("Start Date", "End Date")
    for (col in date_cols) {
      val <- deployment[[col]][i]
      if (!is.na(val) && val != "") {
        if (is.na(as.Date(val, format = "%m/%d/%Y"))) {
          issues <- c(issues, paste("❌ Bad date in", col, "row", i, ":", val, " — should be mm/dd/yyyy"))
        }
      }
    }
    
    # Camera Malfunction Date only required if Camera Functioning == "No"
    cam_func <- deployment$`Camera Functioning`[i]
    if (!is.na(cam_func) && tolower(cam_func) == "no") {
      val <- deployment$`Camera Malfunction Date`[i]
      if (is.na(val) || val == "") {
        issues <- c(issues, paste("❌ Camera Malfunction Date missing in row", i, 
                                  " — required because Camera Functioning = No"))
      } else if (is.na(as.Date(val, format = "%m/%d/%Y"))) {
        issues <- c(issues, paste("❌ Bad date in Camera Malfunction Date row", i, ":", val, 
                                  " — should be mm/dd/yyyy"))
      }
    }
    
    # --- Times ---
    time_cols <- c("Start Time", "End Time")
    for (col in time_cols) {
      val <- deployment[[col]][i]
      if (!is.na(val) && val != "") {
        if (!grepl("^(?:[01]?[0-9]|2[0-3]):[0-5][0-9](:[0-5][0-9])?$", val)) {
          issues <- c(issues, paste("❌ Bad time in", col, "row", i, ":", val, " — should be HH:MM or HH:MM:SS 24h"))
        }
      }
    }
    
    # ---- Numeric checks with suppression ----
    num_cols <- c("Camera Height", "Detection Distance")
    for (col in num_cols) {
      val <- deployment[[col]][i]
      if (!is.na(val) && val != "" && suppressWarnings(is.na(as.numeric(val)))) {
        issues <- c(issues, paste("❌ Non-numeric value in", col, "row", i, ":", val))
      }
      
      if (col == "Longitude" && !is.na(val) && val != "" && as.numeric(val) >= 0) {
        issues <- c(issues, paste("❌ Longitude in row", i, "is not negative:", val,
                                  " — must have a minus sign for western hemisphere"))
      }
    }
    
    # --- Latitude check (must exist and be non-zero) ---
    lat_val <- suppressWarnings(as.numeric(deployment$Latitude[i]))
    if (is.na(lat_val) || lat_val == 0) {
      issues <- c(issues, paste("❌ Latitude missing or invalid (0) in row", i))
    }
    
    # --- Longitude check (must exist and be non-zero; auto-fix sign if positive) ---
    long_val <- suppressWarnings(as.numeric(deployment$Longitude[i]))
    if (is.na(long_val) || long_val == 0) {
      issues <- c(issues, paste("❌ Longitude missing or invalid (0) in row", i))
    } else if (long_val > 0) {
      # Auto-fix longitude sign
      fixed_val <- -abs(long_val)
      deployment$Longitude[i] <- fixed_val
      message("🛠 Fixed Longitude in row ", i, ": ", long_val, " → ", fixed_val)
    }
    
    # --- Camera Orientation ---
    val <- deployment$`Camera Orientation`[i]
    valid_cardinals <- c("N","NE","E","SE","S","SW","W","NW")
    if (!is.na(val) && val != "") {
      val_upper <- toupper(val)
      val_numeric <- suppressWarnings(as.numeric(val))
      if (!(val_upper %in% valid_cardinals | (!is.na(val_numeric) & val_numeric >= 0 & val_numeric <= 359))) {
        issues <- c(issues, paste("❌ Invalid Camera Orientation in row", i, ":", val,
                                  "— must be N/NE/.../NW or 0–359 degrees"))
      }
    }
    
    # --- Camera Functioning ---
    val <- deployment$`Camera Functioning`[i]
    if (!is.na(val) && val != "") {
      val_lower <- tolower(val)
      if (!(val_lower %in% c("yes","no"))) {
        issues <- c(issues, paste("❌ Invalid Camera Functioning value in row", i, ":", val, " — must be Yes or No"))
      }
    }
    
  } # end row loop
  
  # ---- Output ----
  if (length(issues) == 0) {
    message("✅ Deployments file is formatted correctly!")
  } else {
    warning("⚠️ Issues found in deployments file:\n", paste(issues, collapse = "\n"))
  }
  return(deployment)
}

# -------------------------------------------------------------------
# QC: check_images (as before)
# -------------------------------------------------------------------

check_images <- function(images, deployments, survey_year) {
  issues <- c()
  fixes  <- c()
  
  # ---- Required columns ----
  required_cols <- c("Site Name", "Latitude", "Longitude", "Timestamp",
                     "Species", "Sighting Count", "Image URL", "Cluster ID")
  missing_cols <- setdiff(required_cols, names(images))
  if (length(missing_cols) > 0) {
    stop(paste0(
      "❌ Missing required column(s): ",
      paste(missing_cols, collapse = ", "),
      "\nPlease fix the images file and re-run the function."
    ))
  }
  
  # ---- Trim whitespace from character columns ----
  char_cols <- names(images)[sapply(images, is.character)]
  images[char_cols] <- lapply(images[char_cols], trimws)
  
  # ---- Fix Site Names missing leading zero ----
  fix_idx <- which(grepl("^[A-Z]+_[0-9]$", images$`Site Name`))
  if (length(fix_idx) > 0) {
    corrected_names <- sub("^(.*_)([0-9])$", "\\10\\2", images$`Site Name`[fix_idx])
    fixes <- paste(images$`Site Name`[fix_idx], "→", corrected_names)
    images$`Site Name`[fix_idx] <- corrected_names
    message("🛠 Fixed Site Names:\n", paste(unique(fixes), collapse = "\n"))
  }
  
  # ---- Site Name validation ----
  valid_regex <- "^[A-Z]{4}_(0[1-9]|[1-9][0-9])$|^[A-Z]+_(0[1-9]|[1-9][0-9]|10[0-9]|[1-9][0-9]{2})$"
  invalid_idx <- which(!grepl(valid_regex, images$`Site Name`))
  if (length(invalid_idx) > 0) {
    site_counts <- table(images$`Site Name`[invalid_idx])
    for (sn in names(site_counts)) {
      issues <- c(issues, paste0(
        "❌ Invalid Site_Name: ", sn,
        " — found ", site_counts[[sn]], " occurrence(s); must follow 'PARK_##/PARK_###' or 'UNIT_##/UNIT_###'."
      ))
    }
  }
  
  # ---- Latitude / Longitude checks ----
  lat_num <- suppressWarnings(as.numeric(images$Latitude))
  lon_num <- suppressWarnings(as.numeric(images$Longitude))
  
  bad_lat <- which(is.na(lat_num) | lat_num == 0)
  if (length(bad_lat) > 0) {
    issues <- c(issues, paste0("❌ ", length(bad_lat), " image(s) have missing or zero Latitude"))
  }
  
  bad_lon <- which(is.na(lon_num) | lon_num == 0)
  if (length(bad_lon) > 0) {
    issues <- c(issues, paste0("❌ ", length(bad_lon), " image(s) have missing or zero Longitude"))
  }
  
  # Auto-fix positive longitudes
  fix_lon <- which(!is.na(lon_num) & lon_num > 0 & lon_num <= 180)
  if (length(fix_lon) > 0) {
    old_vals <- lon_num[fix_lon]
    new_vals <- -abs(old_vals)
    images$Longitude[fix_lon] <- new_vals
    message("🛠 Fixed Longitude for ", length(fix_lon), " image(s) ")
  }
  
  # ---- Required fields check ----
  for (field in c("Species", "Sighting Count", "Image URL")) {
    missing_rows <- which(is.na(images[[field]]) | images[[field]] == "")
    if (length(missing_rows) > 0) {
      issues <- c(issues, paste0("❌ Missing ", field, " in ", length(missing_rows), " image(s)"))
    }
  }
  
  # ---- Cluster ID numeric check ----
  cluster_num <- suppressWarnings(as.numeric(images$`Cluster ID`))
  bad_cluster <- which(!is.na(images$`Cluster ID`) & images$`Cluster ID` != "" & is.na(cluster_num))
  if (length(bad_cluster) > 0) {
    issues <- c(issues, paste0("❌ ", length(bad_cluster), " image(s) have non-numeric Cluster ID"))
  }
  
  # ---- Cross-check Site Names with deployments ----
  bad_site_match <- which(!images$`Site Name` %in% deployments$`Site Name`)
  if (length(bad_site_match) > 0) {
    issues <- c(issues, paste0("❌ ", length(bad_site_match), " image(s) have Site Names not found in deployments"))
  }
  
  # ---- Timestamp parsing and year check ----
  ts_parsed <- suppressWarnings(as.POSIXct(images$Timestamp,
                                           format = "%Y-%m-%d %H:%M:%S",
                                           tz = "UTC"))
  missing_time <- is.na(ts_parsed) & !is.na(images$Timestamp) & images$Timestamp != ""
  if (any(missing_time)) {
    ts_parsed[missing_time] <- suppressWarnings(
      as.POSIXct(
        paste0(images$Timestamp[missing_time], " 00:00:00"),
        format = "%Y-%m-%d %H:%M:%S",
        tz = "UTC"
      )
    )
  }
  bad_ts <- which(!is.na(images$Timestamp) & images$Timestamp != "" & is.na(ts_parsed))
  if (length(bad_ts) > 0) {
    site_counts <- table(images$`Site Name`[bad_ts])
    for (sn in names(site_counts)) {
      issues <- c(issues, paste0("❌ Bad Timestamp format at site ", sn,
                                 " — ", site_counts[[sn]], " occurrence(s); should be yyyy-mm-dd HH:MM:SS."))
    }
  }
  
  # Year check
  wrong_year_idx <- which(!is.na(ts_parsed) & format(ts_parsed, "%Y") != survey_year)
  if (length(wrong_year_idx) > 0) {
    site_years <- split(format(ts_parsed[wrong_year_idx], "%Y"),
                        images$`Site Name`[wrong_year_idx])
    for (sn in names(site_years)) {
      unique_years <- unique(site_years[[sn]])
      issues <- c(issues, paste0(
        "❌ Timestamp year mismatch at site ", sn,
        " — found year(s): ", paste(unique_years, collapse = ", "),
        "; expected year: ", survey_year, "."
      ))
    }
  }
  
  # ---- First image vs. Start Date check ----
  dt_images <- data.table::data.table(images)
  dt_images[, ts := ts_parsed]
  
  # Get first image per site
  first_img <- dt_images[!is.na(ts), .(first_img = min(ts)), by = `Site Name`]
  
  # Prepare deployment start dates
  dt_deploy <- data.table::data.table(deployments)
  dt_deploy[, dep_start := as.Date(`Start Date`, format = "%m/%d/%Y")]
  
  # Merge and compare
  merged_dates <- merge(first_img,
                        dt_deploy[, .(`Site Name`, dep_start)],
                        by = "Site Name",
                        all.x = TRUE)
  merged_dates[, first_img_date := as.Date(first_img)]
  
  mismatched_sites <- merged_dates[first_img_date != dep_start, `Site Name`]
  if (length(mismatched_sites) > 0) {
    issues <- c(issues, paste0(
      "❌ Date mismatch for ", length(mismatched_sites), " site(s): ",
      paste(mismatched_sites, collapse = ", "),
      " — first image date does not match deployment Start Date."
    ))
  }
  
  # ---- Final output ----
  if (length(issues) == 0) {
    message("✅ Images file is formatted correctly and all Site Names match deployments!")
  } else {
    warning("⚠️ Issues found in images file:\n", paste(issues, collapse = "\n"))
  }
  
  return(images)
}

# -------------------------------------------------------------------
# FORMAT DEPLOYMENTS (trim to 56 days if needed)
# -------------------------------------------------------------------

format_deployments <- function(deployments) {
  # Ensure dplyr is available
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package 'dplyr' is required for format_deployments()", call. = FALSE)
  }
  if (!requireNamespace("lubridate", quietly = TRUE)) {
    stop("Package 'lubridate' is required for format_deployments()", call. = FALSE)
  }
  
  # Parse dates with multiple format attempts (including 2-digit years)
  parse_date_safe <- function(x) {
    if (inherits(x, "Date")) return(x)
    if (inherits(x, "POSIXct")) return(as.Date(x))
    
    # Try lubridate first (handles 2-digit years automatically)
    parsed <- lubridate::parse_date_time(
      x,
      orders = c("mdy", "ymd", "m-d-y", "y-m-d", "mdY", "Ymd"),
      quiet = TRUE
    )
    
    # If still NA, try as.Date with common formats (including 2-digit year)
    if (any(is.na(parsed))) {
      parsed <- as.Date(x, tryFormats = c("%m/%d/%y", "%m/%d/%Y", "%Y-%m-%d", "%m-%d-%Y", "%d/%m/%Y"))
    }
    
    if (inherits(parsed, "POSIXct")) parsed <- as.Date(parsed)
    return(parsed)
  }
  
  deps <- deployments %>%
    dplyr::mutate(
      Start_Date = parse_date_safe(`Start Date`),
      End_Date   = parse_date_safe(`End Date`)
    )
  
  # If deployment > 56 days, trim End Date
  deps <- deps %>%
    dplyr::mutate(
      operational_days = as.numeric(difftime(End_Date, Start_Date, units = "days")),
      End_Date = dplyr::if_else(operational_days > 56, 
                                Start_Date + 56, 
                                End_Date),
      `End Date` = format(End_Date, "%m/%d/%Y")
    ) %>%
    dplyr::select(-Start_Date, -End_Date, -operational_days)
  
  return(deps)
}

# -------------------------------------------------------------------
# TRIM IMAGES TO 56 DAYS
# -------------------------------------------------------------------

trim_images_56days <- function(images) {
  # Ensure Timestamp is POSIXct
  if (!inherits(images$Timestamp, "POSIXt")) {
    images$Timestamp <- suppressWarnings(
      as.POSIXct(images$Timestamp, tz = "UTC")
    )
  }
  
  # Determine first timestamp per site
  first_times <- images %>%
    dplyr::group_by(`Site Name`) %>%
    dplyr::summarize(first_date = min(Timestamp, na.rm = TRUE), .groups = "drop")
  
  # Join to images to calculate days from first image
  images_check <- images %>%
    dplyr::left_join(first_times, by = "Site Name") %>%
    dplyr::mutate(days_from_start = as.numeric(difftime(Timestamp, first_date,
                                                        units = "days")))
  
  # Check if any rows exceed 56 days
  rows_to_remove <- images_check %>% dplyr::filter(days_from_start > 56)
  
  if (nrow(rows_to_remove) == 0) {
    message("✅ No images taken after 56 days. Dataset left unchanged.")
    return(images)
  } else {
    # Count rows removed per site
    removal_counts <- rows_to_remove %>%
      dplyr::group_by(`Site Name`) %>%
      dplyr::summarize(removed = n(), .groups = "drop")
    
    message("⚠️ Images beyond 56 days were removed:")
    for (i in seq_len(nrow(removal_counts))) {
      message(paste0("  Site ", removal_counts$`Site Name`[i], ": ",
                     removal_counts$removed[i], " row(s) removed"))
    }
    
    # Keep only rows within 56 days
    images_trimmed <- images_check %>%
      dplyr::filter(days_from_start <= 56) %>%
      dplyr::select(-first_date, -days_from_start)
    
    return(images_trimmed)
  }
}

# -------------------------------------------------------------------
# SUMMARY HELPERS (deployment, images, deer)
# -------------------------------------------------------------------

summarize_deployments <- function(deployments) {
  deployments %>%
    dplyr::mutate(
      Start_Date = as.Date(`Start Date`, format = "%m/%d/%Y"),
      End_Date   = as.Date(`End Date`,   format = "%m/%d/%Y"),
      operational_days = as.numeric(difftime(End_Date, Start_Date,
                                             units = "days"))
    ) %>%
    dplyr::summarise(
      n_cameras         = dplyr::n(),
      mean_days         = mean(operational_days, na.rm = TRUE),
      total_camera_days = sum(operational_days, na.rm = TRUE),
      .groups = "drop"
    )
}

summarize_images_by_species <- function(images) {
  images %>%
    dplyr::group_by(Species) %>%
    dplyr::summarise(
      total_images     = dplyr::n(),
      total_detections = sum(as.numeric(`Sighting Count`), na.rm = TRUE),
      cameras_detected = dplyr::n_distinct(`Site Name`),
      .groups = "drop"
    ) %>%
    dplyr::arrange(dplyr::desc(total_detections))
}

standardize_deer_species <- function(images) {
  if (any(tolower(images$Species) %in% c("white-tailed deer", "white tailed deer"))) {
    images <- images %>%
      dplyr::mutate(
        Species = dplyr::if_else(
          tolower(Species) %in% c("white-tailed deer", "white tailed deer"),
          "Deer",
          Species
        )
      )
    message("🦌 Standardized all 'white-tailed deer' classifications to 'Deer' for consistency.")
  } else {
    message("✅ No 'white-tailed deer' entries found — all deer species names already standardized.")
  }
  images
}

deer_summary_per_site <- function(images) {
  deer_images <- images %>%
    dplyr::filter(grepl("Deer", Species, ignore.case = TRUE))
  
  deer_images %>%
    dplyr::group_by(`Site Name`) %>%
    dplyr::summarise(
      deer_images     = dplyr::n(),
      deer_detections = sum(as.numeric(`Sighting Count`), na.rm = TRUE),
      .groups = "drop"
    )
}

deer_counts_per_camera <- function(images) {
  images %>%
    dplyr::filter(Species == "Deer") %>%
    dplyr::group_by(`Site Name`, Latitude, Longitude) %>%
    dplyr::summarise(
      total_deer = sum(as.numeric(`Sighting Count`), na.rm = TRUE),
      .groups = "drop"
    )
}

deer_daily_detections <- function(images) {
  images %>%
    dplyr::filter(Species == "Deer") %>%
    dplyr::mutate(Date = as.Date(Timestamp)) %>%
    dplyr::group_by(`Site Name`, Date) %>%
    dplyr::summarise(
      detections = sum(as.numeric(`Sighting Count`), na.rm = TRUE),
      .groups = "drop"
    )
}

