# simulate_registry_full.R
# Full simulation per user spec (1e6 persons, diag_data with pre-AMI comorbs, mortality, etc.)
# Run locally. Uses data.table for performance.
set.seed(20251014)

# Dependencies
if (!requireNamespace("data.table", quietly=TRUE)) install.packages("data.table")
library(data.table)

# PARAMETERS
n_total <- 1e6L
n_AMI_m <- 2500L
n_AMI_f <- 1500L
study_start <- as.IDate("2010-01-01")
study_end   <- as.IDate("2024-12-31")
years <- 2010:2024

# comorb prevalences among AMI patients (as fractions)
comorb_prev_ami <- list(
  HF = c(male = 0.05, female = 0.06),
  Hyperchol = c(male = 0.29, female = 0.31),
  Diabetes = c(male = 0.15, female = 0.15),
  Hypertension = c(male = 0.27, female = 0.35)
)
# non-AMI prevalences = 0.40 * comorb_prev_ami
scale_nonami <- 0.40
comorb_prev_nonami <- lapply(comorb_prev_ami, function(x) x * scale_nonami)

# ICD pools
ICD_AMI <- "DI21"
ICD_HF <- "DI50"
ICD_HYPERCHOL <- "DE78"
ICD_DIABETES <- c("DE10","DE11","DE12","DE13","DE14")
ICD_HYPERT <- c("DI10","DI11","DI12","DI13","DI14","DI15")
other_diag <- c("DX01","DX02","DX03","DX04","DX05","DX06","DX07","DX08","DX09")

# AMI incidence decline weights (2010->2024, 30% lower in 2024 vs 2010)
year_weights <- seq(1.0, 0.7, length.out = length(years))
year_probs <- year_weights / sum(year_weights)

# LOS discrete distribution (median ~7, IQR ~3-11)
los_vals <- 1:30
los_probs <- exp(-abs(los_vals - 7)/3.0)
los_probs <- los_probs / sum(los_probs)
sample_los <- function(n) sample(los_vals, size = n, replace = TRUE, prob = los_probs)

# Utility: random date between two dates (inclusive)
rand_date <- function(n, start = study_start, end = study_end) {
  start_i <- as.integer(start)
  end_i <- as.integer(end)
  as.IDate(start_i + sample(0:(end_i - start_i), size = n, replace = TRUE), origin = "1970-01-01")
}

# CREATE BASIC POPULATION
ids <- seq_len(n_total)
sex <- sample(c(1L,2L), n_total, replace = TRUE) # 1 male, 2 female

# Ages baseline; override AMI ages later
ages <- sample(18:90, n_total, replace = TRUE)

# birth dates (approx) relative to ref date 2017-01-01
ref_date <- as.IDate("2017-01-01")
birth_d <- ref_date - round(ages * 365.25)

basic_dt <- data.table(id = ids, sex = sex, birth_d = birth_d)
basic_dt[, death_d := as.IDate(NA)]

# SELECT AMI PATIENTS BY SEX (exact counts)
male_ids <- basic_dt[sex == 1L, id]
female_ids <- basic_dt[sex == 2L, id]
set.seed(20251014) # ensure reproducible selection
ami_ids_m <- sample(male_ids, n_AMI_m)
ami_ids_f <- sample(female_ids, n_AMI_f)
ami_ids <- c(ami_ids_m, ami_ids_f)
ami_set <- rep(FALSE, n_total); ami_set[ami_ids] <- TRUE

# 6% have recurrent AMI
n_repeat <- ceiling(length(ami_ids) * 0.06)
repeat_ids <- sample(ami_ids, n_repeat)

# Adjust ages for AMI patients to match requested medians/IQRs
sd_male <- (75 - 55) / 1.35
sd_female <- (83 - 63) / 1.35
ages[ami_ids_m] <- pmax(18L, pmin(100L, round(rnorm(length(ami_ids_m), mean = 65, sd = sd_male))))
ages[ami_ids_f] <- pmax(18L, pmin(100L, round(rnorm(length(ami_ids_f), mean = 73, sd = sd_female))))
basic_dt[, age := ages]

# CONTACTS PER PERSON
contact_prob <- 0.92
has_contact <- runif(n_total) < contact_prob
# For those with contacts assign 2-3 contacts (user preference)
n_contacts <- integer(n_total)
n_contacts[has_contact] <- sample(2:3, sum(has_contact), replace = TRUE)
# Ensure AMI patients have at least one contact reserved for AMI main
n_contacts[ami_ids] <- pmax(n_contacts[ami_ids], 1L)
total_contacts_est <- sum(n_contacts)
message("Estimated total contacts: ", format(total_contacts_est, big.mark=","))

# VECTORISED NON-AMI CONTACT GENERATION (fast)
nonami_ids <- ids[!ami_set]
n_nonami <- length(nonami_ids)
# number of non-AMI contacts per non-AMI person = n_contacts (already set)
nonami_contacts <- n_contacts[nonami_ids]

# Expand ids
if (n_nonami > 0) {
  id_rep <- rep(nonami_ids, times = nonami_contacts)
  N_nonami_rows <- length(id_rep)
  nonami_dt <- data.table(id = id_rep)
  # contact_type: 1 inpatient (30%), 2 outpatient (70%)
  nonami_dt[, contact_type := sample(c(1L,2L), .N, replace = TRUE, prob = c(0.3, 0.7))]
  # dates
  nonami_dt[, date_adm := rand_date(.N)]
  # discharge logic: outpatients same day; inpatients get LOS sampled
  nonami_dt[contact_type == 1L, date_dis := date_adm + sample_los(.N)]
  nonami_dt[contact_type == 2L, date_dis := date_adm]
  # primary diag_a: mix of other_diag and comorb codes (more likely others)
  all_choices <- c(other_diag, ICD_HF, ICD_HYPERCHOL, ICD_DIABETES, ICD_HYPERT)
  nonami_dt[, diag_a := sample(all_choices, .N, replace = TRUE)]
  # number of diag_b per row: 1 or 2
  nb <- sample(1:2, N_nonami_rows, replace = TRUE, prob = c(0.85, 0.15))
  # assign diag_b as semicolon-concatenated codes
  nonami_dt[, diag_b := vapply(nb, function(k) paste(sample(all_choices, k, replace = FALSE), collapse = ";"), character(1))]
} else {
  nonami_dt <- data.table(id = integer(), date_adm = as.IDate(character()), date_dis = as.IDate(character()),
                          diag_a = character(), diag_b = character(), contact_type = integer())
}

# AMI MAIN ADMISSIONS (one per AMI patient) — choose index dates with decreasing year probs
n_ami_pat <- length(ami_ids)
ami_index_years <- sample(years, n_ami_pat, replace = TRUE, prob = year_probs)
ami_date_adm <- as.IDate(sapply(ami_index_years, function(y) as.character(sample(seq(as.Date(sprintf("%d-01-01", y)), as.Date(sprintf("%d-12-31", y)), by="day"), 1))))
ami_date_dis <- ami_date_adm + sample_los(n_ami_pat)
ami_main_dt <- data.table(id = ami_ids, date_adm = ami_date_adm, date_dis = ami_date_dis,
                          diag_a = ICD_AMI, diag_b = NA_character_, contact_type = 1L)

# REPEAT AMI (6%): choose gap 29-1825 days (>=28 days after) and within study_end (or adjust)
if (length(repeat_ids) > 0) {
  gap_days <- sample(29:1825, length(repeat_ids), replace = TRUE)
  first_ami_map <- setNames(ami_date_adm, ami_ids)
  repeat_adm <- as.IDate(first_ami_map[as.character(repeat_ids)]) + gap_days
  # clamp to study_end if beyond
  repeat_adm[repeat_adm > study_end] <- study_end - sample(0:30, sum(repeat_adm > study_end), replace = TRUE)
  repeat_dis <- repeat_adm + sample_los(length(repeat_ids))
  repeat_dt <- data.table(id = repeat_ids, date_adm = repeat_adm, date_dis = repeat_dis,
                          diag_a = ICD_AMI, diag_b = NA_character_, contact_type = 1L)
} else {
  repeat_dt <- data.table(id = integer(), date_adm = as.IDate(character()), date_dis = as.IDate(character()),
                          diag_a = character(), diag_b = character(), contact_type = integer())
}

# PRE-AMI COMORBIDITY CONTACTS FOR AMI PATIENTS
# For each AMI patient and each comorb type, decide presence per comorb_prev_ami.
# If present, create 2-3 pre-AMI contacts (evenly distributed within 10 years before AMI).
pre_rows <- list()
pre_idx <- 1L

for (k in seq_along(ami_ids)) {
  pid <- ami_ids[k]
  first_date <- ami_date_adm[k]
  s <- sex[pid]
  gender_label <- ifelse(s == 1L, "male", "female")
  window_start <- pmax(study_start, first_date - 3650) # 10 years before
  window_end <- first_date - 1
  # for each comorb
  # HF
  if (runif(1) < comorb_prev_ami$HF[gender_label]) {
    npre <- sample(2:3, 1)
    adm_dates <- as.IDate(window_start + sample(0:(as.integer(window_end - window_start)), npre, replace = TRUE))
    for (d in adm_dates) {
      pre_rows[[pre_idx]] <- list(id = pid, date_adm = d, date_dis = d + sample_los(1)[1], diag_a = ICD_HF, diag_b = ICD_HF, contact_type = 1L)
      pre_idx <- pre_idx + 1L
    }
  }
  # Hyperchol
  if (runif(1) < comorb_prev_ami$Hyperchol[gender_label]) {
    npre <- sample(2:3, 1)
    adm_dates <- as.IDate(window_start + sample(0:(as.integer(window_end - window_start)), npre, replace = TRUE))
    for (d in adm_dates) {
      pre_rows[[pre_idx]] <- list(id = pid, date_adm = d, date_dis = d + sample_los(1)[1], diag_a = ICD_HYPERCHOL, diag_b = ICD_HYPERCHOL, contact_type = 2L)
      pre_idx <- pre_idx + 1L
    }
  }
  # Diabetes
  if (runif(1) < comorb_prev_ami$Diabetes[gender_label]) {
    npre <- sample(2:3, 1)
    codes <- sample(ICD_DIABETES, npre, replace = TRUE)
    adm_dates <- as.IDate(window_start + sample(0:(as.integer(window_end - window_start)), npre, replace = TRUE))
    for (i2 in seq_len(npre)) {
      d <- adm_dates[i2]; code <- codes[i2]
      pre_rows[[pre_idx]] <- list(id = pid, date_adm = d, date_dis = d + sample_los(1)[1], diag_a = code, diag_b = code, contact_type = 2L)
      pre_idx <- pre_idx + 1L
    }
  }
  # Hypertension
  if (runif(1) < comorb_prev_ami$Hypertension[gender_label]) {
    npre <- sample(2:3, 1)
    codes <- sample(ICD_HYPERT, npre, replace = TRUE)
    adm_dates <- as.IDate(window_start + sample(0:(as.integer(window_end - window_start)), npre, replace = TRUE))
    for (i2 in seq_len(npre)) {
      d <- adm_dates[i2]; code <- codes[i2]
      pre_rows[[pre_idx]] <- list(id = pid, date_adm = d, date_dis = d + sample_los(1)[1], diag_a = code, diag_b = code, contact_type = 2L)
      pre_idx <- pre_idx + 1L
    }
  }
}

if (length(pre_rows) > 0) {
  pre_dt <- rbindlist(pre_rows)
} else {
  pre_dt <- data.table(id = integer(), date_adm = as.IDate(character()), date_dis = as.IDate(character()),
                       diag_a = character(), diag_b = character(), contact_type = integer())
}

# Non-AMI population: add comorbidity-specific contacts with scaled prevalences
# For each non-AMI id, for each comorb type decide presence with scaled prevalence; create 0-2 contacts
nonami_comorb_rows <- list(); ni <- 1L
nonami_ids_vec <- nonami_dt[, unique(id)]
# We'll sample for a subset in vectorized manner to speed up
# For memory reasons we do in chunks
chunk_ids <- split(nonami_ids_vec, ceiling(seq_along(nonami_ids_vec)/50000))
for (ch in seq_along(chunk_ids)) {
  ids_chunk <- chunk_ids[[ch]]
  for (pid in ids_chunk) {
    s <- sex[pid]
    gl <- ifelse(s==1L, "male", "female")
    # HF
    if (runif(1) < comorb_prev_nonami$HF[gl]) {
      npre <- sample(0:2, 1, prob = c(0.5,0.35,0.15))
      if (npre>0) {
        adm_dates <- rand_date(npre)
        for (d in adm_dates) {
          nonami_comorb_rows[[ni]] <- list(id = pid, date_adm = d, date_dis = d + sample_los(1)[1], diag_a = ICD_HF, diag_b = ICD_HF, contact_type = sample(c(1L,2L),1,prob=c(0.2,0.8)))
          ni <- ni + 1L
        }
      }
    }
    # Hyperchol
    if (runif(1) < comorb_prev_nonami$Hyperchol[gl]) {
      npre <- sample(0:2, 1, prob = c(0.5,0.35,0.15))
      if (npre>0) {
        adm_dates <- rand_date(npre)
        for (d in adm_dates) {
          nonami_comorb_rows[[ni]] <- list(id = pid, date_adm = d, date_dis = d + sample_los(1)[1], diag_a = ICD_HYPERCHOL, diag_b = ICD_HYPERCHOL, contact_type = 2L)
          ni <- ni + 1L
        }
      }
    }
    # Diabetes
    if (runif(1) < comorb_prev_nonami$Diabetes[gl]) {
      npre <- sample(0:2, 1, prob = c(0.5,0.35,0.15))
      if (npre>0) {
        codes <- sample(ICD_DIABETES, npre, replace = TRUE)
        adm_dates <- rand_date(npre)
        for (i2 in seq_len(npre)) {
          nonami_comorb_rows[[ni]] <- list(id = pid, date_adm = adm_dates[i2], date_dis = adm_dates[i2] + sample_los(1)[1], diag_a = codes[i2], diag_b = codes[i2], contact_type = 2L)
          ni <- ni + 1L
        }
      }
    }
    # Hypertension
    if (runif(1) < comorb_prev_nonami$Hypertension[gl]) {
      npre <- sample(0:2, 1, prob = c(0.5,0.35,0.15))
      if (npre>0) {
        codes <- sample(ICD_HYPERT, npre, replace = TRUE)
        adm_dates <- rand_date(npre)
        for (i2 in seq_len(npre)) {
          nonami_comorb_rows[[ni]] <- list(id = pid, date_adm = adm_dates[i2], date_dis = adm_dates[i2] + sample_los(1)[1], diag_a = codes[i2], diag_b = codes[i2], contact_type = 2L)
          ni <- ni + 1L
        }
      }
    }
  }
  message("Processed non-AMI comorb chunk ", ch, "/", length(chunk_ids))
}

if (length(nonami_comorb_rows) > 0) nonami_comorb_dt <- rbindlist(nonami_comorb_rows) else nonami_comorb_dt <- data.table()

# COMBINE ALL DIAGNOSIS ROWS
# non-AMI main rows (nonami_dt), comorb rows for non-AMI (nonami_comorb_dt), AMI main (ami_main_dt), pre-AMI comorb (pre_dt), repeat AMI (repeat_dt)
diag_dt <- rbindlist(list(nonami_dt, nonami_comorb_dt, pre_dt, ami_main_dt, repeat_dt), use.names = TRUE, fill = TRUE)
setorder(diag_dt, id, date_adm)

# Ensure diag_b is a semicolon-separated string (for rows where diag_b may be vector)
diag_dt[, diag_b := as.character(diag_b)]
# For rows where diag_b is NA or "NA", fill with diag_a or sensible default
diag_dt[is.na(diag_b) | toupper(diag_b) == "NA", diag_b := diag_a]

# Now: ensure that for each AMI patient, the index AMI row's diag_b includes all comorbidities that were present before AMI
# Build a pre-AMI map per AMI id: look for diagnoses in the 10-year window prior to AMI and gather relevant comorb codes
ami_first <- diag_dt[diag_a == ICD_AMI, .(first_ami = min(date_adm)), by = id]
setkey(ami_first, id)
# For speed, subset diag_dt to rows that could be prior comorbs (within 10 years prior to their AMI)
for (row in seq_len(nrow(ami_first))) {
  pid <- ami_first$id[row]
  first_date <- ami_first$first_ami[row]
  window_start <- pmax(study_start, first_date - 3650)
  prior_codes <- unique(diag_dt[id==pid & date_adm >= window_start & date_adm < first_date, diag_a])
  # Filter to keep only comorb families of interest and map family codes to canonical ones
  # canonical codes to add: ICD_HF, ICD_HYPERCHOL, any ICD_DIABETES, any ICD_HYPERT
  add_codes <- character(0)
  if (ICD_HF %in% prior_codes) add_codes <- c(add_codes, ICD_HF)
  if (ICD_HYPERCHOL %in% prior_codes) add_codes <- c(add_codes, ICD_HYPERCHOL)
  diabs <- intersect(prior_codes, ICD_DIABETES)
  if (length(diabs)>0) add_codes <- c(add_codes, diabs[1]) # add one diabetes code
  hyps <- intersect(prior_codes, ICD_HYPERT)
  if (length(hyps)>0) add_codes <- c(add_codes, hyps[1])
  # Also if a comorb wasn't present but should be according to prevalence, we may add it probabilistically (to maintain prevalence)
  s <- basic_dt[id==pid, sex]
  g <- ifelse(s==1L,"male","female")
  if (!ICD_HF %in% add_codes && runif(1) < comorb_prev_ami$HF[g]) add_codes <- c(add_codes, ICD_HF)
  if (!ICD_HYPERCHOL %in% add_codes && runif(1) < comorb_prev_ami$Hyperchol[g]) add_codes <- c(add_codes, ICD_HYPERCHOL)
  if (!any(add_codes %in% ICD_DIABETES) && runif(1) < comorb_prev_ami$Diabetes[g]) add_codes <- c(add_codes, sample(ICD_DIABETES,1))
  if (!any(add_codes %in% ICD_HYPERT) && runif(1) < comorb_prev_ami$Hypertension[g]) add_codes <- c(add_codes, sample(ICD_HYPERT,1))
  # Ensure at least one diag_b exists
  if (length(add_codes)==0) add_codes <- sample(c(ICD_HYPERCHOL, ICD_HF, sample(ICD_DIABETES,1)),1)
  # locate the index AMI row(s) for this patient and set diag_b (concatenate unique codes)
  diag_dt[id==pid & diag_a==ICD_AMI & date_adm==first_date, diag_b := paste(unique(unlist(strsplit(diag_b,";")), add_codes), collapse=";"), by=.(id)]
}

# Final fixes for date logic:
# coerce types
diag_dt[, date_adm := as.IDate(date_adm)]
diag_dt[, date_dis := as.IDate(date_dis)]
# Outpatient: date_dis == date_adm (ensure)
diag_dt[contact_type==2L, date_dis := date_adm]
# Inpatient: ensure at least 1 day stay
diag_dt[contact_type==1L & date_dis <= date_adm, date_dis := date_adm + 1L]

# Ensure all rows have diag_a and diag_b
diag_dt[is.na(diag_a), diag_a := sample(c(other_diag, ICD_HF, ICD_HYPERCHOL, ICD_DIABETES, ICD_HYPERT), .N, replace = TRUE)]
diag_dt[is.na(diag_b) | diag_b=="" , diag_b := diag_a]

# ---------------------------
# MORTALITY: AMI patients
# ---------------------------
# Target cumulative probabilities for males at days: 0,10,30,182,365 -> p = 0, 0.08, 0.13, 0.17, 0.21
t_points <- c(0L,10L,30L,182L,365L)
p_male <- c(0.0, 0.08, 0.13, 0.17, 0.21)
p_female <- p_male * 0.8
max_day <- 365L * 5L
days_vec <- 0:max_day
male_cdf <- approx(x = t_points, y = p_male, xout = days_vec, rule = 2)$y
female_cdf <- approx(x = t_points, y = p_female, xout = days_vec, rule = 2)$y

sample_surv_days <- function(gender) {
  u <- runif(1)
  cdf <- if (gender==1L) male_cdf else female_cdf
  idx <- which(cdf >= u)[1]
  if (is.na(idx)) return(Inf)
  return(as.integer(idx - 1L))
}

# Map first AMI dates
first_ami <- diag_dt[diag_a==ICD_AMI, .(first_ami = min(date_adm)), by=.(id)]

# For each AMI patient, sample survival days and set death_d in basic_dt if within study_end
basic_dt[, death_d := as.IDate(NA)]
for (r in seq_len(nrow(first_ami))) {
  pid <- first_ami$id[r]
  firstd <- first_ami$first_ami[r]
  g <- basic_dt[id==pid, sex]
  surv <- sample_surv_days(g)
  if (is.finite(surv) && surv > 0L) {
    dd <- firstd + surv
    if (dd <= study_end) basic_dt[id==pid, death_d := dd]
  }
}

# NON-AMI mortality: lower, increasing with age
ann_prob <- function(age) { 0.0008 * exp((age - 40)/30) }
# for non-AMI persons (no AMI first event), sample death year with low prob
non_ami_ids2 <- setdiff(basic_dt$id, first_ami$id)
# vectorize: compute annual probabilities across years and then draw event
for (pid in non_ami_ids2) {
  age_i <- basic_dt[id==pid, as.integer((as.IDate("2017-01-01") - birth_d)/365.25)]
  annual <- ann_prob(age_i)
  year_weights2 <- seq(1.0, 1.12, length.out = length(years))
  probs <- annual * year_weights2
  p_total <- 1 - prod(1 - probs)
  if (runif(1) < p_total) {
    probs_norm <- probs / sum(probs)
    yidx <- sample(seq_along(years), size = 1, prob = probs_norm)
    y <- years[yidx]
    dd <- rand_date(1, start = as.IDate(paste0(y,"-01-01")), end = as.IDate(paste0(y,"-12-31")))
    basic_dt[id==pid, death_d := dd]
  }
}

# Cleanup: drop auxiliary columns if any
basic_dt[, age := NULL]

# Final checks
message("Unique persons in diag: ", diag_dt[, uniqueN(id)])
message("Total diag rows: ", nrow(diag_dt))
message("AMI unique patients: ", diag_dt[diag_a==ICD_AMI, uniqueN(id)])
message("Writting files to data/ ...")


# Correct MI age distribution ---------------------------------------------

library(data.table)

# Identify AMI patients
ICD_AMI <- "DI21"
ami_ids <- diag_dt[diag_a == ICD_AMI, unique(id)]

# Parameters for target distributions
target_params <- list(
  male   = list(mean = 65, sd = (75 - 55) / 1.35),
  female = list(mean = 73, sd = (83 - 63) / 1.35)
)

# Re-sample ages for AMI patients only
for (sex_code in c(1L, 2L)) {
  g <- ifelse(sex_code == 1L, "male", "female")
  ids_sub <- basic_dt[id %in% ami_ids & sex == sex_code, id]
  if (length(ids_sub) > 0) {
    new_ages <- round(rnorm(
      n = length(ids_sub),
      mean = target_params[[g]]$mean,
      sd   = target_params[[g]]$sd
    ))
    # Clip ages to 18–100
    new_ages <- pmax(18L, pmin(100L, new_ages))
    basic_dt[id %in% ids_sub, adj_age := new_ages]
  }
}

# For non-AMI keep original
basic_dt[!id %in% ami_ids, adj_age := as.integer((as.IDate("2017-01-01") - birth_d) / 365.25)]

# Optional: check distributions
basic_dt[id %in% ami_ids, .(
  median = median(adj_age),
  IQR_L  = quantile(adj_age, .25),
  IQR_H  = quantile(adj_age, .75)
), by = sex]

# Update birth dates so ages match 2017 reference
ref_date <- as.IDate("2017-01-01")
basic_dt[, birth_d := ref_date - round(adj_age * 365.25)]

# Clean up and save
basic_dt[, adj_age := NULL]

# Tidy --------------------------------------------------------------------

basic_data <- basic_dt |>
  dplyr::mutate(sex = factor(sex, levels = c(1,2))) |>
  dplyr::mutate(id = factor(id)) |>
  dplyr::mutate(dplyr::across(c(birth_d, death_d), ~as.Date(.x)))

diag_data <- diag_dt |>
  dplyr::mutate(
    dplyr::across(c(id, contact_type), ~ as.factor(.x)),
    dplyr::across(c(date_adm, date_dis), ~ as.Date(.x)),
    dplyr::across(c(diag_a, diag_b), ~ as.character(.x))
  ) |>
  tidyr::separate_rows(diag_b, sep = ";")


basic_data <- basic_data |>
  rename(
    "birth_date" = birth_d,
    "death_date" = death_d
  )

colnames(diag_data)

diag_data <- diag_data |>
  rename(
    "adm_date" = date_adm,
    "dis_date" = date_dis
  )

basic_data <- basic_data |>
  dplyr::mutate(
    across(c(birth_date, death_date), ~ as.Date(.x))
  )

diag_data <- diag_data |>
  dplyr::mutate(
    dplyr::across(c(adm_date, dis_date), ~ as.Date(.x))
  )

saveRDS(basic_data, here::here("data_raw/basic_data.rds"))
saveRDS(diag_data, here::here("data_raw/diag_data.rds"))

message("Done. Files: data/basic_data.csv, data/diag_data.csv (and .rds)")

