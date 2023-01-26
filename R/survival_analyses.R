library(broom)
library(dplyr)
library(ggplot2)
library(here)
library(lubridate)
library(parallel)
library(purrr)
library(readxl)
library(rlang)
library(survival)
library(survminer)
library(wesanderson)
library(writexl)

options(mc.cores = detectCores() - 1)

# Help functions
# https://stackoverflow.com/questions/68084740
catch_warn <- function(expr) {
  res <- suppressWarnings(tryCatch(
    expr = {
      withCallingHandlers(
        expr = expr,
        warning = function(w) {
          parent <- parent.env(environment())
          parent$warning_arg <- w
        }
      )
    }
  ))
  w <- if ("warning_arg" %in% ls()) {
    trimws(as.character(warning_arg))
  } else {
    NULL
  }
  return(list(result = res, warning = w))
}
# https://github.com/cran/RcmdrPlugin.KMggplot2/blob/master/R/geom-stepribbon.r
source(paste0("https://raw.githubusercontent.com/cran/RcmdrPlugin.KMggplot2/",
              "master/R/geom-stepribbon.r"))

# Set project directory
i_am("R/survival_analyses.R")

# Output directory
outdir <- here("results", paste0("survival_analyses_",
                                 format(Sys.Date(), "%Y%m%d")))
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

# Load data
raw_data <- read_xlsx(here(
  "data/Criteres-non-rep-cery-GENERALISE_Pour_Jerome_score_activite_V3.xlsx"))

# Recode dates
raw_data <- raw_data %>%
  mutate(DATE_INIT_TTT = as.Date(DATE_INIT_TTT), DATE_OBS = as.Date(DATE_OBS))

# Check that date of treatment initialization, sex, age, ethnicity and
# pheno_2c19 are well defined
raw_data %>%
  select(IPP, DATE_OBS, DATE_INIT_TTT, Sexe, Age_c, groupe_ethnique,
         PHENO_2C19, PHENO_EXTREME) %>%
  group_by(IPP) %>%
  filter(length(unique(DATE_INIT_TTT)) > 1 | any(is.na(DATE_INIT_TTT)) |
         length(unique(Sexe)) > 1 | any(is.na(Sexe)) |
         length(unique(Age_c)) > 1 | any(is.na(Age_c)) |
         length(unique(groupe_ethnique)) > 1 | any(is.na(groupe_ethnique)) |
         length(unique(PHENO_2C19)) > 1 | any(is.na(PHENO_2C19))) %>%
  as.data.frame()

# Redefine ethnicity
ethnicity <- raw_data %>%
  select(IPP, groupe_ethnique) %>%
  group_by(IPP, groupe_ethnique) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(IPP) %>%
  filter(n == max(n)) %>%
  select(!n)
raw_data <- raw_data %>%
  select(!groupe_ethnique) %>%
  left_join(ethnicity, by = "IPP")
rm(ethnicity)

# Time independent covariates
fixed_covariates <- raw_data %>%
  select(IPP, sex = Sexe, age = Age_c, ethnicity = groupe_ethnique) %>%
  unique()
if (any(duplicated(fixed_covariates$IPP))) stop("dup ipp")

# Check that p_phenotype_2c19, pheno_extreme and escitalopram dosage are well
# defined
raw_data %>%
  select(IPP, DATE_OBS, P_phenotype_2C19, PHENO_EXTREME,
         DOSE_CITESCIT_AJUSTEE) %>%
  group_by(IPP, DATE_OBS) %>%
  filter(length(unique(DOSE_CITESCIT_AJUSTEE)) > 1 |
         length(unique(P_phenotype_2C19)) > 1 | any(is.na(P_phenotype_2C19)) |
         length(unique(PHENO_EXTREME)) > 1 | any(is.na(PHENO_EXTREME))) %>%
  as.data.frame()

# Fix escitalopram dosage
dosage <- raw_data %>%
  select(IPP, DATE_OBS, DOSE_CITESCIT_AJUSTEE) %>%
  group_by(IPP, DATE_OBS) %>%
  summarise(DOSE_CITESCIT_AJUSTEE = DOSE_CITESCIT_AJUSTEE[1], .groups = "drop")
raw_data <- raw_data %>%
  select(!DOSE_CITESCIT_AJUSTEE) %>%
  left_join(dosage, by = c("IPP", "DATE_OBS"))
rm(dosage)

# Time dependent covariates
varying_covariates <- raw_data %>%
  mutate(time0 = as.numeric(DATE_OBS - DATE_INIT_TTT)) %>%
  select(IPP, time0, dosage = DOSE_CITESCIT_AJUSTEE) %>%
  unique()
if (any(duplicated(varying_covariates[c("IPP", "time0")])))
  stop("dup ipp/date")

# Survival analyses
groups <- c("PHENO_2C19", "P_phenotype_2C19", "PHENO_EXTREME")
events <- grep("^CRIT_", names(raw_data), value = TRUE)
min_cens_times <- c(0, 60, 90, 120, 180, 365)
R <- unlist(recursive = FALSE, lapply(groups %>% setNames(., .), function(g) {
  unlist(recursive = FALSE, lapply(events %>% setNames(., .), function(e) {
    lapply(min_cens_times %>% setNames(., .), function(m) {
      list(group = g, event = e, min_cens_time = m)
    })
  }))
}))
R <- mclapply(R, function(r) {
  #print(do.call(paste, r))
  limits <- c(60, 90, 120, 180, 365, Inf) %>% setNames(., .)
  # print(paste(r$group, r$event))
  # Survival data (for KM analysis)
  # * If an event and a non-event are observed simultaneously for the same
  #   IPP/Date, we keep only the event.
  # * The group is defined at the first observation
  surv_data_0 <- raw_data %>%
    group_by(IPP, date_obs = DATE_OBS, date_init_ttt = DATE_INIT_TTT,
             group = !!sym(r$group)) %>%
    summarise(event = sum(!!sym(r$event)) > 0, .groups = "drop") %>%
    arrange(IPP, date_obs) %>%
    group_by(IPP)
  surv_data_1 <- surv_data_0 %>%
    mutate(event = cumsum(event)) %>%
    summarise(
      date_first_obs = date_obs[1],
      group = group[1],
      date_init_ttt = date_init_ttt[1],
      date_before_event = suppressWarnings(max(date_obs[event == 0])),
      date_after_event = suppressWarnings(min(date_obs[event == 1])),
      event = any(event > 0)
    ) %>%
    mutate(
      group = factor(group,
        if (r$group == "PHENO_EXTREME") 0:1 else c("EMIM", "UM", "PM")),
      time_first_obs = as.numeric(date_first_obs - date_init_ttt),
      time_before_event = as.numeric(date_before_event - date_init_ttt),
      time_before_event = if_else(event, time_before_event,
                                  pmax(time_before_event, r$min_cens_time)),
      time_after_event = as.numeric(date_after_event - date_init_ttt),
      observed = !(event & time_first_obs == time_after_event |
                   !event & time_first_obs == time_before_event),
      time_surv = case_when(
        !observed ~ NA_real_,
        event ~ time_after_event,
        TRUE ~ time_before_event
      ),
      surv = Surv(time = time_first_obs, time2 = time_surv, event = event)
    )
  # Survival data by interval (for Cox regression)
  surv_data_2 <- surv_data_0 %>%
    filter(row_number() <= if (any(event)) min(which(event)) else Inf) %>%
    mutate(
      group = if_else(row_number() == 1, group, lag(group)),
      time1 = as.numeric(date_obs - date_init_ttt),
      time0 = if_else(row_number() == 1, time1, lag(time1)),
      time1 = if_else(event | row_number() < n(), time1,
                      pmax(time1, r$min_cens_time))
    ) %>%
    filter(time1 > time0) %>%
    mutate(
      group = factor(group,
        if (r$group == "PHENO_EXTREME") 0:1 else c("EMIM", "UM", "PM")),
      surv = Surv(time0, time1, event)
    ) %>%
    select(IPP, group, time0, time1, event, surv) %>%
    left_join(fixed_covariates, by = "IPP") %>%
    left_join(varying_covariates, by = c("IPP", "time0"))
  # Overview of individual survivals
  sorted_ipps <- surv_data_2 %>%
    filter(row_number() == n()) %>%
    arrange(time1) %>%
    pull(IPP)
  surv_paths <- surv_data_2 %>%
    group_by(IPP) %>%
    mutate(IPP = factor(IPP, sorted_ipps)) %>%
    ggplot(aes(colour = group)) +
    geom_segment(aes(y = IPP, yend = IPP, x = time0, xend = time1)) +
    geom_point(aes(y = IPP, x = time1, shape = event)) +
    labs(x = "days") +
    theme_bw() +
    theme(axis.text.y = element_text(size = 6), legend.position = "bottom")
  for (l in limits) {
    surv_paths = surv_paths +
      geom_vline(xintercept = l, color = "grey", linetype = "dashed")
  }
  # KM analysis
  km_fit <- survfit(surv ~ group, surv_data_1)
  names(km_fit$strata) <- sub("^group=", "", names(km_fit$strata))
  # Log rank tests
  lr_tests_1 <- map_dfr(limits, function(l) {
    d <- surv_data_1 %>%
      select(IPP, group, time_first_obs, time_surv, event) %>%
      filter(time_first_obs < l) %>%
      mutate(
        group = droplevels(group),
        b = time_surv > l,
        time_surv = if_else(b, l, time_surv),
        event = if_else(b, FALSE, event),
        surv = Surv(time = time_first_obs, time2 = time_surv, event = event)
      )
    pv <- catch_warn(summary(coxph(surv ~ group, d))$sctest[["pvalue"]])
    tibble(
      limit = l,
      logrank_pv = pv$result,
      warning = pv$warning
    )
  })
  # KM survival table
  km_surv_tab <- tidy(km_fit)
  # KM curves
  km_figs <- map(limits, function(l) {
    pv <- round(pull(filter(lr_tests_1, limit == l), logrank_pv), 3)
    l <- min(l, max(surv_data_1$time_surv, na.rm = TRUE))
    p <- ggsurvplot(
      km_fit,
      data = surv_data_1,
      #pval = pv,
      conf.int = TRUE,
      risk.table = TRUE,
      risk.table.col = "strata",
      xlim = c(0, l),
      break.time.by = l / 10,
      ggtheme = theme_bw(),
      palette = wes_palette("Darjeeling1", length(km_fit$strata))
    )
    p$plot <- p$plot + annotate("text", x = 0, y = 0, label = paste("p =", pv))
    return(p)
  })
  # Cox regressions (according to the observed treatment duration)
  cox_fits <- map(limits, function(l) {
    d <- surv_data_2 %>%
      filter(time0 < l) %>%
      mutate(
        group = droplevels(group),
        b = time1 > l,
        time1 = if_else(b, l, time1),
        event = if_else(b, FALSE, event),
        surv = Surv(time0, time1, event)
      )
    z <- catch_warn(coxph(surv ~ group, d))
    cox_fit <- z$result
    attr(cox_fit, "warning") <- z$warning
    # Cox survival table
    grps <- levels(d$group)
    cox_surv_tab <- do.call(bind_rows, lapply(grps, function(g) {
      survfit(cox_fit, newdata = data.frame(group = g)) %>%
        tidy() %>%
        mutate(strata = g)
    }))
    # Cox curve
    pv <- summary(cox_fit)$sctest[["pvalue"]]
    cox_fig <- cox_surv_tab %>%
      ggplot(aes(x = time, y = estimate, colour = strata, fill = strata)) +
      geom_step() +
      #geom_stepribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .5,
      #                color = NA) +
      annotate("text", x = 0, y = 0, label = paste("p =", round(pv, 3))) +
      labs(x = "days", y = "survival probability") +
      theme_bw()
    # Graphical Test of Proportional Hazards
    cox_fig_ph_test <- tryCatch(ggcoxzph(cox.zph(cox_fit))[[1]],
                                error = function(e) NULL)
    list(cox_fit = cox_fit, cox_surv_tab = cox_surv_tab, cox_fig = cox_fig,
         cox_fig_ph_test = cox_fig_ph_test)
  })
  # Log rank tests using survival data by interval
  lr_tests_2 <- map_dfr(names(cox_fits), function(l) {
    fit <- cox_fits[[l]]$cox_fit
    pv <- summary(fit)$sctest[["pvalue"]]
    tibble(
      limit = as.numeric(l),
      logrank_pv = pv,
      warning = attr(fit, "warning")
    )
  })
  # Lists of results
  surv_data <- list(
    survival_data              = surv_data_1,
    survival_data_by_intervall = surv_data_2
  )
  fits <- map(cox_fits, ~ .x$cox_fit) %>%
    setNames(., paste0("cox_", names(.))) %>%
    append(list(km = km_fit), after = 0)
  lr_tests <- list(
    fixed_group         = lr_tests_1,
    time_varying_groups = lr_tests_2
  )
  surv_tabs <- map(cox_fits, ~ .x$cox_surv_tab) %>%
    setNames(., paste0("cox_", names(.))) %>%
    append(list(km = km_surv_tab), after = 0)
  surv_figs <- append(
    km_figs %>% setNames(., paste0("km_", names(.))),
    map(cox_fits, ~ .x$cox_fig) %>% setNames(., paste0("cox_", names(.)))
  )
  ph_test_figs <- map(cox_fits, ~ .x$cox_fig_ph_test)
  risk_ratios <- map(cox_fits, ~ {
    tidy(.x$cox_fit, exponentiate = TRUE, conf.int=TRUE) %>%
      select(term, estimate, std.error, conf.low, conf.high, p.value)
  }) %>%
    setNames(., paste0("cox_", names(.)))
  # Differences
  # https://dominicmagirr.github.io/post/2022-01-18-be-careful-with-standard-errors-in-survival-survfit/
  diffs <- lapply(names(surv_tabs), function(z) {
    K <- if (grepl("cox", z)) min(as.numeric(sub("cox_", "", z)), 365) else 365
    do.call(rbind, lapply(limits[limits <= K], function(l) {
      d <- surv_tabs[[z]] %>%
        mutate(var.log.survival = std.error^2) %>%
        group_by(strata) %>%
        filter(time == suppressWarnings(max(time[time <= l])))
      if (nrow(d) < 2) return(NULL)
      s <- d$strata
      e <- d$estimate
      w <- d$var.log.survival
      do.call(rbind, lapply(1:(length(s) - 1), function(i) {
        do.call(rbind, lapply((i + 1):length(s), function(j) {
          z_log <- abs(log(e[i]) - log(e[j])) / sqrt(w[i] + w[j])
          pv <- if (e[i] < 1 & e[j] < 1) 2 * (1 - pnorm(z_log)) else NA
          data.frame(
            day = l,
            strata_1 = s[i],
            strata_2 = s[j],
            surv_1 = e[i],
            surv_2 = e[j],
            p_value = 2 * (1 - pnorm(z_log))
          )
        }))
      }))
    }))
  })
  rownames(diffs) <- NULL
  # All results in a list
  list(surv_data = surv_data, surv_paths = surv_paths, fits = fits,
       lr_tests = lr_tests, surv_tabs = surv_tabs, surv_figs = surv_figs,
       ph_test_figs = ph_test_figs, risk_ratios = risk_ratios, diffs = diffs)
})

# Export results
for (s in names(R)) {
  o <- file.path(outdir, s)
  if (!dir.exists(o)) dir.create(o)
  r <- R[[s]]
  write_xlsx(r$surv_data, file.path(o, "survival_data.xlsx"))
  svg(file.path(o, "surv_paths.svg"), width = 14, height = 35)
  print(r$surv_paths)
  dev.off()
  write_xlsx(r$lr_tests, file.path(o, "logrank_tests.xlsx"))
  write_xlsx(r$surv_tabs, file.path(o, "survival_tables.xlsx"))
  for (l in names(r$surv_figs)) {
    svg(file.path(o, paste0("fig_", l, ".svg")), width = 14, height = 7)
    print(r$surv_figs[[l]])
    dev.off()
  }
  for (l in names(r$ph_test_figs)) {
    if (!is.null(r$ph_test_figs[[l]])) {
      svg(file.path(o, paste0("fig_ph_test_", l, ".svg")),
          width = 14, height = 7)
      print(r$ph_test_figs[[l]])
      dev.off()
    }
  }
  write_xlsx(r$risk_ratios, file.path(o, "risk_ratios.xlsx"))
  write_xlsx(r$diffs, file.path(o, "diffs.xlsx"))
}
rm(s, o, r, l)

# Session info
sink(file.path(outdir, "sessionInfo.txt"))
print(sessionInfo(), locale = FALSE)
sink()
