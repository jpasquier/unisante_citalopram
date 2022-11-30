library(broom)
library(dplyr)
library(ggplot2)
library(haven)
library(parallel)
library(rlang)
library(survival)
library(survminer)
library(wesanderson)
library(writexl)

options(mc.cores = detectCores() - 1)

# Help function
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
    as.character(warning_arg)
  } else {
    NULL
  }
  return(list(result = res, warning = w))
}

# Working directory
setwd('~/Projects/Consultations/Coumau Aude (CYP2C19)')

# Output directory
outdir <- paste0('results/survival_analyses_',
                 format(Sys.Date(), '%Y%m%d'))
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

# Load data
raw_data <- read_dta(
  'data/Criteres-non-rep-cery-GENERALISE_Pour_Jerome_score_activite_V3.dta')

# Survival analyses
groups <- c('PHENO_2C19', 'P_phenotype_2C19', 'PHENO_EXTREME')
events <- c('CRIT_FULL', 'CRIT_QUET', 'CRIT_MIRTA', 'CRIT_ARIPI', 'CRIT_LIT',
            'CRIT_DEUX_AD', 'CRIT_ARRET', 'CRIT_SWITCH', 'CRIT_DOSE_MAX',
            'CRIT_DOSE_MIN', 'CRIT_SWITCH_ARRET_QUET',
            'CRIT_SWITCH_ARRET_DEUX_AD', 'CRIT_SWITCH_ARRET')
R <- unlist(recursive =FALSE, lapply(groups %>% setNames(., .), function(g) {
  lapply(events %>% setNames(., .), function(e) list(group = g, event = e))
}))
R <- mclapply(R, function(r) {
  # print(paste(r$group, r$event))
  # Survival data
  # * If an event and a non-event are observed simultaneously for the same
  #   IPP/Date, we keep only the event.
  # * The group is defined at the first observation
  surv_data <- raw_data %>%
    group_by(IPP, date_obs = DATE_OBS, date_init_ttt = DATE_INIT_TTT,
             group = !!sym(r$group)) %>%
    summarise(event = sum(!!sym(r$event)) > 0, .groups = 'drop') %>%
    arrange(IPP, date_obs) %>%
    group_by(IPP) %>%
    mutate(
      event = cumsum(event)
    ) %>%
    summarise(
      date_first_obs = date_obs[1],
      group = group[1],
      date_init_ttt = date_init_ttt[1],
      date_before_event = suppressWarnings(max(date_obs[event == 0])),
      date_after_event = suppressWarnings(min(date_obs[event == 1])),
      event = any(event > 0)
    ) %>%
    mutate(
      group = factor(group),
      observed = !(event & date_first_obs == date_after_event |
                   !event & date_first_obs == date_before_event),
      time_first_obs = as.numeric(date_first_obs - date_init_ttt),
      obs_365 = observed & time_first_obs <= 365,
      time_before_event = if_else(
        observed, as.numeric(date_before_event - date_init_ttt), NA_real_),
      time_after_event = if_else(
        observed, as.numeric(date_after_event - date_init_ttt), NA_real_),
      time_surv = if_else(event, time_after_event, time_before_event),
      surv = Surv(time = time_first_obs, time2 = time_surv, event = event)
    )
  # Overview of individual survivals
  surv_paths <- surv_data %>%
    filter(observed) %>%
    mutate(IPP = factor(IPP, IPP[order(time_surv)])) %>%
    ggplot(aes(colour = group)) +
    geom_segment(aes(y = IPP, yend = IPP, x = 0, xend = time_first_obs),
                 linetype = 'dotted') +
    geom_segment(aes(y = IPP, yend = IPP, x = time_first_obs,
                     xend = time_surv)) +
    geom_point(aes(y = IPP, x = time_before_event), shape = 1) +
    geom_point(aes(y = IPP, x = time_after_event), shape = 4) +
    geom_vline(xintercept = 365, color = 'grey', linetype = 'dashed') +
    theme(axis.text.y = element_text(size = 6),
          legend.position = 'bottom')
  # KM analyses
  # Log rank test with left truncated data:
  # see https://stat.ethz.ch/pipermail/r-help/2009-August/399999.html
  d <- list(lt     = filter(surv_data, observed),
            lt_365 = filter(surv_data, obs_365))
  km_fits <- list(
    lt     = survfit(surv ~ group, d$lt),
    lt_365 = survfit(surv ~ group, d$lt_365)
  ) %>%
    lapply(function(fit) {
      names(fit$strata) <- sub('^group=', '', names(fit$strata))
      fit
    })
  for (z in names(km_fits)) attr(km_fits[[z]], 'data') <- d[[z]]
  # Cox analyses
  cox_fits <- list(
    lt = catch_warn(coxph(surv ~ group, filter(surv_data, observed))),
    lt_365 = catch_warn(coxph(surv ~ group, filter(surv_data, obs_365)))
  )
  # KM curves
  km_figs <- lapply(names(km_fits) %>% setNames(., .), function(s) {
    fit <- km_fits[[s]]
    data <- attr(fit, 'data')
    pv <- round(summary(cox_fits[[s]]$result)$sctest['pvalue'], 3)
    if (!is.null(cox_fits[[s]]$warning)) paste0(pv, '*')
    if (s == 'lt') {
      p <- ggsurvplot(
        fit,
        data = data,
        pval = pv,
        conf.int = TRUE,
        risk.table = TRUE,
        risk.table.col = "strata",
        ggtheme = theme_bw(),
        palette = wes_palette('Darjeeling1', length(fit$strata))
      )
      p$plot <- p$plot +
        geom_vline(xintercept = 365, color = 'grey', linetype = 'dashed')
    } else {
      p <- ggsurvplot(
        fit,
        data = data,
        pval = pv,
        conf.int = TRUE,
        risk.table = TRUE,
        risk.table.col = "strata",
        ggtheme = theme_bw(),
        xlim = c(0, 365),
        break.time.by = 30,
        palette = wes_palette('Darjeeling1', length(fit$strata))
      )
    }
    if (!is.null(cox_fits[[s]]$warning)) {
      p <- p +
        labs(caption = paste('* Warning (calculation of the p-value):',
                             cox_fits[[s]]$warning))
    }
    p
  })
  # Differences at 1, 2, 3, 6 and 12 months
  # https://dominicmagirr.github.io/post/2022-01-18-be-careful-with-standard-errors-in-survival-survfit/
  diffs <- do.call(rbind, lapply(round(365.2425 / 12 * c(1:3, 6, 12)),
                                 function(d) {
    fit <- km_fits$lt_365
    tidy(fit) %>%
      mutate(
        var.survival = estimate^2 * std.error^2,
        var.log.survival = std.error^2,
      ) %>%
      group_by(strata) %>%
      filter(time == suppressWarnings(max(time[time <= d]))) %>%
      {
        s <- .$strata
        e <- .$estimate
        v <- .$var.survival
        w <- .$var.log.survival
        do.call(rbind, lapply(1:(length(s) - 1), function(i) {
          do.call(rbind, lapply((i + 1):length(s), function(j) {
            z_plain <- abs(e[i] - e[j]) / sqrt(v[i] + v[j])
            z_log <- abs(log(e[i]) - log(e[j])) / sqrt(w[i] + w[j])
            data.frame(
              day = d,
              fit = 'lt_365',
              strata_1 = s[i],
              strata_2 = s[j],
              surv_1 = e[i],
              surv_2 = e[j],
              var_surv_1 = v[i],
              var_surv_2 = v[j],
              log_surv_1 = log(e[i]),
              log_surv_2 = log(e[j]),
              var_log_surv_1 = w[i],
              var_log_surv_2 = w[j],
              p_value = 2 * (1 - pnorm(z_log))
            )
          }))
        }))
      }
  }))
  rownames(diffs) <- NULL
  list(surv_data = surv_data, surv_paths = surv_paths, km_fits = km_fits,
       cox_fits = cox_fits, km_figs = km_figs, diffs = diffs)
})

# Export results
for (s in names(R)) {
  o <- file.path(outdir, s)
  if (!dir.exists(o)) dir.create(o)
  r <- R[[s]]
  write_xlsx(path = file.path(o, 'survival_tables.xlsx'), list(
    all                    = tidy(r$km_fits$lt),
    observed_before_day365 = tidy(r$km_fits$lt_365),
    survival_data          = r$surv_data
  ))
  svg(file.path(o, 'surv_paths.svg'), width = 14, height = 35)
  print(r$surv_paths)
  dev.off()
  for (z in names(r$km_figs)) {
    svg(file.path(o, paste0('km_', z, '.svg')),
        width = 14, height = 7)
    print(r$km_figs[[z]])
    dev.off()
  }
  write_xlsx(r$diffs, file.path(o, 'diffs.xlsx'))
}
rm(s, o, r, z)

# Session info
sink(file.path(outdir, 'sessionInfo.txt'))
print(sessionInfo(), locale = FALSE)
sink()
