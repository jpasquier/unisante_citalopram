library(broom)
library(dplyr)
library(ggplot2)
library(readxl)
library(survival)
library(survminer)
library(wesanderson)
library(writexl)

# Working directory
setwd('~/Projects/Consultations/Coumau Aude (CYP2C19)')

# Criterion
K <- 2

# Output directory
outdir <- paste0('results/survival_analyses_', format(Sys.Date(), '%Y%m%d'),
                 '/criterion', K)
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

# Load data
rawDataFile <- if (K == 1) {
  'data/Criteres-non-rep-cery-GENERALISE_Pour_Jerome.xlsx'
} else {
  'data/Criteres-non-rep-cery-GENERALISE_Pour_Jerome2.xlsx'
}
rawData <- read_xlsx(rawDataFile) %>%
  mutate(DATE_OBS = as.Date(DATE_OBS), DATE_INIT_TTT = as.Date(DATE_INIT_TTT))

# Check that there is only one DATE_INIT_TTT and one PHENO_2C19 per IPP
all(aggregate(DATE_INIT_TTT ~ IPP, rawData,
              function(x) length(unique(x)) == 1)[[2]])
all(aggregate(PHENO_2C19 ~ IPP, rawData,
              function(x) length(unique(x)) == 1)[[2]])

# If an event and a non-event are observed simultaneously for the same
# IPP/Date, we keep only the event.
rawData <- rawData %>%
  group_by(IPP, DATE_OBS, DATE_INIT_TTT, PHENO_2C19) %>%
  summarise(Critere = sum(Critere) > 0, .groups = 'drop')

# Sort data by IPP and DATE_OBS
rawData <- arrange(rawData, IPP, DATE_OBS)

# Survival data
survData <- rawData %>%
  select(IPP, DateInitTTT = DATE_INIT_TTT, Pheno2C19 = PHENO_2C19) %>%
  unique() %>%
  left_join(
    by = 'IPP',
    rawData %>%
      group_by(IPP) %>%
      mutate(
        DATE_OBS = as.Date(DATE_OBS),
        Event = cumsum(Critere)
      ) %>%
      summarise(
        DateFirstObs = min(DATE_OBS),
        DateBeforeEvent = suppressWarnings(max(DATE_OBS[Event == 0])),
        DateAfterEvent = suppressWarnings(min(DATE_OBS[Event == 1])),
        Event = any(Event > 0)
      )
  ) %>%
  mutate(
    Pheno2C19 = factor(Pheno2C19),
    # Pheno2C19 = factor(case_when(
    #               Pheno2C19 %in% c('EMIM', 'UM') ~ 'EMIM/UM',
    #               Pheno2C19 == 'PM' ~ 'PM'
    #             )),
    Observed = !(Event & DateFirstObs == DateAfterEvent |
                 !Event & DateFirstObs == DateBeforeEvent),
    TimeFirstObs = as.numeric(DateFirstObs - DateInitTTT),
    TimeBeforeEvent = if_else(
      Observed, as.numeric(DateBeforeEvent - DateInitTTT), NA_real_),
    TimeAfterEvent = if_else(
      Observed, as.numeric(DateAfterEvent - DateInitTTT), NA_real_),
    TimeSurv = if_else(Event, TimeAfterEvent, TimeBeforeEvent),
    # TimeSurv = if_else(Event, (TimeAfterEvent - TimeBeforeEvent) / 2,
    #                    TimeBeforeEvent),
    SurvNaive = Surv(TimeSurv, Event),
    SurvInterval = Surv(time = TimeBeforeEvent, time2 = TimeAfterEvent,
                        type = 'interval2'),
    SurvLeftTrunc = Surv(time = TimeFirstObs, time2 = TimeSurv, event = Event)
  )

# Overview of individual survivals
surv.paths <- survData %>%
  filter(Observed) %>%
  mutate(IPP = factor(IPP, IPP[order(TimeSurv)])) %>%
  ggplot(aes(colour = Pheno2C19)) +
  geom_segment(aes(y = IPP, yend = IPP, x = 0, xend = TimeFirstObs),
               linetype = 'dotted') +
  geom_segment(aes(y = IPP, yend = IPP, x = TimeFirstObs, xend = TimeSurv)) +
  geom_point(aes(y = IPP, x = TimeBeforeEvent), shape = 1) +
  geom_point(aes(y = IPP, x = TimeAfterEvent), shape = 4) +
  geom_vline(xintercept = 365, color = 'grey', linetype = 'dashed') +
  theme(axis.text.y = element_text(size = 6),
        legend.position = 'bottom')
svg(file.path(outdir, 'surv_paths.svg'), width = 14, height = 35)
print(surv.paths)
dev.off()


# KM analyses
# Log rank test with left truncated data:
# see https://stat.ethz.ch/pipermail/r-help/2009-August/399999.html
km.fit <- list(
  naive    = survfit(SurvNaive ~ Pheno2C19, survData, conf.type = 'plain'),
  interval = survfit(SurvInterval ~ Pheno2C19, survData, conf.type = 'plain'),
  ltrunc   = survfit(SurvLeftTrunc ~ Pheno2C19, survData, conf.type = 'plain')
)

# Log rank test p-values
km.sdf.naive <- survdiff(SurvNaive ~ Pheno2C19, survData)
1 - pchisq(km.sdf.naive$chisq, length(km.sdf.naive$n) - 1)
cox.fit.naive <- coxph(SurvNaive ~ Pheno2C19, survData)
summary(cox.fit.naive)$sctest['pvalue']
cox.fit.ltrunc <- coxph(SurvLeftTrunc ~ Pheno2C19, survData)
summary(cox.fit.ltrunc)$sctest['pvalue']

# KM curves
km.figs <- lapply(names(km.fit) %>% setNames(., .), function(s) {
  fit <- km.fit[[s]]
  lapply(c(full = 1, truncated = 2), function(k) {
    if (k == 1) {
      p <- ggsurvplot(
        fit,
        conf.int = TRUE,
        risk.table = TRUE,
        risk.table.col = "strata",
        ggtheme = theme_bw(),
        palette = wes_palette('Darjeeling1', length(fit$strata))
      )
    } else {
      p <- ggsurvplot(
        fit,
        conf.int = TRUE,
        risk.table = TRUE,
        risk.table.col = "strata",
        ggtheme = theme_bw(),
        palette = wes_palette('Darjeeling1', length(fit$strata)),
        xlim = c(0, 400),
        break.time.by = 30
      )
    }
    p$plot <- p$plot +
      geom_vline(xintercept = 365, color = 'grey', linetype = 'dashed')
    ttl <- c(
      naive    = 'Right censoring',
      interval = 'Interval censoring',
      ltrunc   = 'Left truncation and right censoring'
    )[s]
    p + labs(title = ttl)
  })
})
km.figs <- unlist(km.figs, recursive = FALSE)
for (s in names(km.figs)) {
  fig <- km.figs[[s]]
  svg(file.path(outdir, paste0('km_', sub('\\.', '_', s), '.svg')),
      width = 14, height = 7)
  print(fig)
  dev.off()
}
rm(fig, s)

# Differences at 12 months
# https://dominicmagirr.github.io/post/2022-01-18-be-careful-with-standard-errors-in-survival-survfit/
diff.12m <- do.call(rbind, lapply(c('naive', 'ltrunc'), function(u) {
  fit <- km.fit[[u]]
  tidy(fit) %>%
    mutate(
      strata = sub('^Pheno2C19=', '', strata),
      var.survival = estimate^2 * std.error^2,
      var.log.survival = std.error^2,
    ) %>%
    group_by(strata) %>%
    filter(time == max(time[time <= 365])) %>%
    {
      s <- .$strata
      e <- .$estimate
      v <- .$var.survival
      w <- .$var.log.survival
      do.call(rbind, lapply(1:(length(s) - 1), function(i) {
        do.call(rbind, lapply((i + 1):length(s), function(j) {
          z.plain <- abs(e[i] - e[j]) / sqrt(v[i] + v[j])
          z.log <- abs(log(e[i]) - log(e[j])) / sqrt(w[i] + w[j])
          data.frame(
            fit = u,
            strata1 = s[i],
            strata2 = s[j],
            estimate1 = e[i],
            estimate2 = e[j],
            variance1 = v[i],
            variance2 = v[j],
            p.value.plain = 2 * (1 - pnorm(z.plain)),
            p.value.log = 2 * (1 - pnorm(z.log))
          )
        }))
      }))
    }
}))
rownames(diff.12m) <- NULL
write_xlsx(diff.12m, file.path(outdir, 'diff_12m.xlsx'))

###############################################################################
# survfit(SurvLeftTrunc ~ Pheno2C19, survData, conf.type = 'log') %>%
#   tidy() %>%
#   mutate(ci.upr = exp(log(estimate) + qnorm(0.975) * std.error),
#          ci.lwr = exp(log(estimate) - qnorm(0.975) * std.error)) %>%
#   select(estimate, std.error, conf.high, ci.upr, conf.low, ci.lwr)%>%
#   as.data.frame()
###############################################################################

# Session info
sink(file.path(outdir, 'sessionInfo.txt'))
print(sessionInfo(), locale = FALSE)
sink()
