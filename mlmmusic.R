#'---
#'title: "Prediction in Music: Continuous Topography"
#'author:
#' - Francesco Mantegna
#' - Phillip M. Alday
#'date: "June 2021 (from commit: `r system('git rev-parse --short HEAD', intern = TRUE)`)"
#'---

library("here")
library("tidyverse")  # for plotting and data manipulation
library("lme4")
library("car")
library("effects")
# library("caret") # for multicollinarity check
library("MASS") # for contr.sdif
library("lattice")
# library("R.matlab")
# library("eegUtils") # from https://github.com/craddm/eegUtils/
library("knitr")
options(contrasts = c("contr.Sum","contr.Poly"))

if(is_html_output()) options(width=120)
knitting <- is_html_output() || is_latex_output()

#' For the voltage, we have 59 preprocessed channels:
#' 64 electrodes - 4 oculograms - 1 mastoid
dat <- read_csv(here("mlm_inputmusic.csv"))

dat$condition <- factor(dat$condition)

#' This sets the contrasts for condition as `aug4 < dominant` and `mean(aug4, dominant) < tonic`
ch <- contr.Helmert(c("aug4","dominant","tonic"))
print(ch)
colnames(ch) <- c("dominant>aug4", "tonic>mean(dominant,aug4)")
contrasts(dat$condition) <- ch
print(ch)

#' We omit the multicollinearity check because we have an orthogonal design.

# X <- model.matrix(~ BS + condition * x * y * z,data=dat)
# lcs <- findLinearCombos(X)
# colnames(X)[lcs$remove]

# so that I don't have to remember the different options different optimizers have...
bobyqa <- lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=1e6))
nlopt <- lmerControl(optimizer="nloptwrap",optCtrl=list(maxeval=1e6))

#' # Confirmatory: N5
#
#' These were extracted from a previous model fit. they will speed up future
#' fits if they're right, but not prevent convergence to the correct values if
#' they're wrong.

#+ theta_n5
theta_n5 <- c(0.18648574040049,-0.0443144981555326,0.0307819282071667,
           0.161675512694921,-0.0355269822704287,0.180234962926686,
           0.0625669832864284,-0.00162793523354175,0.0225908831188544,
           0.103424412839264,-0.035597248934163,0.117197823819823)

#' This gives a warning about scaling; however,the model still converges on the
#' original scale, so I would leave it. If you add scaling in, then I would only
#' do it for N5 and BS and not the for scalp topography.

#+ condition_model_n5, cache=TRUE, dependson="theta_n5"
system.time(m_n5 <- lmer(N5 ~ 1 + BS + condition * group * x * y * z +
                        (1 + condition | subject_id) +
                        (1 + condition | scale_id),
                      data=dat,
                      control=nlopt,
                      start=theta_n5,
                      REML=FALSE))

print(summary(m_n5),correlation=FALSE,symbolic.cor=TRUE)

s <- summary(m_n5)
capture.output(s, file = "mlmsummary.txt")

sprintf("Number of fixed-effect correlations > 0.1: %d",sum(as.matrix(vcov(m_n5)) > 0.1))

if(!knitting) png('model_n5.png')
plot(m_n5)
if(!knitting) dev.off()

if(!knitting) png('qqplot_n5.png')
qqmath(m_n5)
if(!knitting) dev.off()

if(!knitting) png('byitem_n5.png')
qqmath(ranef(m_n5,condVar=TRUE))
if(!knitting) dev.off()

if(!knitting) png('dotplot_n5.png')
dotplot(ranef(m_n5,condVar=TRUE))
if(!knitting) dev.off()

#+ diagnostics_n5, cache=TRUE, fig.width=10, fig.height=10
if(!knitting) png('resvsfit_n5.png')
fortify.merMod(m_n5, drop_na(dat)) %>%
  ggplot(aes(.fitted,.resid)) +
  geom_point(aes(color=condition),alpha=0.3) +
  facet_wrap( ~ subject_id) +
  geom_hline(yintercept=0) +
  theme_light() +
  labs(title="Residuals vs. Fitted",
       subtitle="Should be cloud shaped",
       x="Fitted",
       y="Residuals")
if(!knitting) dev.off()

#' For the ANOVA-style display, we focus on terms involving condition and group.
#+ anova_n5, cache = TRUE, dependson="condition_model"
kable(a <- Anova(m_n5))
kable(a[str_detect(rownames(a),"condition|group") & a$`Pr(>Chisq)` < 0.05,])
capture.output(a, file = "ANOVA.txt")

#' The above summary can also be expressed graphically.

# Wald is the fastest method; boot the most accurate and slowest, profile is a
# the middle road
ci <- confint(m_n5,method="Wald")
# prof <- profile(m_n5, parallel="multicore", ncpus=4)
# ci <- confint(prof)

#+ coefplot_n5, cache=TRUE, fig.width=7, fig.height=5
ci.gg <- ci %>%
  as_tibble(rownames = "coef") %>%
  dplyr::filter(substr(coef,1,4) != ".sig") %>% # omit random effects
  mutate(est = fixef(m_n5),
         coef = factor(coef, levels=rev(names(fixef(m_n5))))) %>%
  dplyr::filter(coef != "(Intercept)") %>% # intercept
  dplyr::filter(!str_detect(coef,"BS")) %>% # baseline terms
  mutate(topo=as.logical(str_detect(coef,"x|y|z")))

gg <- ggplot(mapping=aes(x=coef,y=est,ymin=`2.5 %`, ymax=`97.5 %`)) +
  geom_pointrange(fatten=1.5) +
  geom_hline(yintercept=0, linetype="dashed") +
  labs(x="Coefficient",y="Estimate (standardized difference)") +
  coord_flip() +
  theme_light()

if(!knitting) png('coef_global_n5.png')
gg %+% subset(ci.gg, !topo) + ggtitle("Global Effects")
if(!knitting) dev.off()

if(!knitting) png('coef_lateral_n5.png')
gg %+% subset(ci.gg, as.logical(str_detect(coef,"x"))) + ggtitle("Lateral Effects")
if(!knitting) dev.off()

if(!knitting) png('coef_saggital_n5.png')
gg %+% subset(ci.gg, as.logical(str_detect(coef,"y"))) + ggtitle("Saggital Effects")
if(!knitting) dev.off()

roi <- c("F7", "F3", "Fz", "F4", "F8", "FT7", "FC3", "FC4", "FT8")

ref_grid <- dat %>% dplyr::select(condition,x,y,z, channel) %>% unique()

ref_grid_roi <- subset(ref_grid, channel %in% roi)

# we need max x because that's symmetric lateral spread
max_x <- max(ref_grid_roi$x)
# we need max and min y that's saggital spread
max_y <- max(ref_grid_roi$y)
min_y <- min(ref_grid_roi$y)
# we don't need to worry about max/min z because we get that for free

dat_roi <- subset(dat, abs(x) <= max_x & y < max_y  & y > min_y)
dat_roi$BS <- 0

# re.form=~0 gets rid of the random effects, so that we're focused on the
# estimated differences *at the population level*. Not that this also means
# we're essentially doing marginalized fixed effect CIs and not prediction
# intervals for future data
set.seed(42)
dat_roi$sim <- simulate(m_n5, nsim=1000, newdata=dat_roi, re.form=~0)

roi_effs <- dat_roi %>% group_by(condition, group) %>%
  # colMeans marginalizes across electrodes in each ROI within each simulation
  # so you then compute stats across simulations
  summarize(mean=mean(colMeans(sim)),
            lower=quantile(colMeans(sim),probs=0.085),
            upper=quantile(colMeans(sim), probs=0.915))

ggplot(roi_effs, aes(x=condition,
                     color=group,fill=group,
                     y=mean,
                     ymin=lower, ymax=upper)) +
  geom_pointrange(position=position_dodge(width=0.1)) +
  labs(title="Estimated effects in frontal ROI for the N5 time window",
       subtitle="with 83% confidence intervals",
       y="µV",
       x="Condition",
       color="Group",
       fill="Group") +
  theme_light()

#' # Exploratory: P3
#
#' These were extracted from a previous model fit. they will speed up future
#' fits if they're right, but not prevent convergence to the correct values if
#' they're wrong.

#+ theta_p3
theta_p3 <- c(0.18648574040049,-0.0443144981555326,0.0307819282071667,
              0.161675512694921,-0.0355269822704287,0.180234962926686,
              0.0625669832864284,-0.00162793523354175,0.0225908831188544,
              0.103424412839264,-0.035597248934163,0.117197823819823)

#' This gives a warning about scaling; however,the model still converges on the
#' original scale, so I would leave it. If you add scaling in, then I would only
#' do it for P3 and BS and not the for scalp topography.

#+ condition_model_p3, cache=TRUE, dependson="theta_p3"
system.time(m_p3 <- lmer(P3 ~ 1 + BS + condition * group * x * y * z +
                           (1 + condition | subject_id) +
                           (1 + condition | scale_id),
                         data=dat,
                         control=bobyqa,
                         start=theta_p3,
                         REML=FALSE))

print(summary(m_p3),correlation=FALSE,symbolic.cor=TRUE)

s <- summary(m_p3)
capture.output(s, file = "mlmsummary.txt")

sprintf("Number of fixed-effect correlations > 0.1: %d",sum(as.matrix(vcov(m_p3)) > 0.1))

if(!knitting) png('model_p3.png')
plot(m_p3)
if(!knitting) dev.off()

if(!knitting) png('qqplot_p3.png')
qqmath(m_p3)
if(!knitting) dev.off()

if(!knitting) png('byitem_p3.png')
qqmath(ranef(m_p3,condVar=TRUE))
if(!knitting) dev.off()

if(!knitting) png('dotplot_p3.png')
dotplot(ranef(m_p3,condVar=TRUE))
if(!knitting) dev.off()

#+ diagnostics_p3, cache=TRUE, fig.width=10, fig.height=10
if(!knitting) png('resvsfit_p3.png')
fortify.merMod(m_p3, drop_na(dat)) %>%
  ggplot(aes(.fitted,.resid)) +
  geom_point(aes(color=condition),alpha=0.3) +
  facet_wrap( ~ subject_id) +
  geom_hline(yintercept=0) +
  theme_light() +
  labs(title="Residuals vs. Fitted",
       subtitle="Should be cloud shaped",
       x="Fitted",
       y="Residuals")
if(!knitting) dev.off()

#' For the ANOVA-style display, we focus on terms involving condition and group.
#+ anova_p3, cache = TRUE, dependson="condition_model"
kable(a <- Anova(m_p3))
kable(a[str_detect(rownames(a),"condition|group") & a$`Pr(>Chisq)` < 0.05,])
capture.output(a, file = "ANOVA.txt")

#' The above summary can also be expressed graphically.

# Wald is the fastest method; boot the most accurate and slowest, profile is a
# the middle road
ci <- confint(m_p3,method="Wald")
# prof <- profile(m_p3, parallel="multicore", ncpus=4)
# ci <- confint(prof)

#+ coefplot_p3, cache=TRUE, fig.width=7, fig.height=5
ci.gg <- ci %>%
  as_tibble(rownames = "coef") %>%
  dplyr::filter(substr(coef,1,4) != ".sig") %>% # omit random effects
  mutate(est = fixef(m_p3),
         coef = factor(coef, levels=rev(names(fixef(m_p3))))) %>%
  dplyr::filter(coef != "(Intercept)") %>% # intercept
  dplyr::filter(!str_detect(coef,"BS")) %>% # baseline terms
  mutate(topo=as.logical(str_detect(coef,"x|y|z")))

gg <- ggplot(mapping=aes(x=coef,y=est,ymin=`2.5 %`, ymax=`97.5 %`)) +
  geom_pointrange(fatten=1.5) +
  geom_hline(yintercept=0, linetype="dashed") +
  labs(x="Coefficient",y="Estimate (standardized difference)") +
  coord_flip() +
  theme_light()

if(!knitting) png('coef_global_p3.png')
gg %+% subset(ci.gg, !topo) + ggtitle("Global Effects")
if(!knitting) dev.off()

if(!knitting) png('coef_lateral_p3.png')
gg %+% subset(ci.gg, as.logical(str_detect(coef,"x"))) + ggtitle("Lateral Effects")
if(!knitting) dev.off()

if(!knitting) png('coef_saggital_p3.png')
gg %+% subset(ci.gg, as.logical(str_detect(coef,"y"))) + ggtitle("Saggital Effects")
if(!knitting) dev.off()

roi <- c("F7", "F3", "Fz", "F4", "F8", "FT7", "FC3", "FC4", "FT8")

ref_grid <- dat %>% dplyr::select(condition,x,y,z, channel) %>% unique()

ref_grid_roi <- subset(ref_grid, channel %in% roi)

# we need max x because that's symmetric lateral spread
max_x <- max(ref_grid_roi$x)
# we need max and min y that's saggital spread
max_y <- max(ref_grid_roi$y)
min_y <- min(ref_grid_roi$y)
# we don't need to worry about max/min z because we get that for free

dat_roi <- subset(dat, abs(x) <= max_x & y < max_y  & y > min_y)
dat_roi$BS <- 0

# re.form=~0 gets rid of the random effects, so that we're focused on the
# estimated differences *at the population level*. Not that this also means
# we're essentially doing marginalized fixed effect CIs and not prediction
# intervals for future data
set.seed(42)
dat_roi$sim <- simulate(m_p3, nsim=1000, newdata=dat_roi, re.form=~0)

roi_effs <- dat_roi %>% group_by(condition, group) %>%
  # colMeans marginalizes across electrodes in each ROI within each simulation
  # so you then compute stats across simulations
  summarize(mean=mean(colMeans(sim)),
            lower=quantile(colMeans(sim),probs=0.085),
            upper=quantile(colMeans(sim), probs=0.915))

ggplot(roi_effs, aes(x=condition,
                     color=group,fill=group,
                     y=mean,
                     ymin=lower, ymax=upper)) +
  geom_pointrange(position=position_dodge(width=0.1)) +
  labs(title="Estimated effects in frontal ROI for the P3 time window",
       subtitle="with 83% confidence intervals",
       y="µV",
       x="Condition",
       color="Group",
       fill="Group") +
  theme_light()

#' # Exploratory: N5 from P3
#'

#+ theta_n5_from_p3
theta_n5_from_p3 <- theta_n5

#' This gives a warning about scaling; however,the model still converges on the
#' original scale, so I would leave it. If you add scaling in, then I would only
#' do it for P3 and BS and not the for scalp topography.

#+ condition_model_n5_from_p3, cache=TRUE, dependson="theta_n5_from_p3"
system.time(m_n5_from_p3_additive <- lmer(N5 ~ 1 + BS + P3 + condition * group * x * y * z +
                                   (1 + condition | subject_id) +
                                   (1 + condition | scale_id),
                                 data=dat,
                                 control=bobyqa,
                                 start=theta_n5_from_p3,
                                 REML=FALSE))

system.time(m_n5_from_p3_bycond <- lmer(N5 ~ 1 + BS + P3 * condition * group * x * y * z +
                                   (1 + condition | subject_id) +
                                   (1 + condition | scale_id),
                                 data=dat,
                                 control=bobyqa,
                                 start=getME(m_n5_from_p3_additive,"theta"),
                                 REML=FALSE))

anova(m_n5, m_n5_from_p3_additive, m_n5_from_p3_bycond)

m_n5_from_p3 <- m_n5_from_p3_bycond

print(summary(m_n5_from_p3),correlation=FALSE,symbolic.cor=TRUE)

s <- summary(m_n5_from_p3)
capture.output(s, file = "mlmsummary.txt")

sprintf("Number of fixed-effect correlations > 0.1: %d",sum(as.matrix(vcov(m_n5_from_p3)) > 0.1))

if(!knitting) png('model_n5_from_p3.png')
plot(m_n5_from_p3)
if(!knitting) dev.off()

if(!knitting) png('qqplot_n5_from_p3.png')
qqmath(m_n5_from_p3)
if(!knitting) dev.off()

if(!knitting) png('byitem_n5_from_p3.png')
qqmath(ranef(m_n5_from_p3,condVar=TRUE))
if(!knitting) dev.off()

if(!knitting) png('dotplot_n5_from_p3.png')
dotplot(ranef(m_n5_from_p3,condVar=TRUE))
if(!knitting) dev.off()

#+ diagnostics_n5_from_p3, cache=TRUE, fig.width=10, fig.height=10
if(!knitting) png('resvsfit_n5_from_p3.png')
fortify.merMod(m_n5_from_p3, drop_na(dat)) %>%
  ggplot(aes(.fitted,.resid)) +
  geom_point(aes(color=condition),alpha=0.3) +
  facet_wrap( ~ subject_id) +
  geom_hline(yintercept=0) +
  theme_light() +
  labs(title="Residuals vs. Fitted",
       subtitle="Should be cloud shaped",
       x="Fitted",
       y="Residuals")
if(!knitting) dev.off()

#' For the ANOVA-style display, we focus on terms involving condition and group.
#+ anova_n5_from_p3, cache = TRUE, dependson="condition_model"
kable(a <- Anova(m_n5_from_p3))
kable(a[str_detect(rownames(a),"condition|group") & a$`Pr(>Chisq)` < 0.05,])
capture.output(a, file = "ANOVA.txt")

#' The above summary can also be expressed graphically.

# Wald is the fastest method; boot the most accurate and slowest, profile is a
# the middle road
ci <- confint(m_n5_from_p3,method="Wald")
# prof <- profile(m_n5_from_p3, parallel="multicore", ncpus=4)
# ci <- confint(prof)

#+ coefplot_n5_from_p3, cache=TRUE, fig.width=7, fig.height=5
ci.gg <- ci %>%
  as_tibble(rownames = "coef") %>%
  dplyr::filter(substr(coef,1,4) != ".sig") %>% # omit random effects
  mutate(est = fixef(m_n5_from_p3),
         coef = factor(coef, levels=rev(names(fixef(m_n5_from_p3))))) %>%
  dplyr::filter(coef != "(Intercept)") %>% # intercept
  dplyr::filter(!str_detect(coef,"BS")) %>% # baseline terms
  mutate(topo=as.logical(str_detect(coef,"x|y|z")))

gg <- ggplot(mapping=aes(x=coef,y=est,ymin=`2.5 %`, ymax=`97.5 %`)) +
  geom_pointrange(fatten=1.5) +
  geom_hline(yintercept=0, linetype="dashed") +
  labs(x="Coefficient",y="Estimate (standardized difference)") +
  coord_flip() +
  theme_light()

if(!knitting) png('coef_global_n5_from_p3.png')
gg %+% subset(ci.gg, !topo) + ggtitle("Global Effects")
if(!knitting) dev.off()

if(!knitting) png('coef_lateral_n5_from_p3.png')
gg %+% subset(ci.gg, as.logical(str_detect(coef,"x"))) + ggtitle("Lateral Effects")
if(!knitting) dev.off()

if(!knitting) png('coef_saggital_n5_from_p3.png')
gg %+% subset(ci.gg, as.logical(str_detect(coef,"y"))) + ggtitle("Saggital Effects")
if(!knitting) dev.off()

roi <- c("F7", "F3", "Fz", "F4", "F8", "FT7", "FC3", "FC4", "FT8")

ref_grid <- dat %>% dplyr::select(condition,x,y,z, channel) %>% unique()

ref_grid_roi <- subset(ref_grid, channel %in% roi)

# we need max x because that's symmetric lateral spread
max_x <- max(ref_grid_roi$x)
# we need max and min y that's saggital spread
max_y <- max(ref_grid_roi$y)
min_y <- min(ref_grid_roi$y)
# we don't need to worry about max/min z because we get that for free

dat_roi <- subset(dat, abs(x) <= max_x & y < max_y  & y > min_y)
dat_roi$BS <- 0
dat_roi$P3 <- 0

# re.form=~0 gets rid of the random effects, so that we're focused on the
# estimated differences *at the population level*. Not that this also means
# we're essentially doing marginalized fixed effect CIs and not prediction
# intervals for future data
set.seed(42)
dat_roi$sim <- simulate(m_n5_from_p3, nsim=1000, newdata=dat_roi, re.form=~0)

roi_effs <- dat_roi %>% group_by(condition, group) %>%
  # colMeans marginalizes across electrodes in each ROI within each simulation
  # so you then compute stats across simulations
  summarize(mean=mean(colMeans(sim)),
            lower=quantile(colMeans(sim),probs=0.085),
            upper=quantile(colMeans(sim), probs=0.915))

ggplot(roi_effs, aes(x=condition,
                     color=group,fill=group,
                     y=mean,
                     ymin=lower, ymax=upper)) +
  geom_pointrange(position=position_dodge(width=0.1)) +
  labs(title="Estimated effects in frontal ROI for the N5 time window after adjusting for the P3",
       subtitle="with 83% confidence intervals",
       y="µV",
       x="Condition",
       color="Group",
       fill="Group") +
  theme_light()

