#'---
#'title: "Prediction in Music: Averaging within ROI"
#'author:
#' - Phillip M. Alday
#' - Francesco Mantegna
#'date: "June 2021 (from commit: `r system('git rev-parse --short HEAD', intern = TRUE)`)"
#'---

library("here")
library("tidyverse")  # for plotting and data manipulation
library("lme4")
library("car")
library("effects")
library("MASS") # for contr.sdif
library("lattice")
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

#' We restrict further analysis to the ROI seen in studies by Koelsch and colleagues

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
dat_roi <- dat_roi %>%
    group_by(condition, group, subject_id, scale_id) %>%
    summarize_at(vars(N5, P3, BS), mean)

#' We omit the multicollinearity check because we have an orthogonal design.

# X <- model.matrix(~ BS + condition * x * y * z,data=dat)
# lcs <- findLinearCombos(X)
# colnames(X)[lcs$remove]

# so that I don't have to remember the different options different optimizers have...
bobyqa <- lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=1e6))
nlopt <- lmerControl(optimizer="nloptwrap",optCtrl=list(maxeval=1e6))

#' # Confirmatory: N5


#+ condition_model_n5, cache=TRUE
system.time(m_n5 <- lmer(N5 ~ 1 + BS + condition * group  +
                             (1 + condition | subject_id) +
                             (1 + condition | scale_id),
                         data=dat_roi,
                         control=bobyqa,
                         REML=FALSE))

print(summary(m_n5),correlation=FALSE,symbolic.cor=TRUE)
print(summary(rePCA(m_n5)))
print(vif(m_n5))

s <- summary(m_n5)
capture.output(s, file = "mlmsummary_roi.txt")

sprintf("Number of fixed-effect correlations > 0.1: %d",sum(as.matrix(vcov(m_n5)) > 0.1))

plot(m_n5)
qqmath(m_n5)
qqmath(ranef(m_n5,condVar=TRUE))
dotplot(ranef(m_n5,condVar=TRUE))

#' For the ANOVA-style display, we focus on terms involving condition and group.
#+ anova_n5, cache = TRUE
kable(a <- Anova(m_n5))
kable(a[str_detect(rownames(a),"condition|group") & a$`Pr(>Chisq)` < 0.05,])

#' The above summary can also be expressed graphically.

# Wald is the fastest method; boot the most accurate and slowest, profile is a
# the middle road
ci <- confint(m_n5,method="Wald")
# prof <- profile(m_n5, parallel="snow", ncpus=4, devtol=1e-06)
# ci <- confint(prof)

#+ coefplot_n5, cache=TRUE, fig.width=7, fig.height=5
ci.gg <- ci %>%
    as_tibble(rownames = "coef") %>%
    dplyr::filter(substr(coef,1,4) != ".sig") %>% # omit random effects
    mutate(est = fixef(m_n5),
           coef = factor(coef, levels=rev(names(fixef(m_n5))))) %>%
    dplyr::filter(coef != "(Intercept)") %>% # intercept
    dplyr::filter(!str_detect(coef,"BS")) # baseline terms

gg <- ggplot(mapping=aes(x=coef,y=est,ymin=`2.5 %`, ymax=`97.5 %`)) +
    geom_pointrange(fatten=1.5) +
    geom_hline(yintercept=0, linetype="dashed") +
    labs(x="Coefficient",y="Estimate (standardized difference)") +
    coord_flip() +
    theme_light()

gg %+% ci.gg + ggtitle("N5")


#' # Exploratory: P3
#

#+ condition_model_p3, cache=TRUE
system.time(m_p3 <- lmer(P3 ~ 1 + BS + condition * group +
                             (1 + condition | subject_id) +
                             (1 + condition | scale_id),
                         data=dat_roi,
                         control=bobyqa,
                         REML=FALSE))

print(summary(rePCA(m_p3)))
#' Okay, we need to redo the RE
system.time(m_p3 <- lmer(P3 ~ 1 + BS + condition * group +
                             (1 + condition | subject_id) +
                             (1  | scale_id),
                         data=dat_roi,
                         control=bobyqa,
                         REML=FALSE))

print(summary(m_p3),correlation=FALSE,symbolic.cor=TRUE)
print(summary(rePCA(m_p3)))
print(vif(m_p3))


s <- summary(m_p3)
capture.output(s, file = "mlmsummary.txt")

sprintf("Number of fixed-effect correlations > 0.1: %d",sum(as.matrix(vcov(m_p3)) > 0.1))

plot(m_p3)
qqmath(m_p3)
qqmath(ranef(m_p3,condVar=TRUE))
dotplot(ranef(m_p3,condVar=TRUE))

#' For the ANOVA-style display, we focus on terms involving condition and group.
#+ anova_p3, cache = TRUE
kable(a <- Anova(m_p3))
kable(a[str_detect(rownames(a),"condition|group") & a$`Pr(>Chisq)` < 0.05,])

#' The above summary can also be expressed graphically.

# Wald is the fastest method; boot the most accurate and slowest, profile is a
# the middle road
ci <- confint(m_p3, method="Wald")
# prof <- profile(m_p3, parallel="snow", ncpus=4)
# ci <- confint(prof)

#+ coefplot_p3, cache=TRUE, fig.width=7, fig.height=5
ci.gg <- ci %>%
    as_tibble(rownames = "coef") %>%
    dplyr::filter(substr(coef,1,4) != ".sig") %>% # omit random effects
    mutate(est = fixef(m_p3),
           coef = factor(coef, levels=rev(names(fixef(m_n5))))) %>%
    dplyr::filter(coef != "(Intercept)") %>% # intercept
    dplyr::filter(!str_detect(coef,"BS")) # baseline terms

gg <- ggplot(mapping=aes(x=coef,y=est,ymin=`2.5 %`, ymax=`97.5 %`)) +
    geom_pointrange(fatten=1.5) +
    geom_hline(yintercept=0, linetype="dashed") +
    labs(x="Coefficient",y="Estimate") +
    coord_flip() +
    theme_light()

gg %+% ci.gg + ggtitle("P3")

#' # Exploratory: N5 from P3
#'

#+ condition_model_n5_from_p3, cache=TRUE
system.time(m_n5_from_p3_additive <- lmer(N5 ~ 1 + BS + P3 + condition * group +
                                              (1 + condition | subject_id) +
                                              (1 + condition | scale_id),
                                          data=dat_roi,
                                          control=bobyqa,
                                          REML=FALSE))
print(summary(rePCA(m_n5_from_p3_additive)))
print(vif(m_n5_from_p3_additive))

system.time(m_n5_from_p3_bycond <- lmer(N5 ~ 1 + BS + P3 * condition * group +
                                            (1 + condition | subject_id) +
                                            (1 + condition | scale_id),
                                        data=dat_roi,
                                        control=bobyqa,
                                        start=getME(m_n5_from_p3_additive,"theta"),
                                        REML=FALSE))

print(summary(rePCA(m_n5_from_p3_bycond)))
print(vif(m_n5_from_p3_bycond))

anova(m_n5, m_n5_from_p3_additive, m_n5_from_p3_bycond)

m_n5_from_p3 <- m_n5_from_p3_bycond

print(summary(m_n5_from_p3),correlation=FALSE,symbolic.cor=TRUE)

s <- summary(m_n5_from_p3)

sprintf("Number of fixed-effect correlations > 0.1: %d",sum(as.matrix(vcov(m_n5_from_p3)) > 0.1))

plot(m_n5_from_p3)
qqmath(m_n5_from_p3)
qqmath(ranef(m_n5_from_p3,condVar=TRUE))
dotplot(ranef(m_n5_from_p3,condVar=TRUE))

#' For the ANOVA-style display, we focus on terms involving condition and group.
#+ anova_n5_from_p3, cache = TRUE
kable(a <- Anova(m_n5_from_p3))
kable(a[str_detect(rownames(a),"condition|group") & a$`Pr(>Chisq)` < 0.05,])

#' The above summary can also be expressed graphically.

# Wald is the fastest method; boot the most accurate and slowest, profile is a
# the middle road
ci <- confint(m_n5_from_p3,method="Wald")
# prof <- profile(m_n5_from_p3, parallel="snow", ncpus=4)
# ci <- confint(prof)

#+ coefplot_n5_from_p3, cache=TRUE, fig.width=7, fig.height=5
ci.gg <- ci %>%
    as_tibble(rownames = "coef") %>%
    dplyr::filter(substr(coef,1,4) != ".sig") %>% # omit random effects
    mutate(est = fixef(m_n5_from_p3),
           coef = factor(coef, levels=rev(names(fixef(m_n5))))) %>%
    dplyr::filter(coef != "(Intercept)") %>% # intercept
    dplyr::filter(!str_detect(coef,"BS")) # baseline terms

gg <- ggplot(mapping=aes(x=coef,y=est,ymin=`2.5 %`, ymax=`97.5 %`)) +
    geom_pointrange(fatten=1.5) +
    geom_hline(yintercept=0, linetype="dashed") +
    labs(x="Coefficient",y="Estimate (standardized difference)") +
    coord_flip() +
    theme_light()

gg %+% ci.gg + ggtitle("N5 from P3")

