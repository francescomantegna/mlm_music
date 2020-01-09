
library("here")
library("tidyverse")  # for plotting and data manipulation
library("lme4")
library("car")
library("effects")
# library("caret") # for multicollinarity check
library("MASS") # for contr.sdif
library("lattice")
library("R.matlab")
library("eegUtils") # from https://github.com/craddm/eegUtils/
library("dplyr")
options(contrasts = c("contr.Sum","contr.Poly"))

#' For the voltage, we have 59 preprocessed channels:
#' 64 electrodes - 4 oculograms - 1 mastoid
dat <- read_csv(here("mlm_inputmusic.csv"))

dat$condition <- factor(dat$condition)


# this sets the contrasts for condition as aug4 < dominant < tonic
contrasts(dat$condition) <- contr.sdif(c("aug4","dominant","tonic"))

# we omit the multicollinearity check because we have an orthogonal design
# X <- model.matrix(~ BS + condition * x * y * z,data=dat)
# lcs <- findLinearCombos(X)
# colnames(X)[lcs$remove]

# so that I don't have to remember the different options different optimizers have...
bobyqa <- lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=1e6))
nlopt <- lmerControl(optimizer="nloptwrap",optCtrl=list(maxeval=1e6))

# these were extracted from a previous model fit. they will speed up future fits if they're right,
# but not prevent convergence to the correct values if they're wrong
# they will speed up

#+ theta
theta <- c(0.18648574040049,-0.0443144981555326,0.0307819282071667,
           0.161675512694921,-0.0355269822704287,0.180234962926686,
           0.0625669832864284,-0.00162793523354175,0.0225908831188544,
           0.103424412839264,-0.035597248934163,0.117197823819823)

# this gives a warning about scaling; however, the model still converges on the
# original scale, so I would leave it. If you add scaling in, then I would only
# do it for N5 and BS and not the for scalp topography

#+ condition_model, cache=TRUE, depends="theta"
system.time(m <- lmer(N5 ~ 1 + BS + condition * group * x * y * z +
                        (1 + condition | subject_id) +
                        (1 + condition | scale_id),
                      data=dat,
                      control=nlopt,
                      start=theta,
                      REML=FALSE))

print(summary(m),correlation=FALSE,symbolic.cor=TRUE)

s <- summary(m)
capture.output(s, file = "mlmsummary.txt")

sprintf("Number of fixed-effect correlations > 0.1: %d",sum(as.matrix(vcov(m)) > 0.1))

png('model.png')
plot(m)
dev.off()

png('qqplot.png')
qqmath(m)
dev.off()

png('byitem.png')
qqmath(ranef(m,condVar=TRUE))
dev.off()

png('dotplot.png')
dotplot(ranef(m,condVar=TRUE))
dev.off()

#+ diagnostics, cache=TRUE, fig.width=10, fig.height=10
png('resvsfit.png')
fortify.merMod(m, drop_na(dat)) %>%
  ggplot(aes(.fitted,.resid)) +
  geom_point(aes(color=condition),alpha=0.3) +
  facet_wrap( ~ subject_id) +
  geom_hline(yintercept=0) +
  theme_light() +
  labs(title="Residuals vs. Fitted",
       subtitle="Should be cloud shaped",
       x="Fitted",
       y="Residuals")
dev.off()

#' For the ANOVA-style display, we can also omit the baseline interval because
#' we don't actually care about it, even if we have to model it.
#+ anova, cache = TRUE
(a <- Anova(m))
a[!str_detect(rownames(a),"BS"),]
capture.output(a, file = "ANOVA.txt")
#' The above summary can also be expressed graphically.

# Wald is the fastest method; boot the most accurate and slowest, profile is a
# the middle road
ci <- confint(m,method="Wald")

#+ coefplot, cache=TRUE, fig.width=7, fig.height=5
ci.gg <- ci %>%
  as_tibble(rownames = "coef") %>%
  dplyr::filter(substr(coef,1,4) != ".sig") %>% # omit random effects
  mutate(est = fixef(m),
         coef = factor(coef, levels=rev(names(fixef(m))))) %>%
  dplyr::filter(coef != "(Intercept)") %>% # intercept
  dplyr::filter(!str_detect(coef,"BS")) %>% # baseline terms
  mutate(topo=as.logical(str_detect(coef,"x|y|z")))

png('ggplot.png')
gg <- ggplot(mapping=aes(x=coef,y=est,ymin=`2.5 %`, ymax=`97.5 %`)) +
  geom_pointrange(fatten=1.5) +
  geom_hline(yintercept=0, linetype="dashed") +
  labs(x="Coefficient",y="Estimate (standardized difference)") +
  coord_flip() +
  theme_light()

gg %+% subset(ci.gg, !topo) + ggtitle("Global Effects")

gg %+% subset(ci.gg, as.logical(str_detect(coef,"x"))) + ggtitle("Lateral Effects")
gg %+% subset(ci.gg, as.logical(str_detect(coef,"y"))) + ggtitle("Saggital Effects")
dev.off()
