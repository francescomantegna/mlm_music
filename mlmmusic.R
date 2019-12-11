
library("here")
library("tidyverse")  # for plotting and data manipulation
library("lme4")
library("car")
library("effects")
library("caret")
library("lattice")
library("R.matlab")
library("eegUtils") # from https://github.com/craddm/eegUtils/
library("dplyr")
options(contrasts = c("contr.Sum","contr.Poly"))

#' For the voltage, we have 59 preprocessed channels:
#' 64 electrodes - 4 oculograms - 1 mastoid
dat <- read_csv("mlm_inputmusic.csv")

X <- model.matrix(~ BS + condition * x * y * z,data=dat)
lcs <- findLinearCombos(X)
colnames(X)[lcs$remove]

# so that I don't have to remember the different options different optimizers have...
bobyqa <- lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=1e6))
nlopt <- lmerControl(optimizer="nloptwrap",optCtrl=list(maxeval=1e6))

system.time(m <- lmer(scale(N5) ~ 1 + scale(BS) + condition * group * scale(x) * scale(y) * scale(z) +
                        (1 + condition | subject_id) +
                        (1 + condition | scale_id),
                      data=dat,
                      control=bobyqa,
                      REML=FALSE))

print(summary(m),correlation=FALSE,symbolic.cor=TRUE)

s <- summary(m)
capture.output(s, file = "mlmsummary.txt")

sprintf("Number of fixed-effect correlations > 0.1: %d",sum(as.matrix(vcov(m)) > 0.1))

jpeg('model.jpg')
plot(m)
dev.off()

jpeg('qqplot.jpg')
qqmath(m)
dev.off()

jpeg('byitem.jpg')
qqmath(ranef(m,condVar=TRUE))
dev.off()

jpeg('dotplot.jpg')
dotplot(ranef(m,condVar=TRUE))
dev.off()

#+ diagnostics, cache=TRUE, fig.width=10, fig.height=10
jpeg('resvsfit.jpg')
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

jpeg('ggplot.jpg')
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
