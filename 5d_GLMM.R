library(tidyverse)
library(lme4)
library(lmerTest)
library(emmeans)
library(car)
library(MuMIn)
library(agridat)

data("beall.webworms")
d1 <- beall.webworms
?beall.webworms  ## info about the beall.webworms dataset
head(d1) ## view data set

hist(d1$y) ## view histogram of raw data

ggplot(d1, aes(x=spray, y=y, fill=lead)) + geom_violin(bw=.5) + stat_summary(fun=mean, geom="point", color="black", size=3, position=position_dodge(.9))

## FIT GLM
glm1 <- glm(y ~ spray * lead, data=d1, family=poisson)

plot(glm1)
hist(glm1$residuals)
plot(glm1$residuals~r3$fitted.values)  ## residuals should be evenly dispersed around 0 across the range of x's
abline(h=0) 
qqPlot(glm1$residuals)
boxplot(glm1$residuals ~ d1$trt)  ## variances should be homogeneous for each group


## diagnose model fit with DHARMa package
library(DHARMa)
sim_glm1 <- simulateResiduals(fittedModel = glm1, n = 250)
plot(sim_glm1)

plotQQunif(sim_glm1)
plotResiduals(sim_glm1)

## test for overdispersion (ie. more variance in the data than expected)
testDispersion(sim_glm1)

## test for zero-inflation (ie. more zeros in the data than expected)
testZeroInflation(sim_glm1)

##################################################################################
### fit block as a random effect
glmm1 <- glmer(y ~ spray * lead + (1|block), data=d1, family=poisson)
summary(glmm1)

## simulate and check residuals
sim_glmm1 <- simulateResiduals(fittedModel = glmm1, n = 250)
plot(sim_glmm1)

## check for overdispersion using a different method
library(RVAideMemoire)
overdisp.glmer(glmm1) # overdispersion ratio calculator from RVAideMemoire

####################################################################################
## correct for overdispersion using observation-level random effect
d1$obs <- 1:length(d1$y)

glmm2 <-glmer(y ~ spray * lead + (1|block) + (1|obs), data=d1, family=poisson)
summary(glmm2)

sim_glmm2 <- simulateResiduals(fittedModel = glmm2, n = 250)
plot(sim_glmm2)
hist(sim_glmm2)

## additional tests
testUniformity(sim_glmm2)
testDispersion(sim_glmm2)
testZeroInflation(sim_glmm2)


## compare all three models using likelihood ratio test (or AIC)
anova(glmm2,glmm1,glm1)

## Model good? Ok! Look at p-values and means
Anova(glmm2)

##emmeans
emmeans(glmm2, pairwise~spray:lead, type="response")

## compare above means to regular lm. Means from GLM look too low
emmeans(lm(y~spray*lead, data=d1), ~spray:lead, type="response")

##emmeans, correcting for bias (only needed for GLMMs) by calculating the total SD of the random effect variance components
totSD <- VarCorr(glmm2) %>% as.data.frame() %>% summarize(totSD=sum(vcov)) %>% mutate(totSD=sqrt(totSD))

emmeans(glmm2, pairwise~spray:lead, type="response", bias.adj=T, sigma=totSD$totSD)
