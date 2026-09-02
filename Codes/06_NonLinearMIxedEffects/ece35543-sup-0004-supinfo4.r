
# Packages
library(nlme)
library(ggplot2)
library(bbmle)
library(emmeans)

# ---- Data handling --------------------------------------------- 

# Database
lfmc.data <- read.table("data_lfmc.txt", header = TRUE)
str(lfmc.data)

# Data are ordered for temporal analysis
lfmc.data <- lfmc.data[with(lfmc.data,order(site, plot, time)), ] 

# "Plot" is coded as factor
lfmc.data$plot <- as.factor(lfmc.data$plot)
str(lfmc.data)

# The levels of "leaf type" are ordered with GW as the base 
lfmc.data$leaf.type = with(lfmc.data, factor(leaf.type, levels=c("Grass W", "Grass E", "M. spinosum", "S. bracteolactus")))
lfmc.data$plot = with(lfmc.data, factor(plot, levels=c("4", "5", "6", "1", "2", "3")))
str(lfmc.data)

lfmc.data$group <- with(lfmc.data, factor(plot:leaf.type)) # one col ("group") crossing both factors ("plot" x "leaf type") is generated
levels(lfmc.data$group)
lfmc.data$group=with(lfmc.data, factor(group, labels=c("GW plot 4","GW plot 5", "GW plot 6",
                                                       "GE plot 1","SM plot 1", "SS plot 1",
                                                       "GE plot 2","SM plot 2", "SS plot 2",
                                                       "GE plot 3","SM plot 3", "SS plot 3"))) # the factor "group" is relabeled according to the manuscript

lfmc.data.gd <- groupedData(lfmc ~ time | group, data = lfmc.data, order.groups=FALSE) # and a groupedData class object is created 
str(lfmc.data.gd)
lfmc.data.gd$plot = with(lfmc.data.gd, factor(plot, levels=c("4", "5", "6", "1", "2", "3"))) # "plot" is coded as factor in the groupedData object

# ---- Data exploration ------------------------------------------- 

par(mfrow=c(1,1), cex=0.85)
plot(lfmc.data$time, lfmc.data$lfmc, xlab="Time (days)", ylab="LFMC (%)", ylim=c(0,350), main="") # LFMC decreases over time
par(mfrow=c(1,2), cex=0.85)
plot(lfmc.data$leaf.type, lfmc.data$lfmc, xlab="Leaf type", ylab="LFMC (%)", ylim=c(0,350), col="grey", main="") # Shrubs are more heterogeneous
hist(lfmc.data$lfmc, main="", xlab="LFMC (%)", col="grey")

ggplot(data = lfmc.data, aes(x = time, y = lfmc)) + 
       geom_point() + 
       geom_smooth() +
       facet_wrap(~leaf.type) + 
       xlab("Time (days)") + 
       ylab("LFMC (%)") # The LFMC temporal dynamics is plotted by "leaf type" 

par(mfrow=c(1,1))
plot(lfmc.data.gd, xlab="Time (days)", ylab="LFMC (%)", 
     aspect = 0.25, layout=c(3,4), col="black") # The LFMC temporal dynamic of each "leaf type" are plotted by "plot" 


# 1)  NONLINEAR MIXED-EFFECTS MODEL (M1)  #########################################################

source("selfstart.R") # the selfStart function is run 

# 1.1 ==== Nonlinear modeling proccess ============================================================ 

# "Model fitting may stop with errors, or produce obscure warnings, or get stuck at an obviously bad fit,
# depending on subtle changes in the way the model or the starting values are specified" (Bolker et al. 2013).

# Here, some of these problems are illustrated in the first steps. 


# 1.1.1 - Random-effects model (based on a "nlsList" model) ----------------------------------------  
nl.re <- nlme(nlsList(lfmc ~ SSdlf(time, A, w, m, s), data = lfmc.data.gd), control = nlmeControl(maxIter = 5000, msMaxIter = 1500)) # all effects on A, w, m, s are assumed as random. The model variance-covariance is unstructured.   

# .- In this code, the selfStart function (SSdlf) is applied within a nlsList() function. 
# .- nlslList() fits one nls() fixed-effects model for each group ("plot" x "leaf type") as a intermediate step in the fit of the nlme() function.

# The previous model issues a warning. It is important to recognize that the unstructured 
# variance-covariance meaning that 10 (4 + 3 + 2 + 1) parameters are tried to be estimated for
# the random effects. This is very likely to result in an over-parameterized model. 

# Therefore, 
nl.re.1 <- update(nl.re, random = pdDiag(A + w + m + s ~ 1)) # a random-effects model assuming zero correlation among the random-effects is fitted 
# This converges easilty and this means the previous model (nl.re) is over parameterized. 

fxf <- fixef(nl.re.1) # the intercepts of the previous model are obtained to be used in the next step. 


# 1.1.2 - Mixed-effects model ("leaf type" is included as a fixed-effect on A, W, M, S) ----------- 
nl.me <- update(nl.re.1, fixed = list(A + w + m + s ~ leaf.type), 
                start = c(A=fxf[1],0,0,0, w=fxf[2],0,0,0, m=fxf[3],0,0,0, s=fxf[4],0,0,0)) # convergence problems

nl.re.1 # Standar deviation (StdDev) of "s" is small. 

# Therefore, the random-effects on "s" could be removed to achieve convergence 
nl.re.2 <- update(nl.re, random = pdDiag(A + w + m ~ 1))  
fxf <- fixef(nl.re.2)

# Now, the model including the fixed-effects of "leaf type" converges 
nl.me.1 <- update(nl.re.2, fixed = list(A + w + m + s ~ leaf.type),  
                  start = c(A=fxf[1],0,0,0, w=fxf[2],0,0,0, m=fxf[3],0,0,0, s=fxf[4],0,0,0))

# the covariance structure of the random effects is much smaller when the fixed-effects are modeled
nl.me.1
#            A      w       m     Residual
# StdDev:  8.11   0.001   2.47     15.39
nl.re.2
#            A      w       m     Residual
# StdDev: 119.38  10.94   8.93     15.26

# ...suggesting parameters to vary with "leaf type"  
AICtab(nl.re.2, nl.me.1) # which is supported by the AIC comparison

# The small StdDev of "w" suggests to remove the random-effect of plot on this parameter  
nl.me.2 <- update(nl.me.1, random = pdDiag(A + m ~ 1))  
AICtab(nl.me.1, nl.me.2) # and the AIC comparison supports the simpler model

# in addition, if the random-effects on "m" is also removed   
nl.me.3 <- update(nl.me.2, random = pdDiag(A ~ 1)) 
AICtab(nl.me.2, nl.me.3) # the model fit is similar, supporting the simpler model

# Residual checking:
plot(nl.me.3) # heterocedasticity is clear
par(mfrow=c(1,2), cex=0.85)
hist(resid(nl.me.3, type="normalized"), main="", xlab="Residuals", col="grey") 
qqnorm(resid(nl.me.3, type="normalized"))
qqline(resid(nl.me.3, type="normalized"))

# 1.1.3 - Heterocedasticity: variance modeling ---------------------------------------------------- 
nl.me.3.vm <- update(nl.me.3, weights = varIdent(form=~1|leaf.type)) 
AICtab(nl.me.3, nl.me.3.vm) # the fit is clearly improved when variances are modeled

# Residual checking
plot(nl.me.3.vm)
par(mfrow=c(1,2), cex=0.85)
hist(resid(nl.me.3.vm, type="normalized"), main="", xlab="Residuals", col="grey") 
qqnorm(resid(nl.me.3.vm, type="normalized"))
qqline(resid(nl.me.3.vm, type="normalized"))
# the residuals look fine now

# 1.1.4 - Temporal correlation --------------------------------------------------------------------
plot(ACF(nl.me.3.vm), alpha = 0.01) # This plot suggests there is no need for modeling the correlation of the residuals


# 1.1.5 - Evaluation of the fixed-effects ---------------------------------------------------------  
nl.me.ml <- update(nl.me.3.vm, method ="ML") # the model is fitted by Maximum likelihood 

nl.me.ml.1 <- update(nl.me.ml, fixed = list(A + w + m ~ leaf.type, s ~ 1),  
                     start = c(A=fxf[1],0,0,0, w=fxf[2],0,0,0, m=fxf[3],0,0,0, s=fxf[4])) # the model is reduced removing the fixed-effect of "leaf type" on "s" 
AICtab(nl.me.ml, nl.me.ml.1) # the AIC comparison suggests that the more complex model is not justified 

nl.me.ml.2 <- update(nl.me.ml, fixed = list(A + w ~ leaf.type, m + s ~ 1),  
                     start = c(A=fxf[1],0,0,0, w=fxf[2],0,0,0, m=fxf[3], s=fxf[4])) # the model is reduced removing the fixed-effect of "leaf type" on "m" 
AICtab(nl.me.ml.1, nl.me.ml.2) # the more complex model is not justified   

nl.me.ml.3 <- update(nl.me.ml, fixed = list(A ~ leaf.type, w + m + s  ~ 1),  
                     start = c(A=fxf[1],0,0,0, w=fxf[2], m=fxf[3], s=fxf[4])) # the model is reduced removing the fixed-effect of "leaf type" on "W"                    
AICtab(nl.me.ml.2, nl.me.ml.3) # the AIC comparison suggests to "w" to be modeled as a function of "leaf type" 

nl.me.ml.4 <- update(nl.me.ml, fixed = list(w ~ leaf.type, A + m + s  ~ 1),  
                     start = c(A=fxf[1], w=fxf[2],0,0,0, m=fxf[3], s=fxf[4])) # the model is reduced removing the fixed-effect of "leaf type" on "A"                    
AICtab(nl.me.ml.2, nl.me.ml.4) # the AIC comparison suggests to "A" to be modeled as a function of "leaf type" 

                 
# 1.2 ==== Final model =============================================================================  

M1 <- update(nl.me.ml.2, method = "REML")
summary(M1)

# This model assumes:

# .- "A" and "w" vary with leaf type (fixed-effects).
# .- "A" varies with plot (random-effects).   
# .- "m" and "s" are unique for all the plots and leaf types. 
# .- Variance depends on "leaf type"
# .- There is no temporal dependence
# .- Residuals follow a normal distribution

fixef(M1) # Fixed effects 
ranef(M1) # Random effects
intervals(M1) # Confidence intervals


# 1.2.1 - Link between R output and (eq. 9)/Table 2 in manuscript (see Appendix 1) ----------------

fixef(M1)[1] #  "µA0" in table 2 and eq. 9
fixef(M1)[2] #  "A1"  in table 2 and eq. 9 
fixef(M1)[3] #  "A2"  in table 2 and eq. 9
fixef(M1)[4] #  "A3"  in table 2 and eq. 9
fixef(M1)[5] #  "µw0" in table 2 and eq. 9
fixef(M1)[6] #  "w1"  in table 2 and eq. 9 
fixef(M1)[7] #  "w2"  in table 2 and eq. 9
fixef(M1)[8] #  "w3"  in table 2 and eq. 9
fixef(M1)[9] #  "m"   in table 2 and eq. 9
fixef(M1)[10] # "s"   in table 2 and eq. 9


# 1.2.2 - Overall predictions (Table 3) ----------------------------------------------------------- 
Agw <- fixef(M1)[1] # "A" for grasses in the W site (GW): -> "µA0" in table 2 and eq. 9   
Age <- fixef(M1)[1]+fixef(M1)[2] # "A" for grasses in the E site (GE)
Asm <- fixef(M1)[1]+fixef(M1)[3] # "A" for Mullinum spinosum (SM)
Ass <- fixef(M1)[1]+fixef(M1)[4] # "A" for Senecio filaginoides (SS)
Wgw <- fixef(M1)[5] # "w" for grasses in the W site (GW) 
Wge <- fixef(M1)[5]+fixef(M1)[6] # "w" for grasses in the E site (GE) 
Wsm <- fixef(M1)[5]+fixef(M1)[7] # "w" for Mullinum spinosum (SM)
Wss <- fixef(M1)[5]+fixef(M1)[8] # "w" for Senecio filaginoides (SS)
M <- fixef(M1)[9] # "m" for all the leaf types
S <- fixef(M1)[10] # "s" for all the leaf types

GW <- subset(lfmc.data, leaf.type == "Grass W")
GE <- subset(lfmc.data, leaf.type == "Grass E")
SM <- subset(lfmc.data, leaf.type == "M. spinosum")
SS <- subset(lfmc.data, leaf.type == "S. bracteolactus")

par(mfrow=c(2,2), cex=0.75)
# Grasses W
plot(GW$time, GW$lfmc, xlab="", ylab="LFMC (%)", ylim=c(0,100), main="Grasses W site")
curve(((Agw-Wgw)/(1 + exp((M-x)/S))+Wgw), lwd=2, lty=1, add=T, col="red") 
# Grasses E
plot(GE$time, GE$lfmc, xlab="", ylab="LFMC (%)", ylim=c(0,100), main="Grasses E site")
curve(((Age-Wge)/(1 + exp((M-x)/S))+Wge), lwd=2, lty=1, add=T, col="red")
# M. Spinosum
plot(SM$time, SM$lfmc, xlab="Time (days)", ylab="LFMC (%)", ylim=c(0,350), main="M. spinosum", font.main=4)
curve(((Asm-Wsm)/(1 + exp((M-x)/S))+Wsm), lwd=2, lty=1, add=T, col="red") 
# S. filaginoides 
plot(SS$time, SS$lfmc, xlab="Time (days)", ylab="LFMC (%)", ylim=c(0,350), main="S. bracteolactus", font.main=4)
curve(((Ass-Wss)/(1 + exp((M-x)/S))+Wss), lwd=2, lty=1, add=T, col="red") 


# 1.2.3 - Predictions at plot level (Table 3) -----------------------------------------------------
Agw.p4 <- fixef(M1)[1]+ranef(M1)[1,1] # "A" for grasses in the plot 4 (GW) 
Agw.p5 <- fixef(M1)[1]+ranef(M1)[2,1] # "A" for grasses in the plot 5 (GW) 
Agw.p6 <- fixef(M1)[1]+ranef(M1)[3,1] # "A" for grasses in the plot 6 (GW) 
Age.p1 <- fixef(M1)[1]+fixef(M1)[2]+ranef(M1)[4,1] # "A" for grasses in the plot 1 (GW) 
Age.p2 <- fixef(M1)[1]+fixef(M1)[2]+ranef(M1)[7,1] # "A" for grasses in the plot 2 (GW) 
Age.p3 <- fixef(M1)[1]+fixef(M1)[2]+ranef(M1)[10,1] # "A" for grasses in the plot 3 (GW) 
Asm.p1 <- fixef(M1)[1]+fixef(M1)[3]+ranef(M1)[5,1] # "A" for M. spinosum in the plot 1 (SM) 
Asm.p2 <- fixef(M1)[1]+fixef(M1)[3]+ranef(M1)[8,1] # "A" for M. spinosum in the plot 2 (SM)
Asm.p3 <- fixef(M1)[1]+fixef(M1)[3]+ranef(M1)[11,1] # "A" for M. spinosum in the plot 3 (SM)
Ass.p1 <- fixef(M1)[1]+fixef(M1)[4]+ranef(M1)[6,1] # "A" for S. filaginoides in the plot 1 (SM)   
Ass.p2 <- fixef(M1)[1]+fixef(M1)[4]+ranef(M1)[9,1] # "A" for S. filaginoides in the plot 2 (SM) 
Ass.p3 <- fixef(M1)[1]+fixef(M1)[4]+ranef(M1)[12,1] # "A" for S. filaginoides in the plot 3 (SM)
  
par(mfrow=c(2,2), cex=0.75) # (Figure 4)
# Grasses W site 
plot(GW$time, GW$lfmc, xlab="", ylab="LFMC (%)", ylim=c(0,100), main="Grasses W site")
curve(((Agw.p4-Wgw)/(1 + exp((M-x)/S))+Wgw), lwd=1, lty=2, add=T, col="red") # plot 4
curve(((Agw.p5-Wgw)/(1 + exp((M-x)/S))+Wgw), lwd=1, lty=2, add=T, col="red") # plot 5
curve(((Agw.p6-Wgw)/(1 + exp((M-x)/S))+Wgw), lwd=1, lty=2, add=T, col="red") # plot 6
# Grasses E site
plot(GE$time, GE$lfmc, xlab="", ylab="LFMC (%)", ylim=c(0,100), main="Grasses E site")
curve(((Age.p1-Wge)/(1 + exp((M-x)/S))+Wge), lwd=1, lty=2, add=T, col="red") # plot 1
curve(((Age.p2-Wge)/(1 + exp((M-x)/S))+Wge), lwd=1, lty=2, add=T, col="red") # plot 2
curve(((Age.p3-Wge)/(1 + exp((M-x)/S))+Wge), lwd=1, lty=2, add=T, col="red") # plot 3
# M. Spinosum (E site) 
plot(SM$time, SM$lfmc, xlab="Time (days)", ylab="LFMC (%)", ylim=c(0,350), main="M. spinosum", font.main=4)
curve(((Asm.p1-Wsm)/(1 + exp((M-x)/S))+Wsm), lwd=1, lty=2, add=T, col="red") # plot 1
curve(((Asm.p2-Wsm)/(1 + exp((M-x)/S))+Wsm), lwd=1, lty=2, add=T, col="red") # plot 2
curve(((Asm.p3-Wsm)/(1 + exp((M-x)/S))+Wsm), lwd=1, lty=2, add=T, col="red") # plot 3
# S. filaginoides (E site)
plot(SS$time, SS$lfmc, xlab="Time (days)", ylab="LFMC (%)", ylim=c(0,350), main="S. bracteolactus", font.main=4)
curve(((Ass.p1-Wss)/(1 + exp((M-x)/S))+Wss), lwd=1, lty=2, add=T, col="red") # plot 1
curve(((Ass.p2-Wss)/(1 + exp((M-x)/S))+Wss), lwd=1, lty=2, add=T, col="red") # plot 2
curve(((Ass.p3-Wss)/(1 + exp((M-x)/S))+Wss), lwd=1, lty=2, add=T, col="red") # plot 3


# 1.2.4 - Derivatives (drying speed) -------------------------------------------------------------- 
lfmc.GW <- expression((Agw-Wgw)/(1 + exp((M-x)/S))+Wgw) # LFMC dynamics for grasses in the W site
ds.GW <- deriv(lfmc.GW, "x", function.arg = TRUE) # Drying speed for grasses in the W site

lfmc.GE <- expression((Age-Wge)/(1 + exp((M-x)/S))+Wge) # LFMC dynamics for grasses in the E site
ds.GE <- deriv(lfmc.GE, "x", function.arg = TRUE) # Drying speed for grasses in the E site

lfmc.SM <- expression((Asm-Wsm)/(1 + exp((M-x)/S))+Wsm) # LFMC dynamics for M. spinosum
ds.SM <- deriv(lfmc.SM, "x", function.arg = TRUE) # Drying speed for M. spinosum

lfmc.SS <- expression((Ass-Wss)/(1 + exp((M-x)/S))+Wss) # LFMC dynamics for S. filaginoides
ds.SS <- deriv(lfmc.SS, "x", function.arg = TRUE) # # Drying speed for S. filaginoides

xvec <- seq(0, 100, length = 1000)
y.ds.GW <- ds.GW(xvec)
y.ds.GE <- ds.GE(xvec)
y.ds.SM <- ds.SM(xvec)
y.ds.SS <- ds.SS(xvec)

par(mfrow=c(2,2), cex=0.75) # (Figure 4)
plot(xvec, y.ds.GW, xlim=c(0,80), ylim=c(0,4), 
     ylab="Drying speed", xlab="Time (days)", type="n", main="Grasses in the W site")
lines(xvec, -1*(attr(y.ds.GW, "grad")), lwd=2, lty=1, col="red")
plot(xvec, y.ds.GE, xlim=c(0,80), ylim=c(0,4), 
     ylab="Drying speed", xlab="Time (days)", type="n", main="Grasses in the E site")
lines(xvec, -1*(attr(y.ds.GE, "grad")), lwd=2, lty=1, col="red")
plot(xvec, y.ds.SM, xlim=c(0,80), ylim=c(0,4), 
     ylab="Drying speed", xlab="Time (days)", type="n", main="M. spinosum", font.main=4)
lines(xvec, -1*(attr(y.ds.SM, "grad")), lwd=2, lty=1, col="red")
plot(xvec, y.ds.SS, xlim=c(0,80), ylim=c(0,4), 
     ylab="Drying speed", xlab="Time (days)", type="n", main="S. filaginoides", font.main=4)
lines(xvec, -1*(attr(y.ds.SS, "grad")), lwd=2, lty=1, col="red")


# 1.3 ---- Testing differences among leaf types --------------------------------------------------- 

# 1.3.1 - Parameter "A" ----
contrast(emmeans(M1, ~leaf.type, param = "A"), "pairwise")

# .- Maximum LFMC in grasses from the W site is lower than grasses from the E site (Grass W - Grass E)

# .- Maximum LFMC in grasses is lower than shrubs (Grass W - M. spinosum)
#                                                 (Grass W - S. bracteolactus)
#                                                 (Grass E - M. spinosum)
#                                                 (Grass E - S. bracteolactus)

# .- There is no difference between the two shrub species (M. spinosum - S. bracteolactus).


# 1.3.2 - Parameter "w" --------------------------------------------------------------------------- 
contrast(emmeans(M1, ~leaf.type, param = "w"), "pairwise")

# .- Minimum LFMC in grasses from the W site is higher than grasses from the E site (Grass W - Grass E)

# .- Minimum LFMC in grasses is lower than shrubs (Grass W - M. spinosum)
#                                                 (Grass W - S. bracteolactus)
#                                                 (Grass E - M. spinosum)
#                                                 (Grass E - S. bracteolactus)

# .- There is no difference between the two shrub species (M. spinosum - S. bracteolactus).


# 2) NONLINEAR FIXED-EFFECTS MODEL (M2)  ##########################################################

# Response function:  lfmc = (A - w) / (1 + exp((m - time)/s))) + w 


# 2.1 ==== Fit of the nonlinear fixed-effects model =============================================== 

# 2.1.1 - Base nonlinear model --------------------------------------------------------------------
nl.fe.0 <- nls(lfmc ~ ((A - w) / (1 + exp((m - time)/s))) + w, 
               start = c(A=15, m=30, s=-17, w=5),
               data = lfmc.data) # here, time is the only predictor variable 

b <- coef(nl.fe.0) # coefficients of the base model

# 2.1.2 - Leaf type fixed effect ------------------------------------------------------------------
nl.fe.1 <- nls(lfmc ~ ((A[leaf.type] - w[leaf.type])/(1 + exp((m - time)/s))) + w[leaf.type], 
              start = list(A=rep(b[1], 4), m=25, s=-10, w=rep(b[4], 4)),
              data = lfmc.data) # and fit is made using the base model's coefficients as start values

b1 <- coef(nl.fe.1) # coefficients of the model with leaf type and time as predictor variables

# 2.1.3 - Plot fixed effect ----------------------------------------------------------------------- 
nl.fe.2 <- nls(lfmc ~ ((A[group] - w[group])/(1 + exp((m - time)/s))) + w[group], 
               start = list(A=rep(c(b1[1],b1[2],b1[3],b1[4]),3), m=25, s=-10, w=rep(c(b1[7],b1[8],b1[9],b1[10]),3)),
               data = lfmc.data) # the model is fit using "b1" as the start values

AICtab(nl.fe.0, nl.fe.1, nl.fe.2) # including leaf type and plot as predictor variables imrpoves the model fit 

# Residual checking:
plot(nl.fe.2) # heterocedasticity in the residuals
hist(resid(nl.fe.2), main="", xlab="Residuals") # slightly skew distribution   
qqnorm(resid(nl.fe.2))
qqline(resid(nl.fe.2))


# 2.2 ==== Final model ============================================================================= 

M2 <- nl.fe.2
summary(M2)

# This model assumes:

# .- "A" and "w" vary with leaf type and plot (fixed-effects).
# .- "m" and "s" are unique for all the plots and leaf types. 
# .- Variance is homogeneous
# .- Residuals are independent (there is no temporal dependence) 
# .- Residuals following a normal distribution


# 3 - LINEAR MIXED-EFFECTS MODEL (M3) ############################################################## 

# Response function:  lfmc = β0 + β1xTime + β2xGE + β3xSM + β4xSS
#                            β5xGExTime + β6xSMxTime + β7xSSxTime

# where GE, SM, and SS are dummy variables (0 or 1) created to indicate differences between 
# the base level of "leaf type", i.e., grasses from the W site [GW], and the remaining
# levels. 


# 3.1 ==== Modeling process of the linear mixed-effects model ===================================== 

# 3.1.1 - Base linear mixed-effect model ----------------------------------------------------------
l.me <- lme(lfmc ~ time * leaf.type, random = ~1 | plot, data = lfmc.data)

# Residual checking:
plot(l.me) # heterocedasticity in the residuals
plot(ACF(l.me, resType = "normalized", na.action=na.omit), alpha=0.05, grid=TRUE) # temporal correlation is observed at lag 1


# 3.1.2 - Variance modelling ---------------------------------------------------------------------- 
l.me.vm <- update(l.me, weights = varIdent(form=~ 1 | leaf.type))
    
# Residual checking:
plot(l.me.vm) # variances there seem to be well modeled
AIC(l.me, l.me.vm) # modeling the variance improves the model fit 


# 3.1.3 - Temporal correlation modelling ----------------------------------------------------------  
l.me.vm.tm <- update(l.me.vm, correlation=corARMA(p=1,q=1, form=~1))
plot(ACF(l.me.vm.tm, resType = "normalized", na.action=na.omit), alpha=0.05, grid=TRUE) # the temporal dependence is modeled
AIC(l.me.vm, l.me.vm.tm) # modeling the temporal correlation improves the model fit

# Residual checking:
plot(l.me.vm.tm) 
plot(ACF(l.me.vm.tm, resType = "normalized", na.action=na.omit), alpha=0.05, grid=TRUE) 
hist(resid(l.me.vm.tm, type="normalized"), main="", xlab="Residuals") 
qqnorm(resid(l.me.vm.tm, type="normalized"))
qqline(resid(l.me.vm.tm, type="normalized"))
# the residuals look fine


# 3.1.4 - Evaluation of the fixed-effects --------------------------------------------------------- 
l.me.ml <- update(l.me.vm.tm, method ="ML") # the model is fitted by Maximum likelihood 
l.me.ml.1 <- update(l.me.ml, ~.-time:leaf.type) # the model is reduced removing the interaction "time x leaf type" 
AICtab(l.me.ml, l.me.ml.1) # the AIC comparison suggests the interaction term to be important  


# 3.2 ==== Final model ============================================================================ 

M3 <- l.me.vm.tm 
summary(M3)

# This model assumes:

# .- β0 vary with plots as a random-effects (random intercept model)
# .- The residual variance depends on "leaf type"
# .- Temporal correlation
# .- Residuals follow a normal distribution


# 4 - LINEAR REGRESSION (M4) ######################################################################

M4 <- lm(lfmc ~ time * group, data = lfmc.data)
summary(M4)

# This model assumes:

# .- Variance is homogeneous
# .- There is no temporal dependence (residuals are independent) 
# .- Residuals following a normal distribution


# Residual checking:
plot(M4) # a clear pattern in the residuals is observed (model assumptions are not met) 

# 5 - NULL MODEL (M5) #############################################################################

M5 <- lm(lfmc ~ 1, data = lfmc.data)
summary(M5)


# 6 - MODEL COMPARISON ############################################################################

# Models:  

# M1: Nonlinear mixed-effects model
# M2: Nonlinear fixed-effects model
# M3: Linear mixed-effetcs model  
# M4: Lnear fixed-effects model (clasical regression) 
# M5: Null model

# 6.1 ==== AIC comparison (Table 1) ===============================================================  

AICtab(M1, M2, M3, M4, M5,  
       weights = T, delta = TRUE, base=T, sort = TRUE) # the nonlinear mixed-effects model is clearly the best

# 6.2 ==== Residual analysis ======================================================================

par(mfrow=c(2,2), cex = 0.65) # (Figure 5)
v1 <- c(1, 2, 3, 4)
v2 <- c("GW", "GE", "SM", "SS")

# Nonlinear mixed-effects model (M1) 
rM1 <- resid(M1, type = "normalized")
fM1 <- fitted(M1)
plot(fM1, rM1, xlab="Fitted LFMC (%)", ylab="Residuals", ylim=c(-3,3))
abline(a=0, b=0, lwd=2)  
boxplot(rM1 ~ lfmc.data$leaf.type, ylab="Residuals", xlab="Leaf type", xaxt="n", ylim=c(-3,3))
axis(side=1, at=v1, labels=v2, las=1)
plot(rM1 ~ lfmc.data$time, ylab="Residuals", xlab="Time (days)", ylim=c(-3,3))
abline(a=0,b=0, lw=2)
qqnorm(rM1, main="", xlab="Theoretical quantiles", ylab="Sample quantiles")
qqline(rM1)

# Nonlinear fixed-effects model (M2)
rM2 <- resid(M2)
fM2 <- fitted(M2)
plot(fM2, rM2, xlab="Fitted LFMC (%)", ylab="Residuals")
abline(a=0, b=0, lwd=2)  
boxplot(rM2 ~ lfmc.data$leaf.type, ylab="Residuals", xlab="Leaf type", xaxt="n")
axis(side=1, at=v1, labels=v2, las=1)
plot(rM2 ~ lfmc.data$time, ylab="Residuals", xlab="Time (days)")
abline(a=0,b=0, lw=2)
qqnorm(rM2, main="", xlab="Theoretical quantiles", ylab="Sample quantiles")
qqline(rM2)

# Linear mixed-effects model (M3) 
rM3 <- resid(M3, type = "normalized")
fM3 <- fitted(M3)
plot(fM3, rM3, xlab="Fitted LFMC (%)", ylab="Residuals", ylim=c(-3,3))
abline(a=0, b=0, lwd=2)  
boxplot(rM3 ~ lfmc.data$leaf.type, ylab="Residuals", xlab="Leaf type", xaxt="n", ylim=c(-3,3))
axis(side=1, at=v1, labels=v2, las=1)
plot(rM3 ~ lfmc.data$time, ylab="Residuals", xlab="Time (days)", ylim=c(-3,3))
abline(a=0,b=0, lw=2)
qqnorm(rM3, main="", ylab="Sample quantiles", xlab="Theoretical quantiles")
qqline(rM3)

# Clasical regression (M4) 
rM4 <- resid(M4)
fM4 <- fitted(M4)
plot(fM4, rM4, xlab="Fitted LFMC (%)", ylab="Residuals")
abline(a=0, b=0, lwd=2)  
boxplot(rM4 ~ lfmc.data$leaf.type, ylab="Residuals", xlab="", xaxt="n")
axis(side=1, at=v1, labels=v2, las=1)
plot(rM4 ~ lfmc.data$time, ylab="Residuals", xlab="Time (days)")
abline(a=0,b=0, lw=2)
qqnorm(rM4, main="", xlab="Theoretical quantiles", ylab="Sample quantiles")
qqline(rM4)

