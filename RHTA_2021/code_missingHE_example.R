
####################################################################################
#R script for fitting a selection model with missingHE
#guidance for fitting different models in terms of specification, type and missingness assumptions can be found in 
#the online vignettes of the package available at: https://cran.r-project.org/web/packages/missingHE/vignettes/

####################################################################################

##STEP 0: download/install JAGS & preface

#JAGS is a program for analysis of Bayesian models using Markov Chain Monte Carlo (MCMC) methods which is used to 
#fit all models in missingHE. The latest version of \texttt{JAGS} can be dowloaded at \url{http://mcmc-jags.sourceforge.net/} and can be installed on multiple OS.

#R2jags is an R package that allows to fit JAGS models from within R and post-process the output of the model. R2jags is the default 
#package used by missingHE in order to fit the models in JAGS and must therefore be installed and loaded in your own R working directory
#before installing missingHE

#install.packages("R2jags", dependencies = TRUE) #install R2jags (only if not already installed)
#install.packages("missingHE") #install missingHE (only if not already installed)
#devtools::install_github("AnGabrio/missingHE", build_vignettes = TRUE) install dev version of the package 

#We will provide an example on how to perform missing data analysis in missingHE using the in-built CEA dataset related to a pilot RCT: the MenSS study (type help(MenSS) for more info about the study data).


##STEP 1: look at the MenSS data

#The MenSS dataset is directly available after installing and loading the package.
#MenSS is a pilot RCT evaluating the cost-effectiveness of a new digital intervention to reduce the incidence of STI in young men with respect to the SOC.

library(missingHE)

#summarise the data variables 

summary(MenSS)

#see empirical distributions of QALYs and total costs in intervention (t=2) and control (t=1) group

par(mfrow = c(2, 2))
hist(MenSS$e[MenSS$t == 1], main = "", xlab = "QALYs (t=1)", xlim = c(0.5, 1), ylim = c(0, 20))
hist(MenSS$e[MenSS$t == 2], main = "", xlab = "QALYs (t=2)", xlim = c(0.5, 1), ylim = c(0, 20))
hist(MenSS$c[MenSS$t == 1], main = "", xlab = "costs (t=1)", xlim = c(0, 1000), ylim = c(0, 20))
hist(MenSS$c[MenSS$t == 2], main = "", xlab = "costs (t=2)", xlim = c(0, 1000), ylim = c(0, 20))
par(mfrow=c(1, 1))


##STEP 2: fit the model

#We use the function selection  to specify a joint bivariate Normal distribution for the QALYs (e) and total cost (c) variables. 
#We inlcude the baseline utilities (u.0) as covariates in the model of the QALYs to account for the potential imbalance between treatment groups in these variables. 
#We then assume a MNAR structure for the missingness mechanism of the QALYs (me) while keeping a MAR assumption for the missingness model of the total costs (mc) in both treatment groups. 
#In missingHE, the model is implemented usig the following command

model.sel.pre <- selection(data = MenSS, dist_e = "norm", dist_c = "norm",
                       model.eff = e ~ u.0, model.cost = c ~ e,
                       model.me = me ~ e, model.mc = mc ~ 1,
                       type = "MNAR", n.iter = 2000, prior = "default")

#The argumnets of the function have the following interpretations:

#data: must contain the data to analyse, specified in a dataframe format.

#dist_e and dist_c: indicate the assumed effectiveness and total cost distributions, specified as character names
#among a set of predefined choices. Available choices include (but not limited to): Normal ("norm") for
#both outcomes, Beta ("beta") for the effectiveness and Gamma ("gamma") or LogNormal ("lnorm") for the costs.

#model.eff and model.cost: are formulas that specify which variables should be included in the effectiveness and total cost models as covariates.
#A joint bivariate distribution can be assumed by placing e on the right hand side of the
#formula for the total costs. By default, both formulas do not contain any covariate and the model assumes independence between the outcomes.

#model.me and model.mc: are formulas that specify which variables should be included in the effectiveness and total cost models as covariates. 
#It is possible to specify a MNAR mechanism by placing e and c in the formulas for the missing effectiveness and total cost models, respectively.
#By default no covariates are included, implicitly assuming a MAR mechanism.

#type: specifies the type of missing data mechanism, either MAR or MNAR, respectively indicated by their character values.

#n.iter: specifies the number of iterations in each chain of the MCMC algorithm.

#prior: specifies the prior distribtions to be used for the parameters of the model, by default all non-informative. These priors can be overwritten by the user. 

#Other arguments that may be provided are:
#the burnin period to be discarded (n.burnin), 
#the number of the chains (n.chains),
#the thinning interval({n.thin}), 
#the initialised values for the parameters in each chain (inits), 
#the upper and lower bounds of the credible intervals for describing the uncertainty around the imputed values (prob),
#if the model text file should be saved in the current working directory (save_model).

#missingHE allows the user to overwrite the default hyperprior values in the prior distributions for all parameters in the model. 
#For example, we can put informative priors on the parameters determining the association between e and the missingness mechanism as follows: 

my.prior <- list("delta.prior.e" = c(5, 1))

#This changes the default standard normal prior to a normal with mean 5 and sd 1.
#The object my.prior can then be passed as argument to the selection function and run

model.sel.final <- selection(data = MenSS, dist_e = "norm", dist_c = "norm",
                        model.eff = e ~ u.0, model.cost = c ~ e,
                        model.me = me ~ e, model.mc = mc ~ 1,
                        type = "MNAR", n.iter = 2000, prior = my.prior, ppc = TRUE)

#the argument ppc requires to save the posterior predictive results of the model for model assessment (optional)


##STEP 3: summarise the results and assess the model convergence

#print summaries of the posterior results for mean QALYs and total cost parameters from the model:

print(model.sel.final)

#coef is an alternative summary function for all regression parameters

coef(model.sel.final)

#we can also plot the posterior distribution of all imputed outcome data using 

plot(model.sel.final, class = "scatter", outcome = "all")

#class: specifies the type of plot to be displayed (scatter or histogram)

#outcome: specifies for which variable/group to display (all, only effects/costs, only intervention/control)

#We can check summaries of cost-effectiveness using summary:

summary(model.sel.final)

#graphical MCMC diagnostics can be obatined for each model parameter. 
#For example, we can use the diagnostic function to examine posterior traceplots for the mean effects by group:

diagnostic(model.sel.final, type = "traceplot", param = "mu.e")

#or density plots

diagnostic(model.sel.final, type = "denplot", param = "mu.e")

#type: specifies the type of diagnostic to display (see help(diagnostic) for an overview)

#param: specifies the parameter(s) to monitor


##STEP 4: assess model fit to observed data

#Model assessment can be done through posterior predictive checks (absolute fit) and predictive information criteria (relative fit).
#PPCs compare "predictions" from the model with respect to the observed data (not very helpful under MNAR!). 
#For example, we can compare densities of replicated and observed data using the command ppc:

ppc(model.sel.final, type = "dens_overlay", outcome = "all", ndisplay = 25)

#type: specifies the type of check to display

#outcome: specifies the type of variable to display (effects/costs, intervention/control)

#ndisplay: number of replications to show

#PICs are used to assess the relative fit to the observed data of two competing models (again not very helpful under MNAR!).
#Available Bayesian criteria include: DIC, WAIC and LOOIC (see help(pic) for more info).
#For example, we can obtain the DIC from the model by typing:

pic.sel.final<-pic(model.sel.final, criterion = "dic")
pic.sel.final$dic

#criterion: type of PIC to compute


##STEP 5: evaluate cost-effectiveness

#Finally, standardised CEA output (e.g. CE planes, net monetary benefit, etc.) can be obtained by post-processing the results from the model fitted in missingHE using functions from the BCEA package (which must be installed and loaded first).
#For example, CE acceptability curves can be computed by applying the function ceac.plot to the object cea stored inside the model output:

library(BCEA)
ceac.plot(model.sel.final$cea)


###############################################################################

#A comprehensive guide to the use of the package and the interpretation of the output is provided through a series of online vignettes available at: http://127.0.0.1:29713/session/Rvig.26142a03761d.html
#Three main vignettes are currently available:
#Introduction to missingHE: a general guide to the use of the main functions of the package
#Fitting MNAR models in missingHE: how to specify MNAR assumptions for each type of modelling approach available
#Model Customisation in missingHE: how to customise the model in different ways

#Instructions on how to use missingHE to fit and assess different types of models can be accessed by typing help on the different functions of the package 

#A short course with practicals and solutions is available at https://github.com/AnGabrio/short-course

#Package code available at https://github.com/AnGabrio/missingHE

