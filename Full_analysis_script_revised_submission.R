### This script:
###   - Imports data from experimental infection of hamsters and rats
###   - Estimates pP, the probability that each leptospire that reaches the within-host environment persists within the host 
###         --> pP is estimates using data from IP infection experiments
###   - Estimates pc, the probability that each leptospire in the inoculum crosses physical immune barriers and survives to reach the interanl tissue
###         --> pc was estimated using data for each route of inoculation (intact skin, shavbed skin, abraded skin, conjunctiva)
###         --> pc was estimated using a basic model, and using a mixture model assuming pc~Beta(alpha, beta)
###  - Calculates 95% profile CIs for each estimated paramter
###  - Generates all figures in the main text, and supplement
###  - Saves the results
###  - Outputs a Session_Info.txt file to record conditions at the time of last run, prior to submission

## Set working directory, load libraries and functions
rm(list = ls()) # Clear workspace
setwd('~/Dropbox/R/Hamsters/') # Set working directory

## Load packages
library(parallel)
library(ggplot2)
library(binom)
library(RColorBrewer)
library(gridExtra)
library(reshape2)

## Output session info
writeLines(capture.output(sessionInfo()), "code_outputs/sessionInfo.txt")

## Set up a palette that looks good in color, and also shows variation in grayscale for the printed article
cols = colorRampPalette(brewer.pal(4,"YlGnBu"))(25)
cols[26] = 'midnightblue'
plot(1:26, cex = 10, pch = 16, col = cols) # Look at all colors
plot(1:4, cex = 10, pch = 16, col = cols[c(7, 16, 24, 26) ]) # Choose a subset of four that show good divergence in color and in grayscale (print will be in grayscale)
cols = cols[c(7, 16, 24, 26) ]

####################### INPUTS ####################
full.dat = read.csv('raw_data/Hamsterdata.csv', header = T)[,1:4] # Hamster data
rat.data = read.csv('raw_data/Ratdata.csv', header = T)[1:16,] # Rat data

####################### OUTPUTS ###################
## MAIN TEXT FIGURES
## Note - Fig. 1 was a diagram drawn using MS Power Point. See template in "figures" subfolder.
fig2 = 'manuscript_figures/Fig2_par_ests.pdf'
fig3 = 'manuscript_figures/Fig3_model_fits.pdf'
fig4 = 'manuscript_figures/Fig4_prob_inf_given_dose_and_route.pdf'
#likelihood profile of pI estimates, hamsters:
plot1 = 'scratch_figures/pI_profiles.pdf' 
## Likelihood profiles, basic model, pc, hamsters:
plot2 = 'scratch_figures/Profiles_basic_model.pdf' 
figS1 = 'manuscript_figures/Profile_rat_mixture_model.png'
#### CODE OUTPUTS
mixture_prof_output = 'code_outputs/Prof_2D_hamsters.RData'
mixture_prof_rats = 'code_outputs/Prof_2D_rats.RData'
hamdat_summary = 'code_outputs/Hamster_data_table.csv'
ratdat_summary = 'code_outputs/Rat_data_table.csv'
#########################

####################################################################
##########################  SECTION 00 ##############################
####   
####   Define functions necessary to calculate 95% CIs from 1D likelihood profiles
####     These functions are called below
####
####################################################################
## This function calculates the threshold negative log likelihood (nll)
## Candidate paramters values whose profile nll better than this value are contained within the CI.
LR.Threshold = function(NLL_best, df){
  threshold = qchisq(.95, df)/2+NLL_best
  threshold
}

## This function inputs a nll threshold, a vector of profile nlls, and a corresponding vector of paramter grid points (pars.vec)
## This function outputs the highest and lowest parameter values in the grid that fall within the 95% profile CI
LR.CI = function(threshold, nll.vec, pars.vec){
if(any(which(nll.vec > threshold) < which.min(nll.vec))  &  any(which(nll.vec > threshold) > which.min(nll.vec)) ){ #If the minimum is not an end point
  #Find string before and after the min value
  lower = nll.vec[1:which.min(nll.vec)]
  upper = nll.vec[(which.min(nll.vec)-1):length(nll.vec)]
  #Extract low CI from first string, and upper CI from upper string
  CI = c(pars.vec[which.min(abs(lower-threshold))], pars.vec[length(lower)+which.min(abs(upper-threshold))])
}else{
  if(any(which(nll.vec > threshold) < which.min(nll.vec)) == FALSE){
    CI = c(pars.vec[1], pars.vec[which.min(abs(nll.vec-threshold))])
    warning('CI contains first value of profile. Consider adding lower paramter values to profile.')
  }
  if(any(which(nll.vec > threshold) > which.min(nll.vec)) == FALSE){
    CI = c(pars.vec[which.min(abs(nll.vec-threshold))], pars.vec[length(nll.vec)])
    warning('CI contains last value of profile. Consider adding higher paramter values to profile.')
  }
}
CI # Return
}



####################################################################
##########################  SECTION 0 ##############################
####   
####   Reformat data for hamsters and rats
####
####################################################################

####    Format the IP data
####################################################################
# Write a function that creates a summary table with the total number of trials and number of infections (deaths) at a given dose
reformat.route = function(route.dat){
  doses = sort(unique(route.dat$dose))
  tbl = matrix(nrow = length(doses), ncol = 4, dimnames = list(doses, c('n', 'died', 'dose', 'pct.died')))
  for(ii in 1:length(doses)){
    tbl[ii, 'n'] = sum(route.dat[route.dat$dose == doses[ii], 'n'])
    tbl[ii, 'died'] = sum(route.dat[route.dat$dose == doses[ii], 'died'])
    tbl[ii, 'dose'] = doses[ii]
    tbl[ii, 'pct.died'] = tbl[ii,'died']/tbl[ii,'n']
  }
  as.data.frame(tbl)
}
#--------------------- Hamsters: -----------------------
#For each IP dose, create a table of (dose, n total, n died)
IP.raw = subset(full.dat, route == 'IP')[,-1] ## Get rid of the first column, "route", which is known (all IP)

# Visualize the summarized data
IP = reformat.route(IP.raw); IP$route = 'IP'


### --------------------------------------------------------------
####    Format data for all other routes of inoculation, hamsters
####      Repeat the steps used for IP data
####################################################################
#Input route data
intact.raw = subset(full.dat, route == 'Int.skin')
abraded.raw = subset(full.dat, route == 'Epiderm')
shaved.raw = subset(full.dat, route == 'Shv.skin')
conjunctival.raw = subset(full.dat, route == 'Conj')

## Get tables of n and n.died by dose for each route
intact = reformat.route(intact.raw); intact$route = 'intact'
shaved = reformat.route(shaved.raw); shaved$route = 'shaved'
abraded = reformat.route(abraded.raw); abraded$route = 'abraded'
conjunctival = reformat.route(conjunctival.raw); conjunctival$route = 'conjunctival'


####################################################################
####    Format rat data 
####################################################################
#Input route data
# Copy rat data column "infected" to a column called "died" so you can use the same functions to reformat
rat.data$died = rat.data$infected # Note: in the rat model, infection was not fatal, so this is not strictly correct
## But we need to re-naming the rat$infected column to match the column names in the hamster data
## This allows us to feed the rat data into functions originally designed to process hamster data
abraded.rat.raw = subset(rat.data, route == 'abraded')
IP.rat.raw = subset(rat.data, route == 'IP')

## Get tables of n and infected ("n.died") by dose for each route
abraded.rat = reformat.route(abraded.rat.raw); abraded.rat$route = 'abraded.rat'
IP.rat = reformat.route(IP.rat.raw); IP.rat$route = "IP.rat"



####################################################################
##########################  SECTION 1 ##############################
####   
####   1. Use the exponential model to calculate the probability of infection, given dose
####
####################################################################
####    Estimate pI in hamsters
####################################################################
## This function calculates the and outputs the negative log likelihood of the IP data
nll.pP = function(pP, dose, nn, n.infected){
  # Calculate a vector of probabilities of infection at each dose, given dd, and pP
  p_infected = sapply(X = dose, FUN = function(dd) 1-exp(-pP*dd))
  # Find the log likelihood of each row of the data, given pP
  log.lk.vec = mapply(FUN = dbinom, x = n.infected, size = nn, prob = p_infected, log = TRUE)
  -sum(log.lk.vec) # Return the total negative log likelihood
}

## Optimize the likelihood function
pP.opt = optim(par = c(pP = .18), fn = nll.pP, method = 'Brent', dose = IP$dose, n.infected = IP$died, nn = IP$n, lower = 0, upper = 1)
## View the outputs
pP.opt ## $minimum gives the pP estimate, $value gives the negative log likelihood, $counts gives outputs from the optimization algorithm, $convergence = 0 signifies convergence. See ?optim for greater detail

## RENAME THE MAXIMUM LIKELIHOOD VALUE OF pP FOR USE IN EQUATIONS BELOW
pP = pP.opt$par

#Calculate the likelihood profile
pPs = c(seq(0, .1, by = .01), seq(0.1, .3, by = 0.001), seq(.31, 1, by = 0.01)) # Set a grid of candidate pP values
lk.prof = numeric(0) # Evaluate the likelihood at every grid point
for(ii in 1:length(pPs)){
  lk.prof[ii] = nll.pP(pPs[ii], dose = IP$dose, n.infected = IP$died, nn = IP$n)
}

# Plot the likelihood profile, with the red horizontal like indicating the maximum value included in the profile CI
pdf(plot1, height = 4.5)
par(mfrow = c(1,2))
plot(pPs, lk.prof, main = 'Profile, pI, hamsters', xlab = 'pP', ylab = 'negative log likelihood', type = 'l')
threshold = LR.Threshold(pP.opt$value, df = 1)
abline(h = threshold, col = 'red')
CI.pP = LR.CI(threshold, lk.prof, pPs) ## Calculate the CI!!
text(.5, 100, sprintf('Est = %.2g\n CI = %.2g-%.2g', pP, CI.pP[1], CI.pP[2]))
points(pP.opt$minimum, pP.opt$value, col = 'blue', pch = 8, cex = 1.2)
#points(full.model$minimum[1], pP.opt$value, col = 'yellow', pch = 16, cex = .8)
abline(v = pP.opt$minimum, lty = 2)


# View the CI calculated above
CI.pP


##########################  SECTION 1a ##############################
####   
####   Repeat pI estimation with rat data
####
####################################################################

####################################################################
####    Estimate pP using the IP rat data
## Optimize the likelihood function
pP.opt.rat = optim(par = c(pP = .18), fn = nll.pP, method = 'Brent', dose = IP.rat$dose, n.infected = IP.rat$died, nn = IP.rat$n, lower = 0, upper = 1)
## View the outputs
pP.opt.rat

## SAVE THE MAXIMUM LIKELIHOOD VALUE OF pP FOR USE IN EQUATIONS BELOW
pP.rat = pP.opt.rat$par

#Calculate the likelihood profile
pPs.rat = c(seq(0, .1, by = .01), seq(0.1, .3, by = 0.001), seq(.31, 1, by = 0.01))
lk.prof.rat = numeric(0)
for(ii in 1:length(pPs.rat)){
  lk.prof.rat[ii] = nll.pP(pPs.rat[ii], dose = IP.rat$dose, n.infected = IP.rat$died, nn = IP.rat$n)
}

# Plot the likelihood profile, with the red horizontal like indicating the maximum value included in the profile CI
plot(pPs.rat, lk.prof.rat, main = 'Profile, pI, rats', xlab = 'pP', ylab = 'negative log likelihood', type = 'l')
threshold = LR.Threshold(pP.opt.rat$value, df = 1)
abline(h = threshold, col = 'red')
CI.pP.rat = LR.CI(threshold, lk.prof.rat, pPs.rat)
text(.5, 25, sprintf('Est = %.2g\n CI = %.2g-%.2g', pP, CI.pP[1], CI.pP[2]))
points(pP.opt.rat$par, pP.opt.rat$value, col = 'blue', pch = 8, cex = 1.2)
#points(full.model$par[1], pP.opt$value, col = 'yellow', pch = 16, cex = .8)
abline(v = pP.opt.rat$par, lty = 2)
dev.off()

# Output the CI calculated above
CI.pP.rat






####################################################################
###########################   SECTION 2   ##########################
####   
####   2. Use all non-IP data to find the MLE of p_c(route) for each
####      tested route of inoculation. Fit using the BASIC model.
####      
####################################################################
####################################################################

####################################################################
####    Write functions to calculate the likelihood
####################################################################

#### THIS FUNCTION CALCULATES THE NEG LOG LIKELIHOOD OF A FULL NON-IP DATA SET, using the Basic model
nll.pc.pp = function(pc, doses, nn, n.inf, pp){
  ## Calculate the probability of infection at each dose, given pp and pc
  p_inf = -expm1(-doses*pp*pc)
  ## Find the log likelihood of the data at each dose
  log.lk = mapply(FUN = dbinom, x = n.inf, size = nn, prob = p_inf, log = TRUE)
 -sum(log.lk)
}


####################################################################
####   Maximize the likelihood to estimate pc for each site of inoculation
####################################################################
opt.pc.abraded.skin = optimize(f = nll.pc.pp, doses = abraded$dose, nn = abraded$n, n.inf = abraded$died, pp = pP, interval = c(1e-20, .1), tol = 1e-9)
opt.pc.abraded.skin
## Note, we really want to optimize over the interval [0, 1], but a value of exactly 0, or a value that's too high within the range give non-finite values. 
## Instead, set the min to 1e-20 (very close to 0, but not exactly), and the max to 0.1 (low enough to return a finite value)

opt.pc.intact.skin = optimize(f = nll.pc.pp, doses = intact$dose, nn = intact$n, n.inf = intact$died, pp = pP, interval = c(1e-20, 1e-6), tol = 1e-12)
opt.pc.intact.skin


opt.pc.shaved.skin = optimize(f = nll.pc.pp, doses = shaved$dose, nn = shaved$n, n.inf = shaved$died, pp = pP, interval = c(1e-20, 1e-6), tol = 1e-10)
opt.pc.shaved.skin

opt.pc.conjunctival.skin = optimize(f = nll.pc.pp, doses = conjunctival$dose, nn = conjunctival$n, n.inf = conjunctival$died, pp = pP, interval = c(1e-20, 1e-5), tol = 1e-12)
opt.pc.conjunctival.skin



################### Repeat for rat data, abraded ##############################
opt.pc.abraded.rat.skin = optimize(f = nll.pc.pp, doses = abraded.rat$dose, nn = abraded.rat$n, n.inf = abraded.rat$died, pp = pP.rat, interval = c(1e-20, .1), tol = 1e-12)
opt.pc.abraded.rat.skin




####################################################################
####   Calculate 1D likelihood profiles, hamsters
####################################################################
## Set up plots
pdf(plot2)
par(mfrow = c(2,3))

## Intact skin profile
# Define a range of pc values to scan across
pc.is = 10^seq(-25, -3, by = .01)
## Find the nll value at each point in the pc grid
prof.intact = sapply(X = pc.is, FUN = function(pc.in){nll.pc.pp(pc.in, doses = intact$dose, nn = intact$n, n.inf = intact$died, pp = pP)})
## Plot the profile values
plot(log10(pc.is), prof.intact, main = 'pc, hamsters, intact skin', xlab = 'log10(pc)', ylab = 'Neg log Likelihood', type = 'l')
## Plot the best value in blue
points(log10(opt.pc.intact.skin$minimum), opt.pc.intact.skin$objective, pch = 8, col = 'blue', cex = 1.2)
## Calculate the threshold nll value that defines the confidence interval
# Calculate profile CIs using method based on the likelihood ratio of profile likelihood and best likelihood
# See Bolker, Ecological Models and Data in R, 2008 Ch. 6
# pdf at https://ms.mcmaster.ca/~bolker/emdbook/book.pdf
threshold.intact = LR.Threshold(opt.pc.intact.skin$objective, 1)
CI.intact = LR.CI(threshold.intact, prof.intact, pc.is)
abline(v = log10(opt.pc.intact.skin$minimum), lty = 2)
abline(h = threshold.intact, col = 'red')

## Repeate for abraded
pc.as = 10^seq(-4, -.5, by = .01)
## Find the nll value at each point in the pc grid
prof.abraded = sapply(X = pc.as, FUN = function(pc.in){nll.pc.pp(pc.in, doses = abraded$dose, nn = abraded$n, n.inf = abraded$died, pp = pP)})
## Plot the profile values
plot(log10(pc.as), prof.abraded, main = 'pc, hamsters, abraded skin', xlab = 'log10(pc)', ylab = 'Neg log Likelihood', type = 'l')
points(log10(opt.pc.abraded.skin$minimum), opt.pc.abraded.skin$objective, pch = 8, col = 'blue', cex = 1.2)
threshold.abraded = LR.Threshold(opt.pc.abraded.skin$objective, 1)
CI.abraded = LR.CI(threshold.abraded, prof.abraded, pc.as)
abline(v = log10(opt.pc.abraded.skin$minimum), lty = 2)
abline(h = threshold.abraded, col = 'red')


## Repeate for shaved
pc.ss = 10^seq(-20, -5, by = .01)
## Find the nll value at each point in the pc grid
prof.shaved = sapply(X = pc.ss, FUN = function(pc.in){nll.pc.pp(pc.in, doses = shaved$dose, nn = shaved$n, n.inf = shaved$died, pp = pP)})
## Plot the profile values
plot(log10(pc.ss), prof.shaved, main = 'pc, hamsters, shaved skin', xlab = 'log10(pc)', ylab = 'Neg log Likelihood', type = 'l')
points(log10(opt.pc.shaved.skin$minimum), opt.pc.shaved.skin$objective, pch = 8, col = 'blue', cex = 1.2)
threshold.shaved = LR.Threshold(opt.pc.shaved.skin$objective, 1)
CI.shaved = LR.CI(threshold.shaved, prof.shaved, pc.ss)
abline(v = log10(opt.pc.shaved.skin$minimum), lty = 2)
abline(h = threshold.shaved, col = 'red')


## Repeate for conjunctival
pc.cs = 10^seq(-8, -4, by = .01)
## Find the nll value at each point in the pc grid
prof.conjunctival = sapply(X = pc.cs, FUN = function(pc.in){nll.pc.pp(pc.in, doses = conjunctival$dose, nn = conjunctival$n, n.inf = conjunctival$died, pp = pP)})
## Plot the profile values
plot(log10(pc.cs), prof.conjunctival, main = 'pc, hamsters, conjunctival skin', xlab = 'log10(pc)', ylab = 'Neg log Likelihood', type = 'l')
points(log10(opt.pc.conjunctival.skin$minimum), opt.pc.conjunctival.skin$objective, pch = 8, col = 'blue', cex = 1.2)
threshold.conjunctival = LR.Threshold(opt.pc.conjunctival.skin$objective, 1)
CI.conjunctival = LR.CI(threshold.conjunctival, prof.conjunctival, pc.cs)
abline(v = log10(opt.pc.conjunctival.skin$minimum), lty = 2)
abline(h = threshold.conjunctival, col = 'red')






####################################################################
####   Repeat for abraded data, rats
####################################################################
# Repeat for abraded skin, basic model, rats
pc.as.rat = 10^seq(-4, -.5, by = .01)
## Find the nll value at each point in the pc grid
prof.abraded.rat = sapply(X = pc.as.rat, FUN = function(pc.in){nll.pc.pp(pc.in, doses = abraded.rat$dose, nn = abraded.rat$n, n.inf = abraded.rat$died, pp = pP)})
## Plot the profile values
plot(log10(pc.as.rat), prof.abraded.rat, main = 'pc, rat, abraded skin', xlab = 'log10(pc)', ylab = 'Neg log Likelihood', type = 'l')
points(log10(opt.pc.abraded.rat.skin$minimum), opt.pc.abraded.rat.skin$objective, pch = 8, col = 'blue', cex = 1.2)
threshold.abraded.rat = LR.Threshold(opt.pc.abraded.rat.skin$objective, 1)
CI.abraded.rat = LR.CI(threshold.abraded.rat, prof.abraded.rat, pc.as.rat)
abline(v = log10(opt.pc.abraded.rat.skin$minimum), lty = 2)
abline(h = threshold.abraded.rat, col = 'red')
dev.off()






####################################################################
###########################   SECTION 3   ##########################
####   
####   3. Repeat fits to abraded data using the mixture model.
####      
####################################################################
####################################################################
## Function: Output the probability of infection (scalar)
##           Input: dose (dd.in), scalar
p_inf_fun = function(dd.in, pp.in, aa.in, bb.in, species = 'hamsters'){
  ## Function: input scalar pc and dose (dd) values
  ##  Output the value of the integrand in equation ##
  integrand = function(pc){
    exp(-dd.in*pc*pp.in)*dbeta(x = pc, shape1 = aa.in, shape2 = bb.in)
  }
  ## We need to integrate from 0 to 1
  ## This integral is a bit fussy. Calculating two integrals over the range [0,1] and choosing a sensible break point between the intergrals usually ensures convergence, as long as aa and bb are both > 0.1. Integrating over the full range, [0,1] often fails.
  if(species == 'hamsters'){
    ## Move the break point closer to 0 as the dose gets larger
    ## The integral requires a lot of subdivisions to converge near pc=0, and convergence is particularly fussy at high doses when using parameters fitted to hamster data
    ## Moving the break point close to 0 allows the integrator to place more subdivisions in the problematic range, [0, break.pt]
    ## This ultimately prevents the error: "integral is probably divergent"
  break.pt = min(10^-(log10(dd.in)/1.75), .95)
  }else if(species == 'rats'){
    ## Using the best-fit parameters for rats, almost all the density is around pc = 0.02, so set the break point a tiny bit higher than that
    break.pt = 0.025
  }else{
    error('species must be hamsters or rats')
  }
  int_low = integrate(f = integrand, lower = 0, upper = break.pt, rel.tol = 1e-10, stop.on.error = FALSE)
  int_high = integrate(f = integrand, lower = break.pt, upper = 1, rel.tol = 1e-10, stop.on.error = FALSE)
  1-(int_low$value+int_high$value)
}






#### THIS FUNCTION CALCULATES THE NEG LOG LIKELIHOOD OF A FULL NON-IP DATA SET
####    *Follow equation 0.12
####    *Input a named vector of parameters: pars = c('aa' = ??, 'bb' = ??)
####    *Input estimated pp value (scalar)
####    *Input dose - vector of tested doses
####    *Input nn, a vector of the number of experiments at each dose
####    *Input n.inf, a vector of the number infected at each dose
####    *Output the neg log transformed likelihood
nll.aa.BB = function(pars, pp, dose, nn, n.inf, sp = 'hamsters'){
  aa = pars['aa']
  bb = pars['bb']
  
  ## Calculate the probability of infection at each dose
  p_inf = sapply(X = dose, FUN = p_inf_fun, pp.in = pp, aa.in = aa, bb.in = bb, species = sp)
  
  ## Calculate the log likelihood of observing the data at each dose
  log.lk.vec = mapply(FUN = dbinom, x = n.inf, size = nn, prob = p_inf, log = TRUE)
  #print(c(aa, bb, -sum(log.lk.vec)))
  -sum(log.lk.vec) # Output total negative log likelihood
}
## Test function
nll.aa.BB(pars = c(aa = .1, bb = .35), pp = pP, dose = abraded$dose, nn = abraded$n, n.inf = abraded$died)
nll.aa.BB(pars = c(aa = .1, bb = .35), pp = pP, dose = abraded$dose, nn = abraded$n, n.inf = abraded$died, sp = 'rats')

####################################################################
####   Maximize the likelihood to estimate pc for abraded data, hamsters
####################################################################
## Fit to abraded data:
opt.aa.BB.abraded.skin = optim(c(aa = .2, bb = .5), fn = nll.aa.BB, pp = pP, dose = abraded$dose, nn = abraded$n, n.inf = abraded$died)
opt.aa.BB.abraded.skin


xx = seq(0, 1, by = 1e-4)
# Plot the estimated pc distribution using the ML values of alpha and beta
plot(xx, dbeta(x = xx, shape1 = opt.aa.BB.abraded.skin$par[1], shape2 = opt.aa.BB.abraded.skin$par[2]), main = 'pc distribution from mixture model', xlab = 'pc', ylab = 'probability density')


## Find likelihood profiles for alpha and beta
avals = c(seq(.001, .0249, by = .0001), seq(.025, .029, by = .001), seq(.03, 1.75, by = .01))
bvals = c(seq(.001, .0249, by = .0001), seq(.025, .029, by = .001), seq(.03, 8.5, by = .01))
aas = rep(avals, each = length(bvals))
BBs = rep(bvals, length(avals))

# This takes a while to run. Could be parallellized, but I'll just save and re-load later
prof2d = aas*0
for(ii in 1:length(aas)){
  prof2d[ii]= tryCatch({
    nll.aa.BB(pars = c(aa = aas[ii], bb = BBs[ii]), pp = pP, dose = abraded$dose, nn = abraded$n, n.inf = abraded$died)
  }, warning = function(w) {
  }, error = function(e) {
    return(NA)
  })
}
# remove na values
prof2d = data.frame(aa = aas, bb = BBs, prof = prof2d)
save(prof2d, file = mixture_prof_output)
load(mixture_prof_output)



## Plot 2d likelihood profile (surface), all values
ggplot()+
  geom_contour(data = prof2d, mapping = aes(x = aa, y = bb, z = prof, color = stat(level)), binwidth = 2)+
  #scale_fill_gradientn(breaks = 1-c(0, .028, 1), colors = c('red', 'yellow', 'blue'))+  #
  scale_color_distiller(palette = 'RdYlBu')+
  geom_point(aes(x = opt.aa.BB.abraded.skin$par[1], y = opt.aa.BB.abraded.skin$par[2]))
## Repeat, plotting only values within 4 nll units of the best value
ggplot()+
  geom_contour(data = subset(prof2d, prof2d$prof < min(prof2d$prof, na.rm = TRUE)+5), 
               mapping = aes(x = aa, y = bb, z = prof, color = stat(level)), binwidth = 1)+
  #scale_fill_gradientn(breaks = 1-c(0, .028, 1), colors = c('red', 'yellow', 'blue'))+  #
  scale_color_distiller(palette = 'RdYlBu')+
  geom_point(aes(x = opt.aa.BB.abraded.skin$par[1], y = opt.aa.BB.abraded.skin$par[2]))








### Define confidence interval
## Find the LR threshold (the worst likelihood contained in the 2D confidence envelope
threshold = LR.Threshold(NLL_best = opt.aa.BB.abraded.skin$value, df = 2)
## Find the indices of profile values below the threshold
accept = which(prof2d$prof <= threshold)
prof.acc = prof2d[accept, ] # Store accepted profile likelihood values (values within the envelope)
# Find the min and max aa values from in the envelope. These are the marginal CI values.
aa.max.ind = which.max(prof.acc$aa)
aa.min.ind = which.min(prof.acc$aa)
prof.acc[aa.max.ind, ]; threshold ## Check that the profile values are very close to the threshold. If not, you may need to decrease grid size
prof.acc[aa.min.ind, ]; threshold
CI.aa = cbind(c(prof.acc$aa[aa.min.ind], prof.acc$bb[aa.min.ind]), c(prof.acc$aa[aa.max.ind], prof.acc$bb[aa.max.ind]))
colnames(CI.aa) = c('min', 'max')
rownames(CI.aa) = c('aa', 'bb')
## Repeat for min and max bb values. These are the marginal CI values.
bb.max.ind = which.max(prof.acc$bb)
bb.min.ind = which.min(prof.acc$bb)
prof.acc[bb.max.ind,]; threshold
prof.acc[bb.min.ind, ]; threshold
CI.bb = cbind(c(prof.acc$aa[bb.min.ind], prof.acc$bb[bb.min.ind]), c(prof.acc$aa[bb.max.ind], prof.acc$bb[bb.max.ind]))
colnames(CI.bb) = c('min', 'max')
rownames(CI.bb) = c('aa', 'bb')








# ####################################################################
# ####   Repeat for rat data
# ####################################################################
## Fit to abraded data:
opt.aa.BB.abraded.rat.skin = optim(c(aa = 2, bb = 100), fn = nll.aa.BB, pp = pP.rat, dose = abraded.rat$dose, nn = abraded.rat$n, n.inf = abraded.rat$died, sp = 'rats', method = 'L-BFGS-B', lower = .01, upper = 50000)
opt.aa.BB.abraded.rat.skin

# Plot the estimated pc distribution using the ML values of alpha and beta
# Note that all cumulative density is concentrated around the MLE pc estimate from the basic model
# See supplementary material for interpretation
plot(xx, pbeta(xx, shape1 = opt.aa.BB.abraded.rat.skin$par[1], shape2 = opt.aa.BB.abraded.rat.skin$par[2]), main = 'pc distribution from mixture model, rats', xlab = 'pc', ylab = 'probability density')





# ####################################################################
# ####   Calculate 2D likelihood profiles (Mixture model), rats
# ####################################################################
## Find likelihood profiles for alpha and beta
opt.aa.BB.abraded.rat.skin$par
## Set up 2D grid
aas.r = seq(1, 1000, by = 10)
BBs.r = seq(100, 30000, by = 100)
aas.rat = rep(aas.r, each = length(BBs.r))
BBs.rat = rep(BBs.r, times = length(aas.r))

#calculate 2D profile for aa and bb:
prof2d.rat = aas.rat*0
for(ii in 1:length(aas.rat)){
  prof2d.rat[ii] = tryCatch({nll.aa.BB(pars = c(aa = aas.rat[ii], bb = BBs.rat[ii]), pp = pP.rat, dose = abraded.rat$dose, nn = abraded.rat$n, n.inf = abraded.rat$died, sp = 'rats')},
                            warning = function(w) {return(NA)},
                            error = function(e) {return(NA)})
}

save(prof2d.rat, aas.rat, BBs.rat, file = mixture_prof_rats)
load(mixture_prof_rats)

threshold = LR.Threshold(opt.aa.BB.abraded.rat.skin$value, 2) # Define threshold nll for 95% confidence envelope

ggplot()+
  geom_raster(aes(x = aas.rat, y = BBs.rat, fill = prof2d.rat))+
  scale_fill_distiller(palette = 'RdYlBu', limits = min(prof2d.rat, na.rm = TRUE)+c(0, 5))




## Base plot
plotdf = data.frame(alpha = aas.rat, beta = BBs.rat, nll = prof2d.rat) ## 2d profiles
plotdf$line_alphas = opt.pc.abraded.rat.skin$minimum*plotdf$beta/(1-opt.pc.abraded.rat.skin$minimum) ## Add coords for line representing alpha/(alpha+beta) = 0.02
linepts = which(abs(plotdf$alpha/(plotdf$alpha+plotdf$beta) - opt.pc.abraded.rat.skin$minimum) < 5e-5)## extract (alpha, beta) pairs that fall on the line
linepts = linepts[c(1, 15, 34)]
plotdf[linepts, ]
pp = ggplot(subset(plotdf, alpha <=1000))+geom_point(aes(x = alpha, y = beta, color = nll)) +
  scale_color_distiller(palette = "Spectral") +
  geom_contour(aes(x = alpha, y = beta, z = nll), binwidth = 2, color = 'gray20') +
  geom_line(aes(x = line_alphas, y = beta), color = 'white', lwd = 2) +
  geom_point(data = plotdf[linepts[1],], aes(x = alpha , y = beta), shape = 7, size = 2)+
  geom_point(data = plotdf[linepts[2],], aes(x = alpha , y = beta), shape = 8, size = 2)+
  geom_point(data = plotdf[linepts[3],], aes(x = alpha , y = beta), shape = 5, size = 2)+
  geom_text(data = plotdf[linepts[1],], aes(x = alpha +90, y = beta, label = round(nll, 4))) +
  geom_text(data = plotdf[linepts[2],], aes(x = alpha +90, y = beta, label = round(nll, 4))) +
  geom_text(data = plotdf[linepts[3],], aes(x = alpha +90, y = beta, label = round(nll, 4))) +
  ggtitle('2D likelihood surface', subtitle = 'mixture model fitted to data from rats, abraded skin')+
  xlab(expression(alpha)) +
  ylab(expression(beta)) +
  theme_bw()
pp





## Plot beta CDF for the three focal points
xx = seq(0, .05, by = 0.0001)
bdlist = data.frame(xx = xx,
                   pt1 = pbeta(q = xx, shape1 = plotdf[linepts[3],'alpha'], shape2 = plotdf[linepts[3],'beta']),
                   pt2 = pbeta(q = xx, shape1 = plotdf[linepts[2],'alpha'], shape2 = plotdf[linepts[2],'beta']),
                   pt3 = pbeta(q = xx, shape1 = plotdf[linepts[1],'alpha'], shape2 = plotdf[linepts[1],'beta']))
bdlist = melt(bdlist, id.var = 'xx')
bdplot = ggplot(data = bdlist)+
  geom_point(aes(x = xx, y = value, shape = variable, color = variable), size = 2, alpha = .5)+
  geom_line(aes(x = xx, y = value, color = variable), alpha = .7)+
  scale_color_discrete(name = "",
                       breaks=c("pt1", "pt2", "pt3"),
                      labels=c("upper point", "middle point", "lower point")) +
  scale_shape_manual(name = "",
                     values = c(5,8,7),
                     breaks=c("pt1", "pt2", "pt3"),
                     labels=c("upper point", "middle point", "lower point")) +
  xlab('pc')+
  ylab('cumulative density')+
  ggtitle('Fitted Beta(\u03B1,\u03B2) distribution', subtitle = 'selected points in top panel specify (\u03B1, \u03B2) values')+
  theme_bw()
bdplot


outplot = arrangeGrob(pp, bdplot, nrow = 2)
ggsave(filename = figS1, plot = outplot, width = 5, height = 7, dpi = 'print', device = png())


## Calculate delta AIC for fits to abraded skin data
##  Mixture model vs basic model in hamsters
##  Not possible to do this comparison in rats, because there was no real MLE in rat mixture model
AICcalc = function(n.pars, nll){
  2*n.pars+2*nll
}
AIC.basic.ham = AICcalc(n.pars = 2, nll = opt.pc.abraded.skin$objective)
AIC.mixture.ham = AICcalc(n.pars = 3, nll = opt.aa.BB.abraded.skin$value)
del.AIC.hamster = c(AIC.basic.ham, AIC.mixture.ham); del.AIC.hamster = del.AIC.hamster-min(del.AIC.hamster)
del.AIC.hamster



####################################################################
###########################   SECTION 4   ##########################
####   
####   4. Plot results!
####      
####################################################################
####################################################################
## Plot point estimates and CIs of each free paramter, as well as estimated Beta ditribution of pc values
{
pdf(fig2)
par(mgp = c(2, 1, 0), mar = c(3,7,2,2))
layout(mat = matrix(c(1,1,2,3), byrow = T, nrow = 2, ncol = 2))
ham.col = 'dodgerblue3'
rat.col = 'gold'
## Plot1 - point estimates and CIs ----------
xtext = c(-1.9, -8.5, -7.1, -5.4)
ytext = c(1.001, 1.001, 1.001, 1.001)+.0005
#cols = c('firebrick1', 'deepskyblue', 'limegreen', 'chocolate1')
labs = c(expression('p'[c]*' abraded'), 
         expression('p'[c]*' intact'), 
         expression('p'[c]*' shaved'),
         expression('p'[c]*' conjunctival'),
         expression('p'[p]),
         expression('p'[c]*' abraded'),
         expression('p'[p]))
ests = c(opt.pc.abraded.skin$minimum, opt.pc.intact.skin$minimum, opt.pc.shaved.skin$minimum, opt.pc.conjunctival.skin$minimum, pP, opt.pc.abraded.rat.skin$minimum, pP.rat)
CI.lows = c(CI.abraded[1], CI.intact[1], CI.shaved[1], CI.conjunctival[1], CI.pP[1], CI.abraded.rat[1], CI.pP.rat[1])
CI.highs = c(CI.abraded[2], CI.intact[2], CI.shaved[2], CI.conjunctival[2], CI.pP[2], CI.abraded.rat[2], CI.pP.rat[2])
yvals = length(ests):1
plot(log10(ests), yvals, pch = c(rep(21, 5), 22,22), bg = c(rep(ham.col, 5), rat.col, rat.col), col = c(rep('black', 5), 'gray', 'gray'), xlab = 'log10(value)', ylab = '', yaxt = 'n', bty = 'n', xlim = c(-10, 0))
segments(x0 = log10(CI.lows), x1 = log10(CI.highs), y0 = yvals, bty = 'n', col = c(rep(ham.col, 5), rat.col, rat.col))
axis(2, at = yvals, labels = labs, las = 2)
legend(-10, 2.5, c('hamsters', 'rats'), pch = c(21,22), yjust = .5, col = c('black', 'gray'), pt.bg = c(ham.col, rat.col), bty = 'n')
abline(h = 2.5, lty = 2)
mtext(3, text = 'A', at = -12, line = 1, font = 2, xpd = NA)
xlims = par('usr')[1:2]



par(mar = c(3, 4,3,2))
xx = seq(0, 1, by = 0.001)
plot((xx), (dbeta(xx, shape1 = opt.aa.BB.abraded.skin$par['aa'], shape2 = opt.aa.BB.abraded.skin$par['bb'])), lty = 1, pch = 21, type = 'b', lwd = .5, ylab = 'density', xlab = expression('p'[c]), bg = ham.col)
text(x = -.3, 36, 'B', font = 2, xpd = NA)

plot((xx), (dbeta(xx, shape1 = opt.aa.BB.abraded.rat.skin$par['aa'], shape2 = opt.aa.BB.abraded.rat.skin$par['bb'])), lty = 1, bg = rat.col, pch = 22, type = 'b', lwd = .5, ylab = 'density', xlab = expression('p'[c]), col = 'gray')
text(x = -.3, 655, 'C', font = 2, xpd = NA)
dev.off()
}




########################## Plot data vs. model fits ############################
##### Fig. 3:

## Reformat observed data for plotting
dat = reformat.route(intact)
intact.CIs = binom.confint(dat[,'died'], dat[,'n'], conf.level = .95, methods = 'wilson')[,c('lower', 'upper')]
dat = reformat.route(shaved)
shaved.CIs = binom.confint(dat[,'died'], dat[,'n'], conf.level = .95, methods = 'wilson')[,c('lower', 'upper')]
dat = reformat.route(abraded)
abraded.CIs = binom.confint(dat[,'died'], dat[,'n'], conf.level = .95, methods = 'wilson')[,c('lower', 'upper')]
dat = reformat.route(conjunctival)
conjunctival.CIs = binom.confint(dat[,'died'], dat[,'n'], conf.level = .95, methods = 'wilson')[,c('lower', 'upper')]
IP.CIs = binom.confint(IP[,'died'], IP[,'n'], conf.level = .95, methods = 'wilson')[,c('lower', 'upper')]

ratdat = reformat.route(abraded.rat)
rat.CIs = binom.confint(ratdat[,'died'], ratdat[,'n'], conf.level = .95, methods = 'wilson')[,c('lower', 'upper')]
IP.CIs.rat = binom.confint(IP.rat[,'died'], IP.rat[,'n'], conf.level = .95, methods = 'wilson')[,c('lower', 'upper')]


#Define a function to plot probability of infection at different D0s, given parameter estimates
p.inf = function(doses, best.pc, best.pp = pP){
  1-exp(-doses*best.pc*best.pp)
}

## Repeat for prob infection, given  IP inoculation
p.inf.IP = function(dose, best.pp){
  1-exp(-dose*best.pp)
}


## Repeat for prob infection, given mixture model
p.inf.mixture = function(dd.in, pp.in, aa.in, bb.in, sp = 'hamsters'){
  sapply(dd.in, function(DD) {p_inf_fun(DD, pp.in, aa.in, bb.in, sp)})
} 

# Convert color to color + transparency
col2alpha <- function(col, alpha) {
  col_rgb <- col2rgb(col)/255
  rgb(col_rgb[1], col_rgb[2], col_rgb[3], alpha = alpha)
}



{
pdf(fig3, height = 5, width = 7)
  par(mfrow = c(2,3), mgp = c(2,1,0), mar = c(3,3,3,1))

  dd = round(10^seq(0, 9, by = .01))
  plot(dd, p.inf(dd, opt.pc.intact.skin$minimum), ylab = 'Probability of Infection', xlab = 'Dose', type = 'l', lwd = 2, ylim = c(0,1), log = 'x', cex.lab = 1.1, xlim = c(1, 1e8))
  #Plot CI bars
  segments(x0 = intact[,'dose'], y0 = intact.CIs[,1], y1 = intact.CIs[,2])
  segments(x0 = intact[,'dose']*.7, y0 = intact.CIs[,1], x1 = intact[,'dose']*1.3)
  segments(x0 = intact[,'dose']*.7, y0 = intact.CIs[,2], x1 = intact[,'dose']*1.3)
  # data
  points(intact[,'dose'], intact[,'died']/intact[,'n'], pch = 21, cex = 1.5, bg = ham.col)
  text(.1, 1.15, 'A', font = 2, xpd = NA)
  
  
  
  #Shaved
  plot(dd, p.inf(dd, opt.pc.shaved.skin$minimum, pP), ylab = 'Probability of Infection', xlab = 'Dose', type = 'l', lwd = 2, ylim = c(0,1), log = 'x', cex.lab = 1.1, xlim = c(1, 1e8))
  #Plot CI bars
  segments(x0 = shaved[,'dose'], y0 = shaved.CIs[,1], y1 = shaved.CIs[,2])
  segments(x0 = shaved[,'dose']*.7, y0 = shaved.CIs[,1], x1 = shaved[,'dose']*1.3)
  segments(x0 = shaved[,'dose']*.7, y0 = shaved.CIs[,2], x1 = shaved[,'dose']*1.3)
  # data
  points(shaved[,'dose'], shaved[,'died']/shaved[,'n'], pch = 21, cex = 1.5, bg = ham.col)
  text(.1, 1.15, 'B', font = 2, xpd = NA)
  
  #Abraded
  ## Basic model fit and data
  plot(dd, p.inf(dd, opt.pc.abraded.skin$minimum, pP), ylab = 'Probability of Infection', xlab = 'Dose', type = 'l', lwd = 2, ylim = c(0,1), log = 'x',  cex.lab = 1.1, xlim = c(1, 1e8))
  yy.mixture.plot = p.inf.mixture(dd.in = dd, aa.in = opt.aa.BB.abraded.skin$par['aa'], bb.in = opt.aa.BB.abraded.skin$par['bb'], pp.in = pP)
  lines(dd, yy.mixture.plot, lty = 3, lwd = 2)
  #Plot data CI bars
  segments(x0 = abraded[,'dose'], y0 = abraded.CIs[,1], y1 = abraded.CIs[,2])
  segments(x0 = abraded[,'dose']*.7, y0 = abraded.CIs[,1], x1 = abraded[,'dose']*1.3)
  segments(x0 = abraded[,'dose']*.7, y0 = abraded.CIs[,2], x1 = abraded[,'dose']*1.3)
  # data
  points(abraded[,'dose'], abraded[,'died']/abraded[,'n'],  pch = 21, cex = 1.5, bg = ham.col)
  legend('bottomright', legend = c('basic model', 'mixture model'), lty = c(1, 3), bty = 'n', lwd = 2)
  text(.1, 1.15, 'C', font = 2, xpd = NA)

  
  
  #Conjunctval
  plot(dd, p.inf(dd, opt.pc.conjunctival.skin$minimum), ylab = 'Probability of Infection', xlab = 'Dose', type = 'l',  lwd = 2, ylim = c(0,1), log = 'x',  cex.lab = 1.1, xlim = c(1, 1e8))
  #Plot CI bars
  segments(x0 = conjunctival[,'dose'], y0 = conjunctival.CIs[,1], y1 = conjunctival.CIs[,2])
  segments(x0 = conjunctival[,'dose']*.7, y0 = conjunctival.CIs[,1], x1 = conjunctival[,'dose']*1.3)
  segments(x0 = conjunctival[,'dose']*.7, y0 = conjunctival.CIs[,2], x1 = conjunctival[,'dose']*1.3)
  # data
  points(conjunctival[,'dose'], conjunctival[,'died']/conjunctival[,'n'], pch = 21, cex = 1.5, bg = ham.col)
  text(.1, 1.15, 'D', font = 2, xpd = NA)
  
  
  
  
  
  
  
  #Abraded, rats
  plot(dd, p.inf(dd, best.pc = opt.pc.abraded.rat.skin$minimum, best.pp = pP), ylab = 'Probability of Infection', xlab = 'Dose', type = 'l', lwd = 2, ylim = c(0,1), log = 'x', cex.lab = 1.1, xlim = c(1, 1e8), col = 'gray55')

  # Mixture model fits
  yy.mixture.plot.rat = p.inf.mixture(dd.in = dd, aa.in = opt.aa.BB.abraded.rat.skin$par['aa'], bb.in = opt.aa.BB.abraded.rat.skin$par['bb'], pp.in =pP.rat, sp = 'rats')
  lines(dd, yy.mixture.plot.rat, lty = 3, lwd = 2, col = 'gray55')
  
  #Plot data CI bars
  segments(x0 = abraded.rat[,'dose'], y0 = rat.CIs[,1], y1 = rat.CIs[,2], col = 'gray55')
  segments(x0 = abraded.rat[,'dose']*.7, y0 = rat.CIs[,1], x1 = abraded.rat[,'dose']*1.3, col = 'gray55')
  segments(x0 = abraded.rat[,'dose']*.7, y0 = rat.CIs[,2], x1 = abraded.rat[,'dose']*1.3, col = 'gray55')
  
  
  # Data
  points(abraded.rat[,'dose'], abraded.rat[,'died']/abraded.rat[,'n'],  pch = 22, cex = 1.5, bg = rat.col, col = 'gray')
  
  legend('bottomright', legend = c('basic model', 'mixture model'), lty = c(1, 3), bty = 'n', lwd = c(1,2), col = 'gray55')
  text(.1, 1.15, 'E', font = 2, xpd = NA)
  
  
  
  
  
  #IP, hamsters
  plot(dd, p.inf.IP(dose = dd, best.pp = pP), ylab = 'Probability of Infection', xlab = 'Dose', type = 'l', ylim = c(0,1), log = 'x', cex.lab = 1.1, xlim = c(1, 1e8), lwd = 2)
  segments(x0 = IP[,'dose'], y0 = IP.CIs[,1], y1 = IP.CIs[,2], lwd = 1)
  segments(x0 = IP[,'dose']*.7, y0 = IP.CIs[,1], x1 = IP[,'dose']*1.3, lwd = 1)
  segments(x0 = IP[,'dose']*.7, y0 = IP.CIs[,2], x1 = IP[,'dose']*1.3, lwd = 1)
  points(IP[,'dose'], IP[,'died']/IP[,'n'],  pch = 21, cex = 1.7, bg = ham.col)
  

  lines(dd, p.inf.IP(dd, pP.rat), lty = 1, lwd = 1.3, col = 'gray55')
  xx.jitter = 10^(log10(IP.rat[,'dose'])+.1) # Jitter rat CIs so you can see CIs from both species
  segments(x0 = xx.jitter, y0 = IP.CIs.rat[,1], y1 = IP.CIs.rat[,2], col = 'gray50', lwd = .7)
  segments(x0 = IP.rat[,'dose'], y0 = IP.CIs.rat[,1], x1 = IP.rat[,'dose']*1.5, col = 'gray50', lwd = .7)
  segments(x0 = IP.rat[,'dose'], y0 = IP.CIs.rat[,2], x1 = IP.rat[,'dose']*1.5, col = 'gray50', lwd = .7)
  points(xx.jitter, IP.rat[,'died']/IP.rat[,'n'],  pch = 22, cex = 1.2, bg = rat.col, col = 'gray')
  
  legend('bottomright', legend = c('hamsters', 'rats'), lty = 1, bty = 'n', lwd = c(2,2), pch = c(21, 22), pt.bg = c(ham.col, rat.col), col = c('black', 'gray50'))
  text(.1, 1.15, 'F', font = 2, xpd = NA)
  dev.off()
}




















{
pdf(fig4, height = 7)
  # Define a function to make transparent colors
  tns.col = function(col.in, tns.pct){
    rgbvals = col2rgb(col.in)
    tnscol = rgb(rgbvals[1], rgbvals[2], rgbvals[3], tns.pct/100*255,maxColorValue = 255)
    tnscol
  }
  
par(mfrow = c(2,1), mgp = c(2,1,0), mar = c(3,3,3,3))
dd = c(1:8, 10^seq(1.7, 11, by = .15))
p.inf.i =  p.inf(dd, opt.pc.intact.skin$minimum); p.inf.i[which(!is.finite(p.inf.i))] = 1
plot(dd, p.inf.i, ylab = 'P(Infection), hamsters', xlab = 'Dose', type = 'l', col = cols[1], ylim = c(0,1.3), log = 'x', cex.lab = 1.1, lwd = 3, bty= 'n', yaxt = 'n', xaxt = 'n', lty = 3)
axis(1, at = 10^(1:9), labels = 10^(1:9))

lower = p.inf(dd, CI.intact[1]); lower[which(!is.finite(lower))] = 1
upper = p.inf(dd, CI.intact[2]); upper[which(!is.finite(upper))] = 1
polygon(x = c(dd, rev(dd)), y = c(lower, rev(upper)), border = NA, col = col2alpha(cols[1], .5))
axis(side = 2, at = seq(0, 1, by = .1))

## Initialize storage and calculate bounds -- doses at which infection becomes possible, or certain
bounds = matrix(NA, nrow = 2, ncol = 5, dimnames = list(c('min', 'always'), c('intact', 'shaved', 'conj', 'abraded', 'IP')))
bounds[,1] = dd[c(max(which(p.inf.i < 0.1)), min(which(p.inf.i > 0.975)))]

#Shaved
p.inf.s =  p.inf(dd, opt.pc.shaved.skin$minimum); p.inf.s[which(!is.finite(p.inf.s))] = 1
lines(dd, p.inf.s, col = cols[2], lwd = 3, lty = 3)
lower =  p.inf(dd, CI.shaved[1], best.pp = pP); lower[which(!is.finite(lower))] = 1
upper =  p.inf(dd, CI.shaved[2], best.pp = pP); upper[which(!is.finite(upper))] = 1
polygon(x = c(dd, rev(dd)), y = c(lower, rev(upper)), border = NA, col = col2alpha(cols[2], .5))
bounds[,2] = dd[c(max(which(p.inf.s < 0.1)), min(which(p.inf.s > 0.975)))]


#Abraded
yy.mixture.plot = p.inf.mixture(dd.in = dd, pp.in = pP, aa.in = opt.aa.BB.abraded.skin$par[1], bb.in = opt.aa.BB.abraded.skin$par[2])
lines(dd, yy.mixture.plot, col = cols[3], lty = 4, lwd = 3)
############ Mixture model CIs for plotting ###############
# In the mixture model, the confidence envelope is defined by a vague elliptical curve of (alpha, beta) pairs
#  Pairs on this curve have likelihoods near the threshold of the best nll + a threshold determined by the likelihood ratio method
#  Above, we defined prof.acc as the set of all (alpha, beta) pairs whose neg log lilelihoods were under the threshold value
#  Here, we'll downsample those values to extract pairs whose likelihood is very close to the threshold
#  Randomly sample points on or near this curve from the 2D profile
#  Sample 50 points per unit on the BB axis
#  Only sample points whose likelihoods are close to the threshold. THese should give the worst acceptable fits.
aa.coords.reduced = bb.coords.reduced = prof.vals.reduced = NULL
threshold = LR.Threshold(opt.aa.BB.abraded.skin$value,  2)
envelope = NULL
for(ii in 1:8){
  ## The first condition verifies points near the threshold
  ## The second condition verifies that we're sampling relatively evenly from across the full range of beta values
  valid = which(prof.acc$prof>threshold-.01 & prof.acc$bb > (-1+ii) & prof.acc$bb<= (0+ii))
  smp = sample(x = valid, size = min(25, length(valid)), replace = FALSE)
  envelope = rbind(envelope, prof.acc[smp, ])
}
# plot(envelope$aa, envelope$bb) ## The reduced set of points fall on the boundary of the confidence envelope
# ## Now, calculate model-predicted probabilities of infection for each dose, using each pair of (alpha, beta) values in the envelope
wrapper = function(AA, BB){
    tryCatch({
      out = p.inf.mixture(dd.in = dd, pp.in = pP, aa.in =AA, bb.in = BB) },
      error = function(e){
      print(e)
      out = rep(NA, 71)
      },
      warning = function(w){
        print(w)
      },
      finally = function() {return(out)}
      )
}
yy.CIs = mapply(FUN = wrapper, AA = envelope$aa, BB = envelope$bb)



#### Define CI envelope for plotting as the highest and lowest y value resulting from estimates using all sampled points from the CI envelope
yy.low = apply(X = yy.CIs, MARGIN = 1, FUN = min) # Estimates correspond to doses in dd
yy.high = apply(X = yy.CIs, MARGIN = 1, FUN = max)
polygon(x = c(dd, rev(dd)),
        y = c(yy.low, rev(yy.high)),
        col = tns.col(cols[3], 20), border = NA)

bounds[,4] = dd[c(1, min(which(yy.mixture.plot > 0.975)))]; #c(dd[min(1, which.max(yy.high <0.1))], 1e11)


#Conjunctval
p.inf.c =  p.inf(dd, opt.pc.conjunctival.skin$minimum); p.inf.c[which(!is.finite(p.inf.c))] = 1
lines(dd, p.inf.c, col = cols[4], lwd = 3)
lower =  p.inf(dd, CI.conjunctival[1]); lower[which(!is.finite(lower))] = 1
upper =  p.inf(dd, CI.conjunctival[2]); upper[which(!is.finite(upper))] = 1
polygon(x = c(dd, rev(dd)), y = c(lower, rev(upper)), border = NA, col = col2alpha(cols[4], .5))
bounds[,3] = dd[c(max(which(p.inf.c < 0.1)), min(which(p.inf.c > 0.975)))]


## IP inoculation
dd = c(1:9, 10^seq(.95, 11, by = .15))
p.IP =  p.inf.IP(dose = dd, best.pp = pP); p.IP[which(!is.finite(p.IP))] = 1
lines(dd, p.IP, lwd = 3)
lower =  p.inf.IP(dose = dd, best.pp = CI.pP[1])
upper =  p.inf.IP(dose = dd, best.pp = CI.pP[2])
polygon(x = c(dd, rev(dd)), y = c(lower, rev(upper)), border = NA, col = col2alpha('black', .5))
bounds[,5] = dd[c(1, min(which(p.IP > 0.975)))]

#cols = c('firebrick1', 'deepskyblue', 'chocolate1', 'limegreen', 'darkgreen')
ys = seq(1.03, 1.27, length = 5)+.05
segments(x0 = bounds[1,], x1 = bounds[2,], y0 = ys, col = c(cols[c(1,2,4,3)], 'black'), lty = 2, lwd = 2)
arrows(x0 = bounds[2,], x1 = 1e11, y0 = ys, col = c(cols[c(1,2,4,3)], 'black'), lty = 1, length = .1, lwd = 2)
segments(x0 = bounds[1,], y0 = ys, y1 = ys+.02, col = c(cols[c(1,2,4,3)], 'black'), lwd =2)
segments(x0 = bounds[2,], y0 = ys, y1 = ys+.02, col = c(cols[c(1,2,4,3)], 'black'), lwd = 2)

legend(4e9, y = .6, (c('intact', 'shaved', 'conjunctival', 'abraded', 'IP, no barrier')), col = rev(c('black', cols[c(3,4,2,1)])), fill = rev(c('black', cols[c(3,4,2,1)])), cex = .8, border = NA, lty = rev(c(1, 4, 1, 2, 3)), xpd = NA, lwd = rev(c(1,1,2,2,2,2)), bty = 'n')
legend(1e3, y = 1.67, c('infection possible at this dose', 'infection almost certain at this dose'), lty = c(2, 1), cex = .8, xpd = NA, bty = 'n')
text(.1, 1.3, 'A', font = 2, xpd = NA)



#### Repeat for rats
p.inf.a =  p.inf(dd, opt.pc.abraded.rat.skin$minimum, best.pp = pP.rat); p.inf.a[which(!is.finite(p.inf.a))] = 1
plot(dd, p.inf.a, ylab = 'P(Infection), rats', xlab = 'Dose', type = 'l', col = cols[3], ylim = c(0,1.3), log = 'x', cex.lab = 1.1, lwd = 3, bty= 'n', yaxt = 'n', xaxt = 'n', lty = 4)
axis(1, at = 10^(1:9), labels = 10^(1:9))

lower = p.inf(dd, CI.abraded.rat[1]); lower[which(!is.finite(lower))] = 1
upper = p.inf(dd, CI.abraded.rat[2]); upper[which(!is.finite(upper))] = 1
polygon(x = c(dd, rev(dd)), y = c(lower, rev(upper)), border = NA, col = col2alpha(cols[3], .2))
axis(side = 2, at = seq(0, 1, by = .1))

## Initialize storage and calculate bounds -- doses at which infection becomes possible, or certain
bounds = matrix(NA, nrow = 2, ncol = 2, dimnames = list(c('min', 'always'), c('abraded', 'IP')))
bounds[,1] = dd[c(max(which(p.inf.a < 0.1)), min(which(p.inf.a > 0.975)))]


## IP inoculation
lines(dd, p.IP, lwd = 3)
lower =  p.inf.IP(dose = dd, best.pp = CI.pP.rat[1])
upper =  p.inf.IP(dose = dd, best.pp = CI.pP.rat[2])
polygon(x = c(dd, rev(dd)), y = c(lower, rev(upper)), border = NA, col = col2alpha('black', .5))
bounds[,2] = dd[c(1, min(which(p.IP > 0.975)))]

#cols = c('firebrick1', 'deepskyblue', 'chocolate1', 'limegreen', 'darkgreen')
ys = seq(1.05, 1.25, length = 4)+.05
ys = ys[1:2]
segments(x0 = bounds[1,], x1 = bounds[2,], y0 = ys, col = c(cols[3], 'black'), lty = 2, xpd = NA, lwd = 2)
arrows(x0 = bounds[2,], x1 = 1e11, y0 = ys, col = c(cols[3], 'black'), lty = 1, length = .1, xpd = NA, lwd = 2)
segments(x0 = bounds[1,], y0 = ys, y1 = ys+.02, col = c(cols[3], 'black'),xpd = NA, lwd = 2)
segments(x0 = bounds[2,], y0 = ys, y1 = ys+.02, col = c(cols[3], 'black'),xpd = NA, lwd = 2)

legend(4e9, y = .6, (c('abraded', 'IP, no barrier')), col = rev(c('black', cols[3])), fill = rev(c('black', cols[3])), cex = .8, border = NA, lty = (c(4,1)), xpd = NA, lwd = rev(c(1,3)), bty = 'n')
legend(1e3, y = 1.5, c('infection possible at this dose', 'infection almost certain at this dose'), lty = c(2, 1), cex = .8, xpd = NA, bty = 'n')
text(.1, 1.3, 'B', font = 2, xpd = NA)
dev.off()
}







### Write a table template for experimental results
ham.data = rbind(intact, shaved, abraded, conjunctival, IP)
rat.data = rbind(abraded.rat, IP.rat)

## Write output files
write.csv(ham.data, file = hamdat_summary, row.names = FALSE)
write.csv(rat.data, file = ratdat_summary, row.names = FALSE)








