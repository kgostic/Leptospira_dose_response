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

## Set up a cluster so you can run things in parallel lower in the script. Keep this cluster open until all analyses finish running.
library(parallel)
library(ggplot2)
library(binom)
library(TailRank)
library(RColorBrewer)
library(gridExtra)
library(reshape2)
## Output session info
## writeLines(capture.output(sessionInfo()), "code_outputs/sessionInfo.txt")
## Set up a palette that looks good in color, and also shows variation in grayscale for the printed article
cols = colorRampPalette(brewer.pal(4,"YlGnBu"))(25)
cols[26] = 'midnightblue'
plot(1:26, cex = 10, pch = 16, col = cols) # Look at all colors
plot(1:4, cex = 10, pch = 16, col = cols[c(7, 16, 24, 26) ]) # Choose a subset of four that show good divergence in color and in grayscale (print will be in grayscale)
cols = cols[c(7, 16, 24, 26) ]

####################### INPUTS ####################
full.dat = read.csv('raw_data/Hamsterdata.csv', header = T) # Hamster data
rat.data = read.csv('raw_data/Ratdata.csv', header = T) # Rat data

####################### OUTPUTS ###################
## MAIN TEXT FIGURES
## Note - Fig. 1 was a diagram drawn using MS Power Point
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

# Define functions necessary to calculate 95% CIs from 1D likelihood profiles
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
#--------------------- Hamsters: -----------------------
#For each IP dose, create a table of (dose, n total, n died)
IP.dat = subset(full.dat, route == 'IP')[,-1] ## Get rid of the first column, "route", which is known (all IP)
doses = sort(unique(IP.dat$dose), decreasing = T) #vector of doses tested in IP trials
n.IP = n.died.IP = numeric(0) # Count the number of trials at each dose, and the number of infections at each dose (death == infection)
for(ii in 1:length(doses)){
  dose.dat = subset(IP.dat, dose == doses[ii])
  n.IP[ii] = sum(dose.dat$n)
  n.died.IP[ii] = sum(dose.dat$died)
  rm(dose.dat)
}

# Visualize the summarized data
IP.summary = data.frame(doses = doses, n.IP = n.IP, n.died.IP = n.died.IP); IP.summary

# Create a vector that repeats a 1 for every individual infected hamster and a 0 for every hamster that was not infected. 
# Need this format for model input
outcome.pP = rep(c(1,0), times = c(n.died.IP[1], n.died.IP[1]-n.IP[1])) # Initialize using info from the first row of the dose.dat data
for(ii in 2:length(n.IP)){
  outcome.pP = c(outcome.pP, rep(c(1,0), times = c(n.died.IP[ii], n.IP[ii]-n.died.IP[ii])))
}
outcome.pP
doses.pP = rep(doses, n.IP)

# Visualize the reformatted data
data.frame(outcome.pP = outcome.pP, doses.pP = doses.pP)

### --------------------------------------------------------------
####    Format data for all other routes of inoculation, hamsters
####      Repeat the steps used for IP data
####################################################################
#Input route data
intact = subset(full.dat, route == 'Int.skin')
abraded = subset(full.dat, route == 'Epiderm')
shaved = subset(full.dat, route == 'Shv.skin')
conjunctival = subset(full.dat, route == 'Conj')

# Write a function that creates a summary table with the total number of trials and number of infections (deaths) at a given dose
reformat.route = function(route.dat){
  doses = sort(unique(route.dat$dose))
  tbl = matrix(nrow = length(doses), ncol = 3, dimnames = list(doses, c('n', 'died', 'dose')))
  for(ii in 1:length(doses)){
    tbl[ii, 'n'] = sum(route.dat[route.dat$dose == doses[ii], 'n'])
    tbl[ii, 'died'] = sum(route.dat[route.dat$dose == doses[ii], 'died'])
    tbl[ii, 'dose'] = doses[ii]
  }
  tbl
}


## Get tables of n and n.died by dose for each route
intact.tbl = reformat.route(intact)
shaved.tbl = reformat.route(shaved)
abraded.tbl = reformat.route(abraded)
conjunctival.tbl = reformat.route(conjunctival)

## Create vectors that represent the binary outcome of each infection trial in the data
## Create a vector of the same length representing the corresponding dose used in the trial of interest
# Abraded
Doses.abraded = rep(abraded$dose, times = abraded$n)
Outcomes.abraded = rep(c(1, 0), times = c(abraded$died[1], abraded$n[1]-abraded$died[1]))
for(ii in 2:length(abraded$died)){
  temp = rep(c(1, 0), times = c(abraded$died[ii], abraded$n[ii]-abraded$died[ii]))
  Outcomes.abraded = c(Outcomes.abraded, temp)
}

# Intact
Doses.intact = rep(intact$dose, times = intact$n)
Outcomes.intact = numeric(0)
for(ii in 1:length(intact$died)){
  temp = rep(c(1, 0), times = c(intact$died[ii], intact$n[ii]-intact$died[ii]))
  Outcomes.intact = c(Outcomes.intact, temp)
}

# Shaved
Doses.shaved = rep(shaved$dose, times = shaved$n)
Outcomes.shaved = numeric(0)
for(ii in 1:length(shaved$died)){
  temp = rep(c(1, 0), times = c(shaved$died[ii], shaved$n[ii]-shaved$died[ii]))
  Outcomes.shaved = c(Outcomes.shaved, temp)
}

# Conjunctival
Doses.conjunctival = rep(conjunctival$dose, times = conjunctival$n)
Outcomes.conjunctival = numeric(0)
for(ii in 1:length(conjunctival$died)){
  temp = rep(c(1, 0), times = c(conjunctival$died[ii], conjunctival$n[ii]-conjunctival$died[ii]))
  Outcomes.conjunctival = c(Outcomes.conjunctival, temp)
}



####################################################################
####    Format rat data 
####################################################################
#Input route data
# Copy rat data column "infected" to a column called "died" so you can use the same functions to reformat
rat.data$died = rat.data$infected # Note: in the rat model, infection was not fatal, so this is not strictly correct
## But we need to re-naming the rat$infected column to match the column names in the hamster data
## This allows us to feed the rat data into functions originally designed to process hamster data
abraded.rat = subset(rat.data, route == 'abraded')
IP.rat = subset(rat.data, route == 'IP') #

## Get tables of n and infected ("n.died") by dose for each route
abraded.tbl.rat = reformat.route(abraded.rat)
IP.tbl.rat = reformat.route(IP.rat)

# IP
Doses.IP.rat = rep(IP.rat$dose, times = IP.rat$n)
Outcomes.IP.rat = rep(c(1, 0), times = c(IP.rat$infected[1], IP.rat$n[1]-IP.rat$infected[1]))
for(ii in 2:length(IP.rat$infected)){
  temp = rep(c(1, 0), times = c(IP.rat$infected[ii], IP.rat$n[ii]-IP.rat$infected[ii]))
  Outcomes.IP.rat = c(Outcomes.IP.rat, temp)
}

# Abraded
Doses.abraded.rat = rep(abraded.rat$dose, times = abraded.rat$n)
Outcomes.abraded.rat = numeric(0)
for(ii in 1:length(abraded.rat$infected)){
  temp = rep(c(1, 0), times = c(abraded.rat$infected[ii], abraded.rat$n[ii]-abraded.rat$infected[ii]))
  Outcomes.abraded.rat = c(Outcomes.abraded.rat, temp)
}







####################################################################
##########################  SECTION 1 ##############################
####   
####   1. Use the IP data to find the MLE of pI
####
####################################################################
####    Estimate pI in hamsters
####################################################################
## This function calculates the likelihood of the IP data, using equation 0.11 from the methods
##  The output gives the negative log likelihood
nll.pP = function(pP, outcome, dose){
  dens = function(outcome, dose){ifelse(outcome == 1, 1-(1-pP)^dose, (1-pP)^dose)}
  lk.vec = numeric(0)
  for(ii in 1:length(outcome)){
    lk.vec[ii] = dens(outcome[ii], dose[ii])
  }
  -sum(log(lk.vec))
}

## Optimize the likelihood function
pP.opt = optim(par = c(pP = .18), fn = nll.pP, method = 'Brent', outcome = outcome.pP, dose = doses.pP, lower = 0, upper = 1)
## View the outputs
pP.opt ## $par gives the pP estimate, $value gives the negative log likelihood, $counts gives outputs from the optimization algorithm, $convergence = 0 signifies convergence. See ?optim for greater detail

## RENAME THE MAXIMUM LIKELIHOOD VALUE OF pP FOR USE IN EQUATIONS BELOW
pP = pP.opt$par

#Calculate the likelihood profile
pPs = c(seq(0, .1, by = .01), seq(0.1, .3, by = 0.001), seq(.31, 1, by = 0.01)) # Set a grid of candidate pP values
lk.prof = numeric(0) # Evaluate the likelihood at every grid point
for(ii in 1:length(pPs)){
  lk.prof[ii] = nll.pP(pPs[ii], outcome.pP, doses.pP)
}

# Plot the likelihood profile, with the red horizontal like indicating the maximum value included in the profile CI
pdf(plot1, height = 4.5)
par(mfrow = c(1,2))
plot(pPs, lk.prof, main = 'Profile, pI, hamsters', xlab = 'pP', ylab = 'negative log likelihood', type = 'l')
threshold = LR.Threshold(pP.opt$value, df = 1)
abline(h = threshold, col = 'red')
CI.pP = LR.CI(threshold, lk.prof, pPs) ## Calculate the CI!!
text(.5, 100, sprintf('Est = %.2g\n CI = %.2g-%.2g', pP, CI.pP[1], CI.pP[2]))
points(pP.opt$par, pP.opt$value, col = 'blue', pch = 8, cex = 1.2)
#points(full.model$par[1], pP.opt$value, col = 'yellow', pch = 16, cex = .8)
abline(v = pP.opt$par, lty = 2)


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
pP.opt.rat = optim(par = c(pP = .18), fn = nll.pP, method = 'Brent', outcome = Outcomes.IP.rat, dose = Doses.IP.rat, lower = 0, upper = 1)
## View the outputs
pP.opt.rat

## SAVE THE MAXIMUM LIKELIHOOD VALUE OF pP FOR USE IN EQUATIONS BELOW
pP.rat = pP.opt.rat$par

#Calculate the likelihood profile
pPs.rat = c(seq(0, .1, by = .01), seq(0.1, .3, by = 0.001), seq(.31, 1, by = 0.01))
lk.prof.rat = numeric(0)
for(ii in 1:length(pPs.rat)){
  lk.prof.rat[ii] = nll.pP(pPs.rat[ii], Outcomes.IP.rat, Doses.IP.rat)
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
####   2. Use all non-IP data to find the MLE of p_t(route) for each
####      tested route of inoculation. Fit using the BASIC model.
####      
####################################################################
####################################################################

####################################################################
####    Write functions to calculate the likelihood
####################################################################

#### THIS FUNCTION CALCULATES THE PROBABILITY OF OBSERVING ANY PARTICULAR BINARY INFECTION OUTCOME IN THE BASIC MODEL 
####    *This probability is defined in equation 2 of the methods
####    *These probabilities used to calculate the likelihood of the basic model (equation 3)
####    *Note, to facilitate efficient calculation when D0 is large, we use an approximation (Methods, Supporting information)

# ##  For testing of the approximaiton (below), define the exact equation.
# ##  This runs incredibly slowly for D0>10^7
# Prob.inf.given.dose.exact = function(D0, pc, pP){
#   D0 = round(D0) # Make sure D0 is an integer. The binomial is only defined for integers, and so non-integer inputs may return non-finite probabilities.
#   DW = 1:D0 # Wihtin host dose could hypothetically be anywhere from 1 to D0 and still return a non-0 probability of infection
#   sum( ( 1-((1-pP)^DW) )*(dbinom(x = DW, size = D0, prob = pc))) #pc to be estimated
# }

## Use the fast approximation
Prob.inf.given.dose = function(D0, pc, pP){
  D0 = round(D0) # Make sure you've got an integer D0

  ## Solve for the Dw value above which within-host persistence is guaranteed (Prob failure < 1e-20)
  Dw_threshold = min(ceiling(log(1e-30)/log(1-pP)), D0)

  ## Lower partial sum:
  ## For low wihtin-host doses, (Dw<Dw_threshold), the probability of wihtin-host persistence is less than one
  ##  In these cases, we have to calculate the exact probability
  Dw_below = 0:Dw_threshold
  prob_given_below_threshold = sum( (1-((1-pP)^Dw_below)) * dbinom(x = Dw_below, size = D0, prob = pc))

  ## Upper partial sum:
  ## For Dw>Dw_threshold, the probability of within-host persistence is effectively  1, and the expression reduces from:
  ##  sum[Dw_threshold+1:D0](p(DI>0|Dw)p(Dw|D0)) (a)
  ##      to
  ##  sum[Dw_threshold+1:D0](p(Dw|D0))           (b)
  ##  Equation (b) is equivalent to 1-PBB(Dw_threshold+1, D0, alpha, beta), where PBB is the beta-binomial CMF
  if(Dw_threshold<D0){
    prob_given_above_threshold = 1-pbinom(q = Dw_threshold, size = D0, prob = pc)
  }else{
    prob_given_above_threshold = 0
  }
  ## Return total probability
  prob_given_above_threshold+prob_given_below_threshold
}

# ## Test the approximaiton's accuracy
# D0s = runif(1000, 1, 10^6)
# pcs = runif(1000, 0, 1)
# accuracy = function(D0.in, pc.in){
#   Prob.inf.given.dose(D0 = D0.in, pc = pc.in, pP = pP)-Prob.inf.given.dose.exact(D0 = D0.in, pc = pc.in, pP = pP)
# }
# cl = makeCluster(detectCores()-1)
# clusterExport(cl, varlist = c('D0s', 'pcs', 'accuracy', 'pP', 'Prob.inf.given.dose', 'Prob.inf.given.dose.exact'))
# ers = clusterMap(cl, fun = accuracy, D0.in = D0s, pc.in = pcs, SIMPLIFY = TRUE) # calculate error
# stopCluster(cl)
# plot(1:1000, ers) # Plot
# max(ers) # None greater than 3e-12


#### THIS FUNCTION CALCULATES THE NEG LOG LIKELIHOOD OF A FULL NON-IP DATA SET
####    *Follow equation 0.12
####    *Output the neg log transformed likelihood
nll.pc = function(pars, outcome.vec, D0.vec, pP){
  pc = ifelse(pars[1] >=0, pars[1], 0) #Make sure pc > 0
  pc = ifelse(pc > 1, 1, pc) #If pc > 1, reset to 1
  
  p_i = numeric(length(D0.vec)) #Get the probability of infection (success) given D0
  #note: P(failure, or no infection) = 1-p(success)
  for(ii in 1:length(D0.vec)){
    p_i[ii] = Prob.inf.given.dose(D0.vec[ii], pc, pP)
  }
  #Now sum the log probabilities for each observed outcome: 0 or 1
  aa = ifelse(outcome.vec == 1, p_i, 1-p_i)
  -sum(log(aa))
}


####################################################################
####   Maximize the likelihood to estimate pc for each site of inoculation
####################################################################
opt.pc.abraded.skin = optim(par = c(pc = .00005), fn = nll.pc, pP = pP, outcome.vec = Outcomes.abraded, D0.vec = Doses.abraded, method = 'Brent', lower = 0, upper = .2) # Note, upper limit set to 0.2 because higher values return probabilities of infection of exactly 1 at all tested doses using the approximation, and the likelihood becomes degenerate. See profile below for confirmation that optimization returns the correct answer.
opt.pc.abraded.skin

opt.pc.intact.skin = optim(c(pc = 1e-8), fn = nll.pc, pP = pP,  outcome.vec = Outcomes.intact, D0.vec = Doses.intact, method = 'Brent', lower = 0, upper = .000001)
opt.pc.intact.skin
 
opt.pc.shaved.skin = optim(c(pc = 1e-7), fn = nll.pc, pP = pP, outcome.vec = Outcomes.shaved, D0.vec = Doses.shaved, method = 'Brent', lower = 0, upper = 5e-6)
opt.pc.shaved.skin

opt.pc.conjunctival.skin = optim(c(pc = 1e-7), fn = nll.pc, pP = pP, outcome.vec = Outcomes.conjunctival, D0.vec = Doses.conjunctival,  method = 'Brent', upper = 5e-6, lower = 0)
opt.pc.conjunctival.skin

################### Repeat for rat data, abraded ##############################
opt.pc.abraded.skin.rat = optim(c(pc = .00005), fn = nll.pc, pP = pP.rat, method = 'Brent', upper = 1, lower = 0, outcome.vec = Outcomes.abraded.rat, D0.vec = Doses.abraded.rat)
opt.pc.abraded.skin.rat




####################################################################
####   Calculate 1D likelihood profiles, hamsters
####################################################################
## Intact skin
pc.is = c(seq(1e-9, 2e-7, by = 1e-9), seq(1e-7, 1e-4, by = 1e-4))## Define a sequence of pc values across which to calculate the grid of likelihood values
cl = makeCluster(detectCores()-1)
clusterExport(cl = cl, varlist = c('Outcomes.intact', 'Doses.intact', 'pP', 'nll.pc', 'pc.is', 'Prob.inf.given.dose'))
wrapper = function(pc.in){nll.pc(pars = c(pc = pc.in), outcome.vec = Outcomes.intact, D0.vec = Doses.intact, pP = pP)}
prof.intact = parSapply(cl = cl, X = pc.is, FUN = wrapper)
stopCluster(cl)
# plot(pc.is, prof.intact) # Preview likelihood profile
## Add results to scratch .pdf of profiles
neighborhood = max(0, which.min(prof.intact)-100):min(length(pc.is), which.min(prof.intact)+150) # Get values within range of the best value
pdf(plot2)
par(mfrow = c(2,3))
plot(pc.is[neighborhood], prof.intact[neighborhood], main = 'Profile, pc(intact) hamsters', xlab = 'pc', ylab = 'Neg log Likelihood', type = 'l')
threshold.intact = LR.Threshold(opt.pc.intact.skin$value, 1)
CI.intact = LR.CI(threshold.intact, prof.intact, pc.is)
# Calculate profile CIs using method based on the likelihood ratio of profile likelihood and best likelihood
# See Bolker, Ecological Models and Data in R, 2008 Ch. 6
# pdf at https://ms.mcmaster.ca/~bolker/emdbook/book.pdf
# extract the best pc value from MLE
pc.opt.intact = pc.is[which.min(prof.intact)]; pc.opt.intact
#text(6e-8, 15, sprintf('Best est = %.3g\nCI = %.3g - %.3g', pc.opt.intact, CI.intact[1], CI.intact[2]))
points(opt.pc.intact.skin$par, opt.pc.intact.skin$value, pch = 8, col = 'blue', cex = 1.2)
#points(full.model$par[3], opt.pc.intact.skin$value, pch = 16, col = 'green3', cex = .8)
abline(v = opt.pc.intact.skin$par, lty = 2)
abline(h = threshold.intact, col = 'red')





# Repeat for abraded skin, basic model
pc.abraded = seq(0, .2, by = 0.0001)
cl = makeCluster(detectCores()-1)
clusterExport(cl = cl, varlist = c('Outcomes.abraded', 'Doses.abraded', 'pP', 'nll.pc', 'pc.abraded', 'Prob.inf.given.dose'))
wrapper = function(pc.in){nll.pc(pars = c(pc = pc.in), outcome.vec = Outcomes.abraded, D0.vec = Doses.abraded, pP = pP)}
prof.abraded = parSapply(cl = cl, X = pc.abraded, FUN = wrapper)
stopCluster(cl)
# plot(pc.abraded, prof.abraded) # Preview likelihood profile
## Plot
neighborhood = max(0, which.min(prof.abraded)-1000):min(length(pc.abraded), which.min(prof.abraded)+1000)
pc.opt.abraded = opt.pc.abraded.skin$par
plot(pc.abraded[neighborhood], prof.abraded[neighborhood], main = 'Profile, pc(abraded) hamsters', xlab = 'pc', ylab = 'Neg log Likelihood', type = 'l')
threshold.abraded = LR.Threshold(opt.pc.abraded.skin$value,1)
abline(h = threshold.abraded, col = 'red')
CI.abraded = LR.CI(threshold.abraded, prof.abraded, pc.abraded)
#text(.055, 80, sprintf('Best est = %.2g\nCI = %.2g - %.2g', pc.opt.abraded, CI.abraded[1], CI.abraded[2]))
points(opt.pc.abraded.skin$par, opt.pc.abraded.skin$value, pch = 8, col = 'blue', cex = 1.2)
#points(full.model$par[2], opt.pc.abraded.skin$value, pch = 16, col = 'green3', cex = .8)
abline(v = opt.pc.abraded.skin$par, lty = 2)



# Shaved
pc.shaved = c(seq(0, 10^-6, by = 10^-10), seq(10^-6, 1, .001))
cl = makeCluster(detectCores()-1)
clusterExport(cl = cl, varlist = c('Outcomes.shaved', 'Doses.shaved', 'pP', 'nll.pc', 'pc.shaved', 'Prob.inf.given.dose'))
wrapper = function(pc.in){nll.pc(pars = c(pc = pc.in), outcome.vec = Outcomes.shaved, D0.vec = Doses.shaved, pP = pP)}
prof.shaved = parSapply(cl = cl, X = pc.shaved, FUN = wrapper)
stopCluster(cl)
## Plot
neighborhood = max(0, which.min(prof.shaved)-2500):min(length(pc.shaved), which.min(prof.shaved)+7000)
threshold = LR.Threshold(opt.pc.shaved.skin$value,1)
CI.shaved = LR.CI(threshold, prof.shaved, pc.shaved)
pc.opt.shaved = opt.pc.shaved.skin$par
plot(pc.shaved[neighborhood], prof.shaved[neighborhood], main = 'Profile, pc(shaved), hamsters', xlab = 'pc', ylab = 'Neg log Likelihood', type = 'l')
abline(h = threshold, col = 'red')
CI.shaved = LR.CI(threshold, prof.shaved, pc.shaved)
text(3.5e-7, 60, sprintf('Best est = %.3g\nCI = %.3g - %.3g', pc.opt.shaved, CI.shaved[1], CI.shaved[2]))
points(opt.pc.shaved.skin$par, opt.pc.shaved.skin$value, pch = 8, col = 'blue')
abline(v = opt.pc.shaved.skin$par, lty = 2)




# Conjunctival
pc.con = c(seq(1e-8, 1e-5, by = 1e-8))
cl = makeCluster(detectCores()-1)
clusterExport(cl = cl, varlist = c('Outcomes.conjunctival', 'Doses.conjunctival', 'pP', 'nll.pc', 'pc.con', 'Prob.inf.given.dose'))
wrapper = function(pc.in){nll.pc(pars = c(pc = pc.in), outcome.vec = Outcomes.conjunctival, D0.vec = Doses.conjunctival, pP = pP)}
prof.conjunctival = parSapply(cl = cl, X = pc.con, FUN = wrapper)
stopCluster(cl)
#plot(pc.con, prof.conjunctival)
## Plot
pc.opt.conjunctival = opt.pc.conjunctival.skin$par
neighborhood = max(0, which.min(prof.conjunctival)-70):min(length(pc.con), which.min(prof.conjunctival)+170)
threshold = LR.Threshold(opt.pc.conjunctival.skin$value,1)
CI.conjunctival = LR.CI(threshold, prof.conjunctival, pc.con)
plot(pc.con[neighborhood], prof.conjunctival[neighborhood], main = 'Profile, pc(conjunctival), hamsters', xlab = 'pc', ylab = 'Neg log Likelihood', type = 'l')
abline(h = threshold, col = 'red')
CI.conjunctival = LR.CI(threshold, prof.conjunctival, pc.con)
text(4e-6, 100, sprintf('Best est = %.3g\nCI = %.3g - %.3g', pc.opt.conjunctival, CI.conjunctival[1], CI.conjunctival[2]))
points(opt.pc.conjunctival.skin$par, opt.pc.conjunctival.skin$value, pch = 8, col = 'blue')
abline(v = opt.pc.conjunctival.skin$par, lty = 2)




####################################################################
####   Repeat for abraded data, rats
####################################################################
# Repeat for abraded skin, basic model, rats
pc.abraded.rat = seq(0, 1, by = 0.0001)
cl = makeCluster(detectCores()-1)
clusterExport(cl = cl, varlist = c('Outcomes.abraded.rat', 'Doses.abraded.rat', 'pP.rat', 'nll.pc', 'pc.abraded.rat', 'Prob.inf.given.dose'))
wrapper = function(pc.in){nll.pc(pars = c(pc = pc.in), outcome.vec = Outcomes.abraded.rat, D0.vec = Doses.abraded.rat, pP = pP.rat)}
prof.abraded.rat = parSapply(cl = cl, X = pc.abraded.rat, FUN = wrapper)
stopCluster(cl)
#plot(pc.abraded.rat, prof.abraded.rat) # Preview likelihood profile
## Plot
neighborhood = max(0, which.min(prof.abraded.rat)-1000):min(length(pc.abraded.rat), which.min(prof.abraded.rat)+1000)
pc.opt.abraded.rat = opt.pc.abraded.skin.rat$par
plot(pc.abraded.rat[neighborhood], prof.abraded.rat[neighborhood], main = 'Profile, pc(abraded) rats', xlab = 'pc', ylab = 'Neg log Likelihood', type = 'l')
threshold.abraded.rat = LR.Threshold(opt.pc.abraded.skin.rat$value,1)
abline(h = threshold.abraded.rat, col = 'red')
CI.abraded.rat = LR.CI(threshold.abraded.rat, prof.abraded.rat, pc.abraded.rat)
text(.055, 20, sprintf('Best est = %.2g\nCI = %.2g - %.2g', pc.opt.abraded.rat, CI.abraded.rat[1], CI.abraded.rat[2]))
points(opt.pc.abraded.skin.rat$par, opt.pc.abraded.skin.rat$value, pch = 8, col = 'blue', cex = 1.2)
#points(full.model$par[2], opt.pc.abraded.skin$value, pch = 16, col = 'green3', cex = .8)
abline(v = opt.pc.abraded.skin.rat$par, lty = 2)
dev.off()






####################################################################
###########################   SECTION 3   ##########################
####   
####   3. Repeat fits to abraded data using the mixture model.
####      
####################################################################
####################################################################

#### THIS FUNCTION CALCULATES THE PROBABILITY OF OBSERVING ANY PARTICULAR BINARY INFECTION OUTCOME IN THE MIXTURE MODEL 
####    *This probability is defined in equation 4 of the methods
####    *Approximaiton explained in the supplementary materail

Prob.inf.given.dose.mixture = function(D0, alpha, beta, pP){
  D0 = round(D0) # Make sure you've got an integer D0
  ## Solve for the Dw value at which within-host persistence is guaranteed (Prob failure < 1e-20)
  Dw_threshold = min(ceiling(log(1e-20)/log(1-pP)), D0)
  
  ## Now, for low within-host doses (Dw<Dw_threshold), we have to calculate the probability of infection using the exact probability, as the probability of wihtin-host persistence is less than one
  Dw_below = 0:Dw_threshold
  prob_given_below_threshold = sum( (1-((1-pP)^Dw_below)) * TailRank::dbb(x=Dw_below, N=D0, u = alpha, v = beta, log = FALSE) )#pc to be
  ## For Dw>Dw_threshold, the probability of within-host persistence is 1, and the expression reduces from:
  ##  sum[Dw_threshold+1:D0](p(DI>0|Dw)p(Dw|D0)) (a)
  ##      to
  ##  sum[Dw_threshold+1:D0](p(Dw|D0))           (b)
  ##  Equation (b) is equivalent to 1-PBB(Dw_threshold+1, D0, alpha, beta), where PBB is the beta-binomial CMF
  if(Dw_threshold<D0){
    prob_given_above_threshold = 1-TailRank::pbb(Dw_threshold, N = D0, u = alpha, v = beta)
  }else{
    prob_given_above_threshold = 0
  }
  ## Return total probability
  prob_given_above_threshold+prob_given_below_threshold
}

# # ## Test accuracy of approximation
# test_approx = function(D0, alpha, beta, pP){
# D0 = round(D0)
# approximation = Prob.inf.given.dose.mixture(D0, alpha, beta, pP)
# exact_soln = sum((1-((1-pP)^(0:D0)))*TailRank::dbb(x=(0:D0), N=D0, u = alpha, v = beta, log = FALSE) )#pc to be
# ## Return error
# approximation-exact_soln
# }
# n.tests = 1000
# D0.test = runif(n.tests, 1, 10^5) # Randomly select 500 doses, alphas and betas to test the approximation
# aa.test = runif(n.tests, 1, 10^3)
# bb.test = runif(n.tests, 1, 10^5)
# cl = makeCluster(detectCores()-1)
# clusterExport(cl, varlist = c('D0.test', 'aa.test', 'bb.test', 'test_approx', 'Prob.inf.given.dose.mixture', 'pP'))
# ers= clusterMap(cl = cl, fun = test_approx, D0 = D0.test, alpha = aa.test, beta = bb.test, pP = pP, SIMPLIFY = TRUE)
# stopCluster(cl)
# plot(1:n.tests, ers, main = 'errors') # Plot
# max(ers)
# ##### All approximation errors are all < 7e-12 and negligible



#### THIS FUNCTION CALCULATES THE NEG LOG LIKELIHOOD OF A FULL NON-IP DATA SET
####    *Follow equation 0.12
####    *Output the neg log transformed likelihood
nll.aa.BB = function(pars, outcome.vec, D0.vec, pP){
  alpha = pars['aa']
  Beta = pars['BB']
  
  p_i = numeric(length(D0.vec)) #Get the probability of infection (success) given D0
  #note: P(failure, or no infection) = 1-p(success)
  for(ii in 1:length(D0.vec)){
    p_i[ii] = Prob.inf.given.dose.mixture(D0.vec[ii], alpha, Beta, pP)
  }
  # Extract the contribution of each infection trial to the likelihood
  # Include the probability of infection in the product if the outcome was infection
  # And 1-p_i if the outcome was not infection
  lk.vec = ifelse(outcome.vec == 1, p_i, 1-p_i)
  -sum(log(lk.vec)) # Output negative log likelihood
}


####################################################################
####   Maximize the likelihood to estimate pc for abraded data, hamsters
####################################################################
xx = seq(0, 1, by = .01) ## Set a vector of x values spanning 0 to 1 for plotting (below)

## Fit to abraded data:
opt.aa.BB.abraded.skin = optim(c(aa = .3, BB = 6), fn = nll.aa.BB, pP = pP, outcome.vec = Outcomes.abraded, D0.vec = Doses.abraded, method = 'L-BFGS-B', lower = c(1e-6, 1e-6), upper = c(Inf,Inf))
opt.aa.BB.abraded.skin

# Plot the estimated pc distribution using the ML values of alpha and beta
plot(xx, dbeta(x = xx, shape1 = opt.aa.BB.abraded.skin$par[1], shape2 = opt.aa.BB.abraded.skin$par[2]), main = 'pc distribution from mixture model', xlab = 'pc', ylab = 'probability density')


## Find likelihood profiles for alpha and beta
aas = c(.00005, .0001, .0005, .001, .005, seq(0.01, 2, by = 0.01))
BBs = c(.0001, .0005, .001, .005, seq(0.01, 10, by = 0.01))

# # calculate 2D profile for aa and bb:
# cl = makeCluster(detectCores()-1)
# clusterExport(cl = cl, varlist = c('aas', 'BBs', 'Outcomes.abraded', 'Doses.abraded', 'pP', 'nll.aa.BB', 'Prob.inf.given.dose.mixture'))
# one.BB.val = function(BB.in){
# sapply(aas, FUN = function(aa.in){nll.aa.BB(pars = c(aa = aa.in, BB = BB.in), outcome.vec = Outcomes.abraded, D0.vec = Doses.abraded, pP = pP)} )
# }
# prof2d = parSapply(cl, X = BBs, FUN = one.BB.val)
# stopCluster(cl)
# save(prof2d, aas, BBs, file = mixture_prof_output)
load(mixture_prof_output)
aa.mat = matrix(aas, nrow = length(aas), ncol = length(BBs), byrow = F)
bb.mat = matrix(BBs, nrow = length(aas), ncol = length(BBs), byrow = T)
contour(x = (aas), y = (BBs), z = prof2d, levels = seq(2, 10, by = 2)+opt.aa.BB.abraded.skin$value)



### Define confidence interval
## Find the LR threshold (the worst likelihood contained in the 2D confidence envelope
threshold = LR.Threshold(NLL_best = opt.aa.BB.abraded.skin$value, df = 2)
## Find the indices of profile values below the threshold
accept = which(prof2d <= threshold)
prof.acc = prof2d[accept] # Store accepted profile likelihood values (values within the envelope)
aa.coords = as.vector(aa.mat[accept]) # Extract aa and bb pairs from within the accepted indices
bb.coords = as.vector(bb.mat[accept])
# Mind the min and max aa values from in the envelope. These are the marginal CI values.
aa.max.ind = which.max(aa.coords)
aa.min.ind = which.min(aa.coords)
prof.acc[aa.max.ind]; threshold
prof.acc[aa.min.ind]; threshold
CI.aa = cbind(c(aa.coords[aa.min.ind], bb.coords[aa.min.ind]), c(aa.coords[aa.max.ind], bb.coords[aa.max.ind]))
colnames(CI.aa) = c('min', 'max')
rownames(CI.aa) = c('aa', 'bb')
## Repeat for min and max bb values. These are the marginal CI values.
bb.max.ind = which.max(bb.coords)
bb.min.ind = which.min(bb.coords)
prof.acc[bb.max.ind]; threshold
prof.acc[bb.min.ind]; threshold
CI.bb = cbind(c(aa.coords[bb.min.ind], bb.coords[bb.min.ind]), c(aa.coords[bb.max.ind], bb.coords[bb.max.ind]))
colnames(CI.bb) = c('min', 'max')
rownames(CI.bb) = c('aa', 'bb')








# ####################################################################
# ####   Repeat for rat data
# ####################################################################
## Fit to abraded data:
opt.aa.BB.abraded.skin.rat = optim(c(aa = .3, BB = 6), fn = nll.aa.BB, pP = pP.rat, outcome.vec = Outcomes.abraded.rat, D0.vec = Doses.abraded.rat, method = 'L-BFGS-B', lower = c(1e-6, 1e-6), upper = c(Inf,Inf), control = list(maxit = 100), hessian = TRUE)
opt.aa.BB.abraded.skin.rat # Output result

# Plot the estimated pc distribution using the ML values of alpha and beta
# Note that all cumulative density is concentrated around the MLE pc estimate from the basic model
# See supplementary material for interpretation
plot(xx, pbeta(xx, shape1 = opt.aa.BB.abraded.skin.rat$par[1], shape2 = opt.aa.BB.abraded.skin.rat$par[2]), main = 'pc distribution from mixture model, rats', xlab = 'pc', ylab = 'probability density')





# ####################################################################
# ####   Calculate 2D likelihood profiles (Mixture model), rats
# ####################################################################
## Find likelihood profiles for alpha and beta
opt.aa.BB.abraded.skin.rat$par
## Set up 2D grid
aas.rat = seq(25, opt.aa.BB.abraded.skin.rat$par[1]+5000, by = 25)
BBs.rat = seq(1000, opt.aa.BB.abraded.skin.rat$par[2]+50000, by = 500)
# #calculate 2D profile for aa and bb:
# cl = makeCluster(detectCores()-1)
# clusterExport(cl = cl, varlist = c('aas.rat', 'BBs.rat', 'Outcomes.abraded.rat', 'Doses.abraded.rat', 'pP.rat', 'nll.aa.BB', 'Prob.inf.given.dose.mixture'))
# one.BB.val = function(BB.in){
# sapply(aas.rat, FUN = function(aa.in){nll.aa.BB(pars = c(aa = aa.in, BB = BB.in), outcome.vec = Outcomes.abraded.rat, D0.vec = Doses.abraded.rat, pP = pP.rat)} )
# }
# prof2d.rat = parSapply(cl, X = BBs.rat, FUN = one.BB.val)
# stopCluster(cl)
# save(prof2d.rat, aas.rat, BBs.rat, file = mixture_prof_rats)
load(mixture_prof_rats)

threshold = LR.Threshold(opt.aa.BB.abraded.skin.rat$value, 2) # Define threshold nll for 95% confidence envelope
aa.mat = matrix(aas.rat, nrow = length(aas.rat), ncol = length(BBs.rat), byrow = F) # Generate a matrix of aa and bb values, of the same length s the 2d profile
bb.mat = matrix(BBs.rat, nrow = length(aas.rat), ncol = length(BBs.rat), byrow = T)
contour(x = (aas.rat), y = (BBs.rat), z = prof2d.rat, levels = seq(0, threshold, by = .5)+opt.aa.BB.abraded.skin.rat$value) # Plot contours

## Base plot
plotdf = data.frame(alpha = as.vector(aa.mat), beta = as.vector(bb.mat), nll = as.vector(prof2d.rat)) ## 2d profiles
plotdf$line_alphas = opt.pc.abraded.skin.rat$par/0.98*plotdf$beta## Add coords for line representing alpha/(alpha+beta) = 0.02
linepts = which(abs(plotdf$alpha/(plotdf$alpha+plotdf$beta) - opt.pc.abraded.skin.rat$par) < 5e-5)## extract (alpha, beta) pairs that fall on the line
# Choose some arbitrary points to plot on the line and call out their values
linepts = linepts[c(18, 72, 160)]
pp = ggplot(plotdf)+geom_point(aes(x = alpha, y = beta, color = nll)) +
  scale_color_distiller(palette = "Spectral") +
  geom_contour(aes(x = alpha, y = beta, z = nll), binwidth = 1, color = 'gray20') +
  geom_line(aes(x = line_alphas, y = beta), color = 'white', lwd = 2) +
  geom_point(data = plotdf[linepts[1],], aes(x = alpha , y = beta), shape = 7, size = 2)+
  geom_point(data = plotdf[linepts[2],], aes(x = alpha , y = beta), shape = 8, size = 2)+
  geom_point(data = plotdf[linepts[3],], aes(x = alpha , y = beta), shape = 5, size = 2)+
  geom_text(data = plotdf[linepts[1],], aes(x = alpha +900, y = beta, label = round(nll, 4))) +
  geom_text(data = plotdf[linepts[2],], aes(x = alpha +900, y = beta, label = round(nll, 4))) +
  geom_text(data = plotdf[linepts[3],], aes(x = alpha +900, y = beta, label = round(nll, 4))) +
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
AIC.basic.ham = AICcalc(n.pars = 2, nll = opt.pc.abraded.skin$value)
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
par(mfrow = c(2,1), mgp = c(2, 1, 0), mar = c(3,3,2,3))
## Plot1 - point estimates and CIs ----------
xtext = c(-1.9, -8.5, -7.1, -5.4)
ytext = c(1.001, 1.001, 1.001, 1.001)+.0005
#cols = c('firebrick1', 'deepskyblue', 'limegreen', 'chocolate1')
labs = c(expression('p'['C(abraded)']), expression('p'['C(intact)']), expression('p'['C(shaved)']), expression('p'['C(conjunctival)']))
labs = c(expression('p'[c]*' abraded'), 
         expression('p'[c]*' intact'), 
         expression('p'[c]*' shaved'),
         expression('p'[c]*' conjunctival'))
ests = c(opt.pc.abraded.skin$par, opt.pc.intact.skin$par, opt.pc.shaved.skin$par, opt.pc.conjunctival.skin$par)
plot(log10(ests[2]), 1, col = cols[2], xlim = c(-10, 0), ylim = c(.99, 1.01), pch = 16, yaxt = 'n', ylab = '', xlab = expression('value'), xaxt = 'n')
axis(1, at = -10:0, labels = 10^(-10:0))
lines(log10(c(CI.intact[1], CI.intact[2])), c(1,1), col = cols[2], lwd = 1.5)
text(xtext[2], ytext[2], labs[2], col = cols[2], cex = 1)
CI.lows = c(CI.abraded[1], CI.intact[1], CI.shaved[1], CI.conjunctival[1])
CI.highs = c(CI.abraded[2], CI.intact[2], CI.shaved[2], CI.conjunctival[2])
for(ii in c(1, 3, 4)){
  points(log10(ests[ii]), 1, col = cols[ii], pch = 16)
  lines(log10(c(CI.lows[ii], CI.highs[ii])), c(1,1), col = cols[ii], lwd = 1.5)
  text(xtext[ii], ytext[ii], labs[ii], col = cols[ii], cex = 1)
}
# Add Rat estimates
points(log10(opt.pc.abraded.skin.rat$par), .999, col = cols[1])
lines(log10(CI.pP.rat), c(.999, .999), col = cols[1])
# Add pP estimates 
#hamsters
points(log10(pP), 1, col = 'black', pch = 16)
lines(log10(CI.pP), c(1,1), col = 'black')
text(log10(pP)+.1, 1.0015, expression('p'['P']), cex = 1, col = 'black')
#rat CIs
points(log10(pP.opt.rat$par), .998, col = 'black')
lines(log10(CI.pP.rat), c(.998, .998), col = 'black')
legend('topleft', c('hamsters', 'rats'), pch = c(16, 1), bty = 'n')
mtext(3, text = 'A', at = -11, line = 1, font = 2, xpd = NA)
xlims = par('usr')[1:2]


## Plot2 - distribution of pc values ----------
hamcol = 'darkred'
ratcol = 'lightpink3'
xx = c(10^(-10:-2), seq(.01,1, by = .01))
plot(log10(xx), (pbeta(xx, shape1 = opt.aa.BB.abraded.skin$par['aa'], shape2 = opt.aa.BB.abraded.skin$par['BB'])), lty = 1, type = 'b', cex = .7, xlab = expression('p'[c]), ylab = 'cumulative density', xaxt = 'n',  col = hamcol, pch = 16)
lines(log10(xx), pbeta(xx, shape1 = opt.aa.BB.abraded.skin.rat$par['aa'], shape2 = opt.aa.BB.abraded.skin.rat$par['BB']), col = ratcol, pch = 1, lty = 2, type = 'b', cex = .7)
axis(1, at = log10(10^(-10:0)), labels = 10^(-10:0))
mtext(3, text = 'B', at = -11, line = 1, font = 2, xpd = NA)
legend('bottomright', legend = c('hamsters', 'rats'), lty = c(1,2), pch =c(16,1), col = c(hamcol, ratcol), cex = .8)


## Add inset
par(fig = c(0.09+.02,0.4+.02, 0.18, .47), new = T) 
xx = seq(0, 1, by = 0.001)
plot((xx), (dbeta(xx, shape1 = opt.aa.BB.abraded.skin$par['aa'], shape2 = opt.aa.BB.abraded.skin$par['BB'])), lty = 1, col = hamcol, pch = 16, type = 'b', cex = .7, ylab = 'density', xlab = '')
text(x = -.3, 36, 'C', font = 2, xpd = NA)

par(fig = c(0.09+.25,0.4+.25, 0.18, .47), new = T) 
plot((xx), (dbeta(xx, shape1 = opt.aa.BB.abraded.skin.rat$par['aa'], shape2 = opt.aa.BB.abraded.skin.rat$par['BB'])), lty = 2, col = ratcol, pch = 1, type = 'b', cex = .7, ylab = '', xlab = '')
text(x = -.3, 630, 'D', font = 2, xpd = NA)
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
IP.CIs = binom.confint(IP.summary[,'n.died.IP'], IP.summary[,'n.IP'], conf.level = .95, methods = 'wilson')[,c('lower', 'upper')]

ratdat = reformat.route(abraded.rat)
rat.CIs = binom.confint(ratdat[,'died'], ratdat[,'n'], conf.level = .95, methods = 'wilson')[,c('lower', 'upper')]
IP.CIs.rat = binom.confint(IP.rat[,'died'], IP.rat[,'n'], conf.level = .95, methods = 'wilson')[,c('lower', 'upper')]


#Define a function to plot probability of infection at different D0s, given parameter estimates
p.inf = function(Doses, opt.pc, pP.in = pP){
  sapply(X = round(Doses), FUN = Prob.inf.given.dose, pc = opt.pc, pP = pP.in)
}

## Repeat for prob infection, given  IP inoculation
p.inf.IP = function(dose, pP.in){
  1-(1-pP.in)^dose
}

## Repeat for prob infection, given mixture model
p.inf.mixture = function(Doses, alpha.in, beta.in, pP.in = pP){
  sapply(X = round(Doses), FUN = Prob.inf.given.dose.mixture, alpha = alpha.in, beta = beta.in, pP = pP.in)
}





# Convert color to color + transparency
col2alpha <- function(col, alpha) {
  col_rgb <- col2rgb(col)/255
  rgb(col_rgb[1], col_rgb[2], col_rgb[3], alpha = alpha)
}



{
pdf(fig3, height = 5, width = 7)
  par(mfrow = c(2,3), mgp = c(2,1,0), mar = c(3,3,3,1))
  tns.col = function(col.in, tns.pct){
    rgbvals = col2rgb(col.in)
    tnscol = rgb(rgbvals[1], rgbvals[2], rgbvals[3], tns.pct/100*255,maxColorValue = 255)
    tnscol
  }
  dd = round(10^seq(0, 9, by = .1))
  plot(dd, p.inf(dd, opt.pc.intact.skin$par), ylab = 'Probability of Infection', xlab = 'Dose', type = 'l', col = cols[1], lwd = 2, ylim = c(0,1), log = 'x', cex.lab = 1.1, xlim = c(1, 1e8))
  # polygon(x = c(dd, rev(dd)),  # CI shaded
  #         y = c(p.inf(dd, opt.pc = CI.intact[1]), rev(p.inf(dd, opt.pc = CI.intact[2]))),
  #         col = tns.col(cols[1], 40), border = NA)
  #Plot points
  points(intact.tbl[,'dose'], intact.tbl[,'died']/intact.tbl[,'n'], pch = 16, cex = 1.5)
  #Plot CI bars
  segments(x0 = intact.tbl[,'dose'], y0 = intact.CIs[,1], y1 = intact.CIs[,2])
  segments(x0 = intact.tbl[,'dose']*.7, y0 = intact.CIs[,1], x1 = intact.tbl[,'dose']*1.3)
  segments(x0 = intact.tbl[,'dose']*.7, y0 = intact.CIs[,2], x1 = intact.tbl[,'dose']*1.3)
  text(.1, 1.15, 'A', font = 2, xpd = NA)
  
  
  
  #Shaved
  plot(dd, p.inf(dd, opt.pc = opt.pc.shaved.skin$par), ylab = 'Probability of Infection', xlab = 'Dose', type = 'l', col = cols[2], lwd = 2, ylim = c(0,1), log = 'x', cex.lab = 1.1, xlim = c(1, 1e8))
  # polygon(x = c(dd, rev(dd)),
  #         y = c(p.inf(dd, opt.pc = CI.shaved[1]), rev(p.inf(dd, opt.pc = CI.shaved[2]))),
  #         col = tns.col('deepskyblue', 20), border = NA)
  points(shaved.tbl[,'dose'], shaved.tbl[,'died']/shaved.tbl[,'n'], pch = 16, cex = 1.5)
  #Plot CI bars
  segments(x0 = shaved.tbl[,'dose'], y0 = shaved.CIs[,1], y1 = shaved.CIs[,2])
  segments(x0 = shaved.tbl[,'dose']*.7, y0 = shaved.CIs[,1], x1 = shaved.tbl[,'dose']*1.3)
  segments(x0 = shaved.tbl[,'dose']*.7, y0 = shaved.CIs[,2], x1 = shaved.tbl[,'dose']*1.3)
  text(.1, 1.15, 'B', font = 2, xpd = NA)
  
  #Abraded
  ## Basic model fit and data
  plot(dd, p.inf(dd, opt.pc.abraded.skin$par), ylab = 'Probability of Infection', xlab = 'Dose', type = 'l', col = cols[3], lwd = 2, ylim = c(0,1), log = 'x',  cex.lab = 1.1, xlim = c(1, 1e8))
  # polygon(x = c(dd, rev(dd)),  # CI shaded
  #         y = c(p.inf(dd, opt.pc = CI.abraded[1]), rev(p.inf(dd, opt.pc = CI.abraded[2]))),
  #         col = tns.col(cols[3], 40), border = NA)
  points(abraded.tbl[,'dose'], abraded.tbl[,'died']/abraded.tbl[,'n'],  pch = 16, cex = 1.5)
  ## Mixture model best fit
  yy.mixture.plot = p.inf.mixture(dd, opt.aa.BB.abraded.skin$par['aa'], opt.aa.BB.abraded.skin$par['BB'])
  lines(dd, yy.mixture.plot, col = cols[3], lty = 3, lwd = 2)
  #Plot data CI bars
  segments(x0 = abraded.tbl[,'dose'], y0 = abraded.CIs[,1], y1 = abraded.CIs[,2])
  segments(x0 = abraded.tbl[,'dose']*.7, y0 = abraded.CIs[,1], x1 = abraded.tbl[,'dose']*1.3)
  segments(x0 = abraded.tbl[,'dose']*.7, y0 = abraded.CIs[,2], x1 = abraded.tbl[,'dose']*1.3)
  legend('bottomright', legend = c('basic model', 'mixture model'), col = cols[3], lty = c(1, 2), bty = 'n', lwd = 2)
  text(.1, 1.15, 'C', font = 2, xpd = NA)

  
  
  #Conjunctval
  plot(dd, p.inf(dd, opt.pc.conjunctival.skin$par), ylab = 'Probability of Infection', xlab = 'Dose', type = 'l', col = cols[4], lwd = 2, ylim = c(0,1), log = 'x',  cex.lab = 1.1, xlim = c(1, 1e8))
  # polygon(x = c(dd, rev(dd)),  # CI 
  #         y = c(p.inf(dd, opt.pc = CI.conjunctival[1]), rev(p.inf(dd, opt.pc = CI.conjunctival[2]))),
  #         col = tns.col(cols[4], 20), border = NA)
  points(conjunctival.tbl[,'dose'], conjunctival.tbl[,'died']/conjunctival.tbl[,'n'], pch = 16, cex = 1.5)
  #Plot CI bars
  segments(x0 = conjunctival.tbl[,'dose'], y0 = conjunctival.CIs[,1], y1 = conjunctival.CIs[,2])
  segments(x0 = conjunctival.tbl[,'dose']*.7, y0 = conjunctival.CIs[,1], x1 = conjunctival.tbl[,'dose']*1.3)
  segments(x0 = conjunctival.tbl[,'dose']*.7, y0 = conjunctival.CIs[,2], x1 = conjunctival.tbl[,'dose']*1.3)
  text(.1, 1.15, 'D', font = 2, xpd = NA)
  
  
  
  
  
  
  
  #Abraded, rats
  plot(dd, p.inf(dd, opt.pc.abraded.skin.rat$par, pP.in = pP.rat), ylab = 'Probability of Infection', xlab = 'Dose', type = 'l', col = cols[3], lwd = 1, ylim = c(0,1), log = 'x', cex.lab = 1.1, xlim = c(1, 1e8))
  # polygon(x = c(dd, rev(dd)),  # CI shaded
  #         y = c(p.inf(dd, opt.pc = CI.abraded.rat[1]), rev(p.inf(dd, opt.pc = CI.abraded.rat[2]))),
  #         col = tns.col(cols[3], 40), border = NA)
  points(abraded.tbl.rat[,'dose'], abraded.tbl.rat[,'died']/abraded.tbl.rat[,'n'],  pch = 16, cex = 1.5)
  yy.mixture.plot.rat = p.inf.mixture(dd, opt.aa.BB.abraded.skin.rat$par['aa'], opt.aa.BB.abraded.skin.rat$par['BB'], pP.in = pP.rat)
  lines(dd, yy.mixture.plot.rat, col = cols[3], lty = 2, lwd = 3)

  
  #Plot CI bars
  segments(x0 = abraded.tbl.rat[,'dose'], y0 = rat.CIs[,1], y1 = rat.CIs[,2])
  segments(x0 = abraded.tbl.rat[,'dose']*.7, y0 = rat.CIs[,1], x1 = abraded.tbl.rat[,'dose']*1.3)
  segments(x0 = abraded.tbl.rat[,'dose']*.7, y0 = rat.CIs[,2], x1 = abraded.tbl.rat[,'dose']*1.3)
  
  legend('bottomright', legend = c('basic model', 'mixture model'), col = cols[3], lty = c(1, 2), bty = 'n', lwd = c(1,2))
  text(.1, 1.15, 'E', font = 2, xpd = NA)
  
  
  
  
  
  #IP, hamsters
  plot(dd, p.inf.IP(dd, pP), ylab = 'Probability of Infection', xlab = 'Dose', type = 'l', ylim = c(0,1), log = 'x', cex.lab = 1.1, xlim = c(1, 1e8), lwd = 2, lty = 2)
  points(IP.summary[,'doses'], IP.summary[,'n.died.IP']/IP.summary[,'n.IP'],  pch = 15, cex = 1.7)
  segments(x0 = IP.summary[,'doses'], y0 = IP.CIs[,1], y1 = IP.CIs[,2], lwd = 1)
  segments(x0 = IP.summary[,'doses']*.7, y0 = IP.CIs[,1], x1 = IP.summary[,'doses']*1.3, lwd = 1)
  segments(x0 = IP.summary[,'doses']*.7, y0 = IP.CIs[,2], x1 = IP.summary[,'doses']*1.3, lwd = 1)
  
  xx.jitter = 10^(log10(dd)+.1)
  lines(dd, p.inf.IP(dd, pP.rat), lty = 3, lwd = 2, col = 'gray50')
  points(IP.rat[,'dose'], IP.rat[,'died']/IP.rat[,'n'],  pch = 2, cex = 1.5, col = col2alpha(col = 'gray50', alpha = 1))
  xx.jitter = 10^(log10(IP.rat[,'dose'])+.1)
  segments(x0 = xx.jitter, y0 = IP.CIs.rat[,1], y1 = IP.CIs.rat[,2], col = 'gray50', lwd = .7)
  segments(x0 = IP.rat[,'dose'], y0 = IP.CIs.rat[,1], x1 = IP.rat[,'dose']*1.5, col = 'gray50', lwd = .7)
  segments(x0 = IP.rat[,'dose'], y0 = IP.CIs.rat[,2], x1 = IP.rat[,'dose']*1.5, col = 'gray50', lwd = .7)
  
  legend('bottomright', legend = c('hamsters', 'rats'), col = c('black', 'gray50'), lty = c(2, 3), bty = 'n', lwd = c(2,2), pch = c(15, 17))
  text(.1, 1.15, 'F', font = 2, xpd = NA)
  dev.off()
}




















{
pdf(fig4, height = 7)
par(mfrow = c(2,1), mgp = c(2,1,0), mar = c(3,3,3,3))
dd = c(1:9, 10^seq(1, 11, by = .15))
p.inf.i =  p.inf(dd, opt.pc.intact.skin$par); p.inf.i[which(!is.finite(p.inf.i))] = 1
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
p.inf.s =  p.inf(dd, opt.pc.shaved.skin$par); p.inf.s[which(!is.finite(p.inf.s))] = 1
lines(dd, p.inf.s, col = cols[2], lwd = 3, lty = 3)
lower =  p.inf(dd, CI.shaved[1]); lower[which(!is.finite(lower))] = 1
upper =  p.inf(dd, CI.shaved[2]); upper[which(!is.finite(upper))] = 1
polygon(x = c(dd, rev(dd)), y = c(lower, rev(upper)), border = NA, col = col2alpha(cols[2], .5))
bounds[,2] = dd[c(max(which(p.inf.s < 0.1)), min(which(p.inf.s > 0.975)))]


#Abraded
yy.mixture.plot = p.inf.mixture(Doses = dd, alpha.in = opt.aa.BB.abraded.skin$par[1], beta.in = opt.aa.BB.abraded.skin$par[2], pP.in = pP)
lines(dd, yy.mixture.plot, col = cols[3], lty = 4, lwd = 3)
############# Mixture model CIs for plotting ###############
## In the mixture model, the confidence envelope is defined by a vague elliptical curve of (alpha, beta) pairs
##  Pairs on this curve have likelihoods near the threshold of the best nll + a threshold determined by the likelihood ratio method
## Above, we defined aa.coords and bb.coords as the set of all (alpha, beta) pairs whose neg log lilelihoods were under the threshold value
## Here, we'll downsample those values to extract pairs whose likelihood is very close to the threshold
## Randomly sample points on or near this curve from the 2D profile
##  Sample 50 points per unit on the BB axis
##  Only sample points whose likelihoods are close to the threshold. THese should give the worst acceptable fits.
aa.coords.reduced = bb.coords.reduced = prof.vals.reduced = NULL
threshold = LR.Threshold(opt.aa.BB.abraded.skin$value,  2)
for(ii in 1:8){
  valid = which(prof.acc>threshold-.1 & bb.coords > (-1+ii) & bb.coords <= (0+ii))
  smp = sample(x = valid, size = min(50, length(valid)), replace = FALSE)
  aa.coords.reduced = c(aa.coords.reduced, aa.coords[smp])
  bb.coords.reduced = c(bb.coords.reduced, bb.coords[smp])
  prof.vals.reduced = c(prof.vals.reduced, prof.acc[smp])
}
# plot(aa.coords.reduced, bb.coords.reduced) ## The reduced set of points fall on the boundary of the confidence envelope
# ## Now, calculate model-predicted probabilities of infection for each dose, using each pair of (alpha, beta) values in aa.coords.reduced, and bb.coords.reduced
wrapper = function(AA, BB){
  p.inf.mixture(Doses = dd, alpha.in = AA, beta.in = BB, pP.in = pP)
}
yy.CIs = mapply(FUN = wrapper, AA = aa.coords.reduced, BB = bb.coords.reduced)
#### Define CI envelope for plotting as the highest and lowest y value resulting from estimates using all sampled points from the CI envelope
yy.low = apply(X = yy.CIs, MARGIN = 1, FUN = min) # Estimates correspond to doses in dd
yy.high = apply(X = yy.CIs, MARGIN = 1, FUN = max)
polygon(x = c(dd, rev(dd)),
        y = c(yy.low, rev(yy.high)),
        col = tns.col(cols[3], 20), border = NA)

bounds[,4] = dd[c(1, min(which(yy.mixture.plot > 0.975)))]; #c(dd[min(1, which.max(yy.high <0.1))], 1e11)


#Conjunctval
p.inf.c =  p.inf(dd, opt.pc.conjunctival.skin$par); p.inf.c[which(!is.finite(p.inf.c))] = 1
lines(dd, p.inf.c, col = cols[4], lwd = 3)
lower =  p.inf(dd, CI.conjunctival[1]); lower[which(!is.finite(lower))] = 1
upper =  p.inf(dd, CI.conjunctival[2]); upper[which(!is.finite(upper))] = 1
polygon(x = c(dd, rev(dd)), y = c(lower, rev(upper)), border = NA, col = col2alpha(cols[4], .5))
bounds[,3] = dd[c(max(which(p.inf.c < 0.1)), min(which(p.inf.c > 0.975)))]


## IP inoculation
p.IP =  p.inf.IP(dose = dd, pP.in = pP); p.IP[which(!is.finite(p.IP))] = 1
lines(dd, p.IP, lwd = 3)
lower =  p.inf.IP(dose = dd, pP.in = CI.pP[1])
upper =  p.inf.IP(dose = dd, pP.in = CI.pP[2])
polygon(x = c(dd, rev(dd)), y = c(lower, rev(upper)), border = NA, col = col2alpha('black', .5))
bounds[,5] = dd[c(1, min(which(p.IP > 0.975)))]

#cols = c('firebrick1', 'deepskyblue', 'chocolate1', 'limegreen', 'darkgreen')
ys = seq(1.03, 1.27, length = 5)+.05
segments(x0 = bounds[1,], x1 = bounds[2,], y0 = ys, col = c(cols[c(1,2,4,3)], 'black'), lty = 2, lwd = 2)
arrows(x0 = bounds[2,], x1 = 1e11, y0 = ys, col = c(cols[c(1,2,4,3)], 'black'), lty = 1, length = .1, lwd = 2)
segments(x0 = bounds[1,], y0 = ys, y1 = ys+.02, col = c(cols[c(1,2,4,3)], 'black'), lwd =2)
segments(x0 = bounds[2,], y0 = ys, y1 = ys+.02, col = c(cols[c(1,2,4,3)], 'black'), lwd = 2)

legend(4e9, y = .6, (c('intact', 'shaved', 'conjunctival', 'abraded', 'IP, no barrier')), col = rev(c('black', cols[c(3,4,2,1)])), fill = rev(c('black', cols[c(3,4,2,1)])), cex = .8, border = NA, lty = rev(c(1, 4, 1, 2, 3)), xpd = NA, lwd = rev(c(1,1,2,2,2,2)), bty = 'n')
legend(1e3, y = 1.67, c('infection probable at this dose', 'infection almost certain at this dose'), lty = c(2, 1), cex = .8, xpd = NA, bty = 'n')
text(.1, 1.3, 'A', font = 2, xpd = NA)



#### Repeat for rats
p.inf.a =  p.inf(dd, opt.pc.abraded.skin.rat$par, pP.in = pP.rat); p.inf.a[which(!is.finite(p.inf.a))] = 1
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
p.IP =  p.inf.IP(dose = dd, pP.in = pP.rat); p.IP[which(!is.finite(p.IP))] = 1
lines(dd, p.IP, lwd = 3)
lower =  p.inf.IP(dose = dd, pP.in = CI.pP.rat[1])
upper =  p.inf.IP(dose = dd, pP.in = CI.pP.rat[2])
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
legend(1e3, y = 1.5, c('infection probable at this dose', 'infection almost certain at this dose'), lty = c(2, 1), cex = .8, xpd = NA, bty = 'n')
text(.1, 1.3, 'B', font = 2, xpd = NA)
dev.off()
}



### Write a table template for experimental results
colnames(IP.summary) = c('dose', 'n', 'died')
IP.summary = IP.summary[c('n', 'died', 'dose')]
IP.summary = IP.summary[8:1, ]


# Calculate percentages infected
intact.tbl = data.frame(intact.tbl)
intact.tbl$pct = intact.tbl$died/intact.tbl$n
shaved.tbl = data.frame(shaved.tbl)
shaved.tbl$pct = shaved.tbl$died/shaved.tbl$n
conjunctival.tbl = data.frame(conjunctival.tbl)
conjunctival.tbl$pct = conjunctival.tbl$died/conjunctival.tbl$n
abraded.tbl = data.frame(abraded.tbl)
abraded.tbl$pct = abraded.tbl$died/abraded.tbl$n
IP.summary = data.frame(IP.summary)
IP.summary$pct = IP.summary$died/IP.summary$n

ham.data = rbind(c('intact', NA,NA,NA), intact.tbl, c('shaved', NA, NA,NA), shaved.tbl, c('abraded', NA,NA,NA), abraded.tbl, c('conjunctival', NA,NA,NA), conjunctival.tbl, c('IP', NA,NA,NA), IP.summary)


abraded.tbl.rat = data.frame(abraded.tbl.rat)
abraded.tbl.rat$pct = abraded.tbl.rat$died/abraded.tbl.rat$n
IP.tbl.rat = data.frame(IP.tbl.rat)
IP.tbl.rat$pct = IP.tbl.rat$died/IP.tbl.rat$n
rat.data = rbind(c('abraded', NA,NA), abraded.tbl.rat, c('IP', NA,NA), IP.tbl.rat)

## Write output files
write.csv(ham.data, file = hamdat_summary, row.names = FALSE)
write.csv(rat.data, file = ratdat_summary, row.names = FALSE)







