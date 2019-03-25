## This script checks numerically that the probability of infection in the basic model can be written as: P(Infection|Dose) = 1-exp(D0*p_p*p_c)


## Long-form equation:
## Let d represent the expected inoculum dose
## Let D0 represent the actual inoculum size
## Let Dw represent the reaches the within-host environment
## Let DI represent the infecting dose
## P(DI|d) = P(D0|d)P(Dw|D0)P(DI|Dw) ## The probability of a given infectious dose, given the expected inoculum size
## P(I|d) = sum_DI=1:Inf sum_Dw=DI:Inf sum_D0=Dw:Inf (P(D0|d)P(Dw|D0)P(DI|Dw)) ## Long form probability of infection, , given the expected inoculum size

 
## Calculate the summand after algebraic rearrangement into three different poisson densities
compare = function(d, D0, Dw, DI, pp, pc){
  ## This is the long-form equation for P(DI|d), as given in ##
  original = dpois(x = D0, lambda = d)*dbinom(x = Dw, size = D0, prob = pc)*dbinom(x = DI, size = Dw, prob = pp)
  ## This is an alternate representation of the summand, as shown in ##
  rearranged = dpois(x = DI, lambda = d*pc*pp)*dpois(x = D0-Dw, lambda = d*(1-pc))*dpois(x = Dw-DI, lambda = d*pc*(1-pp))
  return(c(original, rearranged))
}

compare(1000, 998, 100, 20, .15, .001)
