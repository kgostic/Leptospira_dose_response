## This script checks numerically that the probability of infection in the basic model can be written as: P(Infection|Dose) = 1-exp(D0*p_p*p_c)


## Long-form equation:
## Let d represent the expected inoculum dose
## Let D0 represent the actual inoculum size
## Let Dw represent the reaches the within-host environment
## Let DI represent the infecting dose
## P(DI|d) = P(D0|d)P(Dw|D0)P(DI|Dw) ## The probability of a given infectious dose, given the expected inoculum size
##    P(Infection|d) = sum_DI=1:Inf sum_Dw=DI:Inf sum_D0=Dw:Inf P(D0|d)P(Dw|D0)P(DI|Dw)) ## Long form probability of infection, , given the expected inoculum size

 
## Calculate the summand after algebraic rearrangement into three different poisson densities
compare = function(d, D0, Dw, DI, pp, pc){
  ## This is the long-form equation for P(DI|d), as given in equation 4 of the main text methods
  original = dpois(x = D0, lambda = d)*dbinom(x = Dw, size = D0, prob = pc)*dbinom(x = DI, size = Dw, prob = pp)
  ## This is an alternate representation of the summand, as shown in equation 5 of the main text methods
  rearranged = dpois(x = DI, lambda = d*pc*pp)*dpois(x = D0-Dw, lambda = d*(1-pc))*dpois(x = Dw-DI, lambda = d*pc*(1-pp))
  return(c(original, rearranged))
}


## Test the comparison function, which will return the value of the summand using the original long-form, and the algebraic rearrangement given in equation 5
compare(1000, 998, 100, 20, .15, .001)


## Test for 1000 randomly chosen values of d, D0, Dw, DI, pp and pc
d.test = round(runif(n = 1000, 1, 10^9))
D0.test = sapply(d.test, function(dd) rpois(n = 1, dd))
Dw.test = round(sapply(D0.test, function(D0) runif(1, 1, D0)))
DI.test = round(sapply(Dw.test, function(Dw) runif(1, 1, Dw)))
pp.test = runif(1000, 0, 1)
pc.test = runif(1000, 0, 1)

## Output comparison
test.outputs = mapply(FUN = compare, d= d.test, D0 = D0.test, Dw = Dw.test, DI = DI.test, pp = pp.test, pc = pc.test)

## Verify that both versions of the summand give exactly the same answer
test.outputs[1,]-test.outputs[2,]

## Q.E.D.

