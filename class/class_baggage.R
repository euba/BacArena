# this .R file contains functions, which are usefull for all classes (therefore not class methods) to perform simple tasks

initArena <- function(specs, #list with species models to analyze
                      #distri, #distribution that bacteria should have on the grid
                      n, m, iter, seed, epsilon, #Arena attributes
                      smax, diffconst, diffmat, gradient, #Substrate attributes
                      model, growth #Bac/Organism attributes
                      ){
  lapply(specs, function(x, n, m, iter, seed, epsilon, growth){ #initialize the bacterial population
    x <- runif(1, 1, n)
    y <- runif(1, 1, m)
    growth <- rnorm()
    return(Bac(x=x, y=y, model=x, growth=1, n=1, m=1, seed=100, iter=100, epsilon=2))
  }, n=n, m=m, iter=iter, seed=seed, epsilon=epsilon, growth=growth)
}
