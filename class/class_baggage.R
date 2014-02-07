# this .R file contains functions, which are usefull for all classes (therefore not class methods) to perform simple tasks

changeObj <- function(mod, newobj){ # changes the objective function of a model according to the given reaction id
            mod@obj_coef <- rep(0, length(mod@obj_coef))
            mod@obj_coef[which(mod@react_id==newobj)] <- 1
            return(mod)
            }
