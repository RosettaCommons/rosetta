#-*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

sort.data.frame <- function(x, by){
  # Author: Kevin Wright
  # with some ideas from Andy Liaw
  # http://tolstoy.newcastle.edu.au/R/help/04/07/1076.html
  
  # x: A data.frame
  # by: A one-sided formula using + for ascending and - for descending
  #     Sorting is left to right in the formula
  
  # Useage is:
  # library(nlme);
  # data(Oats)
  # sort(Oats, by= ~nitro-Variety)
  
  if(by[[1]] != "~")
    stop("Argument 'by' must be a one-sided formula.")
  
  # Make the formula into character and remove spaces
  formc <- as.character(by[2]) 
  formc <- gsub(" ", "", formc) 
  # If the first character is not + or -, add +
  if(!is.element(substring(formc, 1, 1), c("+", "-")))
    formc <- paste("+", formc, sep = "")
  
  # Extract the variables from the formula
  vars <- unlist(strsplit(formc, "[\\+\\-]"))    
  vars <- vars[vars != ""] # Remove any extra "" terms
  
  # Build a list of arguments to pass to "order" function
  calllist <- list()
  pos <- 1 # Position of + or -
  for(i in 1:length(vars)){
    varsign <- substring(formc, pos, pos)
    pos <- pos + 1 + nchar(vars[i])
    if(is.factor(x[, vars[i]])){
      if(varsign == "-") {
        calllist[[i]] <- -rank(x[, vars[i]])
      } else {
        calllist[[i]] <- rank(x[, vars[i]])
      }
    } else {
      if(varsign == "-") {
        calllist[[i]] <- -x[, vars[i]]
      } else {
        calllist[[i]] <- x[,vars[i]]
      }
    }
  }
  return(x[do.call("order", calllist), ])
}

capwords <- function(s, strict = FALSE) {
  #From ?toupper example
  cap <- function(s) paste(toupper(substring(s,1,1)),
    {s <- substring(s,2); if(strict) tolower(s) else s},
    sep = "", collapse = " " )
    sapply(strsplit(s, split = " "), cap, USE.NAMES = !is.null(names(s)))
}