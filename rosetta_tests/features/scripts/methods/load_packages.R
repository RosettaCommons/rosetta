
#######################################################################
#                                                                     #
# Use in scripts to load packages.  If the packages are not present,  #
# it will ask the user if they want to install them from CRAN.        #
#                                                                     #
#######################################################################

load_packages <- function(pkgs, verbose=FALSE){
	failed_pkgs <- c()
	for(pkg in pkgs){
		if(verbose){        
                	success <- require(pkg, quietly=FALSE, warn.conflicts=TRUE, character.only=TRUE)
                } else {
                	success <- require(pkg, quietly=TRUE, warn.conflicts=FALSE, character.only=TRUE)
                }
		if(!success){
			failed_pkgs <- c(failed_pkgs, pkg)
		}
	}
	if(length(failed_pkgs) > 0){
		cat("Failed to load the following packages:\n")
		for(pkg in pkgs){
			cat("    ", pkg, "\n")
		}
		valid_response = FALSE
		while(!valid_response){
			cat("Would you like to try to install them from CRAN? (y/n) \n")
			ans <- "n"
			ans <- readChar(file("stdin"), nchars=1)
			if(substr(ans, 1L, 1L) == "y"){
				valid_response = TRUE

				cat("Checking if lib directory is writable...\n")
				#check if the location to install packages is writable this is
				#similar to what install.packages does, but install.packages
				#checks for interactive().  It appears to not be possible to
				#set interactive() to TRUE. (Note: I've stripped out the
				#windows specific behavior).
				lib=.libPaths()[1L]
				lib_is_writable <- file.info(lib)$isdir && (file.access(lib, 2) == 0)
				if (!lib_is_writable){
					cat(sprintf("'lib = \"%s\"' is not writable\n", lib))
					userdir <- unlist(strsplit(Sys.getenv("R_LIBS_USER"),
            .Platform$path.sep))[1L]
					if (!file.exists(userdir)) {
						msg <- "Would you like to create a personal library\n'%s'\nto install packages into?"
						cat(paste(sprintf(msg, userdir), " (y/n) "))
						ans <- readChar(file("stdin"), nchar=1)
						if (substr(ans, 1L, 1L) == "n")
							stop("unable to install packages")
						if (!dir.create(userdir, recursive = TRUE))
							stop("unable to create ", sQuote(userdir))
					}
					lib <- userdir
				}
				cat(".libPaths()", .libPaths())
				for(pkg in failed_pkgs){
					install.packages(pkg, lib, contriburl=contrib.url("http://cran.fhcrc.org/"))
					success <- require(pkg, quietly=TRUE, warn.conflicts=FALSE, character.only=TRUE)
					
				}
			} else if(substr(ans, 1L, 1L) == "n"){
				valid_response = TRUE
				cat("To proceed please install the following packages from with in R:\n")
				for(pkg in failed_pkgs){
					cat(paste("    install.packages('",pkg,"')\n", sep=""))
				}
				q()
			} else {
				valid_response = FALSE
				cat("Unrecognized response:", ans, "\n")
				cat("please press 'y', or 'n'\n")
			}
		}
	}
}
