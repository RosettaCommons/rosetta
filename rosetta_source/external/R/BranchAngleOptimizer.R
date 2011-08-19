### Read the parameter file into a list of matricies 
read.branch_angle_undefined <- function(filename, numneighbors) {

	# determine how many parameters go with each set
	if (numneighbors == 3) {
		numeach <- 6
	} else if (numneighbors == 4) {
		numeach <- 12
	} else {
		stop("Invalid number of neighbors")
	}

	# read all the parameters into a single long vector
	flat_params <- scan(filename)
	
	if (length(flat_params) %% numeach != 0) {
		stop(paste("Number of parameters not a multiple of", numeach))
	}
	
	paramlist <- list()
	
	# iterate over the long vector and turn it into a list of matricies
	for (i in seq(1, length(flat_params) / numeach)) {
	
		# get a matrix out of the correct range of flat_params
		parammat <- matrix(flat_params[numeach*(i-1)+(1:numeach)], ncol = 2, byrow = TRUE)
		# convert theta0 to radians
		parammat[,2] <- parammat[,2]/180*pi
		# add row and column names
		rownames(parammat) <- c("m1_m2", "m1_b1", "m2_b1", "m1_b2", "m2_b2", "b1_b2")[1:(numeach/2)]
		colnames(parammat) <- c("Ktheta", "theta0")
		
		# append the matrix to the list
		paramlist <- c(paramlist, list(parammat))
	}
	
	paramlist
}

### Make a vector of named parameters from a matrix
namedparams <- function(parammat) {

	params <- as.list(as.vector(parammat))
	names(params) <- paste(rep(colnames(parammat), each = nrow(parammat)), 
	                       rep(rownames(parammat), ncol(parammat)), sep = "_")

	params
}

### Expression for computing the 3 neighbor energy directly from spherical coordinates
paener3 <- expression(Ktheta_m1_m2*(acos(sin(theta_m1)*sin(theta_m2)*cos(phi_m1-phi_m2) + 
                                         cos(theta_m1)*cos(theta_m2))-theta0_m1_m2)^2 +
                      Ktheta_m1_b1*(acos(sin(theta_m1)*sin(theta_b1)*cos(phi_m1-phi_b1) + 
                                         cos(theta_m1)*cos(theta_b1))-theta0_m1_b1)^2 +
                      Ktheta_m2_b1*(acos(sin(theta_m2)*sin(theta_b1)*cos(phi_m2-phi_b1) + 
                                         cos(theta_m2)*cos(theta_b1))-theta0_m2_b1)^2)

### Corresponding derivatives for optimization
paenerderivall3 <- deriv(paener3, c("theta_m2", "theta_b1"))
paenerderivbranch3 <- deriv(paener3, c("theta_b1"))

### Expression for computing the 4 neighbor energy directly from spherical coordinates
paener4 <- expression(Ktheta_m1_m2*(acos(sin(theta_m1)*sin(theta_m2)*cos(phi_m1-phi_m2) + 
                                         cos(theta_m1)*cos(theta_m2))-theta0_m1_m2)^2 +
                      Ktheta_m1_b1*(acos(sin(theta_m1)*sin(theta_b1)*cos(phi_m1-phi_b1) + 
                                         cos(theta_m1)*cos(theta_b1))-theta0_m1_b1)^2 +
                      Ktheta_m2_b1*(acos(sin(theta_m2)*sin(theta_b1)*cos(phi_m2-phi_b1) + 
                                         cos(theta_m2)*cos(theta_b1))-theta0_m2_b1)^2 +
                      Ktheta_m1_b2*(acos(sin(theta_m1)*sin(theta_b2)*cos(phi_m1-phi_b2) + 
                                         cos(theta_m1)*cos(theta_b2))-theta0_m1_b2)^2 +
                      Ktheta_m2_b2*(acos(sin(theta_m2)*sin(theta_b2)*cos(phi_m2-phi_b2) + 
                                         cos(theta_m2)*cos(theta_b2))-theta0_m2_b2)^2 +
                      Ktheta_b1_b2*(acos(sin(theta_b1)*sin(theta_b2)*cos(phi_b1-phi_b2) + 
                                         cos(theta_b1)*cos(theta_b2))-theta0_b1_b2)^2)

### Corresponding derivatives for optimization
paenerderivall4 <- deriv(paener4, c("theta_m2", "phi_b1", "theta_b1", "phi_b2", "theta_b2"))
paenerderivbranch4 <- deriv(paener4, c("phi_b1", "theta_b1", "phi_b2", "theta_b2"))

### Optimize all spherical coordinates for a 4 neighbor system
opt_all_angles3 <- function(parammat, initialangles = NULL) {

	if (is.null(initialangles)) {
		initialangles <- list(theta_m2 = 120, theta_b1 = 120)
		for (i in seq(along = initialangles))
        	initialangles[[i]] <- initialangles[[i]]/180*pi
	}
	
	# phi_b1 enforces planarity
	fixedangles <- list(phi_m1 = 0, theta_m1 = 0, phi_m2 = 0, phi_b1 = pi)
	
	params <- c(namedparams(parammat), fixedangles, initialangles)
    
    fn <- function(x, params) {
        params$theta_m2 <- x[1]
        params$theta_b1 <- x[2]
        eval(paener3, envir = params)
    }
    gr <- function(x, params) {
        params$theta_m2 <- x[1]
        params$theta_b1 <- x[2]
        attr(eval(paenerderivall3, envir = params), "gradient")
    }
    
    optresult <- optim(unlist(initialangles), fn, gr, method = "BFGS", params = params)
    
    lower <- c(0, 0)
    upper <- c(pi, pi)
    if (any(optresult$par < lower || optresult$par > upper))
    	warning("Angles out of range!")
    
    optresult
}

### Optimize branching spherical coordinates for a 3 neighbor system
opt_branch_angles3 <- function(parammat, theta_m2, initialangles = NULL) {

	if (is.null(initialangles)) {
		initialangles <- list(theta_b1 = 120)
		for (i in seq(along = initialangles))
        	initialangles[[i]] <- initialangles[[i]]/180*pi
	}
	
	# phi_b1 enforces planarity
	fixedangles <- list(phi_m1 = 0, theta_m1 = 0, phi_m2 = 0, theta_m2 = theta_m2, phi_b1 = pi)
	
	params <- c(namedparams(parammat), fixedangles, initialangles)
    
    fn <- function(x, params) {
        params$theta_b1 <- x[1]
        eval(paener3, envir = params)
    }
    gr <- function(x, params) {
        params$theta_b1 <- x[1]
        attr(eval(paenerderivbranch3, envir = params), "gradient")
    }
    
    optresult <- optim(unlist(initialangles), fn, gr, method = "BFGS", params = params)
    
    lower <- c(0)
    upper <- c(pi)
    if (any(optresult$par < lower || optresult$par > upper))
    	warning("Angles out of range!")
    
    optresult
}

### Optimize all spherical coordinates for a 4 neighbor system
opt_all_angles4 <- function(parammat, initialangles = NULL) {

	if (is.null(initialangles)) {
		initialangles <- list(theta_m2 = 109.5, phi_b1 = 120, theta_b1 = 109.5, 
		                      phi_b2 = -120, theta_b2 = 109.5)
		for (i in seq(along = initialangles))
        	initialangles[[i]] <- initialangles[[i]]/180*pi
	}
	
	fixedangles <- list(phi_m1 = 0, theta_m1 = 0, phi_m2 = 0)
	
	params <- c(namedparams(parammat), fixedangles, initialangles)
    
    fn <- function(x, params) {
        params$theta_m2 <- x[1]
        params$phi_b1 <- x[2]
        params$theta_b1 <- x[3]
        params$phi_b2 <- x[4]
        params$theta_b2 <- x[5]
        eval(paener4, envir = params)
    }
    gr <- function(x, params) {
        params$theta_m2 <- x[1]
        params$phi_b1 <- x[2]
        params$theta_b1 <- x[3]
        params$phi_b2 <- x[4]
        params$theta_b2 <- x[5]
        attr(eval(paenerderivall4, envir = params), "gradient")
    }
    
    optresult <- optim(unlist(initialangles), fn, gr, method = "BFGS", params = params)
    
    lower <- c(0, -pi, 0, -pi, 0)
    upper <- c(pi, pi, pi, pi, pi)
    if (any(optresult$par < lower || optresult$par > upper))
    	warning("Angles out of range!")
    
    optresult
}

### Optimize branching spherical coordinates for a 4 neighbor system
opt_branch_angles4 <- function(parammat, theta_m2, initialangles = NULL) {

	if (is.null(initialangles)) {
		initialangles <- list(phi_b1 = 119, theta_b1 = 113.5, 
		                      phi_b2 = -119, theta_b2 = 108)
		for (i in seq(along = initialangles))
        	initialangles[[i]] <- initialangles[[i]]/180*pi
	}
	
	fixedangles <- list(phi_m1 = 0, theta_m1 = 0, phi_m2 = 0, theta_m2 = theta_m2)
	
	params <- c(namedparams(parammat), fixedangles, initialangles)
    
    fn <- function(x, params) {
        params$phi_b1 <- x[1]
        params$theta_b1 <- x[2]
        params$phi_b2 <- x[3]
        params$theta_b2 <- x[4]
        eval(paener4, envir = params)
    }
    gr <- function(x, params) {
        params$phi_b1 <- x[1]
        params$theta_b1 <- x[2]
        params$phi_b2 <- x[3]
        params$theta_b2 <- x[4]
        attr(eval(paenerderivbranch4, envir = params), "gradient")
    }
    
    optresult <- optim(unlist(initialangles), fn, gr, method = "BFGS", params = params)
    
    lower <- c(-pi, 0, -pi, 0)
    upper <- c(pi, pi, pi, pi)
    if (any(optresult$par < lower || optresult$par > upper))
    	warning("Angles out of range!")
    
    optresult
}

### Self-starting model for fitting Ktheta given fixed minimum energy and optimal theta
SSparaKtheta <- selfStart(~ energy0 + Ktheta*(x - theta0)^2, function(mCall, data, LHS) {

    value <- 1
    names(value) <- "Ktheta"
    value

}, c("Ktheta"), c("x", "Ktheta", "theta0", "energy0"))

### Convenience function for fitting Ktheta parameter
parafitKtheta <- function(x, y, theta0, energy0, plotpred = FALSE) {

    theta0 <- rep(unname(theta0), length(x))
    energy0 <- rep(unname(energy0), length(x))
    fit <- nls(y ~ SSparaKtheta(x, Ktheta, theta0, energy0), data.frame(x, y, theta0, energy0))
    
    if (plotpred) {
    	plot(x, y, type="l")
        points(x, predict(fit), type = "l", col = "red")
    }

    c(fit$m$getPars(), theta0 = theta0[1], energy0 = energy0[1])
}

### Self-starting model for fitting data to a quadratic function
SSquad <- selfStart(~ a + b*x + c*x^2, function(mCall, data, LHS) {

    xy <- sortedXyData(mCall[["x"]], LHS, data)

    linearfit <- lsfit(xy[,1], xy[,2])

    value <- c(linearfit$coefficients[1], linearfit$coefficients[2], 0)
    names(value) <- mCall[c("a", "b", "c")]
    value

}, c("a", "b", "c"))

### Convenience function for doing quadratic fits
quadfit <- function(x, y, plotpred = FALSE) {

    fit <- nls(y ~ SSquad(x, a, b, c), data.frame(x, y))
    
    if (plotpred)
        points(x, predict(fit), type = "l", col = "red")

    fit$m$getPars()
}

### Calculate coefficients for a 3 neighbor system
calculate_coefs3 <- function(parammat, plot.coefs = FALSE) {

	# first do an initial optimization of all angles to find the overall minimum
	optall <- opt_all_angles3(parammat)
	
	# this will be the matrix everything is stored in 
	coefmat <- matrix(c(optall$value, optall$par), nrow = 1)
	
	# calculate coefficients only for energies less than 20 kcal/mol over the minimum
	energycut <- optall$value + 20
	
	# calculate optimized angles ascending by single degrees from the overall minimum
	previnitial <- as.list(optall$par[-1])
	theta_m2 <- optall$par["theta_m2"] + 1/180*pi
	while((optbranch <- opt_branch_angles3(parammat, theta_m2, previnitial))$value < energycut) {
	
		coefmat <- rbind(coefmat, c(optbranch$value, theta_m2, optbranch$par))
		previnitial <- as.list(optbranch$par)
		theta_m2 <- theta_m2 + 1/180*pi
	}
	
	# calculate optimized angles descending by single degrees from the overall minimum
	previnitial <- as.list(optall$par[-1])
	theta_m2 <- optall$par["theta_m2"] - 1/180*pi
	while((optbranch <- opt_branch_angles3(parammat, theta_m2, previnitial))$value < energycut) {
	
		coefmat <- rbind(c(optbranch$value, theta_m2, optbranch$par), coefmat)
		previnitial <- as.list(optbranch$par)
		theta_m2 <- theta_m2 - 1/180*pi
	}
	
	colnames(coefmat) <- c("energy", names(optall$par))
	
	minindex <- which.min(coefmat[,"energy"])
	
	# The fits for 3 neighbor systems are perfect
	#coef_overall <- parafitKtheta(coefmat[,"theta_m2"], coefmat[,"energy"], optall$par["theta_m2"], 
	#                              optall$value)
	Ktheta <- mean((coefmat[-minindex,"energy"] - coefmat[minindex,"energy"]) /
	          (coefmat[-minindex,"theta_m2"] - coefmat[minindex,"theta_m2"])^2)
	coef_overall <- c(Ktheta = Ktheta, theta0 = coefmat[minindex,"theta_m2"],
	                  energy0 = coefmat[minindex,"energy"])
	
	coef_cols <- colnames(coefmat)[-(1:2)]
	coef_quad <- matrix(nrow = length(coef_cols), ncol = 3)
	colnames(coef_quad) <- c("A", "B", "C")
	rownames(coef_quad) <- coef_cols
	
	for (parname in coef_cols) {
	
		# The fits for 3 neighbor systems are perfect
		#coef_quad[parname,] <- quadfit(coefmat[,"theta_m2"], coefmat[,parname])
		coef_quad[parname,] <- c(lsfit(coefmat[,"theta_m2"], coefmat[,parname])$coef, 0)
	}
	
	if (plot.coefs) {
	
		plot(coefmat[,"theta_m2"]/pi*180, coefmat[,"energy"], pch = 20, 
		     xlab = "M2 Bond Angle (degrees)", ylab = "Energy (kcal/mol)")
		fitenergy <- coef_overall["Ktheta"]*(coefmat[,"theta_m2"] - coef_overall["theta0"])^2 + 
		             coef_overall["energy0"]
		points(coefmat[,"theta_m2"]/pi*180, fitenergy, type = "l", col = "red")
		points(coef_overall["theta0"]/pi*180, coef_overall["energy0"], cex = 2, col = "red")
		#readline()
		
		plot(coefmat[,"theta_m2"]/pi*180, coefmat[,"theta_b1"]/pi*180, pch = 20, 
		     xlab = "M2 Bond Angle (degrees)", ylab = "B1 Bond Angle (degrees)")
		fitangle <- coef_quad["theta_b1","A"] + coef_quad["theta_b1","B"]*coefmat[,"theta_m2"] +
		            coef_quad["theta_b1","C"]*coefmat[,"theta_m2"]^2
		points(coefmat[,"theta_m2"]/pi*180, fitangle/pi*180, type = "l", col = "red")
		points(coefmat[minindex,"theta_m2"]/pi*180, fitangle[minindex]/pi*180, cex = 2, col = "red")
		#readline()
	}
	
	list(overall = coef_overall, quad = coef_quad)
}

### Calculate coefficients for a 4 neighbor system
calculate_coefs4 <- function(parammat, plot.coefs = FALSE) {

	# first do an initial optimization of all angles to find the overall minimum
	optall <- opt_all_angles4(parammat)
	
	# this will be the matrix everything is stored in 
	coefmat <- matrix(c(optall$value, optall$par), nrow = 1)
	
	# calculate coefficients only for energies less than 20 kcal/mol over the minimum
	energycut <- optall$value + 20
	
	# calculate optimized angles ascending by single degrees from the overall minimum
	previnitial <- as.list(optall$par[-1])
	theta_m2 <- optall$par["theta_m2"] + 1/180*pi
	while((optbranch <- opt_branch_angles4(parammat, theta_m2, previnitial))$value < energycut) {
	
		coefmat <- rbind(coefmat, c(optbranch$value, theta_m2, optbranch$par))
		previnitial <- as.list(optbranch$par)
		theta_m2 <- theta_m2 + 1/180*pi
	}
	
	# calculate optimized angles descending by single degrees from the overall minimum
	previnitial <- as.list(optall$par[-1])
	theta_m2 <- optall$par["theta_m2"] - 1/180*pi
	while((optbranch <- opt_branch_angles4(parammat, theta_m2, previnitial))$value < energycut) {
	
		coefmat <- rbind(c(optbranch$value, theta_m2, optbranch$par), coefmat)
		previnitial <- as.list(optbranch$par)
		theta_m2 <- theta_m2 - 1/180*pi
	}
	
	colnames(coefmat) <- c("energy", names(optall$par))
	
	coef_overall <- parafitKtheta(coefmat[,"theta_m2"], coefmat[,"energy"], optall$par["theta_m2"], 
	                              optall$value)
	
	coef_cols <- colnames(coefmat)[-(1:2)]
	coef_quad <- matrix(nrow = length(coef_cols), ncol = 3)
	colnames(coef_quad) <- c("A", "B", "C")
	rownames(coef_quad) <- coef_cols
	
	for (parname in coef_cols) {
	
		coef_quad[parname,] <- quadfit(coefmat[,"theta_m2"], coefmat[,parname])
	}
	
	if (plot.coefs) {
	
		plot(coefmat[,"theta_m2"]/pi*180, coefmat[,"energy"], pch = 20, 
		     xlab = "M2 Bond Angle (degrees)", ylab = "Energy (kcal/mol)")
		fitenergy <- coef_overall["Ktheta"]*(coefmat[,"theta_m2"] - coef_overall["theta0"])^2 + 
		             coef_overall["energy0"]
		points(coefmat[,"theta_m2"]/pi*180, fitenergy, type = "l", col = "red")
		points(coef_overall["theta0"]/pi*180, coef_overall["energy0"], cex = 2, col = "red")
		#readline()
		
		def.par <- par(no.readonly = TRUE)
		layout(matrix(seq_len(4), nrow = 2))
		minindex <- which.min(coefmat[,"energy"])
		
		plot(coefmat[,"theta_m2"]/pi*180, coefmat[,"phi_b1"]/pi*180, pch = 20, cex = .5,
		     xlab = "M2 Bond Angle (degrees)", ylab = "B1 Torsion Offset (degrees)")
		fitangle <- coef_quad["phi_b1","A"] + coef_quad["phi_b1","B"]*coefmat[,"theta_m2"] +
		            coef_quad["phi_b1","C"]*coefmat[,"theta_m2"]^2
		points(coefmat[,"theta_m2"]/pi*180, fitangle/pi*180, type = "l", col = "red")
		points(coefmat[minindex,"theta_m2"]/pi*180, fitangle[minindex]/pi*180, cex = 2, col = "red")
		
		plot(coefmat[,"theta_m2"]/pi*180, coefmat[,"theta_b1"]/pi*180, pch = 20, cex = .5,
		     xlab = "M2 Bond Angle (degrees)", ylab = "B1 Bond Angle (degrees)")
		fitangle <- coef_quad["theta_b1","A"] + coef_quad["theta_b1","B"]*coefmat[,"theta_m2"] +
		            coef_quad["theta_b1","C"]*coefmat[,"theta_m2"]^2
		points(coefmat[,"theta_m2"]/pi*180, fitangle/pi*180, type = "l", col = "red")
		points(coefmat[minindex,"theta_m2"]/pi*180, fitangle[minindex]/pi*180, cex = 2, col = "red")
		
		plot(coefmat[,"theta_m2"]/pi*180, coefmat[,"phi_b2"]/pi*180, pch = 20, cex = .5,
		     xlab = "M2 Bond Angle (degrees)", ylab = "B2 Torsion Offset (degrees)")
		fitangle <- coef_quad["phi_b2","A"] + coef_quad["phi_b2","B"]*coefmat[,"theta_m2"] +
		            coef_quad["phi_b2","C"]*coefmat[,"theta_m2"]^2
		points(coefmat[,"theta_m2"]/pi*180, fitangle/pi*180, type = "l", col = "red")
		points(coefmat[minindex,"theta_m2"]/pi*180, fitangle[minindex]/pi*180, cex = 2, col = "red")
		
		plot(coefmat[,"theta_m2"]/pi*180, coefmat[,"theta_b2"]/pi*180, pch = 20, cex = .5,
		     xlab = "M2 Bond Angle (degrees)", ylab = "B2 Bond Angle (degrees)")
		fitangle <- coef_quad["theta_b2","A"] + coef_quad["theta_b2","B"]*coefmat[,"theta_m2"] +
		            coef_quad["theta_b2","C"]*coefmat[,"theta_m2"]^2
		points(coefmat[,"theta_m2"]/pi*180, fitangle/pi*180, type = "l", col = "red")
		points(coefmat[minindex,"theta_m2"]/pi*180, fitangle[minindex]/pi*180, cex = 2, col = "red")
		par(def.par)
		#readline()
	}
	
	list(overall = coef_overall, quad = coef_quad)
}

# Write out parameters and calculated coefficients
catparams <- function(parammat, paramlist, file = "") {

	for (i in seq_len(nrow(parammat))) {
		cat(as.character(parammat[i,1]), as.character(parammat[i,2]/pi*180), file = file, sep = "\t")
		cat("\n", file = file)
	}
	cat("\n", file = file)

	cat(paramlist$overall[1], paramlist$overall[2]/pi*180, paramlist$overall[3], file = file, 
	    sep = "\t")
	cat("\n\n", file = file)
	
	if (nrow(paramlist$quad) == 1) {
		cat(pi, 0, 0, file = file)
		cat("\n", file = file)
	}
	for (i in seq_len(nrow(paramlist$quad))) {
		cat(paramlist$quad[i,], file = file, sep = "\t")
		cat("\n", file = file)
	}
	cat("\n\n", file = file)
}

process_database <- function(dbpath) {

	undef1path <- file.path(dbpath, "branch_angle", "branch_angle_1_undefined.txt")
	
	if (file.exists(undef1path)) {
	
		parammatlist <- read.branch_angle_undefined(undef1path, 3)
		
		coefs1 <- list()
		for (parammat in parammatlist) {
			
			coefs1 <- c(coefs1, list(calculate_coefs3(parammat)))
		}
		
		user1path <- file.path(dbpath, "branch_angle", "branch_angle_1_user.txt")
		
		user1file <- file(user1path, "a")
		for (i in seq_along(parammatlist)) catparams(parammatlist[[i]], coefs1[[i]], user1file)
		close(user1file)
	
		file.remove(undef1path)
	}

	undef2path <- file.path(dbpath, "branch_angle", "branch_angle_2_undefined.txt")
	
	if (file.exists(undef2path)) {
	
		parammatlist <- read.branch_angle_undefined(undef2path, 4)
		
		coefs2 <- list()
		for (parammat in parammatlist) {
			
			coefs2 <- c(coefs2, list(calculate_coefs4(parammat)))
		}
		
		user2path <- file.path(dbpath, "branch_angle", "branch_angle_2_user.txt")
		
		user2file <- file(user2path, "a")
		for (i in seq_along(parammatlist)) catparams(parammatlist[[i]], coefs2[[i]], user2file)
		close(user2file)
	
		file.remove(undef2path)
	}
}
