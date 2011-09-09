# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.



# vector math

# This is a collection of functions for 3D vector algebra.

# Array math in R:
#  * Use arrays! (otherwise it can be slow)
#  * Read this (http://cran.r-project.org/doc/manuals/R-intro.html#Arrays-and-matrices)
#  * Multi-dimensional arrays are 1D arrays with an extra dimension attrubute see '?dim'
#  * R uses column-major notation (http://en.wikipedia.org/wiki/Row-major_order)

#input:  Two arrays of dimension c(n, d)
#output: An array of dimension c(n, 1) where the i'th element of the output is the
#        the dot product of the i'th rows of the input
vector_dotprod <- function(a, b) rowSums(a*b)

#input:  An array of dimensions c(n, d)
#output: An array of dimension c(n, d) where if the i'th row of the input is (x1,x2, ..., xd)
#        then the i'th of the output is (x1/norm, x2/norm, ..., xd/norm) where
#        norm = sqrt(x1^2 + x2^2 + ... + xd^2)
vector_normalize <- function(m) m / sqrt(rowSums(m^2))

#input:  Two arrays each of dimension c(n, 3)
#output: An array of dimension c(n,3)  where the i'th row of the output is
#        the vector cross product of the i'th rows of the input vectors
vector_crossprod <- function( a, b ) {
	result <- matrix( NA, nrow( a ), 3 )
	result[,1] <- a[,2] * b[,3] - a[,3] * b[,2]
	result[,2] <- a[,3] * b[,1] - a[,1] * b[,3]
	result[,3] <- a[,1] * b[,2] - a[,2] * b[,1]
	result
}

#input:  4 arrays each of dimension c(n, 3)
#output: An array of dimension c(n, 1) which is the dihedral angle about v2 -> v3
#        between v1->v2 to v3->v4 using the right hand rule.
vector_dihedral <- function(v1, v2, v3, v4) {
	w1 <- vector_normalize(v2-v1)
	w2 <- vector_normalize(v3-v2)
	w3 <- vector_normalize(v4-v3)
	x <- vector_dotprod(w1,w2) * vector_dotprod(w2,w3) - vector_dotprod(w1,w3)
	y <- vector_dotprod(w1, vector_crossprod(w2, w3))
	atan2(y,x)
}
