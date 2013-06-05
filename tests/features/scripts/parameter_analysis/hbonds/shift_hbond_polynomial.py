import refit_hbond_polynomials as rhp
import numpy
import blargs

from cvxopt import lapack, matrix

def drange( start, stop, step ) :
    r = start
    while r < stop :
        yield r
        r += step

def read_HBpoly_csv_file( fname ) :
    polys = []
    flines = open( fname ).readlines()
    for line in flines :
        cols = line.split(",")
        deg = int( cols[10] )
        name = cols[1]
        min = float( cols[4] )
        max = float( cols[5] )
        coeffs = []
        for col in xrange(11,11+deg) :
            coeffs.append( float(cols[col]) )
        polys.append( (name,coeffs,min,max) )
    return polys

def closest_to_target( vals, target ) :
    besti = -1
    closest_dist = -1
    for i in xrange(len(vals)) :
        if besti == -1 or abs(vals[i] - target) < closest_dist :
            closest_dist = abs(vals[i]-target)
            besti = i
    return besti

def shift_polynomial( name, coeffs, min, max, shift, count ) :
    xs = [ x for x in drange( min, max, 0.1 ) ]
    newxs = matrix([ x + shift for x in xs ])
    ys = matrix( numpy.polyval( coeffs, xs ) )

    deriv_coeffs = numpy.polyder( coeffs )
    minimum = rhp.find_real_polynomial_roots_near_target( deriv_coeffs )

    print "minimum: ",minimum[ closest_to_target(minimum,1) ],minimum
    x0 = [ minimum[ closest_to_target(minimum,1) ] + shift ]
    y0 = [ -0.5 ]

    print "xmin:", min, "xmax:", max
    x1 = [ min+shift, max+shift ]
    y1 = [ 1.0, 1.0 ]

    newpoly = rhp.fit_polynomial_multiple_constraints( newxs, ys, len(coeffs), x0, y0, x1, y1 )
    rhp.print_polynomial_in_HBPoly_csv_format( name, newpoly, count )

if __name__ == "__main__" :
    with blargs.Parser(locals()) as p :
        p.str("polynomial_file").shorthand("p")
        p.float("shift_value").shorthand("s")
    polys = read_HBpoly_csv_file( polynomial_file )

    count = 0
    for name,coeffs,min,max in polys :
        count += 1
        shift_polynomial( name,coeffs, min, max, shift_value, count )
    
