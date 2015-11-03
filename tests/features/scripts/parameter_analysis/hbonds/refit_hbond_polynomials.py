import blargs

import Tkinter, tkFileDialog, tkSimpleDialog
import numpy
import inspect, warnings
import math
import sys

from cvxopt import lapack, matrix

def fit_polynomial(X, Y, d):
    m = len(X)
    assert(len(Y) == m)

    A = matrix( [[X**k] for k in xrange(d)])

   # make a deep copy of Y
    xls = +Y

    # general least-squares: minimize ||A*x - y||_2
    lapack.gels(+A,xls)
    xls = xls[:d]

    return xls

def fit_polynomial_one_constraint(X, Y, d, x0, y0):
    m = len(X)
    assert(len(Y) == m)

    A = matrix( [[X**k] for k in xrange(d)])

    G = matrix(0.0, (2,d))
    G[0, range(d)] = matrix([[x0**k] for k in xrange(d)])
    G[1, range(1,d)] = matrix([[k*x0**(k-1)] for k in xrange(1,d)])

    # LS fit
    #
    #     minimize    (1/2) * || A*x - Y ||_2^2
    #     subject to  G*x = h
    #
    # Solve as a linear equation
    #
    #     [ A'*A  G' ] [ x ]   [ A'*Y ]
    #     [ G     0  ] [ Y ] = [ 0    ].

    K = matrix(0.0, (d+2,d+2))
    K[:d,:d] = A.T * A
    K[d:,:d] = G
    xls = matrix(0.0, (d+2,1))
    xls[:d] = A.T * Y
    xls[d] = y0
    lapack.sysv(K, xls)
    xls = xls[:d]
    return xls

# x0 is a list of xvals for zero-derivative constraints.
# The fit polynomial must have a derivative of zero
# at these points.  y0 is a list of yvals for x0
# x1 is a list of xvals for constraints where the
# derivative is unconstrained.  y1 is a list of yvals
# for x1.
# x0 must have the same length as y0.
# x1 must have the same length as y1
def fit_polynomial_multiple_constraints(X, Y, d, x0, y0, x1 = None, y1 = None):
    m = len(X)
    assert(len(Y) == m)

    A = matrix( [[X**k] for k in xrange(d)])

    ncst0 = len(x0)
    assert( ncst0 == len(y0))

    ncst1 = 0
    if x1 : ncst1 = len(x1)

    G = matrix(0.0, (2*ncst0+ncst1,d))
    for ii in xrange(ncst0) :
        G[2*ii+0, range(d)]   = matrix([[x0[ii]**k] for k in xrange(d)])
        G[2*ii+1, range(1,d)] = matrix([[k*x0[ii]**(k-1)] for k in xrange(1,d)])
    if x1 :
        for ii in xrange(ncst1) :
            G[ii+2*ncst0,range(d)] = matrix( [[x1[ii]**k] for k in xrange(d)] )

    # LS fit
    #
    #     minimize    (1/2) * || A*x - Y ||_2^2
    #     subject to  G*x = h
    #
    # Solve as a linear equation
    #
    #     [ A'*A  G' ] [ x ]   [ A'*Y ]
    #     [ G     0  ] [ Y ] = [ 0    ].

    K = matrix(0.0, (d+2*ncst0+ncst1,d+2*ncst0+ncst1))
    K[:d,:d] = A.T * A
    K[d:,:d] = G
    xls = matrix(0.0, (d+2*ncst0+ncst1,1))
    xls[:d] = A.T * Y
    for ii in xrange( ncst0 ) :
        xls[d+ii*2] = y0[ii]
    if x1 :
        for ii in xrange( ncst1 ) :
            xls[d+2*ncst0+ii] = y1[ii]
    lapack.sysv(K, xls)
    xls = xls[:d]

    coefs = numpy.array(xls.T)[0]
    coefs = coefs[::-1]
    return coefs


def lineno():
    """Returns the current line number in our program."""
    return inspect.currentframe().f_back.f_lineno

def find_real_polynomial_roots_near_target( coeffs, target = 1.8 ) :
    poly = numpy.poly1d( coeffs )
    lowroot = -1; highroot = -1;
    for root in poly.r :
        #print cols[1], root
        if numpy.abs( numpy.imag(root)) > 1e-2 : continue
        realroot = numpy.real(root)
        if lowroot < 0 :
            if realroot < target :
                lowroot = realroot
        elif realroot < target and target - lowroot > target - realroot :
            lowroot = realroot
        if highroot < 0 :
            if realroot > target :
                highroot = realroot
        elif realroot > target and highroot - target > realroot - target :
            highroot = realroot
    return lowroot, highroot


def find_polynomial_extrema( coeffs, deg, target = 1.8 ) :
    """ where does the polynomial hit the value of +1.0 ? """
    offset = 1.0
    offset_coeffs = list( coeffs )
    offset_coeffs[deg-1] -= offset
    return find_real_polynomial_roots_near_target( offset_coeffs, target )


def refit_polynomial_from_file( fname, degree, hb_weight, verbose ) :
    lines = open( fname ).readlines()
    last_x = -1234;
    xs = []; energies = [];
    x0s = []; y0s = [];
    for line in lines :
        if line[0] == '#' : continue;
        if line == "\n" : continue;
        cols = line.split()
        if cols[0] == "ENERGY_VALUES" :
            assert( len(cols) == 4 )
            xval = float(cols[1])
            #print xval, last_x
            assert( last_x == -1234 or last_x < xval )
            last_x = xval
            xs.append( xval )
            # cols[2] == original hydrogen bond polynomial value
            # cols[3] == new energy term value, e.g. hack_elec
            hb_energy = float(cols[2])
            # for numerical reasons, the hbond energy can sometimes end up slightly negative even though it should
            # be 0; trim away these values by just rounding up to zero
            if hb_energy < 0 and hb_energy > -1e-5 :
                hb_energy = 0.0
            alt_energy = float(cols[3])
            energies.append( (hb_energy, alt_energy) )
        elif cols[0] == "ZERO_DERIVATIVE_CONSTRAINT" :
            assert( len(cols) == 3 )
            x0s.append( float(cols[1]) )
            y0s.append( float(cols[2]) )
        elif cols[0] == "OUTPUT_POLY_NANE" :
            output_name = cols[1]
        else :
            print "ERROR: expected to read 'OUTPUT_POLY_NANE', 'ZERO_DERIVATIVE_CONSTRAINT', or 'ENERGY_VALUES' but read", cols[0]
            sys.exit(1)

    # construct the values for the new polynomial that "absorb" the values of the new term
    ys_start = []; min_y = -1234; min_xval = -1234;
    for x, e in zip( xs, energies ) :
        orig_hbpoly = e[0]
        new_evalue = e[1]
        diff = orig_hbpoly - new_evalue
        ys_start.append( diff )
        if x > 1.1 and x < 2.2 and ( min_y == -1234 or diff < min_y ) :
            min_y = diff
            min_xval = x

    # calculate interpolation values for the hb polynomial.  Shift them up to have a minima at zero, scale them down by 1 / weight,
    # and then shift them back down so their minima is at 0.5
    ys = []
    for y in ys_start :
        ys.append( ( y - min_y ) / hb_weight - 0.5 )

    # xandys = []
    # for x,y in zip(xs,ys) : xandys.append("%f,%f" % (x,y) )
    # print ",".join(xandys)

    xs2, ys2 = tune_control_points( xs, ys, energies, min_xval )
        

    # now fit the polynomial
    x0s.append( min_xval ); y0s.append( -0.5 );

    if verbose :
        print "x0:", min_xval, "y0:", -0.5

    poly = fit_polynomial_multiple_constraints( matrix(xs2), matrix(ys2), degree, x0s, y0s )
    #for coeff in poly : print coeff,
    #print

    print_polynomial_in_HBPoly_csv_format( output_name, poly )

def print_polynomial_in_HBPoly_csv_format( output_name, poly, poly_index=0 ) :
    # ok -- good -- now we should do what?
    # write the output polynomial in the .csv format used by Rosetta
    # a. computing the x-intercepts <-- these really aren't important
    # b. computing the final polynomial values

    low_xint, hi_xint = find_real_polynomial_roots_near_target( poly )
    low_extr, hi_extr = find_polynomial_extrema( poly, poly.size )

    # output in this format:
    # 1,poly_AHdist_g_2,,hbgd_AHdist,1.01702445149,2.96903286337,1.0,1.0,1.473361393366,37.009437349874,6,0.017766989922,-0.797066393700,5.512449246897,-13.205746417374,10.902052999654,-1.393833179673,,,,,

    vals = []
    vals.append( "%d" % poly_index )
    vals.append( output_name )
    vals.append( "" )
    vals.append( "hbgd_AHdist")
    vals.append( "%10.8f,%10.8f" % ( low_extr, hi_extr ) )
    vals.append( "1.0,1.0" )
    vals.append( "%10.8f,%10.8f" % ( low_xint, hi_xint ) )
    vals.append( "%d" % poly.size )
    for coeff in poly :
        vals.append( "%10.8f" % coeff )
    for i in xrange( 11 - poly.size ) : vals.append( "" )

    print ",".join( vals )

def fit_parabola( xs, ys ) :
    x = list(xs)
    M = numpy.matrix( [[ x[0]*x[0], x[0], 1 ], [x[1]*x[1], x[1], 1 ], [ x[2]*x[2], x[2], 1 ]] )
    b = numpy.matrix( ys )
    # print M
    # print b.transpose()
    return numpy.linalg.inv( M ) * b.transpose();

def eval_parabola( coeffs, x ) :
    val = coeffs[0]*x*x + coeffs[1]*x + coeffs[2]
    #print "val", val.item(0,0)
    return val.item(0,0)

def solve_parabola( coeffs, y ) :
    # solve the quadratic equation:
    a = coeffs[0]
    b = coeffs[1]
    c = coeffs[2] - y
    discriminant = b**2 - 4*a*c
    if ( discriminant < 0 ) :
        return [500,-500]
    return ( (-b + math.sqrt( discriminant )) / (2*a),  (-b - math.sqrt( discriminant )) / (2*a) )

def tune_control_points( xs, ys, energies, min_xval ) :
    # trim down to the monotonically increasing range surrounding min_xval
    xs2a,ys2a = [], []
    xs2b,ys2b = [], []
    lasty = -1234
    for x,y in zip(xs,ys) :
        if x >= min_xval :
            if lasty < y :
                xs2a.append(x)
                ys2a.append(y)
            else :
                break
            lasty = y
    lasty = -1234
    for x, y in reversed(zip(xs,ys)) :
        if x < min_xval :
            if lasty < y :
                xs2b.append(x)
                ys2b.append(y)
            else :
                break
            lasty = y
    xs, ys = [],[]
    xs2b.reverse(); ys2b.reverse()
    xs = xs2b + xs2a; ys = ys2b + ys2a
    # print "xs",xs
    # print "ys",ys
    # print


    # lets throw out the xs and ys where ys > 1.0 or that are outside of the range near the min_xval where it is below 1 without
    # crossing above 1
    xs2 = []; ys2 = [];
    hits_1p1_low = False;
    hits_1p1_hi  = False;
    xlow = -1234; xhi = -1234;
    for x,y in zip(xs,ys):
        if x < min_xval and y >= 1.1 and not hits_1p1_low :
            if verbose : print "hits 1p1 low"
            hits_1p1_low = True
            xlow = x
        if x > min_xval and y >= 1.1 and not hits_1p1_hi  :
            if verbose : print "hits 1p1 hi"
            hits_1p1_hi  = True 
            xhi  = x
    if hits_1p1_hi and hits_1p1_low :
        if verbose : print "hits both 1p1 low and hi!"
        for x,y in zip(xs,ys) :
            if x >= xlow and x <= xhi :
                xs2.append(x)
                ys2.append(y)
        return xs2, ys2;

    # otherwise, we need to fudge some control points 

    #print "min_xval", min_xval
    for x,y in zip(xs,ys) :
        if x == min_xval :
            xs2.append(x)
            ys2.append(y)
            break

    # take the good parts of the curve, if any
    if hits_1p1_hi :
        for x,y in zip(xs,ys) :
            if x > min_xval and x <= xhi :
                xs2.append(x)
                ys2.append(y)
    if hits_1p1_low :
        for x, y in zip(xs,ys) :
            if x < min_xval and x >= xlow :
                xs2.append(x)
                ys2.append(y)
    if not hits_1p1_hi :
        # # find the point, if any, at which y stops monotonically increasing
        # xmax = -1234
        # ylast = -1234
        # for x,y in zip(xs,ys) :
        #     if x == min_xval :
        #         ylast = y
        #     elif x > min_xval :
        #         if y < ylast :
        #             xmax = x
        #             break
        #         ylast = y
        # #print "extending hi: xmax", xmax
        # if xmax == -1234 :
        #     # ok -- so, we didn't find a local maximum; this is troubling.
        #     pass
        # else :
        # let's trim backwards from the local maximum and see if we can't
        # draw a line up to +1.1

        xmax = xs[-1] # take the last element as the x maximum

        trim = 0.2
        done = False
        took_halfway_point = False
        while ( not done ) :
            #print "trim", trim
            if xmax - trim < min_xval :
                trim = ( xmax - min_xval ) / 2;
                took_halfway_point = True
            xs3, ys3 = [], []
            for x,y in zip( xs,ys ) :
                if x > min_xval and x <= xmax - trim :
                    #print "adding to xs3/ys3:",x, y
                    xs3.append(x)
                    ys3.append(y)
            # xi,yi     = xs3[-2],ys3[-2]
            # xip1,yip1 = xs3[-1],ys3[-1]
            # slope = ( yip1 - yi ) / ( xip1 - xi )
            coeffs = fit_parabola( xs3[-3:], ys3[-3:] )
            x_of_1p1 = solve_parabola(coeffs,1.1)[0]
            if x_of_1p1 > 3.0 : # 3.0 = MAX_R, the maximum hbond distance considered
                if took_halfway_point :
                    # OH SHIT!
                    # let's punt -- we'll just thow in some x,y pairs 
                    xs3.append( 2.9 ); ys3.append( 1.1 );
                    xs3.append( 2.8 ); ys3.append( 1.0 );
                    break;
                trim += 0.1 # let's try again trimming more of the original points
                continue
            else :
                done = True
        if done :
            x = xmax - trim + 0.1
            while x < x_of_1p1 :
                xs3.append( x ); ys3.append( eval_parabola(coeffs,x) );
                x += 0.1
            # one more for good measure
            xs3.append( x ); ys3.append( eval_parabola(coeffs,x) );
        xs2.extend(xs3); ys2.extend(ys3);

    if not hits_1p1_low :
        # # find the point, if any, at which y stops monotonically increasing
        # xmax = -1234
        # ylast = -1234
        # xsrev,ysrev = list(xs), list(ys); xsrev.reverse(); ysrev.reverse();
        # for x,y in zip(xsrev,ysrev) :
        #     if ylast != -1234:
        #         print "backwards", x, y, ylast
        #     if x == min_xval :
        #         ylast = y
        #     elif x < min_xval :
        #         if y < ylast :
        #             xmax = x
        #             break
        #         ylast = y
        #print "extend low: xmax", xmax
        xmax = xs[0]
        if xmax <= 0 :
            # ok -- so, we didn't find a local maximum -- take everything!
            for x,y in zip(xsrev,ysrev) :
                if x < min_xval :
                    xs2.append(x); ys2.append(y)
        else :
            # let's trim backwards from the local maximum and see if we can't
            # draw a parabola up to +1.1
            trim = 0.2
            done = False
            took_halfway_point = False
            while ( not done ) :
                if xmax + trim > min_xval :
                    trim = ( min_xval - xmax ) / 2;
                    took_halfway_point = True
                #print trim
                xs3, ys3 = [], []
                for x,y in reversed(zip( xs,ys )) :
                    if x < min_xval and x >= xmax + trim :
                        xs3.append(x)
                        ys3.append(y)
                coeffs = fit_parabola( xs3[-3:],ys3[-3:] )
                x_of_1p1 = solve_parabola(coeffs,1.1)[1]
                if coeffs[0] < 0 :
                    if took_halfway_point :
                        # let's punt!
                        xs3.append(0)
                        xs3.append(1.2)
                    trim += 0.1
                    continue
                else :
                    done = True

            x = xmax + trim - 0.1
            while x > 0 :
                nexty = eval_parabola(coeffs,x)
                if nexty > ys3[-1] :
                    xs3.append( x ); ys3.append( nexty );
                else :
                    break
                x -= 0.1
            xs2.extend(xs3); ys2.extend(ys3);

    if verbose :
        print "fitpoints_x = [",
        for x in xs2 : print x,
        print "];"
        print "fitpoints_y = [",
        for y in ys2 : print y,
        print "];"
        print
        print "controlpoints format: "
        ctrpts = []
        for x,y in zip(xs2,ys2) :
            ctrpts.append("%f"%x)
            ctrpts.append("%f"%y)
        print ",".join(ctrpts)
    

    #print xs2, ys2
    return xs2, ys2;

    # crossed_above_1p1 = False
    # for x,y in zip(xs,ys) :
    #     if y <= 1.1 and x >= min_xval and not crossed_above_1p1 :
    #         xs2.append(x)
    #         ys2.append(y)
    #     if y > 1.1 and x >= min_xval : crossed_above_1p1 = True
    # 
    # # try again
    # if not crossed_above_1p1 :
    #     old_max_x = -1234
    #     last_x_val = -1234
    #     last_e_val = -1234
    #     final_e_val = -1234
    #     final_y_val = -1234
    #     for x,e in zip(xs,energies) :
    #         if x >= min_xval and e[0] == 0.0 :
    #             old_max_x = x
    #             #print "old max x:", old_max_x
    #             break;
    #     if old_max_x == -1234 :
    #         print "Error: could not find xval at which the hbond score goes to zero.  Extend hbond energy reporting range"
    #         sys.exit(1)
    # 
    #     # step back two tenths of an angstrom and take the slope from there
    #     for x,e,y in zip(xs,energies,ys) :
    #         if x >= old_max_x - 0.2 :
    #             final_e_val = e[0]
    #             final_x_val = x
    #             final_y_val = y
    #             break
    #         last_e_val = e[0]
    #         last_x_val = x
    # 
    #     xs2 = []; ys2 = [];
    #     if verbose : print "final_e_val",final_e_val,"last_e_val",last_e_val,"final_x_val",final_x_val,"last_x_val",last_x_val
    #     slope = ( final_e_val - last_e_val ) / ( final_x_val - last_x_val )
    #     if verbose : print "slope",slope,"final_y_val",final_y_val
    #     # crossed_above_1p1 = False
    #     for x,y in zip(xs,ys) :
    #         newy = slope * ( x - final_x_val ) + final_y_val # point slope form
    #         #print x, y, newy
    #         if x >= min_xval and x <= final_x_val and not crossed_above_1p1 :
    #             xs2.append(x)
    #             ys2.append(y)
    #         elif x >= min_xval and x > final_x_val :
    #             if newy < 1.1 :
    #                 xs2.append(x)
    #                 ys2.append(newy)
    #             else :
    #                 crossed_above_1p1 = True
    # 
    # crossed_above_1p1 = False
    # save_xs2 = list(xs2); save_ys2 = list(ys2); # keep a deep copy in case we need to take a second pass sweep
    # 
    # revlist = zip(xs,ys,energies)
    # revlist.reverse()
    # xs2.reverse(); ys2.reverse();
    # for x,y,e in revlist :
    #     #print "revlist",x,y,e[0],e[1]
    #     if e[0] <= 1.1 and x < min_xval and not crossed_above_1p1 :
    #         #print "appending",x,y
    #         xs2.append(x)
    #         ys2.append(y)
    #     if e[0] > 1.1 and x < min_xval : crossed_above_1p1 = True
    # xs2.reverse(); ys2.reverse();
    # 
    # # try again
    # if not crossed_above_1p1 :
    #     old_max_x = -1234
    #     last_x_val = -1234
    #     last_e_val = -1234
    #     final_e_val = -1234
    #     final_y_val = -1234
    #     for x,y,e in revlist :
    #         if x <= min_xval and e[0] == 0.0 :
    #             old_max_x = x
    #             #print "old max x:", old_max_x
    #             break;
    #     if old_max_x != -1234 :
    # 
    #         xs2 = list(save_xs2); ys2 = list(save_ys2);
    #         xs2.reverse(); ys2.reverse();
    #         # step back two tenths of an angstrom and take the slope from there
    #         for x,y,e in revlist :
    #             if x <= old_max_x + 0.2 :
    #                 final_e_val = e[0]
    #                 final_x_val = x
    #                 final_y_val = y
    #                 break
    #             last_e_val = e[0]
    #             last_x_val = x
    #         
    #         if verbose : print "final_e_val",final_e_val,"last_e_val",last_e_val,"final_x_val",final_x_val,"last_x_val",last_x_val
    #         slope = ( final_e_val - last_e_val ) / ( final_x_val - last_x_val )
    #         if verbose : print "slope",slope,"final_y_val",final_y_val
    #         for x,y,e in revlist :
    #             newy = slope * ( x - final_x_val ) + final_y_val # point slope form
    #             #print x, y, newy
    #             if x <= min_xval and x >= final_x_val and not crossed_above_1p1 :
    #                 #print "appending original value:",x,y
    #                 xs2.append(x)
    #                 ys2.append(y)
    #             elif x <= min_xval and x < final_x_val :
    #                 if newy < 1.1 :
    #                     #print "appending linearized value:",x,y
    #                     xs2.append(x)
    #                     ys2.append(newy)
    #                 else :
    #                     crossed_above_1p1 = True
    # 
    # xs2.reverse(); ys2.reverse();
    # 
    # if ( verbose ) :
    #     print "fitpoints_x = [",
    #     for x in xs2 : print x,
    #     print "];"
    #     print "fitpoints_y = [",
    #     for y in ys2 : print y,
    #     print "];"
    #     print
    #     print "controlpoints format: "
    #     ctrptstr = ""
    #     for x,y in zip(xs2,ys2) : ctrptstr += ("%f,%f," % (x, y))
    #     print ctrptstr
    # 
    # return xs2,ys2

if __name__ == "__main__" :
    with blargs.Parser(locals()) as p :
        p.str( "polynomial_file" ).shorthand( "f" )
        p.str( "polynomial_list" ).shorthand( "l" )
        p.int( "degree" ).shorthand( "d" )
        p.float( "weight" ).shorthand( "w" ).default( 1.17 )
        p.flag("verbose").shorthand("v")

    fnames = []
    if polynomial_file : fnames.append( polynomial_file )
    if polynomial_list : 
        [ fnames.append( x.strip() ) for x in open( polynomial_list  ).readlines() ]

    if len(fnames) == 0 :
        print "Error: No input files given!"
        sys.exit(1)

    for fname in fnames :
        poly = refit_polynomial_from_file( fname, degree, weight, verbose )
