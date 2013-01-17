#! /usr/bin/env python

"""
Provides functions for verifying results from a test run
"""

import sys, re
from apbs_logger import Logger

error_tolerance = 1e-4
from math import log10, floor



def round_sigfigs( x, sigfigs ):
    """
    Rounds a number to a specified number of significant figures
    """
    return round( x, sigfigs - int( floor( log10( abs( x ) ) ) ) - 1 )



def check_results( computed_result, expected_result, input_file, logger, ocd ):
    """
    Compares computed results to an expected results within some margin of error
    """
    
    # OCD mode requires a match up to 12 significant figures
    if ocd:
        computed_result = round_sigfigs( computed_result, 12 )
        expected_result = round_sigfigs( expected_result, 12 )
   
    # Non-OCD mode only requires a match out to six significan figures
    else:
        computed_result = round_sigfigs( computed_result, 6 )
        expected_result = round_sigfigs( expected_result, 6 )
    
    # Compute the error in the calculation
    error = abs(computed_result - expected_result)/expected_result*100.0
    
    # An exact match after rounding to specifiec precision means the test passed
    if computed_result == expected_result:
        logger.message( "*** PASSED ***" )
        logger.log( "PASSED %.12e" % computed_result )
        
    # Otherwise, test that the error is below error tolerance
    elif error < error_tolerance:
        logger.message( "*** PASSED (with rounding error - see log) ***" )
        logger.log( "PASSED within error (%.12e; expected %.12e; %g%% error)" % ( computed_result, expected_result, error ) )
        
    # If neither is true, the test failed
    else:
        logger.message( "*** FAILED ***" )
        logger.message( "   APBS returned      %.12e" % computed_result )
        logger.message( "   Expected result is %.12e (%g%% error)" % ( expected_result, error ) )
        logger.log( "FAILED (%.12e; expected %.12e; %g%% error)" % ( computed_result, expected_result, error ) )
    
    
    
if __name__ == '__main__':
    print >> sys.stderr, "The python source file %s is a module and not runnable" % sys.argv[ 0 ]
    sys.exit( 1 )
