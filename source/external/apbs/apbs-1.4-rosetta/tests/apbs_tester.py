#! /usr/bin/env python

"""
Provides utility for testing apbs against examples and known results
"""

import sys, os, re, datetime, subprocess, operator
from optparse import OptionParser
from ConfigParser import ConfigParser

# The inputgen utility needs to be accessible, so we add its path
sys.path.insert( 0, "../tools/manip" )
from inputgen import splitInput

from apbs_check_forces import check_forces
from apbs_check_results import check_results
from apbs_logger import Logger

# Matches a floating point number such as -1.23456789E-20
float_pattern = r'([+-]?\d+\.\d+E[+-]\d+)'



def test_binary():
    """
    Ensures that the apbs binary is available
    """
    
    # Attempts to find apbs in the system path first
    try:
        binary = "apbs"
        subprocess.call( [binary, "--version"] )
        return binary
    except OSError:
        pass

    # Next, looks for the apbs binary in the apbs bin directory
    try:
        binary = os.path.abspath( "../bin/apbs" )
        subprocess.call( [binary, "--version"] )
        return binary
    except OSError:
        return ""



def process_serial( binary, input_file ):
    """
    Runs the apbs binary on a given input file
    """
    
    # First extract the name of the input file's base name
    base_name = input_file.split('.')[0]
    
    # The output file should have the same basename
    output_name = '%s.out' % base_name
    
    # Ensure that there are sufficient permissions to write to the output file
    output_file = open( output_name, 'w' )
    
    # Construct the system command and make the call
    command = [ binary, input_file ]
    subprocess.call( command, stdout = output_file )
    
    # Look for the results in the output file
    output_file = open( output_name, 'r' )
    output_text = output_file.read()
    output_pattern = r'Global net (?:ELEC|APOL) energy \= ' + float_pattern
    output_results = [ float( r ) for r in re.findall( output_pattern, output_text ) ]
    
    # Return all the matched results as a list of floating point numbers
    return output_results



def process_parallel( binary, input_file, procs, logger ):
    """
    Performs parallel apbs runs of the input file
    """

    logger.message("Splitting the input file into %d separate files using the inputgen utility" % procs )
    logger.message( "" )

    # Get the base name, and split the input file using inputgen's splitInput
    base_name = input_file.split('.')[0]
    splitInput( input_file )
    
    results = None
    for proc in range( procs ):
    
        # Process each paralle input file and capture the results from each
        proc_input_file = '%s-PE%d.in' % ( base_name, proc )
        proc_results = process_serial( binary, proc_input_file )

        # Log the results from each parallel run
        logger.message( 'Processor %d results:' % proc )
        for proc_result in proc_results:
            logger.message( '  %.12E' % proc_result )
        logger.message( '' )

        # Aggregate the results from each processor
        if results == None:
            results = proc_results
        else:
            results = [ r + p for ( r, p ) in zip( results, proc_results ) ]

    # Return the aggregated results from the parllel run
    return results



def run_test( binary, test_files, test_name, test_directory, logger, ocd ):
    """
    Runs a given test from the test cases file
    """
    
    logger.log( '=' * 80 )
    logger.log( "Test Timestamp: %s" % str( datetime.datetime.now() ) )
    logger.log( "Test Name:      %s" % test_name )
    logger.log( "Test Directory: %s" % test_directory )

    # The net time is initially zero
    net_time = datetime.timedelta( 0 )

    # Change the current working directory to the test directory
    os.chdir( test_directory )

    for ( base_name, expected_results ) in test_files:

        # Get the name of the input file from the base name
        input_file = '%s.in' % base_name

        # If the expected results is 'forces', do a forces test on the input
        if expected_results == 'forces':
            logger.message( '-' * 80 )
            logger.message( 'Testing forces from %s' % input_file )
            logger.message( '' )
            logger.log( 'Testing forces from %s' % input_file )
            start_time = datetime.datetime.now()
            check_forces( input_file, 'polarforces', 'apolarforces', logger )
            
        else:
            logger.message( '-' * 80 )
            logger.message( 'Testing input file %s' % input_file )
            logger.message( '' )
            logger.log( 'Testing %s' % input_file )
            
            # Record the start time before the test runs
            start_time = datetime.datetime.now()

            computed_results = None
            
            # Determine if this is a parallel run
            match = re.search( r'\s*pdime((\s+\d+)+)', open( input_file, 'r' ).read() )
            
            # If it is parallel, get the number of procs and do a parallel run
            if match:
                procs = reduce( operator.mul, [ int(p) for p in match.group( 1 ).split() ] )
                computed_results = process_parallel( binary, input_file, procs, logger )
            
            # Otherwise, just do a serial run
            else:
                computed_results = process_serial( binary, input_file )
                
            # Split the expected results into a list of text values
            expected_results = expected_results.split()
            for i in range( len( expected_results ) ):
            
                # If the expected result is a star, it means ignore that result
                if expected_results[i] == '*':
                    continue
                    
                # Compare the expected to computed results
                computed_result = computed_results[i]
                expected_result = float( expected_results[i] )
                logger.message( "Testing computed result %.12E against expected result %12E" % ( computed_result, expected_result ) )
                check_results( computed_result, expected_result, input_file, logger, ocd )

        # Record the end time after the test
        end_time = datetime.datetime.now()
        elapsed_time = end_time - start_time
        net_time += elapsed_time
        stopwatch = elapsed_time.seconds + elapsed_time.microseconds / 1e6
        
        # Log the elapsed time for this test
        logger.message( "Elapsed time: %f seconds" % stopwatch )
        logger.message( '-' * 80 )

    stopwatch = net_time.seconds + net_time.microseconds / 1e6
    
    # Log the elapsed time for all tests that were run
    logger.message( "Total elapsed time: %f seconds" % stopwatch )
    logger.message( "Test results have been logged" )
    logger.message( '-' * 80 )
    logger.log( "Time:           %d seconds" % stopwatch )

    os.chdir( '../../tests' )



def main():
    """
    Parse command line options and run tests with given options
    """
    
    parser = OptionParser()
    
    #Describes the available options.
    parser.add_option(
        '-c', '--test_config', dest='test_config',
        type='string', default='test_cases.cfg',
        help="Set the test configuration file to FILE", metavar="FILE"
        )
        
    parser.add_option(
        '-t', '--target_test', dest='target_test',
        type='string', action='append', default=[],
        help="Set the test to run to TEST", metavar="TEST"
        )
        
    parser.add_option(
        '-o', '--ocd', action='store_true', dest='ocd',
        help="Run APBS in OCD mode"
        )
        
    parser.add_option(
        '-l', '--log_file', dest='log_file', type='string', default='test.log',
        help="Save the test log to FILE.", metavar="FILE"
        )
    
    # Parse the command line and extract option values
    ( options, args ) = parser.parse_args()

    # Messages will go to stdout, log messages will go to the supplied log file
    message_fd = sys.stdout
    logfile_fd = None
    
    # Verify that the log file is writable
    try:
        logfile_fd = open( options.log_file, 'w' )
    except IOError as err:
        parser.error( "Could't open log_file %s: %s" % ( options.log_file, err.strerror ) )
        
    # Set up the logger with the message and log file descriptor
    logger = Logger( message_fd, logfile_fd )

    # Read the test cases file
    config = ConfigParser()
    config.read( options.test_config )

    # Make sure that the apbs binary can be found
    binary = test_binary()
    if binary == '':
        parser.error( "Coulnd't detect an apbs binary in the path or local bin directory" )

    # Get the names of all the test sections to run.
    test_sections = []
    if 'all' in options.target_test or options.target_test == []:
        test_sections = config.sections()
        print "Testing all sections"
    else:
        test_sections = options.target_test

    print "The following sections will be tested: " + ', '.join( test_sections )
    print '=' * 80

    # Run each test that has been requested
    for test_name in test_sections:
    
        print "Running tests for " + test_name + " section"
        
        # Verify that the test is described in the test cases file
        if test_name not in config.sections():
            print "  " + test_name + " section not found in " + options.test_config
            print "  skipping..."
            continue
            
        # As each test runs, remove the value so it can't be run twice
        test_directory = config.get( test_name, 'input_dir' )
        config.remove_option( test_name, 'input_dir' )
        
        # Run the test!
        run_test( binary, config.items( test_name ), test_name, test_directory, logger, options.ocd )



# If this file is executed as a script, call the main function
if __name__ == '__main__':
    main()







