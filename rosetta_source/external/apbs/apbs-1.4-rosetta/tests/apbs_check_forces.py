#! /usr/bin/env python

"""
Checks computed forces from an apbs run
"""

import sys, re
from apbs_logger import Logger

error_tolerance = 1e-6

class PolarForce:
    """
    Exctracts and compares computations of polar forces
    """
    
    # A crazy regex pattern used to match label/value sets
    pattern=r'\s+(?P<label>[a-zA-Z]+)\s+(?P<x>[+-]?\d\.\d+E[+-]\d+)\s+(?P<y>[+-]?\d\.\d+E[+-]\d+)\s+(?P<z>[+-]?\d\.\d+E[+-]\d+)'

    def __init__( self, label, x, y, z ):
        """
        Constructs a polar force result from supplied values
        """
        self.label = label
        self.x = x
        self.y = y
        self.z = z
        
    def __init__( self, line ):
        """
        Extracts ploar force results from a file at a given line
        """
        m = re.search( self.pattern, line )
        self.label = m.group( 'label' )
        self.x = float( m.group( 'x' ) )
        self.y = float( m.group( 'y' ) )
        self.z = float( m.group( 'z' ) )
        
    def diff( self, other ):
        """
        Compares the value of two polar force field results
        """
        
        diff_dict = {}
        for d in ( 'x', 'y', 'z' ):
            diff_dict[ d ] = abs( getattr( self, d ) - getattr( other, d ) )
        return diff_dict    
        
    def __repr__( self ):
        return "PolarForce{ label:%s, x:%g, y:%g, z:%g}\n" % ( self.label, self.x, self.y, self.z )
        
    def short( self ):
        return self.label
        
        

class ApolarForce( PolarForce ):
    """
    Exctracts and compares computations of apolar forces
    """

    # A crazy regex pattern used to match label/value sets
    pattern=r'\s+(?P<label>[a-zA-Z]+)\s+(?P<atom>\w+)\s+(?P<x>[+-]?\d\.\d+E[+-]\d+)\s+(?P<y>[+-]?\d\.\d+E[+-]\d+)\s+(?P<z>[+-]?\d\.\d+E[+-]\d+)'

    def __init__( self, label, atom, x, y, z ):
        """
        Constructs an apolar force result from supplied values
        """
        super( ApolarForce, self ).__init__( self, x, y, z )
        self.label = label
        
    def __init__( self, line ):
        """
        Extracts aploar force results from a file at a given line
        """
        m = re.search( self.pattern, line )
        self.label = m.group( 'label' )
        self.atom = m.group( 'atom' )
        self.x = float( m.group( 'x' ) )
        self.y = float( m.group( 'y' ) )
        self.z = float( m.group( 'z' ) )
        
    def __rpr__( self ):
        return "ApolarForce{ label:%s, atom:%s, x:%g, y:%g, z:%g}\n" % ( self.label, self.atom, self.x, self.y, self.z )
        
    def short( self ):
        return '%s for %s' % ( self.label, self.atom )
        
        
        
        
def extract_forces( force_class, lines, start_pattern, ):
    """
    Extracts force results
    """
    force_dict = {}
    in_section = False
    in_forces = False
    start_line = -1
    end_line = -1
    for line_number, line_text in enumerate(lines):
        if not in_section:
            if line_text.startswith( start_pattern ):
                in_section = True
        if in_section and not in_forces:
            if re.search( force_class.pattern, line_text ):
                in_forces = True
                start_line = line_number
        if in_section and in_forces:
            if not re.search( force_class.pattern, line_text ):
                end_line = line_number
                break        
    return parse_forces( force_class, lines[ start_line : end_line ] )
        


def parse_forces( force_class, lines ):
    force_dict = {}
    for line in lines:
        force_item = force_class( line )
        force_dict[ force_item.label ] = force_item
    return force_dict
            
        

def compare_force_dicts( test_force_dict, true_force_dict ):
    """
    Compares force dictionaries
    """
    
    for force_key in test_force_dict.keys():
        test_force = test_force_dict[ force_key ]
        true_force = true_force_dict[ force_key ]
        diff_dict = test_force.diff( true_force )

        for ( diff_key, diff_value ) in diff_dict.items():
            test_value = getattr( test_force, diff_key )
            true_value = getattr( true_force, diff_key )
            
            if diff_value == 0.0:
                logger.message( '*** Comparison %s in %s PASSED ***' % ( test_force.short(), diff_key ) )
                logger.log( 'Comparison %s in %s PASSED (%g)' % ( test_force.short(), diff_key, test_value ) )
                
            elif diff_value < error_tolerance:
                logger.message( '*** Comparison %s in %s PASSED (with rounding error - see log)***' % ( test_force.short(), diff_key ) )
                logger.log( 'Comparison %s in %s PASSED within error (%g; expected %g)' %( test_force.short(), diff_key, test_value, true_value ) )
                
            else:
                logger.message( '*** Comparison %s in %s FAILED ***' % ( test_force.short(), diff_key ) )
                logger.message( '   APBS returned %g' % test_value )
                logger.message( '   Expected result is %g (difference of: %g)' % ( true_value, diff_value ) )
                logger.log( 'Comparison %s in %s FAILED (%g; expected %g)' %( test_force.short(), diff_key, test_value, true_value ) )
                
                

def check_forces( input_file, polar_file, apolar_file, logger ):
    
    logger.both( "Checking forces for input file %s" % input_file )
    
    f = None
    try:
        f = open( input_file, 'r' )
    except IOError:
        print >> sys.stderr, "Couldn't read from forces file %s" % input_file
    input_lines = f.readlines()
    
    test_polar_force_dict = extract_forces( PolarForce, input_lines, 'print force' )
    test_apolar_force_dict = extract_forces( ApolarForce, input_lines, 'print APOL force' )
        
    try:
        f = open( polar_file, 'r' )
    except IOError:
        print >> sys.stderr, "Couldn't read from forces file %s" % input_file
    input_lines = f.readlines()
    true_polar_force_dict = parse_forces( PolarForce, input_lines )
    
    try:
        f = open( apolar_file, 'r' )
    except IOError:
        print >> sys.stderr, "Couldn't read from forces file %s" % input_file
    input_lines = f.readlines()
    true_apolar_force_dict = parse_forces( ApolarForce, input_lines )
    
    logger.both( "Checking Polar Forces" )
    compare_force_dicts( test_polar_force_dict, true_polar_force_dict )
    
    logger.both( "Checking Apolar Forces" )
    compare_force_dicts( test_apolar_force_dict, true_apolar_force_dict )
    
    

def test():
    l = open( 'forces.log', 'w' )
    logger = Logger( sys.stderr, l )
    check_forces( 'apbs-forces.out', 'polarforces', 'apolarforces', logger )
        
        
        
if __name__ == '__main__':
    print >> sys.stderr, "The python source file %s is a module and not runnable" % sys.argv[ 0 ]
    sys.exit( 1 )
