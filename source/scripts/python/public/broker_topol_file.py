import sys
import argparse
import top_file

def process_command_line(argv):
    
    parser = argparse.ArgumentParser( description="".join(["This script is for generating topology files to specification.",
                                                            "Using the --strands option to define a list of residue ranges",
                                                            "that are to be used as pairs. A note of caution: put the arguments",
                                                            " of the script in quotes, so the shell doesn't get confused."] ) )

    parser.add_argument("--strands", required=True, type=str, metavar="xx-yy,zz-aa", 
                        help="".join(["The residue ranges (e.g. xx through yy) that are considered strands. The strand's position",
                                      " in the list (starting at zero) is then then used in the --topology option definition for pairings."]) ) 
    parser.add_argument("--topology", required=True, type=str, metavar="w-x:P,y-z:A",
                        help="The ways in which the given strands, defined with the --strands option, are fit together. For example, strand w is paired with strand x in a parallel fashion." )
    parser.add_argument("--len_shift", default=0, type=int, metavar="DELTA",
                        help="Allow strands to shift up to DELTA against each other." )

    args = parser.parse_args(argv[1:])

    strand_defs = []
    for pair in args.strands.split(','):
        try:
            (i_str, j_str) = pair.split('-')
            (i_str, j_str) = ( i_str.lstrip().rstrip(), j_str.rstrip().lstrip() )
            (i, j) = ( int( i_str ), int( j_str ) )
        except:
            raise "Couldn't parse the strand definition: ", pair
            print >> sys.stderr, "Couldn't parse the strand definition: ", pair
            sys.exit(1)
        
        print strand_defs
        strand_defs.append( ( i, j ) )

    strand_pairs = []
    for triplet in [s.rstrip().lstrip() for s in args.topology.split(",")]:
        try:
            (p1, substr) = triplet.split("-")
            (p2, para) = substr.split(":")
            (p1, p2) = (int( p1.lstrip().rstrip() ), int( p2.lstrip().rstrip() ))
            assert( para == "P" or para == "A" )
            strand_pairs.append( (p1, p2, para) )
        except:
            print >> sys.stderr, "Couldn't parse the strand pair definition: ", triplet
            sys.exit(1)

    args.strands = strand_defs
    args.topology = strand_pairs

    if( args.len_shift > 0 ):
        reduced = set()
        for ( s1, s2, para ) in strand_pairs:
            for s in [ s1, s2 ]:
                lindiv = all( [ s1 not in reduced and 
                                s2 not in reduced for ( s1, s2, para ) in strand_pairs
                                if s == s1 or s == s2] )
                if( lindiv ):
                    args.strands[s] = ( args.strands[s][0], 
                                         args.strands[s][1] - args.len_shift )
                    reduced.add( s )
    elif( args.len_shift < 0 ):
        sys.stderr.write( "--len_shift flag values must be >= 0." )
        sys.exit(1)

    return args

def main(argv=None):
    
    args = process_command_line(argv)

    tfile = top_file.TopologyFile( args.strands,
                                   [args.topology] )

    print tfile

    return 1

if __name__ == "__main__":
    status = main(sys.argv)
    sys.exit(status)
