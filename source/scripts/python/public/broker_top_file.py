import itertools

class StrandPairing:
    def __init__( self, s1, s2, parallel, pleat, weight=20 ):
        self.weight = weight
        
        if( parallel == "A" or parallel == "P" ):
            self.parallel = parallel
        else:
            raise ValueError( "Parallelism must be 'A' or 'P'." )

        s1 = ( min(s1), max(s1) )
        s2 = ( min(s2), max(s2) )

        self.s1 = range( min(s1), max(s1)+1 )
        self.s2 = range( min(s2), max(s2)+1 )

        self.pleat = pleat

    def copy( self ):
        return StrandPairing( self.s1, self.s2, self.parallel, self.pleat, self.weight )

    def reverse_parallel( self ):
        if( self.parallel == "A" ):
            self.parallel = "P"
        else:
            self.parallel = "A"

    def reverse_pleating( self ):
        if( self.pleat[1] == 1 ):
            self.pleat = (self.pleat[0], 2)
        else:
            self.pleat = (self.pleat[0], 1)

    def __str__( self ):
        s1 = self.s1
        s2 = self.s2

        pleating = self.pleating_list()

        if( self.parallel == "P" ):
            pair1 = str(min(s1))+"-"+str(min(s2))
            pair2 = str(max(s1))+"-"+str(max(s2))
            register = min(s2) - min(s1)
        elif( self.parallel == "A" ):
            pair1 = str(min(s1))+"-"+str(max(s2))
            pair2 = str(max(s1))+"-"+str(min(s2))
            register = min(s1) + max(s2)

        rendered = " ".join( [ "PAIRSTAT_ENTRY:", "{0:.2f}".format(self.weight), "  .", self.parallel,
                               pair1, "to", pair2, "reg:", str(register), "pleating:",
                               " ".join([str(p) for p in pleating] ) ] )
        return rendered

    def pleating_list( self ):
        (residue, pleat_value ) = self.pleat
        s1 = self.s1
        
        if( residue in s1 ):
            assert( residue >= min(s1) )
            init_pleat = ( ( pleat_value - 1 ) + ( ( residue - min(s1) ) % 2 ) ) % 2
            pleating = [ ( ( i + init_pleat ) % 2 ) + 1 for i in range( len( s1 ) ) ]
            
        else:
            raise IndexError( " ".join(["The residue", str(residue), 
                                        " is not in the reference strand (",
                                        str(min(s1)), "-", str(max(s1)), ")."]))
        
        return pleating

class StrandTopology:
    def __init__(self, index=None):
        self.strand_pairings = []
        self.index=index
    def set_index( self, index ):
        self.index = index
    def add_pairing( self, new_pairing ):
        self.strand_pairings.append( new_pairing )
    def __str__( self ):
        pairings = sorted( self.strand_pairings, key=lambda x: x.s1[0] )

        weight_sum = sum( [ strand.weight for strand in pairings ] )
        n_strands = len( pairings )
        
        header_str = " ".join([ "STRAND_TOPOLOGY", str(n_strands).rjust(2),
                                str(round( weight_sum, 1 )).rjust(4),
                                str(self.index).zfill(4) ])

        line_list = [header_str]
        line_list.extend([ str(pair) for pair in pairings ])
        
        return "\n".join( line_list )
            
class TopologyFile:
    
    def __init__(self, strands, top_specs):
        topologies = []

        for spec in top_specs:
            topologies.append( [ (strands[i], strands[j], p ) for (i,j,p) in spec ] )

        self.topologies = topologies

    def __str__(self):
        
        topols = []

        for topol_spec in self.topologies:
            pair_sets = []
            for pair_spec in topol_spec:
                region_pairing = pair_spec[0], pair_spec[1]
                para_in = pair_spec[2]

                ( small_strand, large_strand ) = sorted( region_pairing, key=lambda x: x[1]-x[0] )
                small_len = small_strand[1] - small_strand[0]
                large_len = large_strand[1] - large_strand[0]
                len_diff = large_len - small_len
                
                shortened_large_strand = (large_strand[0], large_strand[1]-len_diff)
                
                strand_possibilities = []
                for delta in range( 0, len_diff+1 ):
                    tmp_large_strand = [ i+delta for i in shortened_large_strand ]
                    if( small_strand == pair_spec[0] ):
                        (strand1, strand2) = small_strand, tmp_large_strand
                    else:
                        (strand1, strand2) = tmp_large_strand, small_strand

                    pleat = ( delta % 2 ) + 1

                    strand = StrandPairing( s1=strand1, s2=strand2, 
                                            parallel=para_in, pleat=(strand1[0],pleat) )
                    strand_possibilities.append( strand  )
                    
                    strand = strand.copy()
                    strand.reverse_pleating()

                    strand_possibilities.append( strand )
                
                pair_sets.append( strand_possibilities )

            for pair_combo in itertools.product(*pair_sets):
                topo = StrandTopology()
                for strandpair in pair_combo:
                    topo.add_pairing( strandpair )
                topols.append( topo )
                
        index = 1
        for topol in topols:
            topol.set_index( index )
            index += 1

        out_str = ["PAIRING_STATISTICS "+str(index-1)]
        out_str.extend([ str( topol ) for topol in topols ])
        return "\n".join( out_str )
