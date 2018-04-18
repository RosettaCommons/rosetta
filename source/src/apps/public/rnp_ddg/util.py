import string

def read_pdb( filename ):

    coords = {}
    pdb_lines = {}
    sequence = {}

    old_resnum = 0
    old_chain  = ''
    chains = []
    residues = []
    for line in open( filename ):

        if (len(line)>54 and  (line[0:4] == 'ATOM' or line[0:4] == 'HETA' ) ):

            resnum = int( line[22:26] )
            chain = line[21]
            atom_name = line[12:16]
            position = [float(line[30:38]),float(line[38:46]),float(line[46:54])]

            if not ( chain in coords.keys() ):
                coords[chain] = {}
                pdb_lines[chain] = {}
                sequence[ chain ] = {}

            sequence[chain][resnum] = line[17:20]

            if not ( resnum in coords[chain].keys() ):
                coords[chain][resnum] = {}
                pdb_lines[chain][resnum] = {}

            coords[chain][resnum][atom_name] = position
            pdb_lines[chain][resnum][atom_name] = line[:-1]

            if ( len(residues) == 0 or \
                 resnum != old_resnum or chain != old_chain ):
                chains.append( chain )
                residues.append( resnum )
            old_resnum = resnum
            old_chain  = chain

    return ( coords, pdb_lines, sequence, chains, residues )

def make_tag_with_dashes( int_vector, char_vector = 0 ):
    tag = ''
    if not isinstance( char_vector, list ) or len( char_vector ) == 0:
         char_vector = []
         for m in range( len( int_vector ) ): char_vector.append( '' )
    assert( len( char_vector ) == len( int_vector ) )

    start_res = int_vector[0]
    for i in range( 1, len(int_vector)+1 ):
        if i==len( int_vector)  or \
                int_vector[i] != int_vector[i-1]+1 or \
                char_vector[i] != char_vector[i-1] :

            stop_res = int_vector[i-1]
            tag += ' '
            if len( char_vector[i-1] ) > 0:
                assert( len( char_vector[i-1] ) == 1 )
                tag += char_vector[i-1]+':'
            if stop_res > start_res:
                tag += '%d-%d' % (start_res, stop_res )
            else:
                tag += '%d' % (stop_res )

            if ( i < len( int_vector) ): start_res = int_vector[i]

    return tag

def get_surrounding_res( pdbfile, sample_res_list=[], radius=None, verbose=False, keep_res=[] ):
	
	### Parse sample_res_list as a string
	sample_res_list, sample_chain_list = parse_tag( sample_res_list ) 
        
	### Read PDB file
	( coords, pdb_lines, sequence, chains, residues ) = read_pdb( pdbfile )

	### Find surrounding residues if radius != None
	if radius:


                if not keep_res is None:
                        keep_res_list, keep_res_chain_list = parse_tag( keep_res )
                else:
                        keep_res_list, keep_res_chain_list = [], []

		for surr_rsd_idx, surr_rsd in enumerate(residues):
			is_surrounding_res = False		

			for sample_rsd_idx, sample_rsd in enumerate(sample_res_list):
				if is_surrounding_res:	break

				for surr_atom in coords[chains[surr_rsd_idx]][surr_rsd]:
					if is_surrounding_res:	break
					surr_atom_xyz = coords[chains[surr_rsd_idx]][surr_rsd][surr_atom] 

					for sample_atom in coords[sample_chain_list[sample_rsd_idx]][sample_rsd]:
						if is_surrounding_res:	break
						sample_atom_xyz = coords[sample_chain_list[sample_rsd_idx]][sample_rsd][sample_atom]

						distance = math.sqrt( sum( power( subtract( surr_atom_xyz, sample_atom_xyz ) , 2 ) ) )
						if ( math.pow( distance, 2 ) < power( float(radius), 2 ) ):
							if verbose:	print "res "+str(chains[surr_rsd_idx])+':'+str(surr_rsd)+" is a surrounding res, distance() = ", distance
							keep_res_list.append( surr_rsd )
							keep_res_chain_list.append( chains[surr_rsd_idx] )
							is_surrounding_res = True
							break

		if not len( keep_res_list ):

			if not len( sample_res_list ):
				if verbose:	print 'WARNING: len(sample_res_list) == ', len(sample_res_list), ' but radius == ', radius
				if verbose:	print 'Must supply sample_res to find residues within the given radius.' 
			else:
				if verbose:	print 'WARNING: len(keep_res_list) == ', len(surrounding_res_list), ' for radius == ', radius
				if verbose:	print 'Try an expanded radius.'
		
	### All residues are surrounding if radius == None
	else:	
		keep_res_list = residues 
		keep_res_chain_list = chains

	assert( len(keep_res_list) == len(keep_res_chain_list) )

	### Remove sample residues from keep_res_list
	surrounding_res_list = []
	surrounding_res_chain_list = []
	for rsd_idx, rsd in enumerate(keep_res_list):
		chain = keep_res_chain_list[ rsd_idx ]
		if rsd in sample_res_list:
			if chain == sample_chain_list[ sample_res_list.index(rsd) ]:  continue
		surrounding_res_list.append( rsd )
		surrounding_res_chain_list.append( chain )

	assert( len(surrounding_res_list) == len(surrounding_res_chain_list) )

	return surrounding_res_list, surrounding_res_chain_list



def parse_tag( tag, alpha_sort=False ):

    if isinstance( tag, list ):
        tag = string.join( tag, ' ' )
    
    int_vector = []
    char_vector= []
    
    xchar = ''

    tag = tag.replace(',',' ')
    tag = tag.replace(';',' ')
    
    for subtag in tag.split(' '):

        if subtag == '': 
            continue

        if '-' in subtag: # '1-10' or 'A1-10' or 'A:1-10' or 'A:1-A:10'
            ( start, stop ) = subtag.split('-')  
            ( start_idx, start_char ) = parse_tag( start )
            ( stop_idx, stop_char ) = parse_tag( stop )
            assert( ( start_char[0] == stop_char[0] ) or ( stop_char[0] == '' ) )
            if start_char[0] != '': 
                xchar = start_char[0]
            subtag = string.join([xchar+':'+str(x) for x in xrange(start_idx[0],stop_idx[0]+1)], ' ')
            int_vector.extend( parse_tag( subtag )[0] )
            char_vector.extend( parse_tag( subtag )[1] )
            continue
   
        if ':' in subtag: # A:100
            subtag = subtag.split(':')
            xchar = subtag[0]
            xint = int(subtag[-1])
        else: # A100 or 100 or 0100            
            for x in xrange( len( subtag ) ):
                try: # 100
                    xint = int(subtag[x:])
                    break
                except: # A100
                    xchar = subtag[x]
          
        int_vector.append( xint )
        char_vector.append( xchar )

    assert( len(int_vector) == len(char_vector) ) 

    if alpha_sort:
        sorted = zip( char_vector, int_vector )
        sorted.sort()
        [ char_vector, int_vector ] = [ list(l) for l in zip(*sorted) ]
        
    return int_vector, char_vector

