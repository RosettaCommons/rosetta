import networkx as nx
import string
import os
import tempfile
import subprocess
import shutil
import glob

hetatm_map = { '5BU':'  U', ' MG':' MG', 'OMC':'  C', '5MC':'  C', 'CCC':'  C', ' DC':'  C', 'CBR':'  C', 'CBV':'  C', 'CB2':'  C', '2MG':'  G', 'H2U':'H2U', 'PSU':'PSU', '5MU':'  U', 'OMG':'  G', '7MG':'  G', '1MG':'  G', 'GTP':'  G', 'AMP':'  A', ' YG':'  G', '1MA':'  A', 'M2G':'  G', 'YYG':'  G', ' DG':'  G', 'G46':'  G', ' IC':' IC',' IG':' IG', 'ZMP':'ZMP' }

longer_names={'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D',
              'CYS': 'C', 'GLU': 'E', 'GLN': 'Q', 'GLY': 'G',
              'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
              'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S',
              'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
              ' rA': 'a', ' rC': 'c', ' rG': 'g', ' rU': 'u',
              '  A': 'a', '  C': 'c', '  G': 'g', '  U': 'u',
              ' MG': 'Z[MG]',' IC':'c[ICY]',' IG':'g[IGU]',
              'ROS': 'Z[ROS]','HOH':'w[HOH]', 'H2U': 'X[H2U]',
              'PSU': 'X[PSU]', '5MU': 'X[5MU]', 'FME': 'X[FME]',
          'U33': 'X[U33]', '  I': 'X[INO]', 'BRU': 'X[5BU]',
              }

# stem (or stem with single nucleotide bulge) = node
# multiloop = node
# internal loop = node
# hairpin = node

def get_helix_stems( secstruct ):

    def get_bp_list(secstruct):
        left_base_in_bp = []
        bp = []
        for i, char in enumerate(secstruct):
            if char == '(':
                left_base_in_bp.append(i + 1)
            elif char == ')':
                if not left_base_in_bp:
                    raise ValueError("Invalid secstruct!")
                bp.append((left_base_in_bp.pop(), i + 1))
        if left_base_in_bp:
            raise ValueError("Invalid secstruct!")
        for i, char in enumerate(secstruct):
            if char == '[':
                left_base_in_bp.append(i + 1)
            elif char == ']':
                if not left_base_in_bp:
                    raise ValueError("Invalid secstruct!")
                bp.append((left_base_in_bp.pop(), i + 1))
        if left_base_in_bp:
            raise ValueError("Invalid secstruct!")
        for i, char in enumerate(secstruct):
            if char == '{':
                left_base_in_bp.append(i + 1)
            elif char == '}':
                if not left_base_in_bp:
                    raise ValueError("Invalid secstruct!")
                bp.append((left_base_in_bp.pop(), i + 1))
        if left_base_in_bp:
            raise ValueError("Invalid secstruct!")
        for i, char in enumerate(secstruct):
            if char == '<':
                left_base_in_bp.append(i + 1)
            elif char == '>':
                if not left_base_in_bp:
                    raise ValueError("Invalid secstruct!")
                bp.append((left_base_in_bp.pop(), i + 1))
        if left_base_in_bp:
            raise ValueError("Invalid secstruct!")
        return bp

    def get_stems(bp):
        stems = []
        stem = []
        for i, j in bp:
            if (not stem) or (i - 1 == stem[-1][0] and j + 1 == stem[-1][1]):
                stem.append((i, j))
            else:
                stems.append(stem)
                stem = [(i, j)]
        if stem:
            stems.append(stem)
        return stems

    bp = sorted(get_bp_list(secstruct))
    return get_stems(bp)


def graph_from_ss( ss_init, add_closing_pair=True ):

    if add_closing_pair and (ss_init[0] == '.' or ss_init[-1] == '.'):
        ss = '(' + ss_init + ')'
        offset_numbering_for_closing_pair = True
    else:
        ss = ss_init
        offset_numbering_for_closing_pair = False


    elt_assignment = []
    # start with empty assignment
    for res in range( len(ss) ):
        elt_assignment.append( '' )

    # get the stems
    stems = get_helix_stems( ss )

    # assign stems 
    for stem_num, stem in enumerate(stems):
        for pair in stem:
            elt_assignment[ pair[0]-1 ] = 'S-' + str(stem_num) + '-L'
            elt_assignment[ pair[1]-1 ] = 'S-' + str(stem_num) + '-R'


    left_elt = None
    right_elt = None
    # assign loops
    for r_index, res in enumerate(elt_assignment):
        if res.startswith( 'S' ):
            left_elt = res
            continue
        right_elt = None
        for r2 in elt_assignment[r_index:]:
            if r2.startswith('S'):
                right_elt = r2
                break
        if left_elt is None:
            elt_assignment[ r_index ] = 'F'
        elif right_elt is None:
            elt_assignment[ r_index ] = 'T'
        else:
            elt_assignment[ r_index ] = 'L_%s_%s' %(left_elt, right_elt)

    # insert BOGUS loop elements between stems
    new_elt_assignment = []
    new_elt_assignment_orig_resnums = []
    for index, elt in enumerate(elt_assignment):
        if index != 0:
            prev_elt = elt_assignment[index - 1]
            if elt.startswith('S') and prev_elt.startswith('S') and elt != prev_elt:
                new_elt_assignment.append( 'L_%s_%s' %(prev_elt, elt) )
                new_elt_assignment_orig_resnums.append( index -0.5 )
        new_elt_assignment.append( elt )
        new_elt_assignment_orig_resnums.append( index )
    orig_elt_assignment = elt_assignment
    elt_assignment = new_elt_assignment
    elt_assignment_dict = {}

    unassigned_loop_res = []
    # assign hairpins and make list of unassigned loop residues
    for elt_index, elt in enumerate( elt_assignment ):
        if elt.startswith('L'):
            start_stem_num = elt.split('_')[1].split('-')[1]
            stop_stem_num = elt.split('_')[2].split('-')[1]
            if start_stem_num == stop_stem_num:
                elt_assignment[ elt_index ] = 'H' + start_stem_num
                elt_assignment_dict[ elt ] = 'H' + start_stem_num
            else: unassigned_loop_res.append( elt_index )

    loop_num = 0
    HACK_STOP = 0
    while len( unassigned_loop_res ) > 0:
        loop_group = []
        loop_group.append( unassigned_loop_res[0] )
        start_connecting_elt = elt_assignment[ loop_group[0] ].split('_')[1]
        next_connecting_elt = elt_assignment[ loop_group[0] ].split('_')[2]
        finish_connecting_elt = elt_assignment[ loop_group[0] ].split('_')[1]
        if finish_connecting_elt[-1] == 'L':
            finish_connecting_elt = finish_connecting_elt[0:-1] + 'R'
        elif finish_connecting_elt[-1] == 'R':
            finish_connecting_elt = finish_connecting_elt[0:-1] + 'L'
        if next_connecting_elt[-1] == 'L':
            next_connecting_elt = next_connecting_elt[0:-1] + 'R'
        elif next_connecting_elt[-1] == 'R':
            next_connecting_elt = next_connecting_elt[0:-1] + 'L'
        # add all other elements in that specific loop
        for res in unassigned_loop_res:
            if elt_assignment[ res ] == elt_assignment[ loop_group[0] ]:
                loop_group.append( res )

        # and remove those elts from the unassigned_loop_res
        for res in loop_group:
            if res in unassigned_loop_res:
                unassigned_loop_res.remove( res )

        if len( unassigned_loop_res ) == 0: 
            for elt in loop_group:
                elt_assignment_dict[ elt_assignment[elt] ] = 'L' + str(loop_num)
                elt_assignment[ elt ] = 'L' + str(loop_num)
            loop_num += 1
            break

        # check whether there is anything that starts with end_connecting_elt opposite side
        end_connecting_elt = ''
        # or once we've gone in a complete circle
        for res_index, res in enumerate( elt_assignment ):
            # if it's a loop, check whether it belongs here
            if res.startswith( 'L_' ):
                if elt_assignment[ res_index ].split('_')[1] == next_connecting_elt:
                    loop_group.append( res_index )
                    for other_res in unassigned_loop_res:
                        if elt_assignment[ other_res ] == elt_assignment[ res_index ]:
                            loop_group.append( other_res )
                next_connecting_elt = elt_assignment[ loop_group[-1] ].split('_')[2]
                end_connecting_elt = elt_assignment[ loop_group[-1] ].split('_')[2]
                if next_connecting_elt[-1] == 'L':
                    next_connecting_elt = next_connecting_elt[0:-1] + 'R'
                elif next_connecting_elt[-1] == 'R':
                    next_connecting_elt = next_connecting_elt[0:-1] + 'L'
            # otherwise check if we have gone far enough
            if end_connecting_elt == finish_connecting_elt: break
        for res in loop_group:
            if res in unassigned_loop_res:
                unassigned_loop_res.remove( res )

        # ok finished with this element
        for elt in loop_group:
            elt_assignment_dict[ elt_assignment[elt] ] = 'L' + str( loop_num )
            elt_assignment[ elt ] = 'L' + str( loop_num )
        loop_num += 1

    for i in range( len(elt_assignment)):
        if elt_assignment[i].startswith('S'):
            elt_assignment_dict[ elt_assignment[i] ] = 'S' + elt_assignment[i].split('-')[1]
            elt_assignment[i] = 'S' + elt_assignment[i].split('-')[1]

    residues_in_each_elt = {}
    # make a dictionary of residues in each elt assignment group
    for group_name in set( elt_assignment ):
        residues_in_each_elt[ group_name ] = []
        for res_index, res in enumerate(orig_elt_assignment):
            if elt_assignment_dict[ res ] == group_name:
                if offset_numbering_for_closing_pair:
                    residues_in_each_elt[ group_name ].append( res_index-1 )
                else:
                    residues_in_each_elt[ group_name ].append( res_index )
        # if there are no residues in the element, get the "adjacent" residues
        if len( residues_in_each_elt[ group_name ] ) == 0:
            # this must be an element that was added in the new_elt_assignment
            for index, res in enumerate( new_elt_assignment ):
                if res == group_name:
                    if offset_numbering_for_closing_pair:
                        residues_in_each_elt[ group_name ].append( new_elt_assignment_orig_resnums[ index ]-1 )
                    else:
                        residues_in_each_elt[ group_name ].append( new_elt_assignment_orig_resnums[ index ]-1 )


    element_connections = {}
    for res_index, res in enumerate(elt_assignment):
        if res not in element_connections.keys():
            element_connections[ res ] = []
        for r2 in elt_assignment[res_index::-1]:
            if r2 != res:
                if r2 not in element_connections[res]:
                    element_connections[ res ].append( r2 )
                break
        for r2 in elt_assignment[ res_index:]:
            if r2 != res:
                if r2 not in element_connections[res]:
                    element_connections[ res ].append( r2 )
                break

    stem_connections = {}
    # get 4 connections for each stem
    for s_index, stem in enumerate(stems):
        res1_bottom_pair = stem[0][0] # --> C1
        res1_top_pair = stem[-1][0] # --> C2
        res2_top_pair = stem[-1][1] # --> C3
        res2_bottom_pair = stem[0][1] # --> C4
        if res1_bottom_pair - 2 >= 0:
            C1 = elt_assignment[ res1_bottom_pair - 2]
        else: C1 = None
        C2 = elt_assignment[ res1_top_pair ]
        C3 = elt_assignment[ res2_top_pair - 2 ]
        if res2_bottom_pair < len(elt_assignment):
            C4 = elt_assignment[ res2_bottom_pair ]
        else: C4 = None
        stem_connections['S' + str(s_index)] = []
        stem_connections['S' + str(s_index)].append( C1 )
        stem_connections['S' + str(s_index)].append( C2 )
        stem_connections['S' + str(s_index)].append( C3 )
        stem_connections['S' + str(s_index)].append( C4 )


    # make a graph
    G = nx.Graph()
    for elt in element_connections.keys():
        G.add_node( elt, residues = residues_in_each_elt[ elt ] )
        # and add the node with the residues that are in it
        for elt2 in element_connections[ elt ]:
            G.add_edge( elt, elt2 )

    if offset_numbering_for_closing_pair:
        G.remove_node( 'S0' )

    return G, stem_connections

def find_pairs( secstruct, helix_chars=[['(',')'],['[',']'],['{','}'], ['<','>']] ):
    ss = list(secstruct) # turn into a list, so we can do replacements
    all_pairs = []  # list of dictionaries

        # find pairs
    for helix_char_pair in helix_chars:
        for i in range(len(ss)):
            count = 0 
            if ss[i]==helix_char_pair[0]:
                count += 1
                pair = []
                pair.append( i )
                for j in range(i+1,len(ss)):
                    if ss[j]==helix_char_pair[0]: count+=1
                    elif ss[j] ==helix_char_pair[1]: count-=1
                    if count == 0:
                        pair.append( j )
                        all_pairs.append(pair)
                        # get rid of those parentheses now
                        ss[pair[0]] = ''
                        ss[pair[1]] = ''
                        break

    return all_pairs

def residue_graph_from_ss( ss ):
    # make a graph where each node is a residue
    # there are connections between base paired and adjacent residues
    base_pairs = find_pairs( ss )

    G = nx.Graph()
    for res in range( len(ss) ):
        G.add_node( res )
    # add edges for adjacent edges
    for res in range( len(ss) -1 ):
        G.add_edge( res, res + 1 )
    # add edges for base pairs
    for pair in base_pairs:
        G.add_edge( pair[0], pair[1] )

    return G

def pdbslice( pdbfile, subset_or_excise, segment_residues, prefix ):

    if not os.path.exists( pdbfile ):
        print( "ERROR: %s does not exist!" %(pdbfile) )
        exit( 1 )

    use_subset = False
    use_excise = False
    if subset_or_excise == 'subset':
        use_subset = True
    elif subset_or_excise == 'excise':
        use_excise = True

    segment_resnums, segment_chains = parse_tag( segment_residues )

    atomnums = []
    outfile = open( prefix + os.path.basename( pdbfile ), 'w')

    for line in open( pdbfile ):
        if len( line ) > 4  and ( line[:6] == 'ATOM  ' or line[:6] == 'HETATM'):
            currentchain = line[21]
            currentresidue = line[22:26]
            atomnum = int( line[6:11] )
            try:
                currentresidue_val = int(currentresidue)
                if (use_subset and \
                    not matches( currentresidue_val, currentchain, segment_resnums, segment_chains ) ): continue
                if (use_excise and \
                    matches( currentresidue_val, currentchain, segment_resnums, segment_chains ) ): continue
            except:
                continue

            atomnums.append( atomnum )

            outfile.write(line)

        elif len( line ) > 6 and line[:6] == 'CONECT':
            atomnum1 = int( line[6:11] )
            atomnum2 = int( line[11:16] )
            if (atomnum1 in atomnums) and (atomnum2 in atomnums): outfile.write( line )

    outfile.close()

def matches( res, chain, match_res, match_chain ):
    assert( len( match_res ) == len( match_chain ) )
    for n in range( len( match_res ) ):
        if ( res == match_res[n] and ( match_chain[n] == ''  or match_chain[n] == chain ) ):
            return True
    return False

def make_tag( int_vector ):
    tag = ''
    for m in int_vector: tag += ' %d' % m
    return tag

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

def parse_tag( tag, alpha_sort=False ):

    if isinstance( tag, list ):
        tag = ' '.join( tag )
        #tag = string.join( tag, ' ' )

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
            subtag = ' '.join([xchar+':'+str(x) for x in range(start_idx[0],stop_idx[0]+1)])
            #subtag = string.join([xchar+':'+str(x) for x in xrange(start_idx[0],stop_idx[0]+1)], ' ')
            int_vector.extend( parse_tag( subtag )[0] )
            char_vector.extend( parse_tag( subtag )[1] )
            continue
   
        if ':' in subtag: # A:100
            subtag = subtag.split(':')
            xchar = subtag[0]
            xint = int(subtag[-1])
        else: # A100 or 100 or 0100
            for x in range( len( subtag ) ):
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

def get_sequences( pdbname, removechain = 0 ):

    netpdbname = pdbname
    assert( os.path.exists(netpdbname))

    lines = open(netpdbname,'r').readlines()

    oldresnum = '   '
    oldchain = ''
    chain = oldchain
    resnum = oldresnum
    count = 0;
    fasta_line = ''

    sequences = []
    all_chains = []
    all_resnums = []

    sequence = ''
    chains = []
    resnums = []

    for line in lines:
        if (len(line)<20): continue
        line_edit = line
        if (line[0:6] == 'HETATM') & (line[17:20]=='MSE'): #Selenomethionine
            line_edit = 'ATOM  '+line[6:17]+'MET'+line[20:]
            if (line_edit[12:14] == 'SE'):
                line_edit = line_edit[0:12]+' S'+line_edit[14:]
            if len(line_edit)>75:
                if (line_edit[76:78] == 'SE'):
                    line_edit = line_edit[0:76]+' S'+line_edit[78:]
        elif (line[0:6] == 'HETATM') & ( line[17:20] in hetatm_map.keys()):
            line_edit = 'ATOM  '+line[6:17]+ hetatm_map[line[17:20]] + line[20:]

        if (line[0:6] == 'HETATM') & (line[17:20] in longer_names.keys() ):
            line_edit = 'ATOM  '+line[6:]

        if line_edit[0:4] == 'ATOM' or line_edit[0:6]=='HETATM':
            resnum = line_edit[22:26].replace( ' ', '' )
            chain = line_edit[21]

        if ( line[0:3] == 'TER' or ( not chain == oldchain ) ) and len( sequence ) > 0:
            sequences.append( sequence )
            all_chains.append( chains )
            all_resnums.append( resnums )
            sequence = ''
            chains   = []
            resnums  = []
            old_resnum = ''

        if (not (resnum == oldresnum and chain == oldchain) ):
            count = count + 1
            longname = line_edit[17:20]
            if longname in longer_names.keys():
                sequence +=  longer_names[longname]
            else:
                sequence +=  'X'
            resnums.append( int(resnum) )
            chains.append( chain )
            oldresnum = resnum
            oldchain = chain

    if len( sequence ) > 0:
        sequences.append( sequence )
        if ( chain == ' ' ): chain = ''
        all_chains.append( chains )
        all_resnums.append( resnums )

    return ( sequences, all_chains, all_resnums )

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

def read_pdb_with_segid( filename ):

    coords = {}
    pdb_lines = {}
    sequence = {}

    old_resnum = 0
    old_chain  = ''
    old_segid = '    '
    chains = []
    residues = []
    segids = []
    for line in open( filename ):

        if (len(line)>54 and  (line[0:4] == 'ATOM' or line[0:4] == 'HETA' ) ):

            resnum = int( line[22:26] )
            chain = line[21]
            segid = line[72:76]
            atom_name = line[12:16]
            position = [float(line[30:38]),float(line[38:46]),float(line[46:54])]

            if not ( chain in coords.keys() ):
                coords[chain] = {}
                pdb_lines[chain] = {}
                sequence[ chain ] = {}
            if not (segid in coords[chain].keys()):
                coords[chain][segid] = {}
                pdb_lines[chain][segid] = {}
                sequence[ chain ][segid] = {}

            sequence[chain][segid][resnum] = line[17:20]

            if not ( resnum in coords[chain][segid].keys() ):
                coords[chain][segid][resnum] = {}
                pdb_lines[chain][segid][resnum] = {}

            coords[chain][segid][resnum][atom_name] = position
            pdb_lines[chain][segid][resnum][atom_name] = line[:-1]

            if ( len(residues) == 0 or \
                 resnum != old_resnum or chain != old_chain or segid != old_segid ):
                chains.append( chain )
                residues.append( resnum )
                segids.append( segids )
            old_resnum = resnum
            old_chain  = chain
            old_segid  = segid

    return ( coords, pdb_lines, sequence, chains, residues, segids )

def better_reorder_pdb( pdb ):

    assert( pdb[-4:] == '.pdb' )
    outfile = pdb.replace( '.pdb', '.REORDER.pdb' )
    with open( outfile, 'w' ) as fid:
        [ coords_main, lines_main, sequence_main, all_chains, __, all_segids ] = read_pdb_with_segid( pdb )
        chains = []
        segids = []
        for c in all_chains:
                if c not in chains:
                    chains.append( c)

        for chain in chains:
            for segid in lines_main[chain].keys():
                 residues = []
                 for r in lines_main[ chain ][ segid ].keys():
                     residues.append ( r )
                 #residues = lines_main[ chain ][ segid ].keys()
                 residues.sort()
                 for residue in residues:
                     for atom in lines_main[ chain ][segid][ residue ].keys():
                         fid.write(  lines_main[ chain ][ segid ][ residue ][ atom ]+'\n' )



def rna_helix( seq, resnum, output_pdb, rosetta_directory, rosetta_extension, test=False ):
    #Build the helix
    temp = tempfile.NamedTemporaryFile(delete=False)
    cmdline = os.path.expanduser( rosetta_directory ) + '/rna_helix' + rosetta_extension
    if not os.path.exists( cmdline ):
	    print( 'ERROR: ' + cmdline + ' does not exist.' )
	    print( 'Check that you have specified the correct -rosetta_directory and -rosetta_extension' )
	    exit( 1 )
    cmdline += (' -rna::corrected_geo  '+
                '-score:rna_torsion_potential RNA11_based_new ' +
                '-chemical::enlarge_H_lj ')
    cmdline += '-o %s ' % temp.name
    cmdline += '-seq ' + seq + ' '
    if test:
	    cmdline += '-testing:INTEGRATION_TEST '
    p = subprocess.Popen(cmdline.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE )
    out, err = p.communicate()
    if err and len(err):
        print( out )
        print( err )
    assert( not err or not len(err) )

    #Renumber the pdb
    resnums = []
    chains  = []
    segids = []
    for i in resnum.split():  get_resnum_chain(i, resnums, chains, segids )
    renumber_pdb([temp.name], resnums, chains, segids)

    temp.close()
    shutil.move(temp.name, output_pdb)

def get_resnum_chain( input_string, resnums, chains, segids ): # could be of the form A:1-4 or A1-4 or 1-4

    if ( input_string.find( ',' ) > -1 ):
        for substring in string.split( input_string, ',' ): get_resnum_chain( substring, resnums, chains, segids )
        return True
    #assert( len( input_string ) > 0 )
    #if len(segids) == 0: segids = ['    ' for _ in resnums]
    assert( len( resnums ) == len( chains ) )
    assert( len( resnums ) == len( segids ) )
    ok = True

    int_string = input_string
    chain = ''
    segid = '    '

    # extract chain info, if present.
    # colon is clear indicator of chain
    if input_string.count(':') == 2:
        first = input_string.find( ':' )
        chain = input_string[:first]
        second = first+2+input_string[ first+2 : ].find(':')
        segid = input_string[ first+1 : second ]
        int_string = input_string[second+1: ]
    else:
        colon_pos = input_string.find( ':' )
        if colon_pos > -1:
            chain = input_string[:colon_pos]
            int_string = input_string[ colon_pos+1 : ]
        else: # look for non-numeric characters
            pos = 0
            while pos < len( input_string ):
                try:
                    if ( input_string[pos] != '-' ): blah = int( input_string[ pos ] )
                    break
                except:
                    chain += input_string[ pos ]
                pos = pos + 1
            int_string = input_string[ pos: ]

    ints = []
    if len( int_string ) == 0 and len( chain ) == 1: # special case, get the whole chain.
        resnums.append( 'all' ) # means get everything from this chain
        chains.append( chain )
    else:
        ok = get_ints( int_string, ints )
        for i in range( len( ints ) ):
            resnums.append( ints[i] )
            chains.append( chain )
            segids.append( segid )
    assert( len( resnums ) == len( chains ) )
    assert( len( resnums ) == len( segids ) )
    return ok

def renumber_pdb(pdbnames, new_numbers, chains = [], segids = [], retain_atom_num = 0):

    for pdbname in pdbnames:
        lines = open(pdbname,'r').readlines()

        oldresnum = '   '
        count = 0;
        outid  = open( 'temp.txt','w')
        atomnum  = 0
        newchain = ''
        for line in lines:
            line_edit = line
            if line[0:3] == 'TER':
                continue

            if line_edit[0:4] == 'ATOM' or line_edit[0:6] == 'HETATM':
                if not (line[16]==' ' or line[16]=='A'): continue
                atomnum += 1
                oldchain = line_edit[21]
                resnum = line_edit[23:26]
                oldsegid = line_edit[72:76]
                if not resnum == oldresnum:
                    count = count + 1
                    oldresnum = resnum

                if ( count <= len( new_numbers ) ):
                    newnum = '%4d' % new_numbers[ count-1 ]
                    if len( chains ) > 0:  newchain = chains[ count - 1]
                    if len( segids ) > 0:  newsegid = segids[ count - 1]
                    if newchain == '': newchain = line_edit[ 22 ]
                    if newsegid == '': newsegid = line_edit[ 72:76 ]
                else:
                    if len( new_numbers) > 0:
                        print('WARNING! residue number %d is greater than length of input numbering %d' % (count, len( new_numbers) ))
                    newnum = '%4d' % count
                    newchain = oldchain
                    newsegid = oldsegid

                if len(newsegid) == 1:
                    newsegid = newsegid+'   '
                elif len(newsegid) == 2:
                    newsegid = newsegid+'  '
                elif len(newsegid) == 3:
                    newsegid = newsegid+' '

                if retain_atom_num:
                    line_edit = '%s%s%s' % (line_edit[0:22], newnum, line_edit[26:] )
                else:
                    line_edit = '%s%5d%s%s%s%s%s%s' % (line_edit[0:6],atomnum,line[11:21],newchain, newnum, line_edit[26:72], newsegid, line_edit[76:] )


                outid.write(line_edit)

        if ( count < len( new_numbers) ): print('WARNING! number of residues %d is less than length of input numbering %d' % (count, len( new_numbers) ))

        outid.close()

        os.system( 'mv temp.txt '+pdbname )

def get_ints( int_string, value ): # could be of the form 5-8  or -9--5
    if ( len( int_string ) == 0 ): return False

    # find hyphen separator
    dashpos = 0
    if int_string[0] == '-': dashpos = 1  # can't be the first dash
    while dashpos < len( int_string ) and int_string[ dashpos ] != '-': dashpos += 1
    if dashpos == len( int_string ):
        try:
            value.append( int( int_string ) )
        except:
            return False
    else:
        try:
            int_start = int( int_string[          :dashpos] )
            int_stop  = int( int_string[dashpos+1 :       ] )
            assert( int_stop > int_start )
            for m in range( int_start, int_stop+1): value.append( m )
        except:
            return False
    return True

def easy_cat( outfile ):
	
	which_files_to_cat = {}
	
	globfiles = glob.glob( outfile+'/*/*out' )
	
	if len( globfiles ) == 0: globfiles = glob.glob( outfile + '/*out'  )
	
	# Remove any "checkpoint" files from stepwise checkpointing
	globfiles = filter( lambda x: "S_" not in x and "_checkpoint" not in x, globfiles )
	
	# make sure to order 0,1,2,... 10, 11, 12, ... 100, 101, ...
	globfiles_with_length = []
	for f in globfiles:
		globfiles_with_length.append( [len(f), f] )
	
	globfiles_with_length.sort()
	globfiles = []
	for f in globfiles_with_length:
		globfiles.append( f[1] )
	#for x in globfiles: print( x )
	
	for file in globfiles:
	    tag = os.path.basename( file ).replace('.out','')
	    if tag not in which_files_to_cat.keys():
	        which_files_to_cat[tag] = []
	    which_files_to_cat[tag].append( file )
	
	
	for tag in which_files_to_cat.keys():
	    cat_file = tag+".out"
	    cat_outfiles( which_files_to_cat[tag], cat_file )
	
	    lines = os.popen( 'grep SCORE '+cat_file).readlines()
	
	    fid_sc = open( cat_file.replace('.out','.sc'),'w' )
	    for line in lines:
	        fid_sc.write( line )
	    fid_sc.close()

def cat_outfiles( outfiles, output_file ):
	fid = open( output_file, "w" )
	
	sequence_line_found    = 0
	description_line_found = 0
	remark_line_found      = 0
	n_file = -1
	
	for out_f in outfiles:
	    all_files = glob.glob(out_f)
	    for filename in all_files:
	        data = open(filename)
	        n_file += 1
	        for line in data:
	            line = line[:-1]
	            if not line: break
	
	            if line[:9] == 'SEQUENCE:':
	                if sequence_line_found: continue # Should not be any more sequence lines!
	                else: sequence_line_found = 1
	
	            if line.find( 'description' ) > -1:
	                if description_line_found: continue
	                else: description_line_found = 1
	
	            description_index = line.find(' S_')
	            if description_index < 0:
	                description_index = line.find(' F_')
	
	            if description_index >= 0:
	                description_index -= 1 # to get rid of space.
	                tag = line[description_index:]
	                newtag = tag + "_%03d" % n_file
	                line = line[:description_index] + newtag
	
	            if len(line) < 1: continue
	
	            fid.write( line+'\n' )
	
	        data.close()
	
	fid.close()

def extract_lowscore_decoys( silent_file, NSTRUCT, rosetta_directory, rosetta_extension, test=False):
	
	MINI_DIR = rosetta_directory 
	
	tags = []
	
	scoretags = os.popen('head -n 2 '+silent_file).readlines()[1]
	#scoretags = string.split( popen('head -n 2 '+silent_file).readlines()[1] )
	scoretag=''
	
	scorecols  = [-1]
	
	
	binary_silentfile = 0
	remark_lines = os.popen('head -n 7 '+silent_file).readlines()
	for line in remark_lines:
	    if ( len( line ) > 6 and line[:6] == "REMARK" ):
	        remark_tags = line.split()
	        if remark_tags.count('BINARY_SILENTFILE'):
	            binary_silentfile = 1
	        if remark_tags.count('BINARY'):
	            binary_silentfile = 1
	
	coarse = 0
	if os.path.exists( 'remark_tags') and remark_tags.count('COARSE'):
	    coarse = 1
	
	assert(silent_file[-3:] == 'out')
	
	# Make the list of decoys to extract
	lines = os.popen( 'grep SCORE '+silent_file+' | grep -v NATIVE').readlines()
	
	score_plus_lines = []
	for line in lines:
	    cols =  line.split()
	    #cols = string.split( line )
	    score = 0.0
	    try:
	        for scorecol in scorecols: score += float( cols[ abs(scorecol) ] )
	    except:
	        continue
	    score_plus_lines.append( ( score, line ))
	
	score_plus_lines.sort()
	lines = map( lambda x:x[-1], score_plus_lines[:NSTRUCT] )
	
	templist_name = 'temp.%s.list' %(os.path.basename(silent_file))
	
	fid = open(templist_name,'w')
	count = 0
	for line in lines:
	    cols = line.split()
	    #cols = string.split(line)
	    tag = cols[-1]
	    if tag.find('desc') < 0:
	        fid.write(tag+'\n')
	        tags.append(tag)
	        count = count+1
	    if NSTRUCT and count >= NSTRUCT:
	        break
	outfilename = silent_file
	
	fid.close()
	
	#Set up bonds file?
	softlink_bonds_file = 0
	wanted_bonds_file = silent_file+'.bonds'
	wanted_rot_templates_file = silent_file+'.rot_templates'
	bonds_files = glob.glob( '*.bonds')
	if len( bonds_files ) > 0:
	    if not os.path.exists( wanted_bonds_file ):
	        softlink_bonds_file = 1
	        system( 'ln -fs '+bonds_files[0]+' '+wanted_bonds_file )
	        system( 'ln -fs '+bonds_files[0].replace('.bonds','.rot_templates') \
	                +' '+wanted_rot_templates_file )
	
	
	MINI_EXE = MINI_DIR+'extract_pdbs'+rosetta_extension
	
	tag_str = ''
	for tag in tags:
		tag_str += tag + ' '
	command = '%s -load_PDB_components -in:file:silent %s -in:file:tags %s' % ( MINI_EXE, outfilename, tag_str )
	if test:
	        command += ' -testing:INTEGRATION_TEST '
	
#	# Check if this is an RNA run.
#	with open( silent_file ) as fid:
#		line = fid.readline() # Should be the sequence.
#		rna = 0
#		sequence = line.split()[-1]
#		#sequence = string.split(line)[-1]
#		rna = 1
#		for c in sequence:
#		    if not ( c == 'a' or c == 'c' or c == 'u' or c == 'g'):
#		        rna = 0
#		        break
#	if rna: command  += ' -enable_dna -enable_rna '
	
	# Check if this is full atom.
	lines = os.popen('head -n 8 '+outfilename).readlines()
	if len(lines[6].split()) > 10:
	#if len(string.split(lines[6])) > 10:
	    command += ' -fa_input'
	
#	if rna:
#	    command = '%s -load_PDB_components -in::file::silent %s -tags %s  -extract' % \
#	              ( MINI_EXE, outfilename, tag_str )
#	              #( MINI_EXE, outfilename, string.join( tags ) )
#	
#	    if binary_silentfile:
#	        silent_struct_type = 'binary_rna'
#	    else:
#	        silent_struct_type = 'rna'
#	
#	    command = '%s -load_PDB_components -in:file:silent %s -in:file:tags %s -in:file:silent_struct_type %s -mute all ' % \
#	              ( MINI_EXE, outfilename, tag_str, silent_struct_type )
#	
#	    if coarse:
#	        command += " -out:file:residue_type_set coarse_rna "
#	    else:
#	        pass
#	
#	elif ( binary_silentfile ):
	if ( binary_silentfile ):
	
	    command = '%s -load_PDB_components -in:file:silent  %s  -in:file:silent_struct_type binary -in:file:fullatom -in:file:tags %s -mute all' % \
	              ( MINI_EXE, outfilename, tag_str )
	    if test:
	            command += ' -testing:INTEGRATION_TEST '
	    
	    if (scoretags.count('vdw')): command += ' -out:file:residue_type_set centroid '
	
	#print(command)
	os.system(command)
	
	count = 1
	
	for tag in tags:
	    command = 'mv %s.pdb %s.%d.pdb' % (tag,os.path.basename(silent_file),count)
	    #print(command)
	    os.system(command)
	    count += 1
	
	command = 'rm '+templist_name
	#print(command)
	os.system(command)
	
	if (softlink_bonds_file):
	    #system( 'rm '+wanted_bonds_file+' '+wanted_rot_templates_file )
	    print( ' WARNING! WARNING' )
	    print( ' Found a .bonds and .rot_templates file and used it!' )
	
