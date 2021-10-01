# -*- coding: utf-8 -*-

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## Authors: David Thieker

################################################################################################
#
# Generate input or analyze results for Mutation Clusters.
#
# Stabilize_Cluster.py
# -mode is either 'trim' 'build' or 'analyze' [Required]
#     -mode trim: Prepare the PDB file for relax (remove heteroatoms)
#     -mode build: Build the directory tree, RosettaScripts (1 relax and 1 cluster) and jobdefinition file
#     -mode analyze: Analyze the results after design has been performed (collect best models based on total_score and report results such as score per mutation)
# -pdb is starting structure filename. The starting structure should already be relaxed (original rosetta script used AtomTree minimization) [Required]
# -seed is residue positions for cluster seed (mutation aroudn this residue). Format is a comma-delimited list with chain_ID preceding each residue (e.g. "H[100],L[23]). Insertion codes are allowed. [Required if -sym is not selected]
# -fav allows the user to apply a bonus score to native residues, thereby reducing the number of mutations in the output models (i.e. "0,1.5,3") [Optional]
# -lig allows the user to specify ligand residues to keep when mode==trim. Format is a comma-delimited list with chain_ID preceding each residue (i.e. "H[405-408],L[403]") [Optional]
# -sym allows mutations on multiple chains simultaneously (Format is chains:residues (i.e. "ABC[20-30,35,80]") [Required if -seed is not selected]
# -subs is amino acid residues that FastDesign is allowed to substitute during the simulations (if WT residue is not included, it will be forced to mutate) [Optional]
# -keep is amino acid residues that are NOT allowed to mutated during FastDesign (i.e. residues within a binding site) [Optional]
# -keep_sym is the same as -keep but reads the symmetric style of input residues [Optional]
#
# Example of trimming a PDB file:
# python Stabilize_Cluster.py -mode trim -pdb input_file.pdb 
#
# Example of building inputs:
# python Stabilize_Cluster.py -mode build -pdb input_file.pdb -seed H[1,12],L[20-30] -subs ADEFHIKLMNQRSTVWY
#
# Example of analyzing results from a symmetric run:
# python Stabilize_Cluster.py -mode analyze -pdb input_file.pdb -sym ABC[20,35,80] -keep_sym AB[27]
#
################################################################################################


''' ROSIE stabilize_cluster App - tested with Python v3.5 '''

import sys, os, re, shutil
import types, math, numpy, matplotlib
matplotlib.use('Agg') #Command-line only, don't display on screen
import matplotlib.pyplot as plt
import matplotlib.colors as colors

############# Print the RosettaScripts that are required for the Cluster calculations (rs_relax_cluster.xml, relax2.xml, and rs_cluster.xml). #############

def gen_RScripts(symmetric_input):
    xml_relax_cluster = '''\
<ROSETTASCRIPTS>
  <SCOREFXNS>
    <ScoreFunction name="sfxn" weights="ref2015"/>  # Score function without constraints
    <ScoreFunction name="sfxn_cc" weights="ref2015">   # Reweight the score function to include coordinate constraints
      <Reweight scoretype="coordinate_constraint" weight="1.0"/>
    </ScoreFunction>
  </SCOREFXNS>
  <RESIDUE_SELECTORS>
  </RESIDUE_SELECTORS>
  <TASKOPERATIONS>
    <InitializeFromCommandline name="ifcl_to"/> # Accept command-line options
    <IncludeCurrent name="incl_curr_to"/> # Include rotamer from input structure during packing
    <RestrictToRepacking name="repack_only" /> # Turn off design at these positions. Only repack and minimize
    <ExtraRotamersGeneric name="ex12_to" ex1="1" ex2="1" extrachi_cutoff="0"/> # Include additional rotamers during repacking
  </TASKOPERATIONS>
  <CONSTRAINT_GENERATORS>
    <CoordinateConstraintGenerator name="all_bb_ca_cc" sd="2" ca_only="true" native="true" />       # Add backbone constraints (C_alpha-only) for all residues
  </CONSTRAINT_GENERATORS>
  <MOVERS>
    <VirtualRoot name="vr" />
    <AddConstraints name="cc" constraint_generators="all_bb_ca_cc" />
    <FastRelax name="relax_atomTree" scorefxn="sfxn_cc" repeats="4" relaxscript="MonomerRelax2019" task_operations="ifcl_to,incl_curr_to,repack_only,ex12_to" cartesian="false" />
  </MOVERS>
  <PROTOCOLS>
    <Add mover="vr"/>
    <Add mover="cc"/>
    <Add mover="relax_atomTree"/>
  </PROTOCOLS>
  <OUTPUT scorefxn="sfxn"/> # Output score without penalties from constraints (cc/cc2)
</ROSETTASCRIPTS>
'''

# Pulled from protocol_mut_cluster_5-7-2021.xml
    xml_cluster = '''\
<ROSETTASCRIPTS>

  <SCOREFXNS>
    <ScoreFunction name="sfxn_cc" weights="ref2015_cst"/> # Reweight the score function to include coordinate and nativeresidue constraints
    <ScoreFunction name="sfxn" weights="ref2015"/>        # Score function without constraints
  </SCOREFXNS>

  <RESIDUE_SELECTORS>
    <Index name="residue_position_to_explore" resnums="%%seed%%" />                                     # Seed for design
    <Neighborhood name="neighbors" selector="residue_position_to_explore"/>                             # 10 Angstrom sphere around seed residue
    <Not name="not_neighbor" selector="neighbors"/>                                                     # Region outside neighborhood that is not explicitly manipulated
    <Neighborhood name="inner_shell" selector="residue_position_to_explore" distance="7.0"/>            # Inner shell within neighborhood, all residues allowed to mutate
    <Not name="not_mutant_neighbor" selector="inner_shell"/>                                            # Everything that isn't the design shell.
    <And name="outer_shell" selectors="neighbors,not_mutant_neighbor"/>                                 # Outer shell within neighborhood, repack and minimize only (NO mutations)
    <Index name="add_res_NO_design" resnums="%%keep%%" error_on_out_of_bounds_index="false"/>                           # User-selected residues that should never mutate (i.e. an important protein-protein interface)
    <And name="inner_NO_design" selectors="inner_shell,add_res_NO_design"/>                             # Additional user-defined residues that should not mutate
    <Or name="NOT-designable" selectors="outer_shell,inner_NO_design"/>                                 # All residues that should not be mutated
  </RESIDUE_SELECTORS>

  <JUMP_SELECTORS>
    <Interchain name="interchain_jumps"/>
  </JUMP_SELECTORS>

  <MOVE_MAP_FACTORIES>
    <MoveMapFactory name="inside_sphere" chi="false" bb="false"> # Disable all degrees of freedom globally (DOFs)
      <Chi residue_selector="neighbors"/> # Allow sidechain DOFs of neighbors
      <Backbone residue_selector="neighbors"/> # Allow backbone DOFs of neighbors
      <Jumps jump_selector="interchain_jumps"/> # Consider different chains
    </MoveMapFactory>
  </MOVE_MAP_FACTORIES>

  <TASKOPERATIONS>
    <OperateOnResidueSubset name="designable_to" selector="inner_shell">
      <RestrictAbsentCanonicalAASRLT aas="%%aa_subs%%"/> Limit design to these amino acids (all except for Cys)
    </OperateOnResidueSubset>
    <OperateOnResidueSubset name="NO_design_to" selector="NOT-designable" >
      <RestrictToRepackingRLT/> Turn off design (allows repacking and minimization)
    </OperateOnResidueSubset>
    <OperateOnResidueSubset name="NO-pack_design_to" selector="not_neighbor" >
      <PreventRepackingRLT/> Turn off design and repacking (still minimizes)
    </OperateOnResidueSubset>
    <InitializeFromCommandline name="ifcl_to" /> # Accept command-line options
    <IncludeCurrent name="incl_curr_to" /> # Include additional rotamers during repacking
    <ExtraRotamersGeneric name="ex12" ex1="1" ex2="1" extrachi_cutoff="0"/> # Include additional rotamers during repacking
    KeepSequenceSymmetry name="seq_sym_to" setting="true"/> # Only for sequence-symmetric simulations
  </TASKOPERATIONS>

  <CONSTRAINT_GENERATORS>
    <CoordinateConstraintGenerator name="not_neightbor_cc" sd="0.5" ca_only="true" residue_selector="not_neighbor" /> # Add strong backbone constraints (N, C_alpha, C, and O) for residues outside of the sphere around the target residue (residues limited by MoveMap)
    <CoordinateConstraintGenerator name="softsphere_cc" sd="2" ca_only="true" residue_selector="outer_shell" />       # Add weaker backbone constraints (N, C_alpha, C, and O) for residues in outer shell that can repack/minimize but NOT design
  </CONSTRAINT_GENERATORS>
  <MOVERS>
    <VirtualRoot name="vr" />
    <AddConstraints name="cc" constraint_generators="not_neightbor_cc" />
    <AddConstraints name="cc2" constraint_generators="softsphere_cc" />
    <FavorNativeResidue name="favor-native" bonus="%%fav_nat%%"/> # Prefer the WT residue over mutations during the design step
    <FastDesign name="design" relaxscript="MonomerDesign2019" scorefxn="sfxn_cc" task_operations="designable_to,NO_design_to,NO-pack_design_to,ex12,ifcl_to,incl_curr_to" movemap_factory="inside_sphere" />
    FastDesign name="sym_design" relaxscript="MonomerDesign2019" scorefxn="sfxn_cc" task_operations="designable_to,NO_design_to,NO-pack_design_to,ex12,ifcl_to,incl_curr_to,seq_sym_to" movemap_factory="inside_sphere" /> #  Only for sequence-symmetric simulations. Also comment out 'design' mover
    SetupForSequenceSymmetryMover name="setup_sym" sequence_symmetry_behaviour="seq_sym_to"> #  Only for sequence-symmetric simulations
    /SetupForSequenceSymmetryMover>
    <ClearConstraintsMover name="clear-cst"/>
  </MOVERS>

  <PROTOCOLS>
    <Add mover="vr"/>
    <Add mover="cc"/>
    <Add mover="cc2"/>
    <Add mover="favor-native"/>
    <Add mover="design"/>
    Add mover="setup_sym"/> #  Only for sequence-symmetric simulations
    Add mover="sym_design"/> #  Only for sequence-symmetric simulations. Also comment out 'design' mover
    <Add mover="clear-cst"/>
  </PROTOCOLS>

  <OUTPUT scorefxn="sfxn"/> # Output score without penalties from constraints (cc/cc2)

</ROSETTASCRIPTS>
'''

    # Edit the cluster RosettaScript to allow sequence symmetry to be maintained IF user provides symmetric input
    if symmetric_input != 'null':
        initial_groups = symmetric_input.split("[")
        initial_groups[1] = initial_groups[1].replace("]","")
        sym_chains=""
        for chain_id in initial_groups[0] : # Add a residueselector for each chain indicated for symmetric mutations
            xml_cluster=xml_cluster.replace('  </RESIDUE_SELECTORS>','    <Chain name="symmetry_chain_'+chain_id+'" chains="'+chain_id+'"/>\n  </RESIDUE_SELECTORS>')
            sym_chains=sym_chains+",symmetry_chain_"+chain_id
        sym_chains=sym_chains[1:]
        xml_cluster=xml_cluster.replace('KeepSequenceSymmetry','<KeepSequenceSymmetry')
        xml_cluster=xml_cluster.replace('<FastDesign name="design"','FastDesign name="design"')
        xml_cluster=xml_cluster.replace('FastDesign name="sym','<FastDesign name="sym')
        xml_cluster=xml_cluster.replace('SetupForSequenceSymmetryMover n','<SetupForSequenceSymmetryMover n')
        xml_cluster=xml_cluster.replace('SequenceSymmetry residue_selectors="symmetry_chains"/>','<SequenceSymmetry residue_selectors="symmetry_chains"/>')
        xml_cluster=xml_cluster.replace('/SetupForSequence','    <SequenceSymmetry residue_selectors="'+sym_chains+'"/>\n    /SetupForSequence')
        xml_cluster=xml_cluster.replace('/SetupForSequence','</SetupForSequence')
        xml_cluster=xml_cluster.replace('Add mover="setup_sym','<Add mover="setup_sym')
        xml_cluster=xml_cluster.replace('Add mover="sym_design','<Add mover="sym_design')
        xml_cluster=xml_cluster.replace('<Add mover="design','Add mover="design')

    # Initial relax with AtomTree minimization (repeats=1). This step is simply to process the pdb file (add hydrogens) so that -use_truncated_termini can be turned on without improperly treating the real N and C-termini.
    xml_file = 'rs_relax_cluster.xml'
    with open (xml_file, 'w') as xml_file :
        xml_file.write( xml_relax_cluster )
    
    # Third RosettaScript for Site Saturation Mutagenesis (repeats=3). Note that this requires a jobdefinition file
    xml_file2 = 'rs_cluster.xml'
    with open (xml_file2, 'w') as xml_file2:
        xml_file2.write( xml_cluster )

############# This section (~110 lines) is for parsing the PDB file and was extracted from rosetta/tools/pdb_structure.py, written by Andrew Leaver-Fay. #############

# some foundation classes for representing the contents of a PDB
def is_integer( s ) :
    try :
        int(s)
        return True
    except :
        return False

class Residue :
    def __init__( self ) :
        self.resstring = "" # residue index + insertion code -- a string
        self.resname = ""
        self.insertion_code = ""
        self.chain = None #pointer to the containing chain
    # insert this residue into a Chain
    def claim( self, containing_chain ) :
        self.chain = containing_chain
    #return the chain-resstring tuple for this residue that uniquely identifies it
    #Requires that the residue has already been inserted into a chain
    def resid( self ) :
        assert( self.chain )
        return self.chain.chain_name + " " + self.resstring

class Chain :
    def __init__( self ) :
        self.chain_name = ""
        self.residues = []
        self.resmap = {}
    def add_residue( self, residue ):
        self.residues.append( residue )
        self.resmap[ residue.resstring ] = residue
        residue.claim( self )
    def replace_residue( self, newresidue ):
        # in place replacement; keep the original location of the residue in the self.residues array
        copyres = copy.copy( newresidue )
        copyres.claim( self )
        if newresidue.resstring not in self.resmap :
            print("Could not replace residue", newresidue.resstring)
            print(len( self.resmap ))
        assert( newresidue.resstring in self.resmap )
        for i in range(len( self.residues )) :
            if self.residues[ i ].resstring == newresidue.resstring :
                self.residues[ i ] = copyres
        self.resmap[ newresidue.resstring ] = copyres
    def residue( self, resstring ) :
        return self.resmap[ resstring ]

class PDBStructure :
    def __init__( self ):
        self.chains = []
        self.chainmap = {}

    def residue( self, chnm, resstring ) :
        return self.chainmap[ chnm ].residue( resstring )

    def add_chain( self, chain ) :
        self.chains.append( chain )
        self.chainmap[ chain.chain_name ] = chain

    def read_from_lines( self, lines ):
        last_chain = Chain()
        last_residue = Residue()
        for line in lines :
            if line[0:4] == "ATOM" or line[0:6] == "HETATM" :
                chnm = self.chain_name_from_pdbline( line )
                if last_chain.chain_name != "" and chnm != last_chain.chain_name :
                    last_chain.add_residue( last_residue )
                    last_residue = Residue()
                    if last_chain.chain_name not in self.chainmap : self.add_chain( last_chain )
                    if chnm not in self.chainmap :
                        last_chain = Chain()
                    else :
                        # this chain has already been added, but now we have more residues
                        last_chain = self.chainmap[ chnm ]
                if last_chain.chain_name == "":
                    last_chain.chain_name = chnm
                resstring = self.resstring_from_pdbline( line )
                if last_residue.resname != "" and last_residue.resstring != resstring :
                    last_chain.add_residue( last_residue )
                    last_residue = Residue()
                if last_residue.resname == "" :
                    last_residue.resname = self.resname_from_pdbline( line )
                    last_residue.resstring = resstring
                    #print "Read residue", last_residue.resname, last_chain.chain_name, last_residue.resstring
        if last_residue.resname != "":
            last_chain.add_residue( last_residue )
        if last_chain.chain_name != "" and last_chain.chain_name not in self.chainmap :
            self.add_chain( last_chain )
    def pdb_resname_range( self ) :    return ( 17, 20 )
    def pdb_chain_name_range( self ) : return ( 21, 22 )
    def pdb_resstring_range( self ) :  return ( 22, 27 )
    def pdb_occupancy_range( self ) :  return ( 56, 60 )

    def resname_from_pdbline( self, line ):
        return line[ self.pdb_resname_range()[0]:self.pdb_resname_range()[1] ]
    def chain_name_from_pdbline( self, line ):
        return line[ self.pdb_chain_name_range()[0]:self.pdb_chain_name_range()[1] ]
    def resstring_from_pdbline( self, line ):
        return line[ self.pdb_resstring_range()[0]:self.pdb_resstring_range()[1] ].strip()

def pdbstructure_from_file( fname ) :
    pdb = PDBStructure()
    with open( fname ) as f: pdb.read_from_lines( f.readlines() )
    return pdb

############# End of excerpt from rosetta/tools/pdb_structure.py #############

# The user input format was changed from the first version. This function converts the new format (A[2,3-6]) to the old format (A2,A3-A6) to maintain future residue parsing functions
def parse_brackets_input_resid(bracketed_input):
    # Split input "A[257,258-260]" into "A257,A258-A260"
    chains=[]
    # Pull string preceding '['
    for out_brackets in bracketed_input.split('['):
        chain=out_brackets.split(']')[-1]
        chain=chain.replace(",","")
        chains.append(chain)
    del chains[-1]

    positions=[]

    # Pull string within brackets ('[]')
    for in_brackets in bracketed_input.split('['):
        position=in_brackets.split(']')[0]
#        print(position)
        positions.append(position)
    del positions[0]
#    print(len(positions))

    # Apply chain ID to each number inside of brackets
    combined_list=""
    c=0
    while c < len(chains):
        test=(chains[c]+positions[c])
        test=test.replace(",",","+chains[c])
        test=test.replace("-","-"+chains[c])
        combined_list=combined_list+','+test
        c=c+1
#        print(test)
    combined_list=combined_list[1:]
    return combined_list

def parse_input_resid(pdb_file, input_residues):
    '''Convert user-provided numbers to list of individual residues.
    input_residues should be in format "L30,H40-H42A" where first
    char is chain_id and the remainder is residue number. H42A
    refers to an insertion code for residue 42A on chain H.
    Output is a list of tuples with chain_id, PDB_resid, resname,
    and index, which should be a proxy for Pose Numbering
    i.e. [('L', '30', 'TYR', '1'), ('H', '40', 'ALA', '225')]'''

    pdb=pdbstructure_from_file(pdb_file)
# Create a list of every residue within the PDB. Each tuple contains (chain_id, PDB_Number, ResName, Pose_Number) where Pose Numbering is assumed to be index+1 (starts from 0 in python, 1 in Rosetta)
    all_res_list=[]
    for ch in pdb.chains :
        for res in ch.residues :
            index=len(all_res_list)
            pose_num=len(all_res_list)+1
            #pose_num=index+1
            all_res_list.append( (ch.chain_name, res.resstring, res.resname, str(pose_num) ) )
    mutation_list = []
    mut_groups = re.split(",",input_residues)
    for mg in mut_groups:
        '''iterate values over range (Check PDB for insertion codes, i.e. H20-H21A = [H20, H21, H21A]'''
        if "-" in mg :
            bounds = re.split("-",mg)
            start_resid_full=bounds[0]
            end_resid_full=bounds[1]
            chain_id=start_resid_full[0]
            start_resid=start_resid_full[1:]
            end_resid=end_resid_full[1:]
            # Identify the index value of start and end residues
            index_start=[i for i, v in enumerate(all_res_list) if ( (v[0] == chain_id) and (v[1] == start_resid) )].pop() # Should be pose numbering
            index_end=[i for i, v in enumerate(all_res_list) if ( (v[0] == chain_id) and (v[1] == end_resid) )].pop() # Should be pose numbering
            # Create a list of tuples for residues within that range
            mutation_range=[]
            x=index_start
            while (x <= index_end):
                mutation_list.append( ( all_res_list[x] ) )
                x = x + 1
        else : # Parse single residue
            chain_id=mg[0]
            resid=mg[1:]
            index=[i for i, v in enumerate(all_res_list) if ( (v[0] == chain_id) and (v[1] == resid) )].pop() # Should be pose numbering
            resname=( pdb.chainmap[ chain_id ].resmap[ resid ].resname )
            mutation_list.append( ( all_res_list[index] ) )
    num_res=len(all_res_list)
    for residue_tuple in mutation_list:
        for obj in residue_tuple:
            assert obj.isalnum()
    return mutation_list, num_res

### Optional ###
def parse_symmetric_input(symmetric_input):
    resid_list = ""
    initial_groups = symmetric_input.split("[")
    initial_groups[1] = initial_groups[1].replace("]","")
    for chain_id in initial_groups[0] :
        '''parse residues and add chain_id in front'''
        mut_groups = re.split(",",initial_groups[1])
        for mg in mut_groups:
            '''treat a range differently (i.e. 20-22 = [A20-A22]'''
            if "-" in mg:
                bounds = re.split("-",mg)
                lower=chain_id+bounds[0]
                upper=chain_id+bounds[1]
                resid_list=resid_list + ',' + lower + '-' + upper
            else:
                resid_list=resid_list + ',' + chain_id + mg
    return resid_list[1:] # The previous code starts resid_list with a ',' so this returns everything except for that character

''' Create a new list of tuples that matches mutation_list EXCEPT that matching residues from different chains have been combined so that the new tuple is [(chain_ids),(PDB_resid_num),(PDB_resname),(indexes)]'''
def symmetric_grouped_mutations(symmetric_input,mutation_list):
    symmetric_resid_list = []
    initial_groups = symmetric_input.split("[")
    initial_groups[1] = initial_groups[1].replace("]","")
    chains=initial_groups[0]
    # Get a list of residue (PDB numbering)
    all_pdb_num=[]
    for residue_tuple in mutation_list:
        if residue_tuple[0] == chains[0] :
            all_pdb_num.append(residue_tuple[1])
    # Create a new list of tuples that matches mutation_list EXCEPT that matching residues from different chains have been combined so that the new tuple is [(chain_ids),(PDB_resid_num),(PDB_resname),(indexes)]
    all_clustered_indexes=[]
    symmetric_mutation_list=[]
    for each_residue in all_pdb_num :
        clustered_index=""
        resname_three_lett=""
        for residue_tuple in mutation_list:
            if residue_tuple[1] == each_residue :
                clustered_index=clustered_index + ',' + residue_tuple[3]
                resname_three_lett=residue_tuple[2]
        all_clustered_indexes.append( (chains, each_residue, resname_three_lett, clustered_index[1:]) )
    return all_clustered_indexes

def parse_fav_native_input(fav_native):
   return [float(v) for v in fav_native.split(',')]

# Default to only keeping protein residues (ATOM), but if ligands should be kept, check HETATM lines for specific residues. Only the first two elements of each tuple will be used from keep_ligand_list. That list can be created by using parse_input_resid with input_ligand_residues
def trim_pdb(pdb_file, keep_ligand_list):
    with open('trimmed_input_protein.pdb', 'w') as trimmed_pdb:
        with open( pdb_file ) as f:
            for line in ( f.read() ).split('\n'):
                if line[0:4] == "ATOM" :
                     trimmed_pdb.write(line + '\n')
                elif line[0:6] == "HETATM" and ('MSE' in line[17:20]):
                     trimmed_pdb.write(line + '\n')
                elif keep_ligand_list is not None and line[0:6] == "HETATM" :
                    for keep_ligand in keep_ligand_list :
                        chain_id=keep_ligand[0]
                        resid=keep_ligand[1]
                        if line[21] == chain_id and (keep_ligand[1] in line[22:27]):
                            trimmed_pdb.write(line + '\n')

# Convert amino acid code from 1 letter to 3 letter #
dict_321 = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
           'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
           'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
           'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}

dict_123 = {v: k for k,v in list(dict_321.items())}

def aa_3_to_1(x):
    if len(x) % 3 != 0:
        raise ValueError('Input length should be a multiple of three')
    y = ''
    for i in range(len(x)//3):
            y += dict_321[x[3*i:3*i+3]]
    return y

def aa_1_to_3(x):
    if len(x) % 1 != 0:
        raise ValueError('Input length should be a multiple of three')
    y = ''
    for i in range(len(x)//1):
            y += dict_123[x]
    return y

### Build directory tree & create jobdefinition file for site saturation mutagenesis ###

# resname refers to the residues that the amino acid will be replaced with. 
def build_jobdef_directories(list_mutations, aa_subs_option, fav_nat_list, keep_native_list, num_res):
    # Create comma-separated string with POSE residue numbering
    keep_native_string=""
    for start in keep_native_list:
        keep_native_string=keep_native_string+","+start[3]
    keep_native_string=keep_native_string[1:]
    xml_filename = 'jobdef.xml'
    with open (xml_filename, 'w') as xml_file: 
        xml_file.write("<JobDefinitionFile> \n")
    for current_pos in list_mutations:
        chain_id = current_pos[0]
        pdb_resid = current_pos[1]
        pose_resid = current_pos[3]
        position = chain_id + '[' + pdb_resid + ']'
        for fav_value in fav_nat_list:
        # Create target directory if it doesn't exist (seed position with PDB naming) 
             dirName = '{}/{}_bonus'.format(position,fav_value)
             if not os.path.exists(dirName):
                 os.makedirs(dirName)
#                 print(("Directory " , dirName ,  " Created "))
             else:    
                 print(("Directory " , dirName ,  " already exists"))

             jobdef_template = '''\
  <Job>
    <Output>
      <PDB path="{position}/{fav_nat_dir}"/>
    </Output>
    <SecondaryOutput>
      <ScoreFile path="{position}/{fav_nat_dir}" filename="score.sc"/>
    </SecondaryOutput>
    <Options>
      <parser__script_vars value="fav_nat={fav_nat} seed={pose_resid} keep={keep_nat} aa_subs={aa_subs_option}"/>
    </Options>
  </Job>
'''
             # Create jobdefinition XML file for RosettaScripts
             with open (xml_filename, 'a') as xml_file: 
                 xml_file.write( jobdef_template.format(
                     position=position,
                     pose_resid=pose_resid,
                     fav_nat_dir=str(fav_value)+"_bonus", # fav_nat_dir refers to the directory name for fav_nat. It's the same in this case, but is different for the "WT" sim where nothing is mutated.
                     fav_nat=fav_value,
                     keep_nat=keep_native_string,
                     aa_subs_option=aa_subs_option ) )
        ### BLOCK all mutations, this will represent a "WT" simulation where residues repack around the seed
        dirName = '{}/wt_sim'.format(position)
        if not os.path.exists(dirName):
            os.makedirs(dirName)
#            print(("Directory " , dirName ,  " Created "))
        else:
            print(("Directory " , dirName ,  " already exists"))
        # Create jobdefinition XML file for RosettaScripts
        with open (xml_filename, 'a') as xml_file:
            xml_file.write( jobdef_template.format(
                position=position,
                pose_resid=pose_resid,
                fav_nat_dir="wt_sim",
                fav_nat="0",
                keep_nat='1-{}'.format(num_res), # prevent design on all residues
                aa_subs_option=aa_subs_option ) )

    with open (xml_filename, 'a') as xml_file: 
        xml_file.write("</JobDefinitionFile> \n")

#### Analysis functions ####

### Running MPI creates multiple score files so combine into one before further analysis. Assumes score files are in current directory
def cat_scores(out_dir, score_file_prefix):
    # Generate a list of score files that match the score_file_prefix
    filenames = [ f for f in os.listdir(out_dir) if f.startswith(score_file_prefix)]
    # Combine all of them into one file
    with open(out_dir + 'mutant_scores.txt', 'w') as outfile:
        with open (out_dir + filenames[0], 'r') as infile:
            outfile.writelines(infile.readlines()[:2]) # Write the header from the first file (first two lines)
    with open(out_dir + 'mutant_scores.txt', 'a') as outfile:
        for fname in filenames:
            with open(out_dir + fname, 'r') as infile:
                outfile.writelines(infile.readlines()[2:]) # Skip the headers (first two lines)

# Convert numbers to floats for parsing the score file; however some inputs are meant to be text (pdb names) so the input is returned either way
def convertFloat(v):
    try: return float(v)
    except ValueError: return v

def tryFloat(v):
    try: float(v)
    except ValueError: return 0.0

def get_best_model_and_score(score_file):
    ''' select model with lowest total_score '''
    top_n = 1

    score_table = []
    score_table_keys = []
    with open( score_file ) as f:
        lines = ( f.read() ).split('\n')
        score_table_keys = lines[1].split()[1:]  # ignoring first line, ignore first key (it just 'SCORE:')
        if not score_table_keys:
            print('Score table is empty!')
            sys.exit(1)
        score_table_keys[-1] = 'decoy'  # changing description → decoy
        lines = lines[:1] + sorted([l for l in lines[1:] if l and l.split()[0] == 'SCORE:' and l.split()[1] != 'total_score'], key=lambda l: float(l.split()[1]), )   # sorting lines by score, split based on total_score (l.split()[1] but skip the 1st line (Sequence:) and the header line ('total_score')
        score_table = [ dict( list(zip(score_table_keys, list(map(convertFloat, l.split()[1:])))) )  for l in lines[1:]]
        score_table_keys = score_table_keys[-1:] + score_table_keys[:-1]  # making 'decoy' the first column
        best_model = [ m['decoy'] for m in score_table[:top_n] if type(m['decoy'])==str]
    # FROM ROSIEv1:    best_model = [ m['decoy'] for m in score_table[:top_n] if type(m['decoy'])==str  and  '_mut_' in m['decoy'] ]  # we need to check if '_mut_' is in decoy name to avoid rare cases when line is misformed (HPC writing issue?)
        best_score = [ m['total_score'] for m in score_table[:top_n] if type(m['total_score'])==float ]
        return best_model, best_score

def compare_sequences(wt_model,designed_model):

# Create a list of every residue within the WT model (repacked, no design). Each tuple contains (chain_id, PDB_Number, ResName, Pose_Number) where Pose Numbering is assumed to be index+1 (starts from 0 in python, 1 in Rosetta)
    wt_pdb=pdbstructure_from_file(wt_model)
    all_res_list_wt=[]
    for ch in wt_pdb.chains :
        for res in ch.residues :
            index=len(all_res_list_wt)
            pose_num=len(all_res_list_wt)+1
            all_res_list_wt.append( (ch.chain_name, res.resstring, res.resname, str(pose_num) ) )
    # Generate wt_fasta sequence
#    print(start_chain)
    wt_fasta=""
    for aa in all_res_list_wt :
        current_chain=aa[0]
        wt_fasta=wt_fasta+aa_3_to_1(aa[2])
#    print(wt_fasta)

    # Create a list of every residue for the designed model
    design_pdb=pdbstructure_from_file(designed_model)
    all_res_list_design=[]
    for ch in design_pdb.chains :
        for res in ch.residues :
            index=len(all_res_list_design)
            pose_num=len(all_res_list_design)+1
            all_res_list_design.append( (ch.chain_name, res.resstring, res.resname, str(pose_num) ) )
#    print(all_res_list_design)

    designed_fasta=""
    for aa in all_res_list_design :
        current_chain=aa[0]
        designed_fasta=designed_fasta+aa_3_to_1(aa[2])

    fasta_diff=[i for i in range(len(wt_fasta)) if wt_fasta[i] != designed_fasta[i]] # print the index (starts at 0, pose_num = this +1) for each residue that is different within the two FASTA sequences
#    print(fasta_diff)
    if len(fasta_diff) != 0 : # As long as there are mutations between the two files, collect more information for each mutation (chain_id, PDB_number, ResName, Pose_Number)
        new_mut_list_pdb_num=""
        new_mut_list_pose_num=""
        new_mut_list_pymol_pdb_num=""
        num_mut=0

        for pos in fasta_diff : # pos = the index (i.e. pose_num) for each mutation
#            print(all_res_list_wt[pos])
            pos_wt_info=all_res_list_wt[pos]
            pos_des_info=all_res_list_design[pos]
            resname_wt=aa_3_to_1(pos_wt_info[2]) # Convert to single letter aa abbreviation
            chain_id=pos_des_info[0]
            pdb_num=pos_des_info[1]
            resname_des=aa_3_to_1(pos_des_info[2]) # Convert to single letter aa abbreviation
            pose_num=pos_des_info[3]

            new_mut_list_pdb_num=new_mut_list_pdb_num+chain_id+'['+resname_wt+pdb_num+resname_des+']'+' '
            new_mut_list_pose_num=new_mut_list_pose_num+chain_id+'['+resname_wt+pose_num+resname_des+']'+' '
            new_mut_list_pymol_pdb_num=new_mut_list_pymol_pdb_num+pdb_num+'+'
            num_mut=num_mut+1

        # Remove the trailing '_' from constantly appending it in the for loop
        new_mut_list_pdb_num=new_mut_list_pdb_num[:-1]
        new_mut_list_pose_num=new_mut_list_pose_num[:-1]
        new_mut_list_pymol_pdb_num=new_mut_list_pymol_pdb_num[:-1]
        # Create a list that contains 1) comma-separated string of the mutations in pdb_numbering, 2) same in pose_numbering, and 3) the number of mutations between these two models (i.e. ['A_E257D,A_G258E,A_M260D,A_H261E,', 'A_E1D,A_G2E,A_M4D,A_H5E,', '4'])
        # Generate fasta for designed model in fasta format (separate chains)
        start_chain=all_res_list_design[0]
        start_chain=start_chain[0]
        designed_fasta='> Chain '+start_chain+'\n'
        for aa in all_res_list_design :
            current_chain=aa[0]
            if (current_chain == start_chain):
                designed_fasta=designed_fasta+aa_3_to_1(aa[2])
            else:
                start_chain=current_chain
                designed_fasta=designed_fasta+'\n> Chain '+start_chain+'\n'+aa_3_to_1(aa[2])

        new_mut_list=[new_mut_list_pdb_num,new_mut_list_pose_num,str(num_mut),designed_fasta,new_mut_list_pymol_pdb_num]
#        print(new_mut_list)
        return(new_mut_list)
    else :
        new_mut_list=["none","none","0","matches_WT","none"]
#        print(new_mut_list)
        return new_mut_list

### Collect best models and create table of results (Total score of Mutant - WT)
def analysis(list_mutations, fav_nat_list, score_file_prefix):
    # Create directory for PDBs
    results_dir = './results/'
    pdb_models_dir = '{}pdb_models/'.format(results_dir)
    if not os.path.exists(results_dir):
        os.makedirs(results_dir)
        os.makedirs(pdb_models_dir)
#        print(("Directory " , results_dir ,  " Created "))
    else:
        print(("Directory " , results_dir ,  " already exists, placing results in existing directory."))

    open(results_dir + 'fasta.txt', 'w').close() # overwrite if existing so we can append later without worry
    with open(results_dir + 'results.csv', 'w') as results:
        results.write("Seed,Native Bonus,Delta Energy / Num Mutations,Num Mutations,Total Score,Mutations (PDB Numbering), Pymol Selector of Mutations (PDB Numbering), Mutations (Pose Numbering),Model Name")
        ''' 1) determine the score for the WT residue, 2) Subtract each mutant score from the WT score, move the best scoring models to the results_dir, and rename seed$x_nat$y.pdb where $x & $y are the seed_residue and fav_nat_bonus, respectively'''
        for current_pos in list_mutations:
            chain_id = current_pos[0]
            resid = current_pos[1]
            position = chain_id + '[' + resid + ']'

            ''' Determine the score of the WT residue at the given position '''
            out_dir = '{}/wt_sim/'.format(position)
            # MPI generates multiple score files, so combine them into a single file
            cat_scores(out_dir, score_file_prefix) # This generates the file named 'mutant_scores.txt' for WT residue
            score_file = out_dir + 'mutant_scores.txt'
            WT_best_model, WT_best_score = get_best_model_and_score(score_file)
            # Copy the best WT pdb file to the results directory and rename seed$x_wt.pdb
            wt_name='seed_{}_wt.pdb'.format(position)
            best_decoy_wt = out_dir + WT_best_model[0] + '.pdb'
            shutil.copy(best_decoy_wt, pdb_models_dir + wt_name )

            for fav_value in fav_nat_list:
                out_dir = '{}/{}_bonus/'.format(position,fav_value)
                '''Collect the best model, compare to WT, and calculate scores'''
                cat_scores(out_dir, score_file_prefix) # This generates the file named 'mutant_scores.txt'
                score_file = out_dir + 'mutant_scores.txt'
                mut_best_model, mut_best_score = get_best_model_and_score(score_file)
                best_decoy_des = out_dir + mut_best_model[0] + '.pdb'
                delta_mut_score = mut_best_score[0] - WT_best_score[0]
#                print(("Mut_Best_Score"), mut_best_score[0], ("WT_Best_Score"), WT_best_score[0])
                comparison=compare_sequences(best_decoy_wt,best_decoy_des) # Returnes a tuple with a list of mutations in PDB numbering (0), pose numbering (1), the number of mutations (2), the FASTA of the design (3), and residue selectors for PyMol in PDB numbering (4)
                if comparison[2] == '0' : # Skip if there are no mutations in the designed model
                    score_per_mut = 0
                else :
                    score_per_mut = ( delta_mut_score / float(comparison[2]) )
                mutant_name='seed_{}_nat{}.pdb'.format(position, fav_value)
                escaped_mutant_name='seed_'+ chain_id + '\[' + resid + '\]'+ '_nat' + str(fav_value) + '.pdb'
                shutil.copy(best_decoy_des, pdb_models_dir + mutant_name )
                results.write('\n' + position + ',' + str(fav_value) + ',' + str(format(score_per_mut,'.2f')) + ',' + comparison[2] + ',' + str(format(mut_best_score[0],'.2f')) + ',' + comparison[0] + ',' + comparison[4] + ',' + comparison[1] + ',' + escaped_mutant_name)
                with open(results_dir + 'fasta.txt', 'a') as fasta:
                    fasta.write('\n' + '>>> ' + mutant_name + '\n' + comparison[3])


################################# Validators #################################

def input_position_validator(string_input):
    allowed_chars = set('[]ABCDEFGHIJKLMNOPQRSTUVWXYZ-,:1234567890')
    for c in string_input:
        if c not in allowed_chars:
            error = ' illegal characters, should only be a range of positions (i.e. H[1,10-20],L[23]).'
            return error

def upper_alpha_validator(string_input):
    # Must be a string with only capitalized letters
    input_string = string_input
    check_letters = input_string.isalpha()
    check_upper = input_string.isupper()

    if not (check_letters and check_upper):
        error = ' should only contain upper-case letters.'
        return error

def aa_validator(string_input):
    # Must be a string with only amino acid letters
    for c in string_input:
        if c not in 'ACDEFGHIKLMNPQRSTVWY':
            print('The amino acid substitutions should only contain upper-case letters that represent amino acids (ACDEFGHIKLMNPQRSTVWY).')
            sys.exit(1)

def symmetric_validator(symmetric_input):
    try:
        initial_groups = symmetric_input.split("[")
        initial_groups[1] = initial_groups[1].replace("]","")
        error = upper_alpha_validator(initial_groups[0])
        if error:
            error = 'The input for symmetric mutations is incorrect. Outside of the brackets' + error
            return error
        if len(initial_groups[0]) <= 1:
            error = 'The input for symmetric mutations must designate more than one chain.'
            return error
        try:
            symmetric_mutations = parse_symmetric_input(symmetric_input)
        except:
            error = 'Error parsing input for symmetric mutations.' # This is a catch-all for other input errors
            return error
    except:
        error = 'Error with symmetric input. Should be Chain_ID[Residue_Number], i.e. AB[257]'
        return error

def check_path_validator(filename):
    assert '..' not in filename # Avoid malicious redirects #ROSIE-only
    allowed_chars = set('ABCDEFGHIJKLMNOPQRSTUVWXYZ-abcdefghijklmnopqrstuvwxyz_1234567890.')
    for c in filename:
        if c not in allowed_chars:
            print('Invalid path to file.')
            sys.exit(1)

    if os.path.isfile('./' + filename):
        pass
    else:
        error = 'Error locating PDB file: ./' + filename + '\n Is the filename correct?' # This is a catch-all in case the PDB file doesn't exist
        print(error)
        sys.exit(1)

def validate_build_and_analyze(pdb_file, input_resnum, aa_subs_option, fav_native, symmetric_input, keep_native, keep_native_sym):

    check_path_validator(pdb_file)

    if input_resnum != 'null' :
        error = input_position_validator(input_resnum)
        if error:
            error = 'The input positions for mutagenesis contains' + error
            print(error)
            sys.exit(1)
        # Check the PDB file for the specified positions
        try:
            input_resnum=parse_brackets_input_resid(input_resnum)
            list_mutations, num_res=parse_input_resid(pdb_file, input_resnum)
        except:
            error = 'One of the positions was not found in the provided PDB file. Please check your input, you may be using an incorrect numbering scheme.'
            print(error)
            sys.exit(1)
    else :
        input_resnum = ''

    # Optional, user decides which aa substitutions to make
    if aa_subs_option != 'null' :
        aa_validator(aa_subs_option) # Format check
    else :
        aa_subs_option = ''

    if fav_native != 'null':
        allowed_chars = set(',.1234567890')
        for c in fav_native:
            if c not in allowed_chars:
                print('Only ",.1234567890" characters are allowed for the favor_native input.')
                sys.exit(1)
        if (fav_native[-1] == '.') or (fav_native[-1] == ',') or (fav_native[0] == '.') or (fav_native[0] == ','):
            print('-fav_nat cannot begin or end with "," or ".".')
            sys.exit(1)
        try:
             fav_nat_list=parse_fav_native_input(fav_native)
        except:
            error = 'Please check your input for -fav, only floats separated by commas are allowed (i.e. "0,1.5,3")'
            print(error)
            sys.exit(1)
    else:
        fav_nat_list = [ "0.0" ]
    
    # Optional, user decides which aa substitutions are allowed during design
    if keep_native != 'null':
        error = input_position_validator(keep_native)
        if error:
            error = 'The input positions for mutagenesis contains' + error
            print(error)
            sys.exit(1)
        try:
            keep_native=parse_brackets_input_resid(keep_native)
            keep_native_list, num_res = parse_input_resid(pdb_file, keep_native)
        except:
            error = 'One of the positions for -keep was not found in the provided PDB file. Please check your input, you may be using an incorrect numbering scheme.'
            print(error)
            sys.exit(1)
    else:
        keep_native_list = [('999999', '999999', '999999', '999999')] # This residue is assumed to not be in the pdb file so no residues will be skipped with this option
    
    # Optional, symmetric mutations
    if symmetric_input != 'null' :
        error = input_position_validator(symmetric_input)
        if error:
            error = 'The input for symmetric mutations contains' + error
            print(error)
            sys.exit(1)
        error = symmetric_validator(symmetric_input) # Format check
        if error:
            print(error)
            sys.exit(1)
        # Check that the user-provided input is in the proper format AND each residue is present within the PDB file
        sym_mutation_list=parse_symmetric_input(symmetric_input)
        try:
            mutation_list, num_res=parse_input_resid(pdb_file, sym_mutation_list)
        except:
            error = 'One of the positions in the symmetric input was not found in the provided PDB file. Please check your input, you may be using an incorrect numbering scheme.'
            print(error)
            sys.exit(1)
    else :
        symmetric_input = ''

    # Optional, keep_native but with symmetric input (i.e.AB[257-261])
    if keep_native_sym != 'null':
        error = input_position_validator(keep_native_sym)
        if error:
            error = 'The input for Keep Native: Symmetric Input contains' + error
            print(error)
            sys.exit(1)
        error = symmetric_validator(keep_native_sym) # Format check
        if error is not None:
            print(error)
            sys.exit(1)
        keep_sym_mutation_list=parse_symmetric_input(keep_native_sym)
        # Check that the user-provided input is in the proper format AND each residue is present within the PDB file
        try:
            keep_native_list, num_res = parse_input_resid(pdb_file, keep_sym_mutation_list)
        except:
            error = 'One of the positions in the -keep_sym input was not found in the provided PDB file. Please check your input, you may be using an incorrect numbering scheme.'
            print(error)
            sys.exit(1)

def parse_build_analyze_inputs(pdb_file, input_resnum, aa_subs_option, fav_native, symmetric_input, keep_native, keep_native_sym):
    # Optional, user decides which aa substitutions to make
    if aa_subs_option != 'null':
        resname = aa_subs_option
    else:
        # Can exclude mutations for saturation mutagenesis here (i.e. to save speed, skip Pro, etc.)
        resname = "ADEFGHIKLMNPQRSTVWY"
    if fav_native != 'null':
        fav_nat_list=parse_fav_native_input(fav_native)
    else:
        fav_nat_list = [ "0.0" ]
    # Optional, user decides which aa substitutions are allowed during design
    if keep_native != 'null':
        keep_native_list, num_res = parse_input_resid(pdb_file, keep_native)
    else:
        keep_native_list = [('999999', '999999', '999999', '999999')] # This residue is assumed to not be in the pdb file so no residues will be skipped with this option
    # Optional, symmetric mutations
    if symmetric_input != 'null':
#        sym_mutation_list=parse_symmetric_input(symmetric_input) # delete?
#        mutation_list, num_res=parse_input_resid(pdb_file, sym_mutation_list) # delete?
        sym_mutation_string=parse_symmetric_input(symmetric_input) #(i.e. AB:257-261 becomes A257-A261,B257-B261)
        sym_mutation_list, num_res=parse_input_resid(pdb_file, sym_mutation_string) # (i.e. complete tuple list of all residues on all chains(
        list_mutations=symmetric_grouped_mutations(symmetric_input,sym_mutation_list)
    # If symmetric_input is off, check input_mutations for proper formatting AND whether each residue is present within the PDB file
    else:
        input_resnum=parse_brackets_input_resid(input_resnum)
        list_mutations, num_res=parse_input_resid(pdb_file, input_resnum)

    if keep_native_sym != 'null':
        keep_native_sym_mutation_string=parse_symmetric_input(keep_native_sym) #(i.e. AB[257-261] becomes A257-A261,B257-B261)
        keep_native_sym_mutation_list, num_res = parse_input_resid(pdb_file, keep_native_sym_mutation_string) # (i.e. complete tuple list of all residues on all chains(
        keep_native_list=symmetric_grouped_mutations(keep_native_sym,keep_native_sym_mutation_list)

    for residue_tuple in list_mutations:
        for obj in residue_tuple:
            assert obj.replace(',', '').isalnum() 
    for letter in resname:
        assert letter.isalpha()
    return list_mutations, num_res, resname, fav_nat_list, keep_native_list

################################# Primary Functions ########################################

##### Mode 1: Prepare pdb file for relax by removing all hetero atoms except for any specified by ligand_option #####
def trim(pdb_file, ligand_option):
    check_path_validator(pdb_file)
    if ligand_option != 'null' :
        error=input_position_validator(ligand_option)
        if error:
            error = 'The ligand position contains' + error
            print(error)
            sys.exit(1)
        # Check the PDB file for the specified positions
        try:
            ligand_option=parse_brackets_input_resid(ligand_option)
            keep_ligand_list, num_res = parse_input_resid(pdb_file, ligand_option)
        except:
             error = 'One of the LIGAND positions was not found in the provided PDB file. Please check your input, you may be using an incorrect numbering scheme.'
             print(error)
             sys.exit(1)
    else :
        keep_ligand_list = ""

    trim_pdb(pdb_file, keep_ligand_list) # Creates the file 'trimmed_input_protein.pdb' in the current directory

##### Mode 2: Generate RosettaScripts (relax and cluster). The scripts are modified based on the symmetric_inputs so must be performed after the protein is trimmed (it uses Rosetta numbering so removing residues would change relative positions). #####
def RScripts(pdb_file, input_resnum, aa_subs_option, fav_native, symmetric_input, keep_native, keep_native_sym):
    validate_build_and_analyze(pdb_file, input_resnum, aa_subs_option, fav_native, symmetric_input, keep_native, keep_native_sym)
    gen_RScripts(symmetric_input)

##### Mode 3: Build SSM inputs (directory tree and jobdefinition file) #####
def build(pdb_file, input_resnum, aa_subs_option, fav_native, symmetric_input, keep_native, keep_native_sym):
    # To avoid code duplication, a separate function for validating both build and analyze was created. The variables are all returned so they can be passed to run_ssm consistently
    validate_build_and_analyze(pdb_file, input_resnum, aa_subs_option, fav_native, symmetric_input, keep_native, keep_native_sym)
    gen_RScripts(symmetric_input)
    list_mutations, num_res, parsed_aa_subs_option, fav_nat_list, keep_native_list = parse_build_analyze_inputs(pdb_file, input_resnum, aa_subs_option, fav_native, symmetric_input, keep_native, keep_native_sym)
    build_jobdef_directories(list_mutations, parsed_aa_subs_option, fav_nat_list, keep_native_list, num_res)

##### Mode 4: Analyze results (calculate delta_scores, collect best models, and generate heatmap) #####
def analyze(pdb_file, input_resnum, aa_subs_option, fav_native, symmetric_input, keep_native, keep_native_sym):
    score_file_prefix = "score"
    validate_build_and_analyze(pdb_file, input_resnum, aa_subs_option, fav_native, symmetric_input, keep_native, keep_native_sym)
    list_mutations, num_res, resname, fav_nat_list, keep_native_list = parse_build_analyze_inputs(pdb_file, input_resnum, aa_subs_option, fav_native, symmetric_input, keep_native, keep_native_sym)
    analysis(list_mutations, fav_nat_list, score_file_prefix)

################################# Execute Functions ########################################

def register():
    return {
        'stabilize-cluster.trim' : trim,
        'stabilize-cluster.RScripts' : RScripts,
        'stabilize-cluster.build' : build,
        'stabilize-cluster.analyze' : analyze,
    }

