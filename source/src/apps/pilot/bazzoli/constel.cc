// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @brief Extracts constellations of contiguous atoms from a protein structure.
/// @author jk
/// @author Andrea Bazzoli (bazzoli@ku.edu)

///////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                           //
//                                          constel                                          //
//                                                                                           //
// ARGUMENTS:                                                                                //
//  ---> "-s PDBFIL", where PDBFIL is the path to the PDB file containing the protein        //
//                    structure.                                                             //
//  ---> "-ignore_unrecognized_res", to ignore residues of unrecognized amino acid type.     //
//  ---> "-suppress_zero_occ_pdb_output", to print only atoms with nonzero occupancy.        //
//  ---> "-ignore_zero_occupancy false", to keep the original coordinates for atoms with     //
//                                       zero occupancy.                                     //
//                                                                                           //
//  ---> SEARCH OPTION, which must be one of the following options:                          //
//       1. "-pair_target_resnum X -target_chain Y": extracts all constellations formed by   //
//                                                   pairs of residues containing target     //
//                                                   residue X from chain Y.                 //
//       2. "-pair_all_res": extracts all constellations formed by all pairs of residues.    //
//       3. "-pair_target_mutations AB_CD": extracts all constellations that can be obtained //
//                                          with an (A->B, C->D) mutation combination, where //
//                                          A, B, C, and D are one-letter amino acid codes.  //
//       4. "-triple_target_resnum X -target_chain Y": extracts all constellations formed by //
//                                                     triples of residues containing target //
//                                                     residue X from chain Y.               //
//       5. "-triple_all_res": extracts all constellations formed by all triples of residues.//
//                                                                                           //
//       6. "-target_cnl CNLFIL": extracts a single, target constellation as specified in    //
//                                file CNLFIL (see function target_constel() in              //
//                                SearchOptions.cc.                                          //
//                                                                                           //
//  ---> OPTIONAL CONSTRAINTS:                                                               //
//       1. "-cnl_stripped", depriving constellations of the atoms that are closest to the   //
//          backbone, to avoid clash between the rescuing compound and what remains of the   //
//          constellation's residues. Defaults to false. (See function                       //
//          SingResCnlCrea::strip_atoms().)                                                  //
//       2. "-max_atom_sasa X", where X is the maximum allowed SASA (Solvent Accessible      //
//          Surface Area) for an atom in a constellation. (See class FilterBySASA.)          //
//       3. "-interface", extracting only constellations that are shared by two or more      //
//          chains.                                                                          //
//       4. "-aromatic", extracting only constellations that contain at least one aromatic   //
//          ring.                                                                            //
//       5. "-cnl_exclude F", extracting only constellations that do not contain any of the  //
//          residues specified in file F. The format of the file is as follows:              //
//          I1 C1\n                                                                          //
//          ...                                                                              //
//          IN CN\n ,                                                                        //
//          where Ii and Ci are the residue index and the chain ID, respectively, of         //
//          the ith residue to be excluded from constellations (i=1,...,N).                  //
//       6. "-prox_ct_max X", where X is the maximum allowed distance, in Angstroms, from a  //
//          constellation to the N- and C-termini of the chains it belongs to (see class     //
//          FilterByProxTerm).                                                               //
//       7. "-prox_tt_max X", where X is the maximum allowed distance, in Angstroms,         //
//          between the N- and C-termini of any chain that a constellation belongs to.       //
//          This constraint gets activated only upon activation of the "-prox_ct_max"        //
//          constraint (see class FilterByProxTerm).                                         //
//       8. "-prox_nres X", where X is the number of residues forming the N- and C-termini   //
//          in the context of the application of the "-prox_ct_max"  and "-prox_tt_max"      //
//          constraints.                                                                     //
//                                                                                           //
//  ---> TARGET COMPOUND, optional, can be one of the following:                             //
//       1. "-indole_coo", to request constellations that may be rescued by a compound       //
//          containing an indole group and a carboxylic group (see class FilterByIndoleCOO). //
//       2. "-tryptamine", to request constellations that may be rescued by tryptamine (see  //
//          class FilterByTryptamine).                                                       //
//       3. "-amphetamine", to request constellations that may be rescued by amphetamine     //
//          (see class FilterByAmphetamine).                                                 //
//       4. "-histamine", to request constellations that may be rescued by histamine (see    //
//          class FilterByHistamine).                                                        //
//                                                                                           //
// NOTES:                                                                                    //
//  - The names of constellation files produced under all search options have the format     //
//    specified in the notes to functions out_pair_constel() and out_triple_constel()        //
//    (defined in Primitives.cc).                                                            //
//  - The program ignores insertion codes and alternate location identifiers when creating   //
//    output file names. It is therefore recommended to use input PDB files where such       //
//    fields are equal to ' '.                                                               //
//                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////

#include <devel/constel/SingResCnlCrea.hh>
#include <devel/constel/ExcludedFilter.hh>
#include <devel/constel/PairConstelFilters.hh>
#include <devel/constel/FilterByProxTerm.hh>
#include <devel/constel/AromaticFilter.hh>
#include <devel/constel/InterfaceFilter.hh>
#include <devel/constel/FilterBySASA.hh>
#include <devel/constel/MasterFilter.hh>
#include <devel/constel/SearchOptions.hh>
#include <devel/constel/Primitives.hh>
#include <devel/constel/NeighTeller.hh>

#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/chemical/AA.hh>
#include <basic/Tracer.hh>
#include <devel/init.hh>
#include <core/import_pose/import_pose.hh>
#include <basic/options/util.hh>
#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>

#include <string>
#include <cfloat>

static thread_local basic::Tracer TR( "apps.pilot.constel.main" );

////////////////////////////////////////////////////////////////////////////////
//                                    MAIN                                    //
////////////////////////////////////////////////////////////////////////////////

OPT_KEY( Boolean, pair_all_res )
OPT_KEY( Integer, pair_target_resnum )
OPT_KEY( String, pair_target_mutations )
OPT_KEY( Boolean, triple_all_res )
OPT_KEY( Integer, triple_target_resnum )
OPT_KEY( String, target_cnl )

OPT_KEY( Boolean, cnl_stripped )

OPT_KEY( String, target_chain )

OPT_KEY( Real, max_atom_sasa )

OPT_KEY( Boolean, interface )

OPT_KEY( Boolean, aromatic )

OPT_KEY( String, cnl_exclude )

OPT_KEY( Real, prox_ct_max )
OPT_KEY( Real, prox_tt_max )
OPT_KEY( Integer, prox_nres )

OPT_KEY( Boolean, indole_coo )
OPT_KEY( Boolean, tryptamine )
OPT_KEY( Boolean, amphetamine )
OPT_KEY( Boolean, histamine )

#ifdef ROSETTA_FLOAT
	#define DFT_MAX_ATOM_SASA FLT_MAX
#else
	#define DFT_MAX_ATOM_SASA DBL_MAX
#endif

using namespace devel::constel;
using namespace basic::options::OptionKeys;
using basic::options::option;
using core::Size;

int main( int argc, char * argv [] )
{
try {

	NEW_OPT( pair_all_res, "Extracts pair-constellations for all residues", false  );
	NEW_OPT( pair_target_resnum, "Extracts pair-constellations for a target residue", -1 );
	NEW_OPT( pair_target_mutations, "Extracts pair-constellations for a target mutation pair", "**_**" );
	NEW_OPT( triple_all_res, "Extracts triple-constellations for all residues", false );
	NEW_OPT( triple_target_resnum, "Extracts triple-constellations for a target residue", -1 );
	NEW_OPT( target_cnl, "Extracts a single, target constellation as specified in the accompanying file", "");
	NEW_OPT( cnl_stripped, "Deprives constellations of the atoms closest to the backbone", false );
	NEW_OPT( target_chain, "Chain the target residue is in", "A" );
	NEW_OPT( max_atom_sasa, "Maximum allowed SASA for a constellation atom", DFT_MAX_ATOM_SASA );
	NEW_OPT( interface, "Filter to keep only constellations shared by multiple chains", false);
	NEW_OPT( aromatic, "Filter to keep only constellations with at least one aromatic ring", false);
	NEW_OPT( cnl_exclude, "Filter to exclude constellations with given residues", "");
	NEW_OPT( prox_ct_max, "Maximum distance for the proximity of a constellation to a chain terminus", 0 );
	NEW_OPT( prox_tt_max, "Maximum distance for the proximity betweein chain termini" , 10.0 );
	NEW_OPT( prox_nres, "Number of residues forming a chain terminus in the context of proximity evaluation", 10 );
	NEW_OPT( indole_coo, "Filter by indole group plus COO group", false );
	NEW_OPT( tryptamine, "Filter by tryptamine", false );
	NEW_OPT( amphetamine, "Filter by amphetamine", false );
	NEW_OPT( histamine, "Filter by histamine", false);

	devel::init(argc, argv);

	TR << "Remember to set the out::file::suppress_zero_occ_pdb_output flag !!!" << std::endl;
	TR << "You may also want to set -ignore_zero_occupancy false" << std::endl;

	// create native pose from pdb
	core::pose::Pose pose_init;
	std::string const input_pdb_name( basic::options::start_file() );
	core::import_pose::pose_from_pdb( pose_init, input_pdb_name );

	// select constellation-size class
	SingResCnlCrea::init(option[cnl_stripped]);

	// initialize constellation filters
	FilterBySASA::init( option[ max_atom_sasa ], pose_init );
	MasterFilter::addfilt(FilterBySASA::has_low_per_atom_sasa);

	if(option[aromatic])
		MasterFilter::addfilt(has_aromatic);

	if(option[interface])
		MasterFilter::addfilt(at_interface);

	std::string const cnl_excls = option[cnl_exclude];
	if(cnl_excls != "") {
		ExcludedFilter::init(pose_init, cnl_excls);
		MasterFilter::addfilt(ExcludedFilter::hasnt_excluded);
	}

	if(option[prox_ct_max]) {
		FilterByProxTerm::init(pose_init, option[prox_ct_max], option[prox_tt_max],
			option[prox_nres]);
		MasterFilter::addfilt(FilterByProxTerm::sat);
	}

	HBondCommon::init();

	if(option[indole_coo])
		MasterFilter::addfilt(FilterByIndoleCOO::sat);

	if(option[tryptamine])
		MasterFilter::addfilt(FilterByTryptamine::sat);

	if(option[amphetamine])
		MasterFilter::addfilt(FilterByAmphetamine::sat);

	if(option[histamine])
		MasterFilter::addfilt(FilterByHistamine::sat);

	//// perform search requested by command-line option ////
	Size const TOTRES = pose_init.total_residue();

	// option "-pair_all_res"
	if(option[pair_all_res]) {

		using core::conformation::Residue;

		NeighTeller nt(pose_init);
		for(Size i=1; i<TOTRES; ++i ) {
			Residue const& ri = pose_init.residue(i);
			for(Size j=i+1; j<=TOTRES; ++j) {
				Residue const& rj = pose_init.residue(j);
				if(nt.isneigh(ri, rj, pose_init))
					pair_constel_set_idx2(i, j, pose_init);
			}
		}
	}

	// option "-pair_target_resnum"
	int pdbnum = option [ pair_target_resnum ];
	if( pdbnum != -1) {

		std::string const tmp_chain = option[ target_chain ];
		if ( tmp_chain.length() != 1 ) {
			TR << "ERROR!! Chain ID should be one character" << std::endl;
			exit(1);
		}
		char cid = tmp_chain[0];

		pair_constel_set(pdbnum, cid, pose_init);
	}

	// option "-pair_target_mutations"
	std::string tgtmuts = option [ pair_target_mutations ];
	if(tgtmuts[0] != '*') {

		if( ( tgtmuts.length() != 5 ) || ( tgtmuts[2] != '_' ) ) {
			TR << "ERROR!! Format for mutation combination must be \"AB_CD\"" << std::endl;
			exit(1);
		}

		pair_constel_set(tgtmuts, pose_init);
	}

	// option "-triple_all_res"
	if(option[triple_all_res]) {

		utility::vector1<Size> cnl(3);

		Size const UI = TOTRES-2;
		Size const UJ = TOTRES-1;
		Size const UK = TOTRES;

		NeighTeller nt(pose_init);
		for(Size i=1; i<=UI; ++i ) {
			cnl[1] = i;
			for(Size j=i+1; j<=UJ; ++j) {
				cnl[2] = j;
				for(Size k=j+1; k<=UK; ++k) {
					cnl[3] = k;
					if(nt.is_neigh_tree(cnl, pose_init))
						triple_constel_set_idx3(i, j, k, pose_init);
				}
			}
		}
	}

	// option "-triple_target_resnum"
	pdbnum = option [ triple_target_resnum ];
	if( pdbnum != -1) {

		std::string const tmp_chain = option[ target_chain ];
		if ( tmp_chain.length() != 1 ) {
			TR << "ERROR!! Chain ID should be one character" << std::endl;
			exit(1);
		}
		char cid = tmp_chain[0];

		triple_constel_set(pdbnum, cid, pose_init);
	}

	// option "-target_cnl"
	std::string tgtcnl = option[target_cnl];
	if(tgtcnl != "")
		target_constel(tgtcnl, pose_init);


	TR << "TASK COMPLETED" << std::endl;

	return 0;

} // try
catch ( utility::excn::EXCN_Base const & e ) {
	std::cerr << "caught exception " << e.msg() << std::endl;
	return -1;
}

} // main
