// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @brief Extracts constellations from a protein structure.
/// @author jk
/// @author Andrea Bazzoli (ndrbzz@gmail.com)

///////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                           //
//                                          constel                                          //
//                                                                                           //
// APPLICATION-SPECIFIC OPTIONS:                                                             //
//                                                                                           //
//  ---> SEARCH OPTION (required), must be one of the following option combinations:         //
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
//  ---> FILTERING OPTIONS (optional), can be one or more of the following options:         //
//       1. "-cnl_stripped": deprives constellations of the atoms that are closest to the    //
//          backbone, to avoid clash between the rescuing compound and what remains of the   //
//          constellation's residues. Defaults to false. (See function                       //
//          SingResCnlCrea::strip_atoms().)                                                  //
//       2. "-max_atom_sasa X": X is the maximum allowed SASA (Solvent Accessible Surface    //
//          Area) for an atom in a constellation. (See class FilterBySASA.)                  //
//       3. "-chain_interface": extracts only constellations that are shared by two or more  //
//          chains.                                                                          //
//       4. "-aromatic": extracts only constellations that contain at least one aromatic     //
//          ring.                                                                            //
//       5. "-cnl_exclude F": extracts only constellations that do not contain any of the    //
//          residues specified in file F. The format of the file is as follows:              //
//          I1 C1\n                                                                          //
//          ...                                                                              //
//          IN CN\n ,                                                                        //
//          where Ii and Ci are the residue index and the chain ID, respectively, of         //
//          the ith residue to be excluded from constellations (i=1,...,N).                  //
//       6. "-prox_ct_max X": X is the maximum allowed distance, in Angstroms, from a        //
//          constellation to the N- and C-termini of the chains it belongs to (see class     //
//          FilterByProxTerm).                                                               //
//       7. "-prox_tt_max X": X is the maximum allowed distance, in Angstroms, between the   //
//          N- and C-termini of any chain that a constellation belongs to. This constraint   //
//          gets activated only upon activation of the "-prox_ct_max" constraint (see class  //
//          FilterByProxTerm).                                                               //
//       8. "-prox_nres X": X is the number of residues forming the N- and C-termini. This   //
//          option, too, is meant to be applied only in combination with the "-prox_ct_max"  //
//          constraint.                                                                      //
//       9. "-indole_coo": requests constellations that may be rescued by a compound         //
//          containing an indole group and a carboxylic group (see class FilterByIndoleCOO). //
//      10. "-tryptamine": requests constellations that may be rescued by tryptamine (see    //
//          class FilterByTryptamine).                                                       //
//      11. "-amphetamine": requests constellations that may be rescued by amphetamine       //
//          (see class FilterByAmphetamine).                                                 //
//      12. "-histamine": requests constellations that may be rescued by histamine (see      //
//          class FilterByHistamine).                                                        //
//                                                                                           //
// NOTES:                                                                                    //
//  - The names of constellation files produced under all search options have the format     //
//    specified in the notes to functions out_pair_constel() and out_triple_constel()        //
//    (defined in Primitives.cc).                                                            //
//  - The program ignores insertion codes and alternate location identifiers when creating   //
//    output file names.                                                                     //
//                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////

#include <protocols/constel/SingResCnlCrea.hh>
#include <protocols/constel/ExcludedFilter.hh>
#include <protocols/constel/PairConstelFilters.hh>
#include <protocols/constel/FilterByProxTerm.hh>
#include <protocols/constel/AromaticFilter.hh>
#include <protocols/constel/InterfaceFilter.hh>
#include <protocols/constel/FilterBySASA.hh>
#include <protocols/constel/MasterFilter.hh>
#include <protocols/constel/SearchOptions.hh>
#include <protocols/constel/Primitives.hh>
#include <protocols/constel/NeighTeller.hh>

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
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/constel.OptionKeys.gen.hh>

#include <string>
#include <cfloat>

static basic::Tracer TR( "apps.public.constel.main" );

////////////////////////////////////////////////////////////////////////////////
//                                    MAIN                                    //
////////////////////////////////////////////////////////////////////////////////

using namespace protocols::constel;
using namespace basic::options::OptionKeys;
using namespace basic::options::OptionKeys::constel;
using basic::options::option;
using core::Size;

int main( int argc, char * argv [] )
{
	try {

		devel::init(argc, argv);

		// ignore TER records on input so that ATOM serial numbers on output have no gaps
		TR << "setting default value of -no_chainend_ter to true" << std::endl;
		option[out::file::no_chainend_ter].default_value(true);

		// use zero occupancy to discriminate atoms that are not part of a constellation
		TR << "setting default value of -suppress_zero_occ_pdb_output to true" << std::endl;
		option[out::file::suppress_zero_occ_pdb_output].default_value(true);

		// consider atoms whose original occupancy is zero
		TR << "setting default value of -ignore_zero_occupancy to false" << std::endl;
		option[run::ignore_zero_occupancy].default_value(false);

		// ignore unrecognized residues by default
		TR << "setting default value of -ignore_unrecognized_res to true" << std::endl;
		option[in::ignore_unrecognized_res].default_value(true);

		// ignore waters by default
		TR << "setting default value of -ignore_waters to true" << std::endl;
		option[in::ignore_waters].default_value(true);

		// do not load PDB components
		TR << "setting default value of -load_PDB_components to false" << std::endl;
		option[in::file::load_PDB_components].default_value(false);


		// create native pose from pdb
		core::pose::Pose pose_init;
		std::string const input_pdb_name( basic::options::start_file() );
		core::import_pose::pose_from_file( pose_init, input_pdb_name , core::import_pose::PDB_file);

		// select constellation-size class
		SingResCnlCrea::init(option[cnl_stripped]);

		// initialize constellation filters
		FilterBySASA::init( option[ max_atom_sasa ], pose_init );
		MasterFilter::addfilt(FilterBySASA::has_low_per_atom_sasa);

		if ( option[aromatic] ) {
			MasterFilter::addfilt(has_aromatic);
		}

		if ( option[chain_interface] ) {
			MasterFilter::addfilt(at_interface);
		}

		std::string const cnl_excls = option[cnl_exclude];
		if ( cnl_excls != "" ) {
			ExcludedFilter::init(pose_init, cnl_excls);
			MasterFilter::addfilt(ExcludedFilter::hasnt_excluded);
		}

		if ( option[prox_ct_max] ) {
			FilterByProxTerm::init(pose_init, option[prox_ct_max], option[prox_tt_max],
				option[prox_nres]);
			MasterFilter::addfilt(FilterByProxTerm::is_satisfied);
		}

		HBondCommon::init();

		if ( option[indole_coo] ) {
			MasterFilter::addfilt(FilterByIndoleCOO::is_satisfied);
		}

		if ( option[tryptamine] ) {
			MasterFilter::addfilt(FilterByTryptamine::is_satisfied);
		}

		if ( option[amphetamine] ) {
			MasterFilter::addfilt(FilterByAmphetamine::is_satisfied);
		}

		if ( option[histamine] ) {
			MasterFilter::addfilt(FilterByHistamine::is_satisfied);
		}

		//// perform search requested by command-line option ////
		Size const TOTRES = pose_init.size();

		// option "-pair_all_res"
		if ( option[pair_all_res] ) {

			using core::conformation::Residue;

			NeighTeller nt(pose_init);
			for ( Size i=1; i<TOTRES; ++i ) {
				Residue const& ri = pose_init.residue(i);
				for ( Size j=i+1; j<=TOTRES; ++j ) {
					Residue const& rj = pose_init.residue(j);
					if ( nt.isneigh(ri, rj, pose_init) ) {
						pair_constel_set_idx2(i, j, pose_init);
					}
				}
			}
		}

		// option "-pair_target_resnum"
		int pdbnum = option [ pair_target_resnum ];
		if ( pdbnum != -1 ) {

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
		if ( tgtmuts[0] != '*' ) {

			if ( ( tgtmuts.length() != 5 ) || ( tgtmuts[2] != '_' ) ) {
				TR << "ERROR!! Format for mutation combination must be \"AB_CD\"" << std::endl;
				exit(1);
			}

			pair_constel_set(tgtmuts, pose_init);
		}

		// option "-triple_all_res"
		if ( option[triple_all_res] ) {

			utility::vector1<Size> cnl(3);

			Size const UI = TOTRES-2;
			Size const UJ = TOTRES-1;
			Size const UK = TOTRES;

			NeighTeller nt(pose_init);
			for ( Size i=1; i<=UI; ++i ) {
				cnl[1] = i;
				for ( Size j=i+1; j<=UJ; ++j ) {
					cnl[2] = j;
					for ( Size k=j+1; k<=UK; ++k ) {
						cnl[3] = k;
						if ( nt.is_neigh_tree(cnl, pose_init) ) {
							triple_constel_set_idx3(i, j, k, pose_init);
						}
					}
				}
			}
		}

		// option "-triple_target_resnum"
		pdbnum = option [ triple_target_resnum ];
		if ( pdbnum != -1 ) {

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
		if ( tgtcnl != "" ) {
			target_constel(tgtcnl, pose_init);
		}


		TR << "TASK COMPLETED" << std::endl;

		return 0;

	} // try
catch (utility::excn::Exception const & e ) {
	e.display();
	return -1;
}

} // main
