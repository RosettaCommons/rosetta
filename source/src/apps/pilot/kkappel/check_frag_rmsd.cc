// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

// Check the minimimum RMSD of any fragment in the RNA library to our RNA structure of interest
// (i.e. what's the best model we could build given our fragment library?)

/// @file
/// @brief


// libRosetta headers
#include <core/types.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/io/silent/BinarySilentStruct.hh>
#include <core/io/silent/SilentFileOptions.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/rna/RNA_ScoringInfo.hh>
#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <protocols/viewer/viewers.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/init/init.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/import_pose/import_pose.hh>
#include <core/import_pose/pose_stream/PoseInputStream.hh>
#include <core/import_pose/pose_stream/PoseInputStream.fwd.hh>
#include <core/import_pose/pose_stream/PDBPoseInputStream.hh>
#include <core/import_pose/pose_stream/SilentFilePoseInputStream.hh>
#include <basic/database/open.hh>
#include <utility/vector1.hh>
#include <ObjexxFCL/string.functions.hh>
#include <protocols/rna/denovo/setup/RNA_DeNovoPoseInitializer.hh>
#include <protocols/rna/denovo/util.hh>
#include <protocols/rna/denovo/fragments/RNA_MatchType.hh>
#include <protocols/rna/denovo/fragments/RNA_Fragments.hh>
#include <protocols/rna/denovo/fragments/FragmentLibrary.hh>
#include <protocols/rna/denovo/fragments/TorsionSet.hh>
#include <protocols/rna/denovo/libraries/RNA_LibraryManager.hh>
#include <protocols/toolbox/AtomLevelDomainMap.hh>
#include <protocols/stepwise/setup/FullModelInfoSetupFromCommandLine.hh>

#include <core/pose/copydofs/util.hh>
#include <protocols/rna/denovo/options/RNA_FragmentMonteCarloOptions.hh>

// C++ headers
#include <iostream>
#include <string>

// option key includes
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/rna.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/keys/stepwise.OptionKeys.gen.hh>

#include <utility/excn/Exceptions.hh>


using namespace core;
using namespace protocols;
using namespace basic::options::OptionKeys;
using namespace basic::options;
using utility::vector1;


OPT_KEY( String, input_secstruct )


std::string
figure_out_secstruct( pose::Pose & pose ){
	using namespace core::scoring;
	using namespace ObjexxFCL;

	// Need some stuff to figure out which residues are base paired. First score.
	ScoreFunctionOP scorefxn( new ScoreFunction );
	scorefxn->set_weight( rna_base_pair, 1.0 );
	(*scorefxn)( pose );

	std::string secstruct = "";
	FArray1D < bool > edge_is_base_pairing( 3, false );
	char secstruct1( 'X' );
	for ( Size i=1; i <= pose.size() ; ++i ) {
		//TR << i << std::endl;
		protocols::rna::denovo::get_base_pairing_info( pose, i, secstruct1, edge_is_base_pairing );
		secstruct += secstruct1;
	}

	std::cout << "SECSTRUCT: " << secstruct << std::endl;
	//TR << "SECSTRUCT: " << secstruct << std::endl;
	return secstruct;
}

utility::vector1< Real > get_min_frag_rmsd( std::string const & in_file, protocols::rna::denovo::fragments::RNA_Fragments const & fragments, std::string const & exclusion_match_type = "MATCH_EXACT" ) {
	using namespace core::pose;
	using namespace core::pose::copydofs;
	using namespace core::scoring;
	using namespace core::import_pose;
	using namespace core::chemical;
	using namespace core::id;
	using namespace protocols::rna::denovo::options;
	using namespace protocols::rna::denovo::fragments;
	using namespace protocols::stepwise::setup;
	using namespace protocols::toolbox;
	using namespace utility::file;
	using namespace basic::options;

	Size type( MATCH_EXACT );
	if ( exclusion_match_type == "MATCH_EXACT" ) {
		type = MATCH_EXACT;
	} else if ( exclusion_match_type == "MATCH_YR" ) {
		type = MATCH_YR;
	} else if ( exclusion_match_type == "MATCH_ALL" ) {
		type = MATCH_ALL;
	} else {
		utility_exit_with_message( "Illegal value provided for option -exclusion_match_type. Must be MATCH_EXACT (default), MATCH_YR, or MATCH_ALL.");
	}

	ResidueTypeSetCAP rsd_set = ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD );

	PoseOP pose_op = get_pdb_and_cleanup( in_file, rsd_set );
	Pose & pose_input = *pose_op;

	std::string full_secstruct;
	if ( option[ input_secstruct ].user() ) {
		full_secstruct = option[ input_secstruct ]();
	} else {
		full_secstruct = figure_out_secstruct( pose_input );
	}

	// std::string const full_secstruct = figure_out_secstruct( pose_input );
	// std::string const full_secstruct = "XXXXXXXXXX";

	utility::vector1< std::pair< Size, std::string > > secstruct_segments; // position, secstruct
	// obtain_secstruct_segments( secstruct_segments, secstruct );
	// secstruct_segments = split_segments_longer_than_6mers( secstruct_segments );

	Real min_rmsd;
	Size num_below_thresh;
	// Real threshold_rmsd = 0.2;
	Real threshold_rmsd = 0.5;
	utility::vector1< Real > min_rmsds;
	utility::vector1< Size > num_fragments_rmsd_below_thresh;
	for ( core::Size i=1; i<=full_secstruct.size(); ++i ) {
		min_rmsd = 1000.;
		num_below_thresh = 0;
		// search for all fragments that might match sequence
		Size const position = i;
		std::string const secstruct = full_secstruct.substr(i-1, 1);
		//std::string secstruct = full_secstruct[i-1];
		//TR << position << " " << secstruct << std::endl;
		Size const frag_length( secstruct.size());
		std::string sequence = pose_input.sequence().substr( position - 1, 1 ); // uucg loop
		//std::string sequence = pose_input.sequence().substr( position - 1, frag_length ); // uucg loop


		// Fix up sequence (noncanonicals)
		for ( Size ii = 0; ii < sequence.size(); ++ii ) {
			if ( pose_input.residue_type( ii + position ).na_analogue() == chemical::na_rad ) {
				sequence[ ii ] = 'a';
			} else if ( pose_input.residue_type( ii + position ).na_analogue() == chemical::na_rcy ) {
				sequence[ ii ] = 'c';
			} else if ( pose_input.residue_type( ii + position ).na_analogue() == chemical::na_rgu ) {
				sequence[ ii ] = 'g';
			} else if ( pose_input.residue_type( ii + position ).na_analogue() == chemical::na_ura ) {
				sequence[ ii ] = 'u';
			}
		}

		// This seems curious. Of course, we could maybe pass `this`, and speed
		// up our searches. Maybe.
		FragmentLibrary const & library = *( fragments.get_fragment_library_pointer( sequence, secstruct, nullptr, utility::vector1< SYN_ANTI_RESTRICTION >(), type ) );
		//TR << "Found library for " << sequence << " with number of matches: " << library.get_align_depth() << std::endl;

		// Now try all the fragments one-by-one in a "scratch pose"
		// make the scratch pose.
		Pose pose;
		make_pose_from_sequence( pose, sequence, ResidueTypeSetCOP( rsd_set ), false );
		cleanup( pose );

		// get ready to compute rmsd's.
		Real rmsd( 0.0 );
		std::map< Size, Size > res_map;
		for ( Size i = 1; i <= frag_length; i++ ) res_map[ i ] = i + position - 1;
		std::map< AtomID, AtomID > atom_id_map;
		setup_atom_id_map_match_atom_names( atom_id_map, res_map, pose, pose_input );

		// Let's do the loop
		AtomLevelDomainMapOP atom_level_domain_map( new AtomLevelDomainMap( pose ) );
		for ( Size n = 1; n <= library.get_align_depth(); n++ ) {
			TorsionSet torsion_set = library.get_fragment_torsion_set( n );
			fragments.insert_fragment( pose, 1, torsion_set, atom_level_domain_map );
			rmsd = superimpose_pose( pose, pose_input, atom_id_map );
			if ( rmsd < min_rmsd ) {
				min_rmsd = rmsd;
				// should definitely do this at the end, but just want to do something easy now
				std::ostringstream pose_name;
				pose_name << "min_rmsd_fragment_" << position << ".pdb";
				pose.dump_pdb( pose_name.str() );
			}
			if ( rmsd < threshold_rmsd ) {
				num_below_thresh += 1;
			}
		}
		min_rmsds.push_back( min_rmsd );
		num_fragments_rmsd_below_thresh.push_back( num_below_thresh );
	}
	std::cout << num_fragments_rmsd_below_thresh << std::endl;
	return min_rmsds;
}
///////////////////////////////////////////////////////////////////////////////


// Just try to run the homology stuff
void test() {
	using namespace protocols::rna::denovo::fragments;
	using namespace protocols::rna::denovo;
	using namespace protocols::rna::denovo::libraries;

	// get all the RNA fragments
	RNA_Fragments const & all_rna_fragments( RNA_LibraryManager::get_instance()->rna_fragment_library( basic::database::full_name("sampling/rna/RICHARDSON_RNA09.torsions") ) );
	utility::vector1< Real > test = get_min_frag_rmsd( option[ OptionKeys::in::file::native ].value(), all_rna_fragments, "MATCH_EXACT" );
	for ( core::Size i = 1; i <= test.size(); ++i ) {
		std::cout << "Minimum RMSD position " << i << ": " << test[i] << std::endl;
	}
	//std::cout << test << std::endl;

}


///////////////////////////////////////////////////////////////////////////////
void*
my_main( void* )
{
	test();

	protocols::viewer::clear_conformation_viewers();
	exit( 0 );
}


///////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{
	try {
		using namespace basic::options;

		std::cout << std::endl << "Basic usage:  " << argv[0] << "  -s <pdb file> " << std::endl;
		std::cout              << "              " << argv[0] << "  -in:file:silent <silent file> " << std::endl;
		std::cout << std::endl << " Type -help for full slate of options." << std::endl << std::endl;

		utility::vector1< Size > blank_size_vector;
		utility::vector1< std::string > blank_string_vector;
		option.add_relevant( score::weights );
		option.add_relevant( in::file::s );
		option.add_relevant( in::file::silent );
		option.add_relevant( in::file::tags );
		option.add_relevant( in::file::native );
		option.add_relevant( OptionKeys::rna::denovo::fragment_homology_rmsd );
		NEW_OPT( input_secstruct, "secondary structure", "" );

		////////////////////////////////////////////////////////////////////////////
		// setup
		////////////////////////////////////////////////////////////////////////////
		core::init::init(argc, argv);

		////////////////////////////////////////////////////////////////////////////
		// end of setup
		////////////////////////////////////////////////////////////////////////////
		protocols::viewer::viewer_main( my_main );
	} catch ( utility::excn::Exception const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
}
