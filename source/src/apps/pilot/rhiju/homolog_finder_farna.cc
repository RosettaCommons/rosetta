// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file swa_monte_carlo.cc
/// @author Rhiju Das (rhiju@stanford.edu)

// libRosetta headers
#include <core/types.hh>
#include <core/chemical/util.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/pose/copydofs/util.hh>
#include <core/scoring/rms_util.hh>
#include <devel/init.hh>
#include <core/id/AtomID.hh>
#include <core/pose/Pose.hh>
#include <core/pose/annotated_sequence.hh>
#include <protocols/stepwise/setup/FullModelInfoSetupFromCommandLine.hh>
#include <protocols/farna/options/RNA_FragmentMonteCarloOptions.hh>
#include <protocols/farna/fragments/FullAtomRNA_Fragments.hh>
#include <protocols/toolbox/AtomLevelDomainMap.hh>
#include <protocols/viewer/viewers.hh>

//////////////////////////////////////////////////
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh> // for option[ out::file::silent  ] and etc.
#include <basic/options/keys/in.OptionKeys.gen.hh> // for option[ in::file::tags ] and etc.
#include <basic/options/option_macros.hh>

#include <basic/Tracer.hh>
#include <utility/file/FileName.hh>


// C++ headers
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <list>

using namespace core;
using namespace protocols;
using namespace basic::options;
using namespace basic::options::OptionKeys;
using utility::vector1;

static basic::Tracer TR( "apps.pilot.rhiju.homolog_finder_farna" );

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
homolog_finder()
{
  using namespace core::pose;
  using namespace core::pose::copydofs;
  using namespace core::scoring;
  using namespace core::chemical;
  using namespace core::id;
  using namespace protocols::farna::options;
  using namespace protocols::farna::fragments;
  using namespace protocols::stepwise::setup;
  using namespace protocols::toolbox;
  using namespace utility::file;

	// Following could be generalized to fa_standard, after recent unification, but
	// probably should wait for on-the-fly residue type generation.
	ResidueTypeSetCAP rsd_set = ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD );

	// Following could go to a FullModelSetup class.
	// read starting pose(s) from disk
	utility::vector1< std::string > const & input_files = option[ in::file::s ]();
	PoseOP pose_op = get_pdb_and_cleanup( input_files[1], rsd_set );
	Pose & pose_input = *pose_op;

	// set up Fragments
	RNA_FragmentMonteCarloOptions options;
	options.initialize_from_command_line();
	FullAtomRNA_Fragments fragments( options.all_rna_fragments_file() );

	// search for all fragments that might match sequence
	Size const frag_length( 4 );
	Size const position( 3 );
	std::string const sequence = pose_input.sequence().substr( position - 1, frag_length ); // uucg loop
	std::string const secstruct = "XXXX";
	Size type( MATCH_EXACT );

	FragmentLibrary const & library = *( fragments.get_fragment_library_pointer( sequence, secstruct, type ) );
	TR << "Found library for " << sequence << " with number of matches: " << library.get_align_depth();

	// Now try all the fragments one-by-one in a "scratch pose"
	// make the scratch pose.
	Pose pose;
	make_pose_from_sequence( pose, sequence, ResidueTypeSetCOP( rsd_set ), false );
	cleanup( pose );
	protocols::viewer::add_conformation_viewer( pose.conformation(), "current", 400, 400 );

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
		TR << rmsd << std::endl;
	}

	std::string const outfile = option[ out::file::o ]();
	//	if ( outfile.size() > 0 ) pose.dump_pdb( outfile );

}

///////////////////////////////////////////////////////////////
void*
my_main( void* )
{
	clock_t const my_main_time_start( clock() );
	homolog_finder();
	protocols::viewer::clear_conformation_viewers();
	std::cout << "Total time to run " << static_cast<core::Real>( clock() - my_main_time_start ) / CLOCKS_PER_SEC << " seconds." << std::endl;
  exit( 0 );
}


///////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{
	try {

		utility::vector1< core::Size > blank_size_vector;
		devel::init(argc, argv);
		protocols::viewer::viewer_main( my_main );

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
}


