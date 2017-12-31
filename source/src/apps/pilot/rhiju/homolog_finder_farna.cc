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
#include <core/scoring/ScoreFunction.hh>
#include <devel/init.hh>
#include <core/id/AtomID.hh>
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/import_pose/options/RNA_FragmentMonteCarloOptions.hh>
#include <core/fragment/rna/FullAtomRNA_Fragments.hh>
#include <core/fragment/rna/FragmentLibrary.hh>
#include <core/fragment/rna/TorsionSet.hh>
#include <core/pose/rna/util.hh>
#include <protocols/rna/denovo/util.hh>
#include <core/pose/toolbox/AtomLevelDomainMap.hh>
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
using namespace core::import_pose;
using namespace protocols;
using namespace basic::options;
using namespace basic::options::OptionKeys;
using utility::vector1;

static basic::Tracer TR( "apps.pilot.rhiju.homolog_finder_farna" );

utility::vector1< std::pair< Size, std::string > >
split_segments_longer_than_6mers(
	utility::vector1< std::pair< Size, std::string > > const & secstruct_segments
) {
	utility::vector1< std::pair< Size, std::string > > new_segments;

	for ( auto const & segment : secstruct_segments ) {
		if ( segment.second.size() <= 6 ) new_segments.push_back( segment );
		else {
			for ( Size ii = 1; ii <= segment.second.size() - 6; ++ii ) {
				std::pair< Size, std::string > newseg;
				newseg.first  = segment.first + ii - 1;
				newseg.second = segment.second.substr( ii - 1, 6 );
				new_segments.push_back( newseg );
			}
		}
	}

	return new_segments;
}

void
obtain_secstruct_segments(
	utility::vector1< std::pair< Size, std::string > > & secstruct_segments,
	std::string const & secstruct
) {
	// NOTE: X is defined as != H
	bool accumulate_next = false;
	std::string temp = "";
	Size pos = 0;
	for ( Size ii = 1; ii <= secstruct.size(); ++ii ) {
		TR << ii << " " << temp << std::endl;
		if ( secstruct[ ii-1 ] == 'H' ) {
			accumulate_next = false;
			if ( pos != 0 ) secstruct_segments.emplace_back( pos, temp );
			temp = "";
			pos = 0;
		}
		if ( secstruct[ ii-1 ] != 'H' && ! accumulate_next ) {
			accumulate_next = true;
			pos = ii;
		}
		if ( accumulate_next ) {
			if ( secstruct[ ii-1 ] != 'H' ) temp += 'X'; //secstruct[ ii ];
		}
	}
}

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
		TR << i << std::endl;
		core::pose::rna::get_base_pairing_info( pose, i, secstruct1, edge_is_base_pairing );
		secstruct += secstruct1;
	}

	TR << "SECSTRUCT: " << secstruct << std::endl;
	return secstruct;
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
homolog_finder()
{
	using namespace core::pose;
	using namespace core::pose::copydofs;
	using namespace core::scoring;
	using namespace core::chemical;
	using namespace core::id;
	using namespace core::import_pose::options;
	using namespace core::fragment::rna;
	using namespace core::pose::toolbox;
	using namespace utility::file;

	// Following could be generalized to fa_standard, after recent unification, but
	// probably should wait for on-the-fly residue type generation.
	ResidueTypeSetCOP rsd_set = ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD );

	// Following could go to a FullModelSetup class.
	// read starting pose(s) from disk
	utility::vector1< std::string > const & input_files = option[ in::file::s ]();
	PoseOP pose_op = get_pdb_and_cleanup( input_files[1], rsd_set );
	Pose & pose_input = *pose_op;

	// set up Fragments
	RNA_FragmentMonteCarloOptions options;
	options.initialize_from_command_line();
	FullAtomRNA_Fragments fragments( options.all_rna_fragments_file() );

	// Get the secstruct
	std::string const secstruct = figure_out_secstruct( pose_input );

	utility::vector1< std::pair< Size, std::string > > secstruct_segments; // position, secstruct
	obtain_secstruct_segments( secstruct_segments, secstruct );
	secstruct_segments = split_segments_longer_than_6mers( secstruct_segments );

	for ( auto const & elem : secstruct_segments ) {
		// search for all fragments that might match sequence
		Size const position = elem.first;
		std::string const & secstruct = elem.second;
		TR << position << " " << secstruct << std::endl;
		Size const frag_length( secstruct.size());
		std::string const sequence = pose_input.sequence().substr( position - 1, frag_length ); // uucg loop
		Size type( MATCH_EXACT );

		FragmentLibrary const & library = *( fragments.get_fragment_library_pointer( sequence, secstruct, nullptr, utility::vector1< SYN_ANTI_RESTRICTION >(), type ) );
		TR << "Found library for " << sequence << " with number of matches: " << library.get_align_depth() << std::endl;

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

		utility::vector1< Size > foo;
		// Let's do the loop
		AtomLevelDomainMapOP atom_level_domain_map( new AtomLevelDomainMap( pose ) );
		for ( Size n = 1; n <= library.get_align_depth(); n++ ) {
			TorsionSet torsion_set = library.get_fragment_torsion_set( n );
			fragments.insert_fragment( pose, 1, torsion_set, atom_level_domain_map );
			rmsd = superimpose_pose( pose, pose_input, atom_id_map );
			TR << rmsd << std::endl;
			if ( rmsd < 1.0 ) {
				foo.push_back( torsion_set.get_index_in_vall() );
			}
		}
		for ( auto const elem : foo ) {
			std::cout << "elem: " << elem << std::endl;
		}
	}

	std::string const outfile = option[ out::file::o ]();
	// if ( outfile.size() > 0 ) pose.dump_pdb( outfile );
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

	} catch (utility::excn::Exception const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
}
