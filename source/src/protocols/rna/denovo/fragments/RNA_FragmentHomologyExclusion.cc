// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.


// Rosetta Headers
#include <protocols/rna/denovo/fragments/RNA_FragmentHomologyExclusion.hh>
#include <protocols/rna/denovo/fragments/FragmentLibrary.hh>
#include <protocols/rna/denovo/fragments/TorsionSet.hh>

#include <core/pose/Pose.hh>

#include <core/chemical/util.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueType.hh>
#include <core/pose/copydofs/util.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/id/AtomID.hh>
#include <core/pose/Pose.hh>
#include <core/pose/annotated_sequence.hh>
#include <protocols/stepwise/setup/FullModelInfoSetupFromCommandLine.hh>
#include <protocols/rna/denovo/options/RNA_FragmentMonteCarloOptions.hh>
#include <protocols/rna/denovo/fragments/FullAtomRNA_Fragments.hh>
#include <protocols/rna/denovo/util.hh>
#include <protocols/toolbox/AtomLevelDomainMap.hh>

#include <core/types.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/rna.OptionKeys.gen.hh>
#include <set>


namespace protocols {
namespace rna {
namespace denovo {
namespace fragments {

using namespace core;
using namespace basic::options;

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
				newseg.first = segment.first + ii - 1;
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
		//TR << ii << " " << temp << std::endl;
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
		//TR << i << std::endl;
		protocols::rna::denovo::get_base_pairing_info( pose, i, secstruct1, edge_is_base_pairing );
		secstruct += secstruct1;
	}

	//TR << "SECSTRUCT: " << secstruct << std::endl;
	return secstruct;
}

utility::vector1< core::Size > analyze_for_homology( std::string const & in_file, RNA_Fragments const & fragments, std::string const & exclusion_match_type = "MATCH_EXACT" ) {
	using namespace core::pose;
	using namespace core::pose::copydofs;
	using namespace core::scoring;
	using namespace core::chemical;
	using namespace core::id;
	using namespace protocols::rna::denovo::options;
	using namespace protocols::rna::denovo::fragments;
	using namespace protocols::stepwise::setup;
	using namespace protocols::toolbox;
	using namespace utility::file;

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

	std::string const secstruct = figure_out_secstruct( pose_input );

	utility::vector1< std::pair< Size, std::string > > secstruct_segments; // position, secstruct
	obtain_secstruct_segments( secstruct_segments, secstruct );
	secstruct_segments = split_segments_longer_than_6mers( secstruct_segments );

	utility::vector1< Size > foo;
	for ( auto const & elem : secstruct_segments ) {
		// search for all fragments that might match sequence
		Size const position = elem.first;
		std::string const & secstruct = elem.second;
		//TR << position << " " << secstruct << std::endl;
		Size const frag_length( secstruct.size());
		std::string sequence = pose_input.sequence().substr( position - 1, frag_length ); // uucg loop

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
			//TR << rmsd << std::endl;
			if ( rmsd < option[ OptionKeys::rna::denovo::fragment_homology_rmsd ]() ) {
				// For e.g. a 4-mer overlap, this just excludes the first one.
				//foo.push_back( torsion_set.get_index_in_vall() );

				for ( Size ii = 1; ii <= secstruct.size(); ++ii ) {
					foo.push_back( torsion_set.get_index_in_vall() + ii - 1 );
				}
			}
		}
	}
	return foo;
}

RNA_FragmentHomologyExclusion::RNA_FragmentHomologyExclusion( RNA_Fragments const & all_rna_fragments ) {
	// AMW: Insert any vall lines explicitly specified from the command line
	fragment_lines_.insert(
		option[ OptionKeys::rna::denovo::exclude_fragments ]().begin(),
		option[ OptionKeys::rna::denovo::exclude_fragments ]().end() );

	if ( option[ OptionKeys::rna::denovo::exclude_native_fragments ].value() ) {
		auto homologous_lines = analyze_for_homology( option[ OptionKeys::in::file::native ].value(), all_rna_fragments,
			option[ OptionKeys::rna::denovo::exclusion_match_type ]() );
		fragment_lines_.insert( homologous_lines.begin(), homologous_lines.end() );
	}
	if ( option[ OptionKeys::rna::denovo::exclude_fragment_files ].user() ) {
		// Add vall lines homologous to the loops from these poses
		for ( auto const & file : option[ OptionKeys::rna::denovo::exclude_fragment_files ]() ) {
			auto homologous_lines = analyze_for_homology( file, all_rna_fragments, option[ OptionKeys::rna::denovo::exclusion_match_type ]() );
			fragment_lines_.insert( homologous_lines.begin(), homologous_lines.end() );
		}
	}
}

} //fragments
} //denovo
} //rna
} //protocols
