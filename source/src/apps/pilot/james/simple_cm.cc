// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file prof_align.cc
/// @brief
/// @author James Thompson

#include <core/types.hh>
#include <devel/init.hh>
#include <basic/prof.hh>
#include <basic/Tracer.hh>
#include <core/chemical/AA.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ChemicalManager.hh>

#include <basic/options/option.hh>

#include <core/sequence/util.hh>
#include <core/sequence/Sequence.hh>
#include <core/sequence/SequenceProfile.hh>
#include <core/sequence/SequenceAlignment.hh>
#include <core/sequence/NWAligner.hh>
#include <core/sequence/SWAligner.hh>
#include <core/sequence/L1ScoringScheme.hh>
#include <core/sequence/MatrixScoringScheme.hh>
#include <core/sequence/SimpleScoringScheme.hh>
#include <core/sequence/ScoringScheme.fwd.hh>
#include <core/sequence/ScoringSchemeFactory.hh>

#include <core/scoring/rms_util.hh>

#include <core/pose/Pose.hh>

#include <core/io/pdb/pose_io.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/SilentStruct.hh>
#include <core/io/silent/SilentStructFactory.hh>

#include <protocols/comparative_modeling/ThreadingMover.hh>
#include <protocols/comparative_modeling/AlignmentSet.hh>
#include <protocols/jobdist/standard_mains.hh>
#include <protocols/jobdist/not_universal_main.hh>

#include <utility/vector1.hh>
#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>
#include <utility/file/FileName.hh>
#include <utility/file/file_sys_util.hh>

#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/char.functions.hh>
#include <ObjexxFCL/string.functions.hh>

#include <numeric/random/random.hh>

// C++ headers
#include <map>
#include <fstream>
#include <iostream>
#include <string>


// option key includes

#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/frags.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>

#include <utility/excn/Exceptions.hh>



///////////////////////////////////////////////////////////////////////////////

using std::map;
using std::string;
using core::Size;
using core::Real;
using utility::vector1;
using utility::file::FileName;
using ObjexxFCL::format::A;
using ObjexxFCL::format::F;
using ObjexxFCL::format::I;
using ObjexxFCL::string_of;
using namespace basic;
using namespace core::sequence;
using namespace basic::options;
using namespace basic::options::OptionKeys;



int
main( int argc, char* argv [] ) {
	try {

	using core::Real;
	using core::Size;

	// options, random initialization
	devel::init( argc, argv );

	// query and template profiles
	FileName fn1( option[ in::file::pssm ]()[1] );
	FileName fn2( option[ in::file::pssm ]()[2] );

	SequenceProfileOP prof1( new SequenceProfile );
	prof1->read_from_file( fn1 );
	prof1->convert_profile_to_probs( 1.0 ); // was previously implicit in read_from_file()

	SequenceProfileOP prof2( new SequenceProfile );
	prof2->read_from_file( fn2 );
	prof2->convert_profile_to_probs( 1.0 ); // was previously implicit in read_from_file()

	// eliminate leading paths from prof1 and prof2
	prof1->id( FileName( prof1->id() ).base() );
	prof2->id( FileName( prof2->id() ).base() );

	// pdbs
	using namespace core::chemical;
	ResidueTypeSetCAP rsd_set = ChemicalManager::get_instance()->residue_type_set(
		option[ in::file::residue_type_set ]()
	);

	core::pose::Pose template_pose;
	if ( option[ in::file::template_pdb ].user() ) {
		core::import_pose::pose_from_pdb(
			template_pose,
			*rsd_set,
			option[ in::file::template_pdb ]()[1]
		);
	}

	// scoring scheme for aligning profiles
	std::string const scoring_scheme_type( option[ frags::scoring::profile_score ]() );
	ScoringSchemeFactory ssf;
	ScoringSchemeOP ss( ssf.get_scoring_scheme( scoring_scheme_type ) );

	// for the ProfSim scoring scheme, the optimal opening and extension
	// penalties were 2 and 0.2, with a scoring shift of -0.45 applied to
	// all ungapped aligned pairs.
	Real const max_gap_open    (  -1   );
	Real const min_gap_open    (  -5   );
	Real const max_gap_extend  (  -0.5 );
	Real const min_gap_extend  (  -4   );
	Real const open_step_size  (   1   );
	Real const extend_step_size(   0.5 );

	// construct alignments
	NWAligner nw_aligner;
	SWAligner sw_aligner;

	std::string output_fn = option[ out::file::alignment ]();
	utility::io::ozstream output( output_fn );

	SequenceAlignment true_aln;
	if ( option[ in::file::alignment ].user() ) {
		std::string align_fn  = option[ in::file::alignment ]()[1];
		true_aln.read_from_file( align_fn );
		//std::cout << true_aln << std::endl;
		//output << "true_aln" << std::endl << true_aln << std::endl;
		true_aln.comment( "mammoth_alignment" );

		if ( option[ in::file::template_pdb ].user() ) {
			protocols::comparative_modeling::ThreadingMover m(
				true_aln,
				template_pose
			);
			protocols::jobdist::not_universal_main( m );
		}
	}

	using protocols::comparative_modeling::AlignmentSet;
	AlignmentSet align_set;

	for ( Real g_open = min_gap_open; g_open <= max_gap_open;
			g_open += open_step_size
	) {
		for ( Real g_extend = min_gap_extend; g_extend <= max_gap_extend;
				g_extend += extend_step_size
		) {
			ss->gap_open  ( g_open   );
			ss->gap_extend( g_extend );

			SequenceAlignment local_align  = sw_aligner.align( prof1, prof2, ss );
			SequenceAlignment global_align = nw_aligner.align( prof1, prof2, ss );

			if ( option[ in::file::alignment ].user() ) {
				Size n_correct_global
					= n_correctly_aligned_positions( global_align, true_aln );
				Size n_correct_local
					= n_correctly_aligned_positions( local_align, true_aln );

				// create a sub-optimal alignment based on true_aln
				SequenceAlignment temp_aln;
				SequenceProfileOP temp_seq1
					= SequenceProfileOP( new SequenceProfile(*prof1) );
				SequenceProfileOP temp_seq2
					= SequenceProfileOP( new SequenceProfile(*prof2) );

				// it's nice that getting the right alignment is only
				// a matter of putting the gaps into the right places.
				Size insert_pos1( true_aln.sequence(1)->start() ),
					insert_pos2( true_aln.sequence(2)->start() );
				for ( Size i = 1; i <= true_aln.length(); ++i ) {
					if ( true_aln.sequence(1)->is_gap(i) )
						temp_seq1->insert_gap( insert_pos1 );
					else ++insert_pos1;
					if ( true_aln.sequence(2)->is_gap(i) )
						temp_seq2->insert_gap( insert_pos2 );
					else ++insert_pos2;
				}

				temp_aln.add_sequence( temp_seq1 );
				temp_aln.add_sequence( temp_seq2 );
				Real true_score = temp_aln.calculate_score_sum_of_pairs( ss );

				output << "local_align(gap_open = " << g_open << ", "
					" g_extend = " << g_extend << ")";
				output << ", score of optimal aln = " << true_score;
				output << ", n_correct = " << n_correct_local
					<< " n_possible = " << true_aln.length() << std::endl;
				output << local_align;
				output << std::endl;

				output << "global_align (gap_open = " << g_open << ", "
					<< " g_extend = " << g_extend << ")";
				output << ", score of optimal aln = " << true_score;
				output << ", n_correct = " << n_correct_global
					<< " n_possible = " << true_aln.length() << std::endl;
				output << global_align;
				output << std::endl;
			} else {
				output << "local_align (gap_open = " << g_open << ", "
					<< "g_extend = " << g_extend << ")" << std::endl
					<< local_align << std::endl << "--" << std::endl;
				output << "global_align (gap_open = " << g_open << ", "
					<< "g_extend = " << g_extend << ")" << std::endl
					<< global_align << std::endl << "--" << std::endl;
			}

			global_align.comment(
				"nwalign_" + string_of(-1 * g_open) + "_" + string_of(-1 *g_extend)
			);
			local_align.comment(
				"swalign_" + string_of(-1 * g_open) + "_" + string_of(-1 *g_extend)
			);

			// skip alignments that have more than 40% gaps
			if ( local_align.max_gap_percentage() < 0.4 ) {
				std::cout << "local_align.max_gap_percentage() = "
					<< local_align.max_gap_percentage() << std::endl;
				align_set.insert( local_align  );
			}

			if ( global_align.max_gap_percentage() < 0.4 ) {
				std::cout << "global_align.max_gap_percentage() = "
					<< global_align.max_gap_percentage() << std::endl;
				align_set.insert( global_align );
			}

		} // extend
	} // open

	using utility::vector1;
	using core::sequence::SequenceAlignment;
	if ( option[ in::file::template_pdb ].user() ) {
		vector1< SequenceAlignment > aligns = align_set.alignments();
		for ( vector1< SequenceAlignment >::iterator it = aligns.begin(),
				end = aligns.end();
				it != end; ++it
		) {
			// build a threading model of the query pose given the template
			protocols::comparative_modeling::ThreadingMover mover(
				*it,
				template_pose
			);
			protocols::jobdist::not_universal_main( mover );

		} // for ( it in aligns )
	} // if option[ in::file::template_pdb ].user()

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
	}
	
	return 0;

} // int main( int argc, char * argv [] )
