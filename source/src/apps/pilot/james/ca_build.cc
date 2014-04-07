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

#include <core/id/NamedAtomID.hh>

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
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <core/pose/Pose.hh>

#include <core/io/pdb/pose_io.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/SilentStruct.hh>
#include <core/io/silent/SilentStructFactory.hh>

#include <protocols/comparative_modeling/ThreadingMover.hh>
#include <protocols/comparative_modeling/StealLigandMover.hh>
#include <protocols/jobdist/standard_mains.hh>
#include <protocols/jobdist/not_universal_main.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/moves/CompositionMover.hh>
#include <protocols/simple_moves/MinMover.hh>

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

	for ( Real g_open = min_gap_open; g_open <= max_gap_open; g_open += open_step_size ) {
		for ( Real g_extend = min_gap_extend; g_extend <= max_gap_extend; g_extend += extend_step_size ) {
			std::cout << "evaluating " << g_open << "," << g_extend << std::endl;
			ss->gap_open  ( g_open   );
			ss->gap_extend( g_extend );

			SequenceAlignment local_align  = sw_aligner.align( prof1, prof2, ss );
			SequenceAlignment global_align = nw_aligner.align( prof1, prof2, ss );

			output << "local_align (gap_open = " << g_open << ", "
				<< "g_extend = " << g_extend << ")" << std::endl
				<< local_align << "--" << std::endl;
			output << "global_align (gap_open = " << g_open << ", "
				<< "g_extend = " << g_extend << ")" << std::endl
				<< global_align << "--" << std::endl;

			global_align.comment(
				"nwalign_" + string_of(-1 * g_open) + "_" + string_of(-1 *g_extend)
			);
			local_align.comment(
				"swalign_" + string_of(-1 * g_open) + "_" + string_of(-1 *g_extend)
			);

			if ( option[ in::file::template_pdb ].user() ) {
				// build a threading model of the query pose given the template
				protocols::moves::CompositionMover l_container;

				l_container.add_mover(
					new protocols::comparative_modeling::ThreadingMover(
						local_align,
						template_pose
					)
				);

				// mapping from template -> query
				core::id::SequenceMapping map = local_align.sequence_mapping(2,1);
				core::Size t_resi( 59 ); // ALA59 in 2a62
				core::Size q_resi( map[ t_resi ] );
				core::id::NamedAtomID query_anchor( "CA", q_resi );
				core::id::NamedAtomID template_anchor( "CA", t_resi );
				utility::vector1< core::id::NamedAtomID > ligand_indices;
				//ligand_indices.push_back( core::id::NamedAtomID( "CA", 320 ) );
				//ligand_indices.push_back( core::id::NamedAtomID( "CA", 321 ) );
				//ligand_indices.push_back( core::id::NamedAtomID( "CA", 322 ) );
				ligand_indices.push_back( core::id::NamedAtomID( "CA", 211 ) );
				ligand_indices.push_back( core::id::NamedAtomID( "CA", 212 ) );
				ligand_indices.push_back( core::id::NamedAtomID( "CA", 213 ) );

				if ( q_resi == 0 ) {
					std::cerr << "Warning: no mapping for template residue " << t_resi << "." << std::endl;
					std::cerr << "Mapping:" << std::endl << map << std::endl;
				} else {
					l_container.add_mover(
						new protocols::comparative_modeling::StealLigandMover(
							template_pose,
							query_anchor,
							template_anchor,
							ligand_indices
						)
					);
				}

				//core::scoring::ScoreFunctionOP scorefxn = core::scoring::getScoreFunction();
				//l_container.add_mover(
				//	new protocols::simple_moves::MinMover(
				//		new core::kinematics::MoveMap, scorefxn, "dfpmin", 1e-5, true
				//	)
				//);

				//protocols::comparative_modeling::ThreadingMover g_mover(
				//	global_align,
				//	template_pose
				//);

				// skip alignments that have more than 40% gaps
				if ( local_align.max_gap_percentage() < 0.4 ) {
					//std::cout << "local_align.max_gap_percentage() = "
					//	<< local_align.max_gap_percentage() << std::endl;
					protocols::jobdist::not_universal_main( l_container );
				} else {
					std::cout << "rejected alignment with "
						<< local_align.max_gap_percentage() << " percent gaps." << std::endl;
				}
			}
		} // gap_extend
	} // gap_open

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;
} // int main( int argc, char * argv [] )
