// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @brief stupid test file for visual studio c++
/// @details

#include <core/types.hh>
#include <core/pose/Pose.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <core/sequence/util.hh>
#include <core/sequence/Sequence.hh>

#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/SilentStruct.hh>
#include <core/io/silent/ProteinSilentStruct.hh>


#include <core/chemical/ChemicalManager.hh>

#include <core/kinematics/MoveMap.hh>

#include <basic/options/option.hh>

#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/TrialMover.hh>
#include <protocols/simple_moves/BackboneMover.hh>

#include <utility/io/ozstream.hh>

#include <string>


// option key includes

#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/james.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

#include <core/io/silent/ProteinSilentStruct.tmpl.hh>
#include <utility/vector1.hh>
#include <numeric/random/random.fwd.hh>


using namespace basic::options;
using namespace basic::options::OptionKeys;
using namespace ObjexxFCL::format;

namespace protocols {
namespace abinitio {

// the following code leaks memory on Windows if compiled with Visual Studio 2005. The fix is to
// upgrade to Service Pack 1. The important bit is the F() statement, which uses a C++ iostream
// to do its formatted output. See here for more information:
// http://connect.microsoft.com/VisualStudio/feedback/ViewFeedback.aspx?FeedbackID=98861
// tex 2008-09-09
int run_boinc_debug() {
	std::string sequence(
		core::sequence::read_fasta_file(
			option[ in::file::fasta ]()[1]
		)[1]->sequence()
	);
	core::pose::Pose fold_pose;
	core::pose::make_pose_from_sequence(
		fold_pose,
		sequence,
		*( core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::CENTROID ) )
 	);

	core::scoring::ScoreFunctionOP scorefxn( core::scoring::get_score_function() );

	utility::io::ozstream log( "debug_log.txt" );
	// fragments
	//std::string frag_small_file  = option[ in::file::frag3 ]();
	//ConstangLengthFragSetOP fragset_small = new ConstantLengthFragSet;
	//fragset_small->read_fragment_file( frag_small_file, option[ OptionKeys::abinitio::number_3mer_frags ]() );

	// make a MoveMap
	core::kinematics::MoveMapOP movemap( new core::kinematics::MoveMap );
	movemap->set_bb( true );

	// nstruct iterations
	core::Size iterations = option[ out::nstruct ]();

	using core::Size;
	using namespace protocols::moves;
	MonteCarloOP mc( new MonteCarlo( fold_pose, *scorefxn, 1.0 ) );
	simple_moves::SmallMoverOP small_mover( new simple_moves::SmallMover( movemap, 1.0, 1 ) );
	moves::TrialMoverOP smooth_trials( new moves::TrialMover( small_mover, mc ) );

	if ( option[ james::debug ]() ) {
		core::io::silent::SilentFileData sfd;

		using numeric::random::gaussian;
		for ( core::Size i = 1; i <= iterations; ++i ) {
			// wiggle all of the residues a bit.
			for ( core::Size pos = 1; pos <= fold_pose.total_residue(); pos++ ) {
				fold_pose.set_phi  ( pos, fold_pose.phi( pos )   + gaussian() );
				fold_pose.set_psi  ( pos, fold_pose.psi( pos )   + gaussian() );
				fold_pose.set_omega( pos, fold_pose.omega( pos ) + gaussian() );
			}

			(*scorefxn)(fold_pose);
			//			core::Real score = (*scorefxn)(fold_pose);
			//mc.boltzmann( fold_pose );

			log << "iteration " << i << std::endl;

			log << "initializing ProteinSilentStruct!" << std::endl;
			core::io::silent::ProteinSilentStruct pss( fold_pose );
			log << "calling write_silent_struct!" << std::endl;
			sfd.write_silent_struct( pss, "test.out" );
		} // for i
	} else {
		for ( core::Size i = 1; i <= iterations * 100; ++i ) {
			log << F( 8, 3, 3.14159 ) << std::endl;
		}
	} // else

	log.close();

	return 0;
}


} // abinitio
} // protocols
