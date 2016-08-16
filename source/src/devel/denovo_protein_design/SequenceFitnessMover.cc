// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file devel/DenovoProteinDesign/SequenceFitnessMover.cc
/// @brief SequenceFitnessMover methods implemented - this is a fast way of evaluating a designed sequences
/// @brief preference for the native topologoy if you have relax/abinitio decoys from the wildtype
/// @brief - turns out that mini abinitio is much faster than R++, this may not be much of an improvement
/// @brief - given that you lose a lot of searching
/// @author Grant


// Unit Headers
#include <devel/denovo_protein_design/SequenceFitnessMover.hh>
#include <devel/denovo_protein_design/SequenceFitnessMover.fwd.hh>

// Package Headers
#include <protocols/moves/Mover.hh>
#include <protocols/moves/Mover.fwd.hh> //MoverOP

// Project Headers
#include <core/pose/Pose.hh>


#include <core/import_pose/import_pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreType.hh>


#include <basic/options/option.hh>
#include <core/scoring/packstat/compute_sasa.hh>

#include <core/scoring/rms_util.hh>

// Utility Headers
#include <basic/Tracer.hh>
#include <core/types.hh>
#include <utility/exit.hh>

// option key includes

#include <basic/options/keys/in.OptionKeys.gen.hh>

#include <protocols/relax/FastRelax.hh>
#include <utility/vector1.hh>


// C++ Headers

using basic::T;
using basic::Error;
using basic::Warning;
using namespace core;
static THREAD_LOCAL basic::Tracer TR( "devel.DenovoProteinDesign.SequenceFitnessMover" );

namespace devel {
namespace denovo_protein_design {

/// @details
void SequenceFitnessMover::apply( core::pose::Pose & pose ){
	using namespace core::scoring::packstat;
	core::scoring::ScoreFunctionOP fullfxn( ( core::scoring::get_score_function_legacy( core::scoring::PRE_TALARIS_2013_STANDARD_WTS ) ));

	core::pose::Pose testpose;
	core::import_pose::pose_from_file( testpose, basic::options::option[ basic::options::OptionKeys::in::file::native ]().name() , core::import_pose::PDB_file);

	core::pose::Pose opose = testpose;

	// set phi, psi, omega of pose to test pose

	for ( core::Size seqpos = 1; seqpos <= testpose.n_residue(); ++seqpos ) {
		testpose.set_phi   ( seqpos, pose.phi(seqpos) );
		testpose.set_psi   ( seqpos, pose.psi(seqpos) );
		testpose.set_omega ( seqpos, pose.omega(seqpos) );
	}

	protocols::relax::FastRelax fast_relax( fullfxn );
	fast_relax.apply( testpose );

	float RMSD_testpose_opose = core::scoring::native_CA_rmsd ( testpose, opose );
	//  float RMSD_pose_opose = protocols::evaluation::native_CA_rmsd ( pose, opose );


	PosePackData pd = pose_to_pack_data( pose, false );

	core::Real packing_score = compute_packing_score( pd, 0 );


	Energy score = (*fullfxn)( testpose );
	/// Now handled automatically.  fullfxn->accumulate_residue_total_energies( testpose );

	std::cout << " Grant Scores "<< score <<" "<< RMSD_testpose_opose << " " << packing_score << std::endl;


	//  protocols::analysis::PackStatMover packstats();
	//  packstats.apply( testpose  );


}//apply

std::string
SequenceFitnessMover::get_name() const {
	return "SequenceFitnessMover";
}

/// @brief
SequenceFitnessMover::SequenceFitnessMover(
) : Mover()
{
	Mover::type( "SequenceFitnessMover" );
}

SequenceFitnessMover::~SequenceFitnessMover(){}

}//DenovoProteinDesign
}//devel

