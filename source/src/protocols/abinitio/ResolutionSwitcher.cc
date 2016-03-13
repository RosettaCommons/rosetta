// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file KinematicTaskCenter
/// @brief  this class will be handled to a SampleProtocol as a control instance
/// @details responsibilities:
///           know which chainbreaks to penalize and close
///           know which jumps to use during sampling, which (if any) to keep after loop-closing
///           supply a JumpMover if jumps should be moved
///           supply a MoveMap
///           supply a "StrictMoveMap": the protocol should not move anything that is dissallowed in strict_movemap(),
///                      it should try to move just stuff in movemap()
/// should this class also know how to ramp score terms ?
/// handle the titration of constraints ?
/// @author Oliver Lange


// Unit Headers
#include <protocols/abinitio/ResolutionSwitcher.hh>

// Package Headers
#include <protocols/moves/Mover.hh>

// Project Headers
#include <core/chemical/VariantType.hh>


#include <core/conformation/Conformation.hh>

#include <core/pose/Pose.hh>

#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/constraints/ConstraintSet.hh>

#include <protocols/simple_moves/PackRotamersMover.hh>

#include <core/kinematics/MoveMap.hh>

#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>
// Utility Headers
#include <basic/Tracer.hh>

// option key includes


#include <core/chemical/ChemicalManager.fwd.hh>
#include <core/id/SequenceMapping.hh>
#include <core/pose/util.hh>
#include <core/util/SwitchResidueTypeSet.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>


static THREAD_LOCAL basic::Tracer tr( "protocols.general_abinitio", basic::t_info );

namespace protocols {
namespace abinitio {

using namespace core;


//@brief utility function
bool copy_side_chains(
	core::pose::Pose& pose,
	utility::vector1< bool >& needToRepack,
	core::pose::Pose const& fa_input_pose
) {
	// copy sidechain torsions from input pose
	if ( pose.total_residue() != fa_input_pose.total_residue() ) {
		utility_exit_with_message("Mismatch of pose lenght in copy_side_chains(..): " );
	}
	tr.Debug << "copy side chains for residues with * / missing density residues with - ";
	for ( core::Size i=1; i<=pose.total_residue(); ++i ) {
		// if there is missing density in the sidechain, then we need to repack
		// check this my making sure that no SC atom is more than 20A (?) away from CA
		//needToRepack[i] = true;
		numeric::xyzVector< core::Real> ca_pos = fa_input_pose.residue(i).atom("CA").xyz();
		for ( Size j=1; j<=fa_input_pose.residue(i).natoms(); ++j ) {
			if ( (ca_pos - fa_input_pose.residue(i).atom(j).xyz()).length() > 20 ) {
				tr.Debug << "-" << i << " ";
				needToRepack[ i ] = true;
				break; //one bad atom is enough
			}
		}
		//copy sidechains only for non-loop regions
		if ( !needToRepack[ i ] ) {
			tr.Debug <<  "*" << i << " ";
			bool const lower_cut ( pose.residue( i ).has_variant_type( chemical::CUTPOINT_LOWER ) );
			bool const upper_cut ( pose.residue( i ).has_variant_type( chemical::CUTPOINT_UPPER ) );
			pose.replace_residue(i, fa_input_pose.residue(i), true /*orient backbone*/ );
			if ( lower_cut ) core::pose::add_variant_type_to_pose_residue( pose, chemical::CUTPOINT_LOWER, i );
			if ( upper_cut ) core::pose::add_variant_type_to_pose_residue( pose, chemical::CUTPOINT_UPPER, i );
		}
	} //for 1..total_residue
	tr.Debug << " that have not moved from template" << std::endl;
	return true;
} //copy_side_chains

ResolutionSwitcher::ResolutionSwitcher(
	core::pose::Pose const& pose,
	bool /*input pose*/,
	bool start_centroid,
	bool apply_to_centroid
) : apply_to_centroid_( apply_to_centroid ),
	init_pose_( pose ),
	init_fa_( pose.is_fullatom() ),
	start_centroid_( start_centroid ),
	scorefxn_fa_( /* NULL */ ),
	repack_buffer_( 0 ),
	map_cst_from_centroid_to_fa_( true )
{}

//@brief return the pose to start the simulation  -- based on the flags it will be centroid or full-atom
core::pose::Pose ResolutionSwitcher::start_pose() const {
	if ( start_centroid_ ) {
		if ( init_fa_ ) {
			pose::Pose pose( init_pose_ );
			core::util::switch_to_residue_type_set( pose, chemical::CENTROID );
			return pose;
		} else return init_pose_;
	} else if ( init_fa_ ) {
		return init_pose_;
	} else {
		utility_exit_with_message("don't have full-atom pose to start");
		return init_pose_; //to make compiler happy.
		//or switch to full-atom?
	}
}

//@brief switch the pose from its end-of-sampling-resolution to full-atom
// copy sidechains from initial pose if torsions have not moved
void ResolutionSwitcher::apply( pose::Pose &pose ) {
	tr.Debug << "resolution_switch" << std::endl;
	bool const bCopySideChains( start_centroid_ || apply_to_centroid_ ); //we assume that start and apply pose are fullatom that side-chains are fine!
	utility::vector1< bool > needToRepack( pose.total_residue() , false );
	if ( apply_to_centroid_ ) {
		//puts full-atom sidechains on loop regions
		tr.Debug <<" change to full-atom pose " << std::endl;
		core::util::switch_to_residue_type_set( pose, core::chemical::FA_STANDARD );
		pose.conformation().detect_bonds();//apl fix this !
	}

	if ( init_fa_ ) { //find residues that have moved --- test for missing density is in copy_side_chains!
		for ( Size i=1; i<=pose.total_residue(); ++i ) {
			if ( std::abs(init_pose_.phi( i ) - pose.phi( i ) ) > 10
					|| std::abs(init_pose_.psi( i ) - pose.psi( i ) ) > 10 ) {
				tr.Trace << "residue " << i << " has moved " << std::endl;
				for ( Size j = std::max( 1, (int)i-(int) repack_buffer_); j <= pose.total_residue() && j <= i + repack_buffer_; j++ ) {
					needToRepack[ j ] = true;
				} //saftey buffer for repacking: 1 residues on each side of moved stuff
			}
		}
	}

	if ( bCopySideChains && init_fa_ ) {
		copy_side_chains( pose, needToRepack, init_pose_ );
	}

	if ( !scorefxn_fa_ ) scorefxn_fa_ = core::scoring::get_score_function();
	scorefxn_fa_->set_weight(  scoring::coordinate_constraint , 1.0 );
	runtime_assert( pose.is_fullatom() );
	// repack loop + missing-density residues
	core::pack::task::PackerTaskOP taskstd = core::pack::task::TaskFactory::create_packer_task( pose );
	taskstd->restrict_to_repacking();
	taskstd->or_include_current(true);
	if ( init_fa_ ) {
		taskstd->restrict_to_residues( needToRepack );
	}
	if ( !init_pose().is_fullatom() ) {
		if ( map_cst_from_centroid_to_fa_ ) {
			//&& basic::options::option[ basic::options::OptionKeys::chemical::patch_selectors ].user()
			//  && basic::options::option[ basic::options::OptionKeys::chemical::patch_selectors ]()[ 1 ] == "CENTROID_HA" ) {
			//  tr.Warning << "[ WARNING ] ResolutionSwitcher cannot copy constraints if CENTROID_HA is selected... HA has other atom-number... " << std::endl
			//        << "clean way would be to make the constraint set independent of numbers as in Templates.cc " << std::endl;
			//  utility_exit_with_message("[ WARNING ] ResolutionSwitcher cannot copy constraints if CENTROID_HA is selected.");
			pose.constraint_set( init_pose().constraint_set()->remapped_clone( init_pose(), pose ) );
		}
	} else {
		pose.constraint_set( init_pose().constraint_set()->clone() );
	}
	protocols::simple_moves::PackRotamersMover pack1( scorefxn_fa_ , taskstd );
	pack1.apply( pose );

	// quick SC minimization
	core::optimization::AtomTreeMinimizer mzr;
	core::optimization::MinimizerOptions options( "lbfgs_armijo_nonmonotone", 1e-5, true, false );
	core::kinematics::MoveMap mm;
	mm.set_bb( false );
	mm.set_chi( true );
	mzr.run( pose, mm, *scorefxn_fa_, options );

}

std::string
ResolutionSwitcher::get_name() const {
	return "ResolutionSwitcher";
}


}
}
