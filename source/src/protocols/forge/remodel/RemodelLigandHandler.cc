// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/forge/remodel/RemodelLoopMover.cc
/// @brief  Loop modeling protocol based on routines from Remodel and EpiGraft
///         packages in Rosetta++.
/// @author Possu Huang (possu@u.washington.edu)

// unit headers
#include <protocols/forge/remodel/RemodelLigandHandler.hh>
#include <protocols/forge/remodel/RemodelRotamerLinks.hh>
#include <protocols/forge/methods/util.hh>

// package headers

// project headers
#include <basic/Tracer.hh>
#include <basic/MetricValue.hh>
#include <core/conformation/Residue.hh>
#include <core/scoring/disulfides/DisulfideMatchingPotential.hh>
#include <core/util/disulfide_util.hh>
#include <core/conformation/util.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/symmetry/util.hh>
#include <core/io/Remarks.hh>
#include <core/pose/util.hh> // for pdbinfo
#include <core/id/AtomID.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/pose/metrics/CalculatorFactory.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/pack/rotamer_set/RotamerLinks.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/remodel.OptionKeys.gen.hh>
#include <core/pack/task/TaskFactory.hh>
#include <protocols/toolbox/pose_metric_calculators/NeighborhoodByDistanceCalculator.hh>
#include <core/pack/task/ResfileReader.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <protocols/simple_moves/symmetry/SetupForSymmetryMover.hh>
#include <protocols/simple_moves/MinMover.hh>
#include <protocols/constraints_additional/BindingSiteConstraint.hh>

// numeric headers
#include <numeric/random/random.hh>
#include <numeric/xyzVector.hh>
#include <numeric/xyzMatrix.hh>

#include <utility/vector0.hh>
#include <utility/vector1.hh>

//external
#include <Eigen/Dense>

// boost headers

// C++ headers
#include <iostream>
#include <math.h>

using namespace basic::options;

namespace protocols {
namespace forge {
namespace remodel {

// Tracer instance for this file
// Named after the original location of this code
static THREAD_LOCAL basic::Tracer TR( "protocols.forge.remodel.RemodelLigandHandler" );

// RNG
//static numeric::random::RandomGenerator RG( 342342 ); // magic number, don't change


// @brief default constructor
RemodelLigandHandler::RemodelLigandHandler()
{
}


/// @brief copy constructor

/// @brief default destructor
RemodelLigandHandler::~RemodelLigandHandler(){}

/// @brief clone this object
protocols::moves::MoverOP
RemodelLigandHandler::clone() const {
	return RemodelLigandHandler::MoverOP( new RemodelLigandHandler( *this ) );
}

/// @brief create this type of object
protocols::moves::MoverOP
RemodelLigandHandler::fresh_instance() const {
	return RemodelLigandHandler::MoverOP( new RemodelLigandHandler() );
}

void RemodelLigandHandler::apply( core::pose::Pose & pose ){
	minimize(pose);
}

void RemodelLigandHandler::minimize( core::pose::Pose & pose )
{
	using namespace basic::options;
	using namespace core::scoring;
	using namespace core;
	using namespace protocols;

	core::kinematics::MoveMapOP movemap( new core::kinematics::MoveMap );

	//handle constraints
	cst_sfx_ = core::scoring::ScoreFunctionFactory::create_score_function(option[OptionKeys::remodel::cen_sfxn] );
	//turn everything off except constraints
	cst_sfx_->set_weight(vdw,0);
	cst_sfx_->set_weight(rg,0);
	cst_sfx_->set_weight(rama,0);
	cst_sfx_->set_weight(hbond_lr_bb,0);
	cst_sfx_->set_weight(hbond_sr_bb,0);
	cst_sfx_->set_weight(omega,0);
	cst_sfx_->set_weight(atom_pair_constraint, 1.0);
	cst_sfx_->set_weight(coordinate_constraint, 1.0);
	cst_sfx_->set_weight(dihedral_constraint, 1.0);

	//enable score functions; initialize to default and then turn on the
	//constraint weights
	fullatom_sfx_ = scoring::get_score_function();
	fullatom_sfx_->set_weight(atom_pair_constraint, 1.0);
	fullatom_sfx_->set_weight(coordinate_constraint, 1.0);
	fullatom_sfx_->set_weight(dihedral_constraint, 1.0);


	//setup the necessary jumps to manipulate
	/*
	Jump tempJump = repeat_pose.jump(i);

	tempJump.set_rotation( Rot );
	tempJump.set_translation( Trx );
	repeat_pose.conformation().set_jump( i, tempJump );
	*/

	//debug check foldtree

	TR << "TREE with Ligand:" << pose.fold_tree() << std::endl;

	//assume ligand is attached to the last jump. won't work in symmetry mode
	Size jump_id = pose.num_jump();
	TR << "pose jumps count:" << jump_id << std::endl;
	movemap->set_jump(jump_id, true);

	//minimize cst_only
	simple_moves::MinMoverOP minMover( new simple_moves::MinMover( movemap , cst_sfx_, "lbfgs_armijo", 0.01, true) );
	minMover->apply(pose);
	//minimize full atom
	simple_moves::MinMoverOP fullminMover( new simple_moves::MinMover( movemap , fullatom_sfx_, "lbfgs_armijo", 0.01, true) );
	fullminMover->apply(pose);

}

std::string
RemodelLigandHandler::get_name() const {
	return "RemodelLigandHandler";
}



} // remodel
} // forge
} // protocol
