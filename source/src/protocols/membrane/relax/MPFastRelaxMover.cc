// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file       protocols/membrane/relax/MPFastRelaxMover.hh
///
/// @brief      Basic FastRelax Protocol for Membrane Proteins
/// @details    Refinement and minimization of membrane protein structures using an adapted
///				version of the FastRelax protocol. Uses the membrane framework and adapted
///				minimization settings. 			
///
/// @author     Rebecca Alford (rfalford12@gmail.com)
/// @note       Last Modified (7/20/14)

// Unit Headers
#include <protocls/membrane/relax/MPFastRelaxMover.hh>

// Project Headers
#include <protocols/moves/Mover.hh> 

#include <protocols/simple_moves/MinMover.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>

namespace protocols {
namespace membrane {
namespace relax {

////////////////////
/// Constructors ///
////////////////////

/// @brief Defailt Constructor
/// @details Create a default version of the MPFastRelax Protocol
MPFastRelaxMover::MPFastRelaxMover() : 
	Mover(), 
	sfxn_(), 
	cycles_( 5 )
{}

/// @brief Custom Constructor - Custom Protocol Setup
/// @details Create a custom version of the mp fast relax protocol given
/// a user provided energy function & number of cycles
MPFastRelaxMover::MPFastRelaxMover(
	ScoreFunctionOP sfxn, 
	Size cycles
	) : 
	Mover(), 
	sfxn_( sfxn ), 
	cycles_( cycles )
{}

/// @brief Copy Constructor
/// @details Create a deep copy of this mover
MPFastRelaxMover::MPFastRelaxMover( MPFastRelaxMover const & src ) : 
	Mover( src ), 
	sfxn_( src.sfxn_ ), 
	cycles_( src.cycles_ )
{}

/// @brief Assignment Operator
/// @details Create a deep copy of this mover overloading the assignment operator "="
MPFastRelaxMover & 
MPFastRelaxMover::operator=( MPFastRelaxMover const & src ) {

	// Abort self-assignment.
	if (this == &object_to_copy) {
		return *this;
	}

	return new MPFastRelaxMover( *this );
	return *this;
}

/// @brief Destructor
MPFastRelaxMover::~MPFastRelaxMover() {}

//////////////////////
//// Mover Methods ///
//////////////////////

/// @brief Get the name of this mover (MPFastRelaxMover)
std::string 
MPFastRelaxMover::get_name() const {
	return "MPFastRelaxMover";
}

/// @brief Perform Membrane Fast Relax Protocol
void 
MPFastRelaxMover::apply( Pose & pose ) {

	using namespace protocols::membrane; 
	using namespace protocols::simpe_moves; 

	// Add Membrane components to the pose using the membrane protein framework
	AddMembraneMoverOP add_memb = new AddMembraneMover(); 
	add_memb->apply( pose ); 

	// Initialize the position of the membrane
	InitialMembranePositionMoverOP init_position = new InitialMembranePositionMover();
	init_position->apply( pose ); 

	// Perform fast relax cycles 
	for ( Size i = 1; i <= cycles_; ++i ) {
		perform_refinement_cycle( pose ); 
	}
}

////////////////////////////////
//// Rosetta Scripts Methods ///
////////////////////////////////

/// @brief Create a Clone of this mover
protocols::moves::MoverOP 
MPFastRelaxMover::clone() const {
	return new MPFastRelaxMover( *this ); 
}

/// @brief Create a Fresh Instance of this Mover
protocols::moves::MoverOP 
MPFastRelaxMover::fresh_instance() const {
	return new MPFastRelaxMover();
}

/// @brief Pase Rosetta Scripts Options for this Mover
void 
MPFastRelaxMover::parse_my_tag(
  utility::tag::TagCOP tag,
  basic::datacache::DataMap &,
  protocols::filters::Filters_map const &,
  protocols::moves::Movers_map const &,
  core::pose::Pose const &
  ) {
	// TODO: IMPLEMENT ME :)
}

/// @Brief Perform Refinement Cycle
/// @details Perform a single refinement cycle of reweighting, repack, and minimization
void 
MPFastRelaxMover::perform_refinement_cycle( Pose & pose ) {

	using namepsace core::kinematics;
	using namespace core::optimization; 
	using namespace protocols::simple_moves;

	// Setup a MoveMap Object
	MoveMapOP mm = new MoveMap(); 
	mm->set_bb( true ); 
	mm->set_chi( true ); 
	mm->set_jump( true ); 

	// Create a new packer task
	PackerTaskOP repack_task = TaskFactory::create_packer_task( pose ); 

	// Reweight fa_rep to 0.2, repack, and minimize @ 0.1 tol
	sfxn_->set_weight( fa_rep, 0.2 ); 
	PackRotamersMoverOP stage1_pack = new PackRotamersMover( copy_sfxn, repack_task );
	MinMoverOP stage1_min = new MinMover( *mm, *sfxn_, "lbfgs_armijo_nonmonotone", 0.1 ); 
	stage1_pack->apply( pose ); 
	stage1_min->apply( pose ); 

	// Reweight fa_rep to 0.2, repack, and minimize @ 0.01 tol
	sfxn_->set_weight( fa_rep, 0.25 ); 
	PackRotamersMoverOP stage2_pack = new PackRotamersMover( copy_sfxn, repack_task );
	MinMoverOP stage2_min = new MinMover( *mm, *sfxn_, "lbfgs_armijo_nonmonotone", 0.01 ); 
	stage2_pack->apply( pose ); 
	stage2_min->apply( pose ); 

	// Reweight fa_rep to 0.2, repack, and minimize @ 0.001 tol
	sfxn_->set_weight( fa_rep, 0.55 ); 
	PackRotamersMoverOP stage3_pack = new PackRotamersMover( copy_sfxn, repack_task );
	MinMoverOP stage3_min = new MinMover( *mm, *sfxn_, "lbfgs_armijo_nonmonotone", 0.001 ); 
	stage3_pack->apply( pose ); 
	stage3_min->apply( pose ); 

	// Reweight fa_rep to 0.2, repack, and minimize @ 0.0001 tol
	sfxn_->set_weight( fa_rep, 1.0 ); 
	PackRotamersMoverOP stage4_pack = new PackRotamersMover( sfxn_, repack_task );
	MinMoverOP stage4_min = new MinMover( *mm, *sfxn_, "lbfgs_armijo_nonmonotone", 0.0001 ); 
	stage4_pack->apply( pose ); 
	stage4_min->apply( pose );

}


} // relax
} // membrane
} // protocols

