// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/loops/loop_closure/jacobi/JacobiLoopClosureMover.cc
/// @brief Loop closure mover that uses linearization of the kinematics to adjust all unlocked backbone torsion angles
/// of a loop at the same time
/// @author teunhoevenaars (teunhoevenaars@gmail.com)

// Unit headers
#include <protocols/loops/loop_closure/jacobi/JacobiLoopClosureMover.hh>
#include <protocols/loops/loop_closure/jacobi/JacobiLoopClosureMoverCreator.hh>

// Core headers
#include <core/chemical/ResidueConnection.hh>
#include <core/conformation/Residue.hh>
#include <core/kinematics/jacobian/JacobianStructure.hh>
#include <core/kinematics/jacobian/SeriesJacobians.hh>
#include <core/kinematics/jacobian/ModuleType1.hh>
#include <core/kinematics/jacobian/util.hh>
#include <core/kinematics/AtomTree.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/optimization/MinimizerMap.hh>
#include <core/optimization/types.hh>
#include <core/pose/Pose.hh>

// Protocol headers
#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/loops/loops_main.hh>
#include <numeric/MathMatrix_operations.hh>
#include <numeric/MathVector_operations.hh>

// Basic/Utility headers
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>
#include <utility/pointer/memory.hh>

// XSD Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>
#include <protocols/loop_modeling/utilities/rosetta_scripts.hh>

#include <numeric/numeric.functions.hh> // AUTO IWYU For min
#include <utility/fixedsizearray1.hh> // AUTO IWYU For fixedsizearray1

static basic::Tracer TR( "protocols.loops.loop_closure.jacobi.JacobiLoopClosureMover" );

namespace protocols {
namespace loops {
namespace loop_closure {
namespace jacobi {

/////////////////////
/// Constructors  ///
/////////////////////

/// @brief Default constructor
JacobiLoopClosureMover::JacobiLoopClosureMover():
	protocols::moves::Mover( JacobiLoopClosureMover::mover_name() )
{

}

////////////////////////////////////////////////////////////////////////////////
/// @brief Constructor based on a loop
JacobiLoopClosureMover::JacobiLoopClosureMover( protocols::loops::Loop const & loop ){

	if ( loop.start() >= loop.stop() ) {
		utility_exit_with_message("Loop is incorrectly defined: loop.start() >= loop.stop()");
	}

	// create default MoveMap
	core::kinematics::MoveMapOP mm = create_default_movemap(loop);

	// run the constructor initialization
	init_constructor( loop, mm );
}

////////////////////////////////////////////////////////////////////////////////
/// @brief Constructor based on a loop and a MoveMap
JacobiLoopClosureMover::JacobiLoopClosureMover( protocols::loops::Loop const & loop, core::kinematics::MoveMapCOP const & mm ){

	if ( loop.start() >= loop.stop() ) {
		utility_exit_with_message("Loop is incorrectly defined: loop.start() >= loop.stop()");
	}

	// check that omega angles are fixed
	for ( core::uint i = loop.start(); i <= loop.stop(); ++i ) {
		if ( mm->get_bb( i, 3 ) ) {
			if ( verbose_ ) {
				TR.Info << "MoveMap contains free omega angle in residue " << i << " which will be ignored"
					<< std::endl;
			}
		}
	}

	// run the constructor initialization
	init_constructor( loop, mm );
}

////////////////////////////////////////////////////////////////////////////////
/// @brief Copy constructor
JacobiLoopClosureMover::JacobiLoopClosureMover( JacobiLoopClosureMover const & ) = default;

////////////////////////////////////////////////////////////////////////////////
/// @brief Destructor (important for properly forward-declaring smart-pointer members)
JacobiLoopClosureMover::~JacobiLoopClosureMover(){}


////////////////////////////////////////////////////////////////////////////////
/// Mover Methods ///
/////////////////////

/// @brief Set the loop to be closed
void
JacobiLoopClosureMover::set_loop( protocols::loops::Loop const & new_loop )
{
	// copy current MoveMap
	core::kinematics::MoveMapCOP mm = movemap_;
	// if current MoveMap is not set, create a default one
	if ( movemap_ == nullptr ) {
		mm = create_default_movemap(new_loop);
	}
	// re-run init_constructor
	init_constructor( new_loop, mm  );
}

/// @brief Set the MoveMap
void
JacobiLoopClosureMover::set_movemap( core::kinematics::MoveMapOP const & new_mm )
{
	// copy current loop
	protocols::loops::Loop loop = *loop_;
	// re-run init_constructor
	init_constructor( loop, new_mm );
}

/// @brief Set norm-values of allowed rotational [deg] and linear error vectors [Ang] that define successful closure
/// @details error_norm_rot should be supplied in degrees.
void
JacobiLoopClosureMover::set_error_norms(core::Real const error_norm_rot, core::Real const error_norm_lin)
{
	// check that provided norms are positive
	if ( error_norm_rot <= 0 || error_norm_lin <= 0 ) {
		TR.Error << "One of (or both) the provided norms are <= 0, while they need to be positive. Norms kept to current value" <<
			std::endl;
	} else {
		// convert provided norms to radians for use in the protocol
		err_rot_allowed_ = error_norm_rot * numeric::constants::d::degrees_to_radians;
		err_lin_allowed_ = error_norm_lin;
	}
}

/// @brief Apply the mover
void
JacobiLoopClosureMover::apply( core::pose::Pose & pose) {
	using namespace core::kinematics::jacobian;

	// Rest closure success to false
	closure_success_ = false;

	// Initialize the Jacobian modules
	init_apply( pose );

	// number of modules in the chain
	core::Size N_modules = jacobian_chain_->modules_.size();

	// Store torsion angles of input pose, so they can be restored in case loop closure does not succeed
	core::optimization::MinimizerMap minmap;
	minmap.setup( pose, *movemap_ );
	core::optimization::Multivec recover_vars(minmap.nangles());
	minmap.copy_dofs_from_pose(pose, recover_vars);

	// Create current position and Euler angles, expressed in ref. frame connected to N-atom of residue 1 of loop
	core::kinematics::RT current = update_current(pose);

	// calculate error twist, which is a 0-vector [dphi_x, dphi_y, dphi_z, dx, dy, dz], where dphi_i is the projection of the
	// update Euler angle errors on the Cartesian axes
	ModuleType1::Screw dT_as_twist = calculate_error_twist( current, target_ );

	// also create empty dT as 0-vector for use in loop
	numeric::MathVector<core::Real> dT(6);

	// extract separate errors from the twist
	core::Real err_rot = dT_as_twist.first.norm();
	core::Real err_lin = dT_as_twist.second.norm();

	// initialize iteration counter
	core::Size cycles_counter = 0;

	// START CLOSURE CYCLES
	while ( (err_rot > err_rot_allowed_ || err_lin > err_lin_allowed_) && cycles_counter < max_cycles_ ) {

		// translate dT into MathVector for linear calculations
		for ( core::Size i=0; i < 3; ++i ) {
			dT(i) = dT_as_twist.first[i];
			dT(3+i) = dT_as_twist.second[i];
		}

		// get structs with Jacobian matrices for the different modules
		utility::vector1< ModuleType1::jacobian_struct > jacobian_matrices_chain (N_modules);
		for ( core::Size i=1; i <= N_modules; ++i ) {
			jacobian_matrices_chain[i] = jacobian_chain_->modules_[i]->get_jacobian_matrices(pose.conformation());
		}

		// Initially pretend that all modules have 6-DoFs
		utility::vector1< numeric::MathVector<core::Real> > dq_all_6dof(N_modules);

		// Calculate required change in dihedral angles if each Jacobian module would contribute maximally and span 6-Dof space
		for ( core::Size i=1; i <= N_modules; ++i ) {
			dq_all_6dof[i] = jacobian_matrices_chain[i].inv_all * dT;
		}
		// weigh the contributions
		weigh_columns_inversely_squared(dq_all_6dof);

		// create empty correction vectors and final vector
		utility::vector1< numeric::MathVector<core::Real> > dq_w( N_modules );

		if ( jacobian_chain_->modules_.back()->get_number_dofs() == 6 ) {
			// the final vector of vectors is simply a copy of the earlier calculated vectors
			dq_w = dq_all_6dof;
		} else { // if the last module has <6 DoF, correct for that
			utility::vector1< numeric::MathVector<core::Real> > dq_all_correct( N_modules-1 );

			// error is the forward projection of the part attributed to the residues that are not free
			numeric::MathVector<core::Real> dT_err = jacobian_matrices_chain.back().fw_cons * dq_all_6dof[N_modules];
			// divide error over all other modules (if they exist)
			for ( core::Size i = 1; i <= N_modules - 1; ++i ) {
				dq_all_correct[i] = jacobian_matrices_chain[i].inv_dofs * dT_err;
			}
			// get weighted distribution of the correction torsion angles
			weigh_columns_inversely_squared(dq_all_correct);
			// and store the negative twist as correction in the last module (so that sum of dq_all_6dof_w + dq_all_correct_w
			// for last entry is now zero
			dq_all_correct.push_back(jacobian_matrices_chain.back().inv_cons * -dT_err);
			// add contributions for each entry in the module vector
			for ( core::Size i = 1; i <= N_modules; ++i ) {
				dq_w[i] = dq_all_6dof[i] + dq_all_correct[i];
			}
		}

		// apply delta values to pose. dq_w is a matrix (with 0 index), while the modules are 1-indexed
		for ( core::Size i=1; i <= N_modules; ++i ) { // cycle through Jacobian modules
			jacobian_chain_->modules_[i]->apply_delta_torsion_values(pose.conformation(), dq_w[i]);
		}

		// update current position and Euler angles, expressed in ref. frame connected to N-atom of residue 1 of loop
		current = update_current(pose);

		// calculate error twist, which is a 0-vector [dphi_x, dphi_y, dphi_z, dx, dy, dz], where dphi_i is the projection of the
		// update Euler angle errors on the Cartesian axes
		dT_as_twist = calculate_error_twist( current, target_ );

		// update errors from the twist
		err_rot = dT_as_twist.first.norm();
		err_lin = dT_as_twist.second.norm();

		// increment counter
		++cycles_counter;
	}

	// store closure cycles
	closure_cycles_ = cycles_counter;

	// report on outcome of closure
	if ( cycles_counter < max_cycles_ ) { // closure successful
		// report the success in the mover
		closure_success_ = true;
		if ( verbose_ ) {
			TR.Info << "Jacobi loop closure successfully completed in " << closure_cycles_ << " cycles" << std::endl;
		}
	} else { // closure NOT successful
		// restore original torsion angles
		minmap.copy_dofs_to_pose(pose, recover_vars);
		TR.Debug << "Loop closure NOT successful, original torsion angles have been restored."<<std::endl;
	}
}

////////////////////////////////////////////////////////////////////////////////
/// @brief Show the contents of the Mover
void
JacobiLoopClosureMover::show(std::ostream & output) const
{
	protocols::moves::Mover::show(output);
}

////////////////////////////////////////////////////////////////////////////////
/// Rosetta Scripts Support ///
///////////////////////////////

/// @brief parse XML tag (to use this Mover in Rosetta Scripts)
void
JacobiLoopClosureMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap&
) {
	max_cycles_ = tag->getOption< core::Size >( "max_cycles" );
	err_rot_allowed_ = tag->getOption< core::Real >( "err_rot_allowed" ) * numeric::constants::d::degrees_to_radians;
	err_lin_allowed_ = tag->getOption< core::Real >( "err_lin_allowed" );
	verbose_ = tag->getOption< bool >( "verbose" );

	// extract loop and initialize closure mover
	protocols::loops::LoopsOP parsed_loops(protocols::loop_modeling::utilities::parse_loops_from_tag(tag));
	// confirm loop is provided
	if ( parsed_loops == nullptr ) { utility_exit_with_message("Input does not contain a loop: please add a Loop subtag"); }
	// confirm single loop has been provided
	if ( parsed_loops->size() != 1 ) { utility_exit_with_message("Input does not contain exactly one loop"); }

	// copy loop and generate MoveMap based on loop
	protocols::loops::Loop loop = *parsed_loops->begin();
	core::kinematics::MoveMapOP mm = create_default_movemap(loop);

	// run the constructor initialization
	init_constructor( loop, mm );
}

void JacobiLoopClosureMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;

	// here you should write code to describe the XML Schema for the class.  If it has only attributes, simply fill the provided AttributeList.
	// First, add loop to attribute list
	utility::tag::XMLSchemaSimpleSubelementList subelements;
	loop_modeling::utilities::append_subelement_and_attributes_for_parse_loops_from_tag( xsd, subelements, attlist );

	attlist
		+ XMLSchemaAttribute::attribute_w_default("max_cycles", xsct_non_negative_integer, "Maximum number of cycles per run. Quit the run even if not converged at this point.","200")
		+ XMLSchemaAttribute::attribute_w_default("err_rot_allowed", xsct_real, "Maximum norm of rotational error [deg] for closure to be accepted as successful.", "5")
		+ XMLSchemaAttribute::attribute_w_default("err_lin_allowed", xsct_real, "Maximum norm of linear error [Ang] for closure to be accepted as successful.", "0.1")
		+ XMLSchemaAttribute::attribute_w_default("verbose", xsct_rosetta_bool, "Setting that determines whether or not to output non-critical information.", "false");

	protocols::moves::xsd_type_definition_w_attributes_and_repeatable_subelements( xsd, mover_name(), "Performs loop closure on a single loop using the Jacobi algorithm", attlist, subelements );
}


////////////////////////////////////////////////////////////////////////////////
/// @brief required in the context of the parser/scripting scheme
protocols::moves::MoverOP
JacobiLoopClosureMover::fresh_instance() const
{
	return utility::pointer::make_shared< JacobiLoopClosureMover >();
}

/// @brief required in the context of the parser/scripting scheme
protocols::moves::MoverOP
JacobiLoopClosureMover::clone() const
{
	return utility::pointer::make_shared< JacobiLoopClosureMover >( *this );
}

std::string JacobiLoopClosureMover::get_name() const {
	return mover_name();
}

std::string JacobiLoopClosureMover::mover_name() {
	return "JacobiLoopClosureMover";
}



/////////////// Creator ///////////////

protocols::moves::MoverOP
JacobiLoopClosureMoverCreator::create_mover() const
{
	return utility::pointer::make_shared< JacobiLoopClosureMover >();
}

std::string
JacobiLoopClosureMoverCreator::keyname() const
{
	return JacobiLoopClosureMover::mover_name();
}

void JacobiLoopClosureMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	JacobiLoopClosureMover::provide_xml_schema( xsd );
}

////////////////////////////////////////////////////////////////////////////////
/// private methods ///
///////////////////////


std::ostream &
operator<<( std::ostream & os, JacobiLoopClosureMover const & mover )
{
	mover.show(os);
	return os;
}

/// @brief Create default MoveMap
core::kinematics::MoveMapOP
JacobiLoopClosureMover::create_default_movemap( protocols::loops::Loop loop) {
	core::kinematics::MoveMapOP mm(utility::pointer::make_shared<core::kinematics::MoveMap>());
	// with all backbone atoms in the loop moveable
	mm->set_bb_true_range(loop.start(), loop.stop());
	// except the omega angles
	for ( core::uint i = loop.start(); i <= loop.stop(); ++i ) {
		mm->set(core::id::TorsionID(i, core::id::BB, 3), false);
	}

	return mm;
}

/// @brief Initializes data members from constructor input arguments
void
JacobiLoopClosureMover::init_constructor( protocols::loops::Loop const & loop, core::kinematics::MoveMapCOP const & mm ) {
	// create copy of loop object that is managed by the closure mover, but can be passed around to upstream movers
	loop_ = utility::pointer::make_shared< protocols::loops::Loop >(loop);
	// check that cutpoint is defined at loop.stop and if not, change that and give a warning
	if ( loop_->cut() != loop_->stop() ) {
		loop_->set_cut(loop_->stop());
		TR.Info << "Jacobi protocols require loop.cut()=loop.stop(), input loop has been corrected." << std::endl;
	}

	// store pointer to the mm;
	movemap_ = mm;

	// create vector with freely moveable residues from loop and movemap. Make sure free_residues_ is empty to start
	free_residues_.clear();
	for ( core::Size i = loop_->start(); i <= loop_->stop(); ++i ) {
		// if phi and psi backbone torsions are moveable, amend to vector
		if ( movemap_->get_bb(i,1) && movemap_->get_bb(i,2) ) {
			free_residues_.push_back(i);
		}
	}

	// reset the jacobian module because it is no longer valid
	jacobian_chain_ = nullptr;
}

/// @brief Initializes private parameters of the Jacobi mover that depend on the pose.
void
JacobiLoopClosureMover::init_apply( core::pose::Pose & pose ) {
	// Check for amount of free residues in the mover and give warning / error
	if ( free_residues_.size() < 3 ) {
		// if this is the case, first check if the vector is empty
		if ( free_residues_.empty() ) {
			utility_exit_with_message("The loop combined with the MoveMap leaves 0 free residues, please correct loop and/or MoveMap");
		}
		// Otherwise, give warning
		TR.Warning
			<< "The loop combined with the MoveMap has < 3 free residues: likely the loop cannot be closed" << std::endl;
	}

	// Check that assist molecules can be defined in case a loop with less than three residues is defined
	core::kinematics::FoldTree const & ft( pose.fold_tree() );
	if ( free_residues_.size() < 3 &&
			ft.get_residue_edge(free_residues_[1]).start() > (free_residues_[1] - free_residues_.size() - 3) ) {
		utility_exit_with_message(
			"Because loop has < 3 free residues, the start of the edge that includes the first residue of the loop must be > 3 residues before loop.start(): correct FoldTree");
	}

	// prepare FoldTree and move cutpoint if necessary
	prepare_foldtree(pose);

	// Store the internal coordinates (most notably bond angles) that describe the position and orientation of the end of the loop.
	// If 'current-as-target' was called before, the target has already been initialized. Otherwise, use residue defaults
	// default values are stored from the residue_types. The function sets icoor_target_bb_
	if ( !target_initialized_ ) {
		// initialize the atom_ids and icoors
		default_target_icoors(pose);
		// calculate the target
		store_target(pose);
	}

	// Initialize Jacobian object in case it has not been initialized
	if ( jacobian_chain_ == nullptr ) {
		// hardcoded with omega fixed. If you want to include omega as a variable in Jacobi loop-closure, it will also
		// be necessary to revise the 'current_closed_as_target' and 'update_current' definitions, because they now assume that
		// omega remains unchanged throughout the loop closure
		jacobian_chain_ = core::kinematics::jacobian::JacobianStructure(pose.conformation(), free_residues_,
			ref_atom_id_).get_single_chain(1);
	}
}

/// @brief Function that checks and adjusts the fold tree as part of the init_apply method
void
JacobiLoopClosureMover::prepare_foldtree(core::pose::Pose & pose) {
	// check that there are no existing jumps inside the loop section of the pose, because that is incompatible with required cutpoint for this loop closure mover
	core::kinematics::FoldTree ft = pose.fold_tree();
	bool root_in_loop = false;
	core::Size current_cutpoint = loop_->stop()+1; // start with value outside loop.

	// confirm that residue at stop()+1 is the next residue in the polymer chain
	if ( pose.residue(loop_->stop()).connected_residue_at_upper() != loop_->stop()+1 ) {
		utility_exit_with_message("residue at loop.stop()+1 is not polymerically connected to residue at loop.stop(): correct input pose");
	}

	// scan loop for jump-point, root and cutpoint
	core::Size cur_res(loop_->start());
	while ( cur_res != loop_->stop() ) {
		// get next connected residue
		cur_res = pose.residue(cur_res).connected_residue_at_upper();
		// perform checks on that
		if ( ft.is_jump_point(cur_res) ) {
			utility_exit_with_message("jump point exists within the loop: fix the fold-tree or input loop.");
		} else if ( ft.is_root(cur_res) ) {
			root_in_loop = true;
		} else if ( ft.is_cutpoint(cur_res) ) {
			current_cutpoint = numeric::min(cur_res, current_cutpoint); // in case of multiple cut-points, store, the lowest
		}
	}
	if ( root_in_loop ) {
		if ( verbose_ ) {
			TR.Info << "root of fold-tree detected inside loop: moving root to start of loop, i.e. residue "
				<< loop_->start() << std::endl;
		}
		ft.reorder( loop_->start() );
		pose.fold_tree( ft );
	}

	// if no cut-point was found, then give an error message that it is pointless to close a loop without a cutpoint
	if ( current_cutpoint == loop_->stop()+1 ) {
		TR.Error << "You are closing a loop on an input pose that does not have a cut-point within the loop so there is nothing to close." << std::endl;
	}

	// if there is one existing cutpoint in the loop at loop_->stop(), but it does not have cutpoint variant types, exit with a message
	if ( ft.is_cutpoint(loop_->stop()) && ( !pose.residue(loop_->stop()).has_variant_type(core::chemical::CUTPOINT_LOWER) ||
			!pose.residue(pose.residue(loop_->stop()).connected_residue_at_upper()).has_variant_type(core::chemical::CUTPOINT_UPPER) ) ) {
		utility_exit_with_message("cut-point defined at loop.stop() but this residue has no variant-type CUTPOINT_LOWER "
			"and/or the subsequent residue has no variant-type CUTPOINT_UPPER: correct FoldTree.");
	}

	// if a cutpoint exists in the loop, but it is currently not at loop_->stop() then move it to the end of the loop
	if ( current_cutpoint < loop_->cut() ) {
		// define types of the residues that formed the original loop connection
		core::chemical::ResidueTypeCOP const res_type_pre_cut (pose.residue_type(current_cutpoint).get_base_type_cop());
		core::chemical::ResidueTypeCOP const res_type_post_cut (pose.residue_type(pose.residue(current_cutpoint).connected_residue_at_upper()).get_base_type_cop());
		// define atom IDs for first and second atoms of residue after the cut
		core::id::AtomID const fix_atom1(res_type_post_cut->mainchain_atom(1), pose.residue(current_cutpoint).connected_residue_at_upper());
		core::id::AtomID const fix_atom2(res_type_post_cut->mainchain_atom(2), pose.residue(current_cutpoint).connected_residue_at_upper());
		core::id::AtomID const fix_atom3(res_type_post_cut->mainchain_atom(3), pose.residue(current_cutpoint).connected_residue_at_upper());

		// assumption is that dihedral angles of the residues are better guess than the residue defaults. Therefore store for recovery
		core::Real const psi_cutpoint = pose.psi(current_cutpoint) * numeric::constants::d::degrees_to_radians;
		core::Real const omega_cutpoint = pose.omega(current_cutpoint) * numeric::constants::d::degrees_to_radians;
		core::Real const phi_cutpoint = pose.phi(pose.residue(current_cutpoint).connected_residue_at_upper()) * numeric::constants::d::degrees_to_radians;

		// call existing protocol for moving cutpoint
		protocols::loops::set_loop_cutpoint_in_pose_fold_tree(loop_->cut(), pose, loop_->start(), loop_->stop());

		// reset theta and d values to residue defaults and restore dihedrals, because this is not fixed automatically by above function
		// given the fact that residues are now occurring in ascending order (as opposed to before move of the cut, when upstream atoms where in descending order in AtomTree

		//pose.set_dof(core::id::DOF_ID(fix_atom1, core::id::PHI), - res_type_pre_cut -> upper_connect().icoor().phi() );
		pose.set_dof(core::id::DOF_ID(fix_atom1, core::id::PHI), psi_cutpoint );
		pose.set_dof(core::id::DOF_ID(fix_atom1, core::id::THETA), res_type_pre_cut->upper_connect().icoor().theta() );
		pose.set_dof(core::id::DOF_ID(fix_atom1, core::id::D), res_type_pre_cut->upper_connect().icoor().d() );

		//pose.set_dof(core::id::DOF_ID(fix_atom2, core::id::PHI), numeric::constants::d::pi );
		pose.set_dof(core::id::DOF_ID(fix_atom2, core::id::PHI), omega_cutpoint );
		pose.set_dof(core::id::DOF_ID(fix_atom2, core::id::THETA), res_type_post_cut->lower_connect().icoor().theta() );

		//pose.set_dof(core::id::DOF_ID(fix_atom3, core::id::PHI), - res_type_post_cut->lower_connect().icoor().phi() );
		pose.set_dof(core::id::DOF_ID(fix_atom3, core::id::PHI), phi_cutpoint );

		// report success of moving cutpoint
		if ( verbose_ ) {
			TR.Info << "Cut-point moved from " << current_cutpoint << " to " << loop_->stop() << " and the values at " <<
				current_cutpoint + 1 << " were reset to the residue type defaults." << std::endl;
		}
	}
	// check that all residues in the loop are in ascending order, which is requirement for Jacobian method
	runtime_assert_string_msg(pose.fold_tree().get_residue_edge(loop_->stop()).start() < loop_->start(), "failed to reorganize FoldTree so that residues in loop appear in ascending order");
}

/// @brief extract internal coordinates of the upper loop connection (between the final residue of the loop and the first residue that follows)
/// @details this function is called at initialization, and can later be manually overwritten using set_icoors
void
JacobiLoopClosureMover::default_target_icoors(core::pose::Pose const & pose){
	// STEP 1: define key atoms
	// store the atomID of the reference atom, which is the CA-atom of the first residue in the loop
	ref_atom_id_ = core::id::AtomID( pose.residue(loop_->start()).mainchain_atom(2), loop_->start() );

	// define types of the residues that form the loop connection at the end of the loop. Can use loop_->stop() + 1
	// because was already checked at initialization whether this corresponds to connected_residue_at_upper
	core::chemical::ResidueTypeCOP res_type_loop_end (pose.residue_type(loop_->stop()).get_base_type_cop());
	core::chemical::ResidueTypeCOP res_type_target (pose.residue_type(loop_->stop() + 1).get_base_type_cop());
	// define target atom numbers
	utility::fixedsizearray1< core::Size, 6 > target_bb_atoms;
	utility::vector1<core::id::AtomID > target_bb_atom_ids{6};
	utility::vector1< std::string > atom_names{6};

	for ( core::Size i=1; i <= 3; ++i ) {
		// in entries 1-3, store the backbone atoms of loop-end residue (le)
		core::Size target_le_atom_num = res_type_loop_end->mainchain_atom(i);
		target_bb_atoms[i] = target_le_atom_num;
		target_bb_atom_ids[i] = core::id::AtomID(target_le_atom_num, loop_->stop());
		atom_names[i] = pose.residue(loop_->stop()).atom_name(target_le_atom_num);
		// in entries 4-6, store the backbone atoms of the first backbone residue (bb))
		core::Size target_bb_atom_num = res_type_target->mainchain_atom(i);
		target_bb_atoms[3+i] = target_bb_atom_num;
		target_bb_atom_ids[3+i] = core::id::AtomID(target_bb_atom_num, loop_->stop() + 1);
		atom_names[3+i] = pose.residue(loop_->stop() + 1).atom_name(target_bb_atom_num);
	}

	// STEP 2: get internal coordinates.
	// define vector with phi angles, theta angles and distances for connections C=N', N'-CA', and CA'-C' (in this order),
	// where ['] indicates atoms in the residue in the backbone after the loop.
	// N.B. phi means the DIHEDRAL angle (i.e. phi, psi, or omega) associated to the stub in question.
	utility::fixedsizearray1<core::Real, 3> targets_phi;
	utility::fixedsizearray1<core::Real, 3> targets_theta;
	utility::fixedsizearray1<core::Real, 3> targets_d;

	//targets_phi[1]  = - res_type_loop_end -> upper_connect().icoor().phi(); // [rad] inverted because of Rosetta right-hand-rule convention
	targets_phi[1]  = pose.psi(loop_->stop()) * numeric::constants::d::degrees_to_radians; // [rad] psi angle of preceding residue is the phi angle of the stub of atom N
	targets_theta[1]= res_type_loop_end -> upper_connect().icoor().theta(); // [rad]
	targets_d[1]    = res_type_loop_end -> upper_connect().icoor().d(); // [Ang]

	//targets_phi[2]  = numeric::constants::d::pi; // [rad] set omega standard to 180 degrees;
	targets_phi[2]  = pose.omega(loop_->stop()) * numeric::constants::d::degrees_to_radians; // [rad] set omega standard to 180 degrees;
	targets_theta[2]= res_type_target -> lower_connect().icoor().theta(); // [rad]
	targets_d[2]    = pose.atom_tree().atom(target_bb_atom_ids[5]).distance(pose.atom_tree().atom(target_bb_atom_ids[4])); // [Ang]

	//targets_phi[3]  = - res_type_target->lower_connect().icoor().phi(); // [rad] inverted because of Rosetta right-hand-rule convention
	targets_phi[3]  = pose.phi(loop_->stop()+1) * numeric::constants::d::degrees_to_radians; // [rad] set omega standard to 180 degrees;
	targets_theta[3]= pose.dof(core::id::DOF_ID(target_bb_atom_ids[6], core::id::DOF_Type(2))); // [rad] NB in previous version was: targets_theta[3]= res_type_target -> icoor(target_bb_atoms[6]).theta(); // [rad]
	targets_d[3]    = pose.dof(core::id::DOF_ID(target_bb_atom_ids[6], core::id::DOF_Type(3))); // [Ang] NB in previous version was: targets_d[3]    = res_type_target -> icoor(target_bb_atoms[6]).d(); // [Ang]

	// If the connecting residue is not a jump point, then information can be stored from the current connection between the
	// last residue of the loop and the first residue of the subsequent backbone.
	if ( ! pose.fold_tree().is_cutpoint(loop_->stop()) and ! pose.fold_tree().is_jump_point(loop_->stop() + 1) ) {

		// if the loop end is currently not a cut-point, then in all expected cases the cut-point is 'upstream', meaning
		// that residues appear in descending order in the atom-tree. Because at the end of initialization, the cutpoint
		// is moved to the final residue of the loop, the residues will appear in inverted order in loop closure. This
		// also alters the atoms that define the phi, theta and d values, which is reflected in the expressions below
		if ( pose.fold_tree().get_residue_edge(loop_->stop()).start() > loop_->stop() ) {
			for ( core::Size i=1; i <= 3; ++i ) {
				targets_theta[i] = pose.atom_tree().atom(target_bb_atom_ids[i+1]).dof(core::id::THETA); // [rad]
				targets_d[i] = pose.atom_tree().atom(target_bb_atom_ids[i+2]).dof(core::id::D); // [Ang]
			}
		} else if ( pose.fold_tree().get_residue_edge(loop_->stop()).start() < loop_->stop() ) {
			// residues appear in ascending order
			for ( core::Size i=1; i <= 3; ++i ) {
				// starting with the first backbone atom of the residue after the loop
				targets_theta[i] = pose.atom_tree().atom(target_bb_atom_ids[3+i]).dof(core::id::THETA); // [rad]
				targets_d[i] = pose.atom_tree().atom(target_bb_atom_ids[3+i]).dof(core::id::D); // [Ang]
			}
		} else {
			// use default values. In this case information about the connection is lost
			TR.Warning << "Unexpected behavior. The values for the connecting phi, theta and d are set to default,"
				" while the target atom is not a jump." << std::endl;
		}
	}

	// step 4: store target internal coordinates as icoors, invert phi angles to be consistent with Icoor definitions.
	// construction vector with residue numbers to facilitate loop
	utility::vector1< core::chemical::AtomICoor > icoor_targets_temp{3};
	for ( core::Size i=1; i <= 3; ++i ) {
		icoor_targets_temp[i] = core::chemical::AtomICoor(atom_names[3 + i], -targets_phi[i], targets_theta[i],
			targets_d[i], atom_names[2 + i], atom_names[1 + i], atom_names[i],
			*res_type_target);
	}
	// store atom IDs and internal coordinates
	target_bb_atom_ids_ = target_bb_atom_ids;
	icoor_targets_ = icoor_targets_temp;
}

/// @brief Calculate the target position and Euler angles of the final C-atom in the loop w.r.t the reference frame connected to ref_atom_id_
/// @details This function calculates the pose (position and orientation) of the last backbone atom of the loop (C-atom [C] of the end residue [e] of the loop [l], 'Cle').
/// For that, the stored internal coordinates of the connection to the backbone are used to relate Cl to the first CA-atom of the backbone (CAbb) after the loop
/// This CAbb atom in term is fixed w.r.t. the CA atom of the start residue of the loop (CA1).
void
JacobiLoopClosureMover::store_target(core::pose::Pose const & pose) {
	// Initialize the homogeneous transformation using the differential stub from the reference atom to the final
	// target atom. In this way the loop target is not influenced by any motion outside the loop.

	// stub of the CA atom of the first residue in the loop
	core::kinematics::Stub const loop_start_stub = pose.atom_tree().atom(ref_atom_id_).get_stub();

	// stub of the C-atom of the target residue for closure (i.e. first residue after the loop). Stub is explicitly defined
	// because when obtained from the atom_tree() the stub can be defined in two ways, depending on the direction of the edge
	core::kinematics::Stub target_end_stub;
	target_end_stub.from_four_points(pose.atom_tree().atom(target_bb_atom_ids_[6]).xyz(),
		pose.atom_tree().atom(target_bb_atom_ids_[6]).xyz(),
		pose.atom_tree().atom(target_bb_atom_ids_[5]).xyz(),
		pose.atom_tree().atom(target_bb_atom_ids_[4]).xyz());
	// get the relative position and orientation from the start and end stubs
	core::kinematics::RT const RT_CA1_Cbb(loop_start_stub, target_end_stub);

	// store as an RT object
	target_ = RT_CA1_Cbb;
	target_initialized_ = true;
}

void
JacobiLoopClosureMover::current_closed_as_target( core::pose::Pose const & pose ){
	using core::chemical::AtomICoor;

	// confirm that there is no cutpoint in place across the loop end
	runtime_assert_msg(!pose.fold_tree().is_cutpoint(loop_->stop()), "Cannot initialize current as target when cut-point is already in place at loop end");

	// run default target
	default_target_icoors(pose);

	// extract current Euler angles and position
	store_target(pose);
}

/// @brief Calculate the current position and Euler angles of the CA-atom following the loop w.r.t the stub connected to ref_atom_id_
core::kinematics::RT
JacobiLoopClosureMover::update_current(core::pose::Pose const & pose) {
	//Get current position and Euler angles of loop end point using stub
	core::kinematics::RT const RT_CA1_Cle(pose.atom_tree().atom(ref_atom_id_).get_stub(),
		pose.atom_tree().atom(target_bb_atom_ids_[3]).get_stub());
	numeric::HomogeneousTransform<core::Real> const H_CA1_Cle(RT_CA1_Cle.get_rotation(), RT_CA1_Cle.get_translation());

	// Homogeneous matrix from final C-atom of the loop to the N-atom of the target (first residue post-loop)
	core::Real const psi_loop_end = pose.psi( loop_->stop() );
	numeric::xyzMatrix<core::Real> const R_Cle_Nbb = numeric::x_rotation_matrix_degrees(psi_loop_end) *
		numeric::z_rotation_matrix_radians(icoor_targets_[1].theta());
	numeric::HomogeneousTransform<core::Real> const H_Cle_Nbb(R_Cle_Nbb,
		icoor_targets_[1].d() * R_Cle_Nbb.col_x() );

	// Homogeneous matrix from N-atom to CA-atom of the target (first residue post-loop)
	numeric::xyzMatrix<core::Real> const R_Nbb_Cabb = numeric::x_rotation_matrix_radians(-icoor_targets_[2].phi()) *
		numeric::z_rotation_matrix_radians(icoor_targets_[2].theta());
	numeric::HomogeneousTransform<core::Real> const H_Nbb_Cabb(R_Nbb_Cabb,
		icoor_targets_[2].d() * R_Nbb_Cabb.col_x() );

	// Homogeneous matrix from CA-atom to C-atom of the target (first residue post-loop)
	numeric::xyzMatrix<core::Real> const R_Cabb_Cbb = numeric::x_rotation_matrix_radians(-icoor_targets_[3].phi()) *
		numeric::z_rotation_matrix_radians(icoor_targets_[3].theta());
	numeric::HomogeneousTransform<core::Real> const H_Cabb_Cbb(R_Cabb_Cbb,
		icoor_targets_[3].d() * R_Cabb_Cbb.col_x() );

	numeric::HomogeneousTransform<core::Real> const H_current(H_CA1_Cle * H_Cle_Nbb * H_Nbb_Cabb * H_Cabb_Cbb);

	// extract current Euler angles and position
	return core::kinematics::RT(H_current.rotation_matrix(), H_current.point());
}

} //jacobi
} //loop_closure
} //loops
} //protocols
