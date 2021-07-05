// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/loops/loop_closure/jacobi/JacobiLoopClosureMover.hh
/// @brief Class definition and method declarations for JacobiLoopClosureMover
/// @author Teun Hoevenaars (teunhoevenaars@gmail.com)

#ifndef INCLUDED_protocols_loops_loop_closure_jacobi_JacobiLoopClosureMover_HH
#define INCLUDED_protocols_loops_loop_closure_jacobi_JacobiLoopClosureMover_HH

// Unit headers
#include <protocols/loops/loop_closure/jacobi/JacobiLoopClosureMover.fwd.hh>
#include <protocols/moves/Mover.hh>

// Protocol headers
#include <protocols/loops/Loop.fwd.hh>

// Core headers
#include <core/chemical/AtomICoor.hh>
#include <core/id/AtomID.hh>
#include <core/kinematics/MoveMap.fwd.hh>
#include <core/kinematics/RT.hh>
#include <core/kinematics/jacobian/SeriesJacobians.fwd.hh>
#include <numeric/constants.hh>

// Basic/Utility headers
//#include <basic/datacache/DataMap.fwd.hh>
//#include <utility/tag/XMLSchemaGeneration.fwd.hh> //transcluded from Mover

namespace protocols {
namespace loops {
namespace loop_closure {
namespace jacobi {

///@brief Iterative loop closure mover, based on a linearization of the kinematics, which adjusts all free backbone torsion angles
///simultaneously in each iteration.
class JacobiLoopClosureMover : public protocols::moves::Mover {

public:

	/////////////////////
	/// Constructors  ///
	/////////////////////

	/// @brief Default constructor
	JacobiLoopClosureMover();

	/// @brief Copy constructor (not needed unless you need deep copies)
	JacobiLoopClosureMover( JacobiLoopClosureMover const & src );

	/// @brief Constructor based on a loop
	JacobiLoopClosureMover( protocols::loops::Loop const & loop );

	/// @brief Constructor based on a loop and a MoveMap
	JacobiLoopClosureMover( protocols::loops::Loop const & loop, core::kinematics::MoveMapCOP const & mm );

	/// @brief Destructor (important for properly forward-declaring smart-pointer members)
	~JacobiLoopClosureMover() override;

public:

	/////////////////////
	/// Mover Methods ///
	/////////////////////

	/// @brief Set the loop to be closed
	void
	set_loop( protocols::loops::Loop const & new_loop );

	/// @brief Get the loop to be closed
	protocols::loops::LoopCOP
	get_loop() { return loop_; };

	/// @brief Set the MoveMap
	void
	set_movemap( core::kinematics::MoveMapOP const & new_mm );

	/// @brief Get the MoveMap
	core::kinematics::MoveMapCOP
	get_movemap() { return movemap_; };

	/// @brief Get pointer to SeriesJacobian object
	core::kinematics::jacobian::SeriesJacobiansOP
	get_jacobian_chain() {return jacobian_chain_; };

	/// @brief Initializes private parameters of the Jacobi mover that depend on the pose
	void
	init_apply( core::pose::Pose & pose);

	/// @brief Execute the loop closure
	void
	apply( core::pose::Pose & pose ) override;

	void
	show( std::ostream & output = std::cout ) const override;

	/// @brief Get boolean indicating whether last closure was successfully completed
	bool
	last_closure_success(){ return closure_success_; };

	/// @brief Returns the number of cycles that were needed for the last closure
	core::Size
	last_closure_cycles() { return closure_cycles_; };

	/// @brief Set maximum number of iterations to find a loop closure solution
	void
	set_max_cycles(core::Size max_cycles){ max_cycles_ = max_cycles; }

	/// @brief Get value for maximum number of iterations to find a loop closure solution
	core::Size
	get_max_cycles() { return max_cycles_; }

	/// @brief Set norm-values of allowed rotational [deg] and linear error vectors [Ang] that define successful closure
	/// @details error_norm_rot should be supplied in degrees.
	void
	set_error_norms(core::Real error_norm_rot, core::Real error_norm_lin);

	/// @brief Get value for norm of allowed rotational error vector [deg] that defines successful closure
	core::Real
	get_allowed_norm_rot() { return err_rot_allowed_ * numeric::constants::d::radians_to_degrees; }

	/// @brief Get value for norm of allowed linear error vector [Ang] that defines successful closure
	core::Real
	get_allowed_norm_lin() { return err_lin_allowed_; }

	/// @brief Set whether or not to output non-critical information.
	/// @details Default is false
	void set_verbose( bool setting ) { verbose_ = setting; }

	/// @brief Set current loop connection with subsequent residue as target connection.
	void
	current_closed_as_target( core::pose::Pose const & pose );

public:

	///////////////////////////////
	/// Rosetta Scripts Support ///
	///////////////////////////////

	/// @brief parse XML tag (to use this Mover in Rosetta Scripts)
	void
	parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap &
	) override;

	/// @brief required in the context of the parser/scripting scheme
	protocols::moves::MoverOP
	fresh_instance() const override;

	/// @brief required in the context of the parser/scripting scheme
	protocols::moves::MoverOP
	clone() const override;

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

private: // METHODS

	/// @brief Create default MoveMap
	core::kinematics::MoveMapOP
	create_default_movemap( protocols::loops::Loop loop);

	/// @brief Initializes data members from constructor input arguments
	void
	init_constructor(protocols::loops::Loop const & loop, core::kinematics::MoveMapCOP const & mm);

	/// @brief Function that checks and adjusts the fold tree as part of the init_apply method
	void
	prepare_foldtree(core::pose::Pose & pose);

	/// @brief extract internal coordinates of the upper loop connection (between the final residue of the loop and the first residue that follows)
	void
	default_target_icoors( core::pose::Pose const & pose );

	/// @brief Set the target based on residue defaults
	void
	store_target( core::pose::Pose const & pose);

	/// @brief Calculate the current position and Euler angles of the final C-atom in the loop w.r.t the reference frame connected to ref_atom_id_
	core::kinematics::RT
	update_current(core::pose::Pose const & pose);

private: // USER-MANAGED PARAMETERS
	/// @brief number of iterations before closure is aborted
	core::Size max_cycles_ = 200; // Simples closure take ~10 cycles.

	/// @brief Maximum norms of rotational [rad] and linear error [Ang] for closure to be accepted as successful
	core::Real err_rot_allowed_ = 5 * numeric::constants::d::degrees_to_radians; // [rad]
	core::Real err_lin_allowed_ = 0.1; // [Ang]

	/// @brief setting that determines whether non-critical information is printed to tracer. Errors/warnings are always printed
	bool verbose_ = false;

private: // INTERNALLY-MANAGED OBJECTS AND VARIABLES
	/// @brief copy of the input loop
	protocols::loops::LoopOP loop_{nullptr};

	/// @brief copy of the pointer to the movemap
	core::kinematics::MoveMapCOP movemap_{nullptr};

	/// @brief vector1 with the numbers of the residues whose phi and psi angles are allowed to be adjusted by the mover
	utility::vector1<core::Size> free_residues_;

	/// @brief pointer to vector with Jacobian objects
	core::kinematics::jacobian::SeriesJacobiansOP jacobian_chain_{nullptr};

	/// @brief counter of number of cycles in the loop closure while-loop
	core::Size closure_cycles_{0};

	/// @brief BOOl to communicate with user about success/failure of closure
	bool closure_success_{false};

	/// @brief atom id of the first CA atom in the loop, which is fixed to the preceding backbone residue and whose stub is the
	/// reference frame for all vectors and matrices
	/// @details this member variable is set when the apply() function is called because it relies on the conformation
	core::id::AtomID ref_atom_id_{0,0};

	/// @brief atom ids of the backbone atoms of the last residue of the loop and the first residue following
	/// @details this member variable is set when the apply() function is called because it relies on the conformation
	utility::vector1<core::id::AtomID> target_bb_atom_ids_{6};

	/// @ brief target vectors for the loop closure. Are either initialized based on residuetype defaults, or can be set
	/// via current_closed_as_target, before calling apply()
	core::kinematics::RT target_;
	bool target_initialized_ = false;

	/// @brief internal coordinates of the target residue, which determine how the loop must be connected to the backbone.
	/// Stored in this form, and not as homogeneous matrix, because maybe a user in the future would like to change individual values manually
	/// @details this member variable is set when the apply() function is called because it relies on the conformation
	utility::vector1< core::chemical::AtomICoor > icoor_targets_;
};

std::ostream &
operator<<( std::ostream & os, JacobiLoopClosureMover const & mover );

} //jacobi
} //loop_closure
} //loops
} //protocols

#endif //protocols_loops_loop_closure_jacobi_JacobiLoopClosureMover_HH
