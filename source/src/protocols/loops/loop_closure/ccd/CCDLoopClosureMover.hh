// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    protocols/loops/loop_closure/ccd/CCDLoopClosureMover.hh
/// @brief   Class definition and method declarations for CCDLoopClosureMover
/// @author  Phil Bradley
/// @author  Oliver Lange
/// @author  Brian Weitzner
/// @author  Labonte <JWLabonte@jhu.edu>
/// @details Currently contains a classic ab initio implementation of CCD closure, according to Oliver....
/// @note    This file is the result of a refactor of code written by Phil and later wrapped in a Mover by Oliver.

#ifndef INCLUDED_protocols_loops_loop_closure_ccd_CCDLoopClosureMover_HH
#define INCLUDED_protocols_loops_loop_closure_ccd_CCDLoopClosureMover_HH

// Unit Header
#include <protocols/loops/loop_closure/ccd/CCDLoopClosureMover.fwd.hh>
#include <protocols/moves/Mover.hh>

// Package Headers
// TODO: Make loop_ an OP so we can use the .fwd.hh here. ~Labonte
#include <protocols/loops/Loop.hh>
#include <protocols/loops/loop_closure/ccd/RamaCheck.fwd.hh>

// Project Headers
#include <core/types.hh>
#include <core/id/TorsionID.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/kinematics/MoveMap.fwd.hh>

// C++ header
#include <string>

namespace protocols {
namespace loops {
namespace loop_closure {
namespace ccd {

enum ChainDirection {
	forward = 1,
	backward,
	n_chain_directions = backward
};

// Surely this exists elsewhere in Rosetta? ~Labonte
// TODO: Look at Andrew's branch.
enum SecondaryStructureType {
	helix = 1,
	strand,
	coil,
	n_secondary_structure_types = coil
};

/// @brief Close a loop (protein or otherwise) using bi-directional Cyclic Coordinate Descent
///
/// @details This class was initially a wrapper for fast_ccd_loop_closure(), an exposed method in ccd_closure.cc.
/// Before using this mover:
///  1) Set a proper foldtree - usually a regular loop foldtree is used. (see protocols::loops::foldtree_from_loops)
///  2) Add cutpoint variants to the loop cutpoint
///
/// The original description of the CCD algorithm in relation to protein loop closure can be found here:
///  Canutescu AA, Dunbrack RL Jr: Cyclic coordinate descent: a robotics algorithm for protein loop closure. 
///  Protein Sci 2003, 12:963-972.
///
class CCDLoopClosureMover : public moves::Mover {

public:  // Standard methods //////////////////////////////////////////////////
	/// @brief  Empty constructor
	CCDLoopClosureMover();

	/// @brief  Copy constructor
	CCDLoopClosureMover( CCDLoopClosureMover const & object_to_copy );

	/// @brief  Constructor with Loop input option
	CCDLoopClosureMover( protocols::loops::Loop const & loop );

	/// @brief  Constructor with Loop and MoveMap input options
	CCDLoopClosureMover( protocols::loops::Loop const & loop, core::kinematics::MoveMapCOP mm );

	// Destructor
	virtual ~CCDLoopClosureMover();

	// Assignment operator
	CCDLoopClosureMover & operator=( CCDLoopClosureMover const & object_to_copy );

public:  // Standard Rosetta methods //////////////////////////////////////////
	// General methods
	/// @brief  Register options with the option system.
	static void register_options();

	/// @brief  Generate a string representation of CCDLoopClosureMover for debugging purposes.
	virtual void show( std::ostream & output=std::cout ) const;

	// Mover methods
	/// @brief  Return the name of the Mover.
	virtual std::string get_name() const;

	virtual protocols::moves::MoverOP clone() const;

	virtual protocols::moves::MoverOP fresh_instance() const;

	/// @brief  Apply the corresponding move to <pose>.
	virtual void apply( core::pose::Pose & pose );

	/// @brief Called by MoverFactory when constructing new Movers
	virtual
	void parse_my_tag( TagCOP tag, basic::datacache::DataMap &, Filters_map const &,
		moves::Movers_map const &,
		Pose const & );


public:  // Accessors/Mutators ////////////////////////////////////////////////
	/// @brief  Get the Loop to be closed.
	protocols::loops::Loop loop() const;

	/// @brief  Set the Loop to be closed.
	void loop( protocols::loops::Loop new_loop );

	/// @brief  Get the current MoveMap.
	core::kinematics::MoveMap movemap() const;

	/// @brief  Set the MoveMap.
	///  Positions outside of the defined loop will be ignored.
	void movemap( core::kinematics::MoveMapCOP new_movemap );


	/// @brief   Set the maximum change in torsion angle a residue is allowed per closure move for each of the three
	/// secondary structure types.
	/// @param   Order of parameters is helix, strand, coil. Angles are in degrees.
	void
	max_per_move_torsion_delta_per_residue(
			core::Angle input_max_delta_helix,
			core::Angle input_max_delta_strand,
			core::Angle input_max_delta_coil )
	{
		max_per_move_torsion_delta_[ helix ] = input_max_delta_helix;
		max_per_move_torsion_delta_[ strand ] = input_max_delta_strand;
		max_per_move_torsion_delta_[ coil ] = input_max_delta_coil;
	}

	/// @brief   Get the maximum change in torsion angle a residue is allowed for a single move for one of the
	/// three secondary structure types.
	/// @return  angle in degrees
	core::Angle
	max_per_move_torsion_delta_per_residue( SecondaryStructureType secstruct ) const
	{
		return max_per_move_torsion_delta_[ secstruct ];
	}

	/// @brief   Get the maximum allowed torsion angle deviation for the entire closure run for one of the three
	/// secondary structure types ( H, E, L).
	core::Angle max_per_move_torsion_delta_per_residue( char secstruct ) const;

	/// @brief   Get the maximum change in torsion angle a residue is allowed for the entire closure run for one of the
	/// three secondary structure types.
	/// @return  angle in degrees
	core::Angle
	max_total_torsion_delta_per_residue( SecondaryStructureType secstruct ) const
	{
		return max_total_torsion_delta_[ secstruct ];
	}

	/// @brief   Get the maximum allowed torsion angle deviation for the entire closure run for one of the three
	/// secondary structure types ( H, E, L).
	core::Angle max_total_torsion_delta_per_residue( char secstruct ) const;

	/// @brief   Set the maximum change in torsion angle a residue is allowed for the entire closure run for each of
	/// the three secondary structure types.
	/// @param Order of parameters is helix, strand, coil.  Angles are in degrees.
	void
	max_total_torsion_delta_per_residue(
			core::Angle input_max_delta_helix,
			core::Angle input_max_delta_strand,
			core::Angle input_max_delta_coil )
	{
		max_total_torsion_delta_[ helix ] = input_max_delta_helix;
		max_total_torsion_delta_[ strand ] = input_max_delta_strand;
		max_total_torsion_delta_[ coil ] = input_max_delta_coil;
	}

	/// @brief   Get the tolerance for loop closure in Angstroms.
	/// @details A forward and backward splice of RMS over N, CA, and C (or target atoms) must be less than the tolerance for an early
	/// return; otherwise, the algorithm will go through the loop the requested number of cycles.
	/// RMS Deviation is calculation as sqrt( dev / number_of_atoms )
	core::Real tolerance() const { return tolerance_; }

	/// @brief   Set the tolerance for loop closure in Angstroms.
	/// @details A forward and backward splice of RMS deviation over N, CA, and C (or target atoms)must be less than the tolerance for
	/// an early return; otherwise, the algorithm will go through the loop the requested number of cycles.
	/// RMS Deviation is calculation as sqrt( dev / number_of_atoms )
	void tolerance( core::Real input_tolerance ) { tolerance_ =  input_tolerance; }

	/// @brief   Get the maximum number of cycles to attempt if the tolerance is not satisfied.
	/// @details Each cycle includes both forward and backward directions.
	core::Size max_cycles() const { return max_cycles_; }

	/// @brief   Set the maximum number of cycles to attempt if the tolerance is not satisfied.
	/// @details Each cycle includes both forward and backward directions.
	void max_cycles( core::Size input_ccd_cycles ) { max_cycles_ = input_ccd_cycles; }

	// FIXME: This function name is misleading! Nothing is being checked when this is called.
	/// @brief  Are closure moves checked using the rama scores and the Metropolis criterion?
	bool check_rama_scores() const { return check_rama_scores_; }

	/// @brief  Set whether or not closure moves are checked using the rama scores and the Metropolis criterion.
	void check_rama_scores( bool setting ) { check_rama_scores_ = setting; }

	/// @brief  Are two-body (neighbor-dependent) Ramachandran maps being used?
	bool use_rama_2B() const { return use_rama_2b_; }

	/// @brief  Set whether or not two-body (neighbor-dependent) Ramachandran maps should be used.
	void use_rama_2B( bool setting ) { use_rama_2b_ = setting; }

	/// @brief  Get the RMS deviation of the target atoms after completion of loop closure.
	/// @details sqrt( dev / number_of_atoms )
	core::Real deviation() const { return deviation_; }

	/// @brief  Get the average change in the main-chain torsion angles of the loop after completion of loop closure.
	/// @note   Change this name to match the datum. ~Labonte
	core::Angle	torsion_delta() const	{ return average_change_in_torsion_angle_; }

	/// @brief  Get the average change in rama score for the residues in the loop after completion of loop closure.
	/// @note   This value will only be meaningful if the option for checking the rama score is turned on.
	/// @note   Change this name to match the datum. ~Labonte
	core::Real rama_delta() const { return average_change_in_rama_score_; }

	/// @brief  Get the number of cycles used to close the loop.
	core::Size actual_cycles() const { return actual_cycles_; }

	/// @brief  Get a pointer to the RamaCheck instance being used.
	RamaCheckBaseOP rama() const;

public:  // Other Public Methods //////////////////////////////////////////////
	/// @brief Return true if the forward and backward RMS deviations are each lower than the tolerance value.
	bool success() const;


private:  // Private methods //////////////////////////////////////////////////
	// Initialize data members from arguments.
	void init( protocols::loops::Loop const & loop, core::kinematics::MoveMapCOP movemap );

	// Initialize data members from option system.
	void init_options();

	// Copy all data members from <from> to <to>.
	void copy_data( CCDLoopClosureMover & to, CCDLoopClosureMover const & from ) const;

	// Return a the coordinates of the atoms to be overlapped for a given residue.
	utility::vector1< core::PointPosition > get_anchors( core::conformation::Residue const & residue ) const;

	// Adjust the residue number and atom number if necessary when determining the connectivity across a torsion 
	void index_pair_in_range( core::uint & pos, core::uint & atom, core::Size const n_mainchain_atoms ) const;

	void get_torsion_axis(
			core::pose::Pose const & pose,
			core::uint const seqpos,
			core::uint const torsion_num,
			core::Vector & axis_atom_coords,
			core::Vector & axis_unit_vector ) const;

	core::Angle calculate_ccd_angle(
			core::pose::Pose const & pose,
			core::uint const pos,
			core::uint const torsion,
			ChainDirection const direction );

	// This method can be overridden to change the way maximum deviations are determined
	virtual void get_maximum_torsion_deltas_for_residue(
			core::pose::Pose const & pose,
			core::uint const seqpos,
			core::Real & per_move_allowed_delta,
			core::Real & total_allowed_delta ) const;

	// This method should be reimplemented in derived classes that may want to use alternate formulations of CCD
	virtual void adjust_residue_to_minimize_deviation(
			core::pose::Pose & pose,
			core::pose::Pose const & starting_pose,
			core::uint const seqpos,
			ChainDirection const direction,
			core::Real const max_per_move_torsion_delta,
			core::Real const max_total_torsion_delta );

	// This method performs most of the work of the CCD closure.
	void close_loop_in_single_direction(
			core::pose::Pose & pose,
			core::pose::Pose const & starting_pose,
			ChainDirection const direction );

	void compute_closure_metrics( core::pose::Pose const & pose, core::pose::Pose const & starting_pose );


private:  // Private data /////////////////////////////////////////////////////
	protocols::loops::Loop loop_;
	core::kinematics::MoveMapCOP movemap_;

	utility::vector1< core::Angle > max_per_move_torsion_delta_;  // indexed by SecondaryStructureType
	utility::vector1< core::Angle > max_total_torsion_delta_;  // indexed by SecondaryStructureType
	core::Distance tolerance_;  // in Angstroms

	core::Size max_cycles_;
	bool check_rama_scores_;
	bool use_rama_2b_;

	core::Real deviation_;  // RMSD in Angstroms
	core::Angle average_change_in_torsion_angle_;

	mutable RamaCheckBaseOP rama_;  // Ramachandran object for checking rama scores
	utility::vector1< core::Real > starting_rama_scores_;
	core::Real average_change_in_rama_score_;

	Size actual_cycles_;


private:  // Constants ////////////////////////////////////////////////////////
	static utility::vector1< core::SSize > const STEP_SIZE;  // indexed by ChainDirection
	static core::Distance const MALARKEY;
	static core::Real const BAD_SCORE;

private: // Private Static Methods ////////////////////////////////////////////
	static std::map< char, SecondaryStructureType > const & sec_struc_char_to_enum_map();

};  // class CCDLoopClosureMover

// Insertion operator (overloaded so that CCDLoopClosureMover can be "printed" in PyRosetta).
std::ostream & operator<< ( std::ostream & os, CCDLoopClosureMover const & mover );

} // namespace ccd
} // namespace loop_closure
} // namespace loops
} // namespace protocols

#endif  // INCLUDED_protocols_loops_loop_closure_ccd_CCDLoopClosureMover_HH
