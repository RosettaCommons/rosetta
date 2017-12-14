// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is protocolsoped by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    MeanField.hh

/// @brief   Declarations and simple accessor/mutator definitions for MeanField
/// @author  Aliza Rubenstein (aliza.rubenstein@gmail.com)

#ifndef INCLUDED_protocols_mean_field_MeanField_HH
#define INCLUDED_protocols_mean_field_MeanField_HH

// Unit header
#include <protocols/mean_field/MeanField.fwd.hh>

// Package headers
#include <protocols/mean_field/RotMatrix.fwd.hh>

// Project headers
#include <core/pack/interaction_graph/FixedBBInteractionGraph.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/pack/task/PackerTask.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>
#include <protocols/mean_field/jagged_array.hh>
#include <protocols/mean_field/jagged_array.functions.hh>
#include <core/types.hh>

// Numeric headers


// C++ headers
#include <iostream>

namespace protocols {
namespace mean_field {

/// @details calculator which conducts mean-field algorithm  on either a single, non-designable pose or used as base class of other MeanFields
class MeanField : public utility::pointer::ReferenceCount {
public:
	// Standard methods ////////////////////////////////////////////////////////
	// @brief Standard constructor, no default constructor, initializes lambda_memory, tolerance, temp, and threshold to standard values
	MeanField( core::Size const option,
		core::pose::PoseOPs & poses,
		utility::vector1< core::pack::task::PackerTaskOP > const & tasks,
		core::scoring::ScoreFunctionOP scfxn );

	// @brief Optional constructor, initializes lambda_mem, tolerance, temp, and threshold to input values
	MeanField( core::Size const option,
		core::pose::PoseOPs & poses,
		utility::vector1< core::pack::task::PackerTaskOP > const & tasks,
		core::scoring::ScoreFunctionOP scfxn,
		core::Real lambda_mem,
		core::Real tolerance,
		core::Real temp,
		core::Real threshold );

	/// @brief Destructor
	virtual ~MeanField();


	// Standard Rosetta methods ////////////////////////////////////////////////

	/// @brief  Generate string representation of MeanField for debugging purposes.
	void show( std::ostream & output=std::cout ) const;

	/// @brief Insertion operator (overloaded so that MeanField can be "printed" in PyRosetta).
	friend std::ostream & operator<<( std::ostream & output, MeanField const & object_to_output );

	/// @brief tries to converge to a final RotMatrix
	void converge();

	/// @brief simply calls the converge method to find a converged RotMatrix
	virtual
	void process();

	/// @brief build a RotMatrix for a specific pose in poses_
	void build_rot_matrix_for_pose( core::Size pose_ind );

	/// @brief delete a pose
	virtual
	void delete_pose( core::Size pose_ind );

	// Accessors/Mutators

	/// @brief returns const pointer to RotMatrix
	virtual
	RotMatrixCOP rot_matrix() const;

	/// @brief returns const pointer to alternate RotMatrix
	inline
	RotMatrixCOP alt_rot_matrix() const
	{
		return alt_rot_matrix_;
	}

	/// @brief returns const reference to energies_matrix
	inline
	jagged_array < core::Real > const & energies_matrix() const
	{
		return energies_matrix_;
	}

	/// @brief get lambda memory, used in convergence process
	inline
	core::Real lambda_memory() const
	{
		return lambda_memory_;
	}

	/// @brief get tolerance, used to determine if converged
	inline
	core::Real tolerance() const
	{
		return tolerance_;
	}

	/// @brief get temperature (kT)
	inline
	core::Real temperature() const
	{
		return temperature_;
	}

	/// @brief get threshold
	inline
	core::Real threshold() const
	{
		return threshold_;
	}

	/// @brief get init option (used to initialize RotMatrix)
	inline
	core::Real init_option() const
	{
		return init_option_;
	}

	/// @brief get number of poses
	/// @details while this is trivial for MeanField (poses_.size() == 1), useful for derived classes with FlexBB
	inline
	core::Size num_poses() const
	{
		return poses_.size();
	}

	/// @brief get number of residues that are packable
	core::Size num_packed() const;

private:
	// Private methods /////////////////////////////////////////////////////////
	//no default constructor, uncopyable
	MeanField();
	MeanField( MeanField const & object_to_copy );
	MeanField & operator=( MeanField const & object_to_copy );

	///used by public methods of algorithm
	/// @brief calculate alt_rot_matrix_ by calculating an energies_matrix_ and then converting energies to Boltzmann probs
	void calc_alt_rot_matrix ( );

	/// @brief calculates energies_matrix_ based on the energies in pig_ and the weights in rot_matrix_
	void calc_energies_matrix ( );

	/// @brief converts the energies_matrix_ to the alt_rot_matrix_ by calculating the Boltzmann probabilities
	void convert_energies_to_alt_rm ( );

	/// @brief returns a bool, based on whether the maximum difference after rmsd-ing rot_matrix_ with alt_rot_matrix_ is < tolerance_
	bool has_converged () const;

	// Private data ////////////////////////////////////////////////////////////
	RotMatrixOP rot_matrix_;
	RotMatrixOP alt_rot_matrix_;
	jagged_array < core::Real > energies_matrix_;

	core::Real lambda_memory_;
	core::Real tolerance_; //used to determine if tolerance reached
	core::Real temperature_; //kT
	core::Real init_option_; //currently 0 or 1
	core::Real threshold_; //truncate energies that are higher than threshold

	core::pack::interaction_graph::FixedBBInteractionGraphOP pig_; //store a OP to the IG since only one IG can be used in each MeanField
	core::pose::PoseOPs & poses_;
	utility::vector1 < core::pack::task::PackerTaskOP > tasks_;
	core::scoring::ScoreFunctionOP scfxn_;

};  // class MeanField
}  // namespace mean_field
}  // namespace protocols

#endif  // INCLUDED_protocols_mean_field_MeanField_HH
