// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is protocolsoped by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    FlexBBMeanField.hh

/// @brief   Declarations and simple accessor/mutator definitions for FlexBBMeanField.
/// @author  arubenstein

#ifndef INCLUDED_protocols_mean_field_FlexBBMeanField_HH
#define INCLUDED_protocols_mean_field_FlexBBMeanField_HH

// Unit header
#include <protocols/mean_field/FlexBBMeanField.fwd.hh>

// Package headers
#include <protocols/mean_field/MeanField.hh>

// Project headers


// Utility headers

// Numeric headers


// C++ headers
#include <iostream>

namespace protocols {
namespace mean_field {
/// @details calculator which conducts mean-field algorithm  on several, non-designable poses to create a single, averaged RotMatrix
class FlexBBMeanField : public protocols::mean_field::MeanField {
public:
	// Standard methods ////////////////////////////////////////////////////////
	// @brief Standard constructor, no default constructor, initializes lambda_memory, tolerance, temp, and threshold to standard values
	FlexBBMeanField( core::Size const option, core::pose::PoseOPs & poses, utility::vector1 < core::pack::task::PackerTaskOP > tasks, core::scoring::ScoreFunctionOP scfxn );

	// @brief Optional constructor, initializes lambda_mem, tolerance, temp, and threshold to input values
	FlexBBMeanField( core::Size const option, core::pose::PoseOPs & poses, utility::vector1 < core::pack::task::PackerTaskOP > tasks, core::scoring::ScoreFunctionOP scfxn,
		core::Real lambda_mem, core::Real tolerance, core::Real temp, core::Real threshold );

	/// @brief Destructor
	~FlexBBMeanField();


	// Standard Rosetta methods ////////////////////////////////////////////////

	/// @brief  Generate string representation of FlexBBMeanField for debugging purposes.
	void show( std::ostream & output=std::cout ) const;

	/// @brief Insertion operator (overloaded so that FlexBBMeanField can be "printed" in PyRosetta).
	friend std::ostream & operator<<( std::ostream & output, FlexBBMeanField const & object_to_output );

	/// @brief overrides process() method in MeanField to add convergence step for each pose and averaging step for all rot_matrices_
	virtual
	void process();

	/// @brief deletes a pose from the list of poses
	virtual
	void delete_pose( core::Size pose_ind );

	/// @brief returns the expected (averaged) RotMatrix (exp_rot_matrix_)
	virtual
	RotMatrixCOP rot_matrix() const;

	/// @brief returns the expected (averaged) EnergyMatrix (exp_energy_matrix_)
	//added 6/23/15
	inline
	jagged_array < core::Real > energy_matrix() const { return exp_energy_matrix_; }

	/// @brief renumbers rotamers with a backbone-independent rotamer numbering scheme
	void renumber_rotamers();

	/// @brief calculate Boltzmann weights of each backbone
	virtual
	void calc_bb_boltz_probs();

	/// @brief calculate expected (averaged) RotMatrix
	void calc_exp_value_rot_matrix();

	// Accessors/Mutators

	/// @brief used by derived classes for write-access to bb_boltz_probs_
	/// @remarks would prefer for the method to be protected, but disallowed by Rosetta Coding Conventions
	inline
	jagged_array < core::Real > & bb_boltz_probs() { return bb_boltz_probs_; }

	inline
	utility::vector1 < jagged_array < core::Real > > const & bb_boltz_probs_per_aa() { return bb_boltz_probs_per_aa_; }

	/// @brief returns const references to the utility vector of rot_matrices representing each pose in the backbone ensemble
	inline
	utility::vector1 < RotMatrix > const & rot_matrices() const { return rot_matrices_; }

	/// @brief returns const references to the utility vector of energy_matrices representing each pose in the backbone ensemble
	inline
	utility::vector1 < jagged_array < core::Real > > const & energy_matrices() const { return energy_matrices_; }

protected:

	//jagged_array of Boltzmann probabilities of each backbone position. bb_boltz_probs_[1][2] is prob of backbone 2 at pos 1
	jagged_array < core::Real > bb_boltz_probs_;

	//vector of Boltzmann probabilities of each aa in each backbone at each position.  bb_boltz_probs_per_aa_[1][1][2] is prob of aa 2
	//at position 1 in backbone 1.
	utility::vector1 < jagged_array < core::Real > > bb_boltz_probs_per_aa_;


private:
	// Private methods /////////////////////////////////////////////////////////

	/// no default constructor, uncopyable
	FlexBBMeanField();
	FlexBBMeanField( FlexBBMeanField const & object_to_copy );
	FlexBBMeanField & operator=( FlexBBMeanField const & object_to_copy );

	// Private data ////////////////////////////////////////////////////////////
	utility::vector1 < RotMatrix > rot_matrices_;
	RotMatrixOP exp_rot_matrix_;
	utility::vector1 < jagged_array < core::Real > > energy_matrices_;
	//added 6/23/15
	jagged_array < core::Real > exp_energy_matrix_;

	//jagged_array of Boltzmann weights (not yet divided by partition function) of each backbone position.
	//bb_boltz_weight[1][2] is weight of backbone 2 at pos 1
	jagged_array < core::Real > bb_boltz_weights_;

	utility::vector1 < core::Size > nrot_per_pos_; //total numbers of rot at each position, across backbones

};  // class FlexBBMeanField

}  // namespace mean_field
}  // namespace protocols

#endif  // INCLUDED_protocols_mean_field_FlexBBMeanField_HH
