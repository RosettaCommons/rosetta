// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is protocolsoped by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    FlexBBDesignMeanField.hh

/// @brief   Declarations and simple accessor/mutator definitions for FlexBBDesignMeanField.
/// @author  arubenstein

#ifndef INCLUDED_protocols_mean_field_FlexBBDesignMeanField_HH
#define INCLUDED_protocols_mean_field_FlexBBDesignMeanField_HH

// Unit header
#include <protocols/mean_field/FlexBBDesignMeanField.fwd.hh>

// Package headers
#include <protocols/mean_field/FlexBBMeanField.hh>
#include <protocols/mean_field/AAMatrix.hh>

// Project headers


// Utility headers
#include <utility/pointer/ReferenceCount.hh>

// Numeric headers


// C++ headers
#include <iostream>

namespace protocols {
namespace mean_field {
/// @details
/// @remarks
class FlexBBDesignMeanField : public protocols::mean_field::FlexBBMeanField {
public:
	// Standard methods ////////////////////////////////////////////////////////
	// @brief Standard constructor, no default constructor, initializes lambda_memory, tolerance, temp, and threshold to standard values
	FlexBBDesignMeanField( core::Size const option,
		core::pose::PoseOPs & poses,
		utility::vector1 < core::pack::task::PackerTaskOP > tasks,
		core::scoring::ScoreFunctionOP scfxn
	);

	// @brief Optional constructor, initializes lambda_mem, tolerance, temp, and threshold to input values
	FlexBBDesignMeanField( core::Size const option,
		core::pose::PoseOPs & poses,
		utility::vector1 < core::pack::task::PackerTaskOP > tasks,
		core::scoring::ScoreFunctionOP scfxn,
		core::Real lambda_mem,
		core::Real tolerance,
		core::Real temp,
		core::Real threshold
	);


	/// @brief Destructor
	~FlexBBDesignMeanField();


	// Standard Rosetta methods ////////////////////////////////////////////////

	/// @brief  Generate string representation of FlexBBDesignMeanField for debugging purposes.
	void show( std::ostream & output=std::cout ) const;

	/// @brief Insertion operator (overloaded so that FlexBBDesignMeanField can be "printed" in PyRosetta).
	friend std::ostream & operator<<( std::ostream & output, FlexBBDesignMeanField const & object_to_output );

	/// @brief overrides process() method in FlexBBMeanField to init AAMatrix for each rot_matrix and average them together
	virtual
	void process();

	/// @brief deletes a pose from the list of poses
	virtual
	void delete_pose( core::Size pose_ind );

	/// @brief overrides FlexBBMeanField::calc_bb_boltz_probs to filter out bb_boltz_probs that aren't designable
	virtual
	void calc_bb_boltz_probs();

	/// @brief returns averaged AAMatrix
	AAMatrixCOP aa_matrix() const;

	/// @brief calculates expected value (averaged) AAMatrix
	void calc_exp_value_aa_matrix();

	// Accessors/Mutators

	/// @brief returns a const vector of the AAMatrices constructed from each pose
	utility::vector1 < AAMatrix > const & aa_matrices() const
	{
		return aa_matrices_;
	}

private:
	// Private methods /////////////////////////////////////////////////////////
	/// no default constructor, uncopyable
	FlexBBDesignMeanField();
	FlexBBDesignMeanField( FlexBBDesignMeanField const & object_to_copy );
	FlexBBDesignMeanField & operator=( FlexBBDesignMeanField const & object_to_copy );

	// Private data ////////////////////////////////////////////////////////////

	utility::vector1 < AAMatrix > aa_matrices_;
	AAMatrixOP exp_aa_matrix_;

};  // class FlexBBDesignMeanField

}  // namespace mean_field
}  // namespace protocols

#endif  // INCLUDED_protocols_mean_field_FlexBBDesignMeanField_HH
