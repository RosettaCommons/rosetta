// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is protocolsoped by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    DesignMeanField.hh

/// @brief   Declarations and simple accessor/mutator definitions for DesignMeanField.
/// @author  arubenstein

#ifndef INCLUDED_protocols_mean_field_DesignMeanField_HH
#define INCLUDED_protocols_mean_field_DesignMeanField_HH

// Unit header
#include <protocols/mean_field/DesignMeanField.fwd.hh>

// Package headers
#include <protocols/mean_field/AAMatrix.fwd.hh>
#include <protocols/mean_field/MeanField.hh>


// Project headers

// Utility headers

// Numeric headers


// C++ headers
#include <iostream>

namespace protocols {
namespace mean_field {
/// @details calculator which conducts mean-field algorithm  on one, designable pose to create a RotMatrix and AAMatrix
class DesignMeanField : public protocols::mean_field::MeanField {
public:
	// Standard methods ////////////////////////////////////////////////////////
	// @brief Standard constructor, no default constructor, initializes lambda_memory, tolerance, temp, and threshold to standard values
	DesignMeanField( core::Size option,
		core::pose::PoseOPs & poses,
		utility::vector1 < core::pack::task::PackerTaskOP > tasks,
		core::scoring::ScoreFunctionOP scfxn
	);

	// @brief Optional constructor, initializes lambda_mem, tolerance, temp, and threshold to input values
	DesignMeanField( core::Size option,
		core::pose::PoseOPs & poses,
		utility::vector1 < core::pack::task::PackerTaskOP > tasks,
		core::scoring::ScoreFunctionOP scfxn,
		core::Real lambda_mem,
		core::Real tolerance,
		core::Real temp,
		core::Real threshold
	);

	/// @brief Destructor
	~DesignMeanField();


	// Standard Rosetta methods ////////////////////////////////////////////////

	/// @brief  Generate string representation of DesignMeanField for debugging purposes.
	void show( std::ostream & output=std::cout ) const;

	/// @brief Insertion operator (overloaded so that FlexBBMeanField can be "printed" in PyRosetta).
	friend std::ostream & operator<<( std::ostream & output, DesignMeanField const & object_to_output );

	/// @brief overrides process() method in MeanField to initialize AAMatrix from RotMatrix
	virtual
	void process();

	// Accessors/Mutators
	/// @brief returns const reference to AAMatrix
	AAMatrixCOP aa_matrix() const;

private:
	// Private methods /////////////////////////////////////////////////////////

	/// no default constructor, uncopyable
	DesignMeanField();
	DesignMeanField( DesignMeanField const & object_to_copy );
	DesignMeanField & operator=( DesignMeanField const & object_to_copy );


	/// @brief initialize the aa_matrix.
	/// @remarks This can theoretically be called at any time but should only be called after the RotMatrix has converged
	void init_aa_matrix();

	// Private data ////////////////////////////////////////////////////////////
	AAMatrixOP aa_matrix_;

};  // class DesignMeanField

}  // namespace mean_field
}  // namespace protocols

#endif  // INCLUDED_protocols_mean_field_DesignMeanField_HH
