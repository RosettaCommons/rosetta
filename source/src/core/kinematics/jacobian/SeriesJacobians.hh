// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/kinematics/jacobian/SeriesJacobians.hh
/// @brief class that defines series of Jacobian modules
/// @author teunhoevenaars (teunhoevenaars@gmail.com)


#ifndef INCLUDED_core_kinematics_jacobian_SeriesJacobians_hh
#define INCLUDED_core_kinematics_jacobian_SeriesJacobians_hh

#include <core/kinematics/jacobian/SeriesJacobians.fwd.hh>

// Utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/VirtualBase.hh>
#include <utility/vector1.hh>

// Class headers
#include <core/id/AtomID.hh>
#include <core/conformation/Conformation.fwd.hh>
#include <core/types.hh>
#include <core/kinematics/jacobian/ModuleType1.fwd.hh>

namespace core {
namespace kinematics {
namespace jacobian {

///@brief List of supported residue types in the loop that is to be represented by a series Jacobian
enum class SeriesJacobianTypeEnum {
	CANONICAL_AA = 1,
	NONCANONICAL_AA, // probably supported, but not yet. Already added to allow enum
	num_types = NONCANONICAL_AA
};

/// @brief The SeriesJacobians class is the mid-level of the Jacobian analysis of a protein's kinematics relations.
/// @author teunhoevenaars (teunhoevenaars@gmail.com)
class SeriesJacobians : public utility::VirtualBase {

public:

	/// @brief typedef a vector containing a series of residue numbers
	typedef utility::vector1<core::Size> residue_series;

	/// @brief No default constructor because SeriesJacobian is not used on its own, but always as part of a JacobianStructure
	SeriesJacobians() = delete;

	/// @brief Constructor based on a vector containing a series of residues
	SeriesJacobians(core::conformation::Conformation const & conformation, residue_series const & residue_set, core::id::AtomID const & ref_atom);

	/// @brief Copy constructor.
	SeriesJacobians(SeriesJacobians const & src);

	/// @brief Destructor.
	~SeriesJacobians() override;

	/// @brief Clone operation: make a copy of this object, and return an owning pointer to the copy.
	SeriesJacobiansOP clone() const;

public: //FUNCTIONS

	///@brief get residues that make up the chain
	utility::vector1< core::Size >
	get_residues(){
		return residue_set_;
	}

	///@brief update all Jacobian matrices in the chain
	void
	update_Jacobian_matrices(core::conformation::Conformation const & conformation);

private: //FUNCTIONS

	///@brief function that initializes the low-level Jacobian modules for a series of canonical amino acids
	utility::vector1< core::kinematics::jacobian::ModuleType1OP >
	init_modules_amino_acids();

	/// @brief function to determine the type of the residues that make up the series
	SeriesJacobianTypeEnum
	determine_residue_series_type( core::conformation::Conformation const & conformation, residue_series const & res_numbers );

public: // VARIABLES

	/// @brief vector with Jacobian modules
	utility::vector1< core::kinematics::jacobian::ModuleType1OP > modules_;

private: //VARIABLES
	/// @brief vector with residue series whose internal DoFs are free to move
	utility::vector1< core::Size > residue_set_;

	/// @brief number of dofs of the serial chain
	core::Size number_dofs_;

	/// @brief reference atom in which all vectors are expressed
	core::id::AtomID ref_atom_;
};

} //jacobian
} //kinematics
} //core

#endif //INCLUDED_core_kinematics_jacobian_SeriesJacobians_hh
