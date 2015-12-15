// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/conformation/membrane/LipidAccInfo.hh
///
/// @brief      Membrane Lipid Accessibility Data
/// @details    Object for storing per-residue lipid exposed and buried surface
///    area values. Predicted from sequence, transmembrane spans, and psiblast
///    prediction using server called from the run_lips.pl script.
///    Last Modified: 7/7/14
///
/// @author  Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_core_conformation_membrane_LipidAccInfo_hh
#define INCLUDED_core_conformation_membrane_LipidAccInfo_hh

// Unit headers
#include <core/conformation/membrane/LipidAccInfo.fwd.hh>

// Package Headers
#include <core/types.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION


namespace core {
namespace conformation {
namespace membrane {

/// @brief      Membrane Lipid Accessibility Data
/// @details    Stores lipid accessibility data derived from OCTOPUS spanning file
///             and psiblast search using run_lips.pl script
class LipidAccInfo : public utility::pointer::ReferenceCount {

public: // constructors

	/// @brief Constructor
	/// @details Create a blank copy of the lipid accessibility data object
	LipidAccInfo();

	/// @brief Custom Constructor
	/// @brief Construct from user-provided lipid Acc Info File
	LipidAccInfo( std::string lipsfile );

	/// @brief Conpy Constructor
	/// @details Create a deep copy of this object
	LipidAccInfo( LipidAccInfo const & src );

	/// @brief Assignment Operator
	/// @details Create a deep copy of this object, overloading the assignment operator
	LipidAccInfo &
	operator=( LipidAccInfo const & src );

	/// @brief Destructor
	~LipidAccInfo();

public: // data access

	/// @brief Access Lipid exposed surface area per-residue
	utility::vector1< core::Real > lipid_exposure();

	/// @details Access Lipid buried surface area per-residue
	utility::vector1< core::Real > lipid_burial();

private: // helper methods

	/// @brief Copy Data
	void copy_data( LipidAccInfo src, LipidAccInfo copy );

private: // data

	// Lipid burial and exposure
	utility::vector1< core::Real > lipid_exposure_;
	utility::vector1< core::Real > lipid_burial_;

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

}; // class LipidAccInfo

} // membrane
} // conformation
} // core

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_conformation_membrane_LipidAccInfo )
#endif // SERIALIZATION


#endif // INCLUDED_core_conformation_membrane_LipidAccInfo_hh

