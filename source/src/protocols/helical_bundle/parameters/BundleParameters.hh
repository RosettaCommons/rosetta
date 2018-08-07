// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/helical_bundle/parameters/BundleParameters.hh
/// @brief  Prototypes and method declarations for the BundleParameters class, a class for holding parameters for helical bundle backbone generation.
/// @author Vikram K. Mulligan (vmullig@uw.edu)


#ifndef INCLUDED_protocols_helical_bundle_parameters_BundleParameters_hh
#define INCLUDED_protocols_helical_bundle_parameters_BundleParameters_hh

// Unit headers
#include <protocols/helical_bundle/parameters/BundleParameters.fwd.hh>

// Package headers
#include <core/conformation/parametric/Parameters.fwd.hh>
#include <core/conformation/parametric/ParametersSet.fwd.hh>
#include <core/conformation/parametric/Parameters.hh>
#include <core/conformation/parametric/ParametersSet.hh>

// Project headers
#include <core/types.hh>

// Utility headers
#include <utility/pointer/access_ptr.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>

// Numeric headers

// C++ headers


#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace protocols {
namespace helical_bundle {
namespace parameters {

/// @brief  Parameters class, used to store sets of parameters for parametric backbone generation.
///
class BundleParameters : public core::conformation::parametric::Parameters
{

public: //Typedefs:

	typedef core::conformation::parametric::Parameters Parameters;
	typedef core::conformation::parametric::ParametersOP ParametersOP;
	typedef core::conformation::parametric::ParametersSet ParametersSet;
	typedef core::conformation::parametric::ParametersSetOP ParametersSetOP;

public:

	/// @brief constructors
	///
	BundleParameters();

	BundleParameters( BundleParameters const & src );

	~BundleParameters();

	/// @brief Copy this residue( allocate actual memory for it )
	///
	ParametersOP clone() const;

	///////////////////
	//Output:
	///////////////////

	/// @brief Get a summary of this ParametersSet object, for output to remark lines of a PDB file.
	///
	virtual void get_pdb_remark(std::stringstream &remark) const;

private:

	/********************************************************************************
	PRIVATE DATA
	*********************************************************************************/

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

}; //class BundleParameters

} // namespace parametric
} // namespace conformation
} // namespace core

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( protocols_helical_bundle_parameters_BundleParameters )
#endif // SERIALIZATION


#endif
