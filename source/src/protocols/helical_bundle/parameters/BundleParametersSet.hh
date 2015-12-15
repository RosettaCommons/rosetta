// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/helical_bundle/parameters/BundleParametersSet.hh
/// @brief  Prototypes and method declarations for the BundleParametersSet class, a class for holding sets of parameters for parametric helical bundle generation.
/// @author Vikram K. Mulligan (vmullig@uw.edu)


#ifndef INCLUDED_protocols_helical_bundle_parameters_BundleParametersSet_hh
#define INCLUDED_protocols_helical_bundle_parameters_BundleParametersSet_hh


// Unit headers
#include <protocols/helical_bundle/parameters/BundleParametersSet.fwd.hh>

// Package headers
#include <core/conformation/parametric/Parameters.fwd.hh>
#include <core/conformation/parametric/ParametersSet.fwd.hh>
#include <core/conformation/parametric/Parameters.hh>
#include <core/conformation/parametric/ParametersSet.hh>
#include <protocols/helical_bundle/parameters/BundleParameters.fwd.hh>
#include <protocols/helical_bundle/parameters/BundleParameters.hh>

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

/// @brief  BundleParametersSet class, used to store sets of parameters for parametric helical bundle generation.
///
class BundleParametersSet : public core::conformation::parametric::ParametersSet
{

public: //Typedefs:

	typedef core::conformation::parametric::Parameters Parameters;
	typedef core::conformation::parametric::ParametersOP ParametersOP;
	typedef core::conformation::parametric::ParametersSet ParametersSet;
	typedef core::conformation::parametric::ParametersSetOP ParametersSetOP;

public:

	/// @brief constructors
	///
	BundleParametersSet();

	BundleParametersSet( BundleParametersSet const & src );

	~BundleParametersSet();

	/// @brief Copy this residue( allocate actual memory for it )
	///
	ParametersSetOP clone() const;

public: //Getters

	/// @brief Returns the symmetry of the bundle created.
	/// @details A value of 0 or 1 indicates no symmetry.  Larger values indicate n-fold radial symmetry
	/// (for example, 3 means threefold radial symmetry about the bundle axis, and each helix defined will
	/// be replicated a total of three times).
	core::Size bundle_symmetry() const { return bundle_symmetry_; }

	/// @brief Returns the number of symmetry copies to generate.
	/// @details A value of 0 means to generate all copies.  Higher values mean to generate only the first N
	/// copies.  For example, if the symmetry were 16 but bundle_symmetry_copies_ were set to 4, only the
	/// first 4 symmetry repeats would be generated.
	core::Size bundle_symmetry_copies() const { return bundle_symmetry_copies_; }

	/// @brief Get the number of helices defined in each symmetry copy of this bundle.
	///
	core::Size n_helices() const { return n_helices_; }


public: //Setters

	/// @brief Sets the symmetry of the bundle created.
	/// @details A value of 0 or 1 indicates no symmetry.  Larger values indicate n-fold radial symmetry
	/// (for example, 3 means threefold radial symmetry about the bundle axis, and each helix defined will
	/// be replicated a total of three times).
	void set_bundle_symmetry( core::Size const val ) { bundle_symmetry_ = val; return; }

	/// @brief Sets the number of symmetry copies to generate.
	/// @details A value of 0 means to generate all copies.  Higher values mean to generate only the first N
	/// copies.  For example, if the symmetry were 16 but bundle_symmetry_copies_ were set to 4, only the
	/// first 4 symmetry repeats would be generated.
	void set_bundle_symmetry_copies( core::Size const val ) { bundle_symmetry_copies_ = val; return; }

	/// @brief Set the number of helices defined in each symmetry copy of this bundle.
	///
	void set_n_helices( core::Size const val ) { n_helices_=val; return; }

	/// @brief Get a summary of this ParametersSet object, for output to remark lines of a PDB file.
	/// @details Default function can be overridden by derived classes.  This version actually outputs
	/// Crick parameter information.
	virtual void get_pdb_remark(std::stringstream &remark) const;

private:

	/********************************************************************************
	PRIVATE DATA
	*********************************************************************************/

	/// @brief The symmetry of the bundle created.
	/// @details A value of 0 or 1 indicates no symmetry.  Larger values indicate n-fold radial symmetry
	/// (for example, 3 means threefold radial symmetry about the bundle axis, and each helix defined will
	/// be replicated a total of three times).
	core::Size bundle_symmetry_;

	/// @brief The symmetry copies to generate.
	/// @details A value of 0 means to generate all copies.  Higher values mean to generate only the first N
	/// copies.  For example, if the symmetry were 16 but bundle_symmetry_copies_ were set to 4, only the
	/// first 4 symmetry repeats would be generated.
	core::Size bundle_symmetry_copies_;

	/// @brief The number of helices defined in each symmetry copy of this bundle.
	///
	core::Size n_helices_;


#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

}; //class BundleParametersSet

} // namespace parameters
} // namespace helical_bundle
} // namespace protocols

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( protocols_helical_bundle_parameters_BundleParametersSet )
#endif // SERIALIZATION


#endif
