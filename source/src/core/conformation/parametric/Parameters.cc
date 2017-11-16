// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/core/conformation/parametric/Parameters.cc
/// @brief  A class for holding sets of parameters for parametric backbone generation.
/// @author Vikram K. Mulligan (vmullig@uw.edu)

// Unit header
#include <core/conformation/parametric/Parameters.hh>

// Package headers

// Project headers

// Basic headers
#include <basic/Tracer.hh>

// Numeric headers

// Utility Headers
#include <utility/vector1.hh>


#ifdef SERIALIZATION
// Utility serialization headers
#include <utility/vector1.srlz.hh>
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/base_class.hpp>
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION


namespace core {
namespace conformation {
namespace parametric {

static basic::Tracer TR( "core.conformation.parametric.Parameters" );

/// @brief Constructor.
///
Parameters::Parameters() :
	residue_list_()
{
}

Parameters::Parameters( Parameters const & src ) :
	utility::pointer::ReferenceCount(),
	utility::pointer::enable_shared_from_this< Parameters >()
{
	residue_list_.clear();
	if ( src.residue_list_.size()>0 ) {
		for ( core::Size i=1, imax=src.residue_list_.size(); i<=imax; ++i ) {
			residue_list_.push_back( src.residue_list_[i]->clone() ); //This copies the residue that was being pointed at.
			//Note that when copying a Conformation, I need to add logic that will ensure that the Parameters objects that result have owning pointers to the residues in the Conformation,
			//rather than to residues that only exist in the Parameters object.
		}
	}
}

Parameters::~Parameters() {}


/// @brief make a copy of this residue( allocate actual memory for it )
///
ParametersOP
Parameters::clone() const
{
	return ParametersOP( new Parameters( *this ) );
}

} // namespace parametric
} // namespace conformation
} // namespace core


#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
core::conformation::parametric::Parameters::save( Archive & arc ) const {
	arc( CEREAL_NVP( residue_list_ ) ); // utility::vector1<core::conformation::ResidueCOP>
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::conformation::parametric::Parameters::load( Archive & arc ) {
	utility::vector1< std::shared_ptr< core::conformation::Residue > > local_residue_list;
	arc( local_residue_list ); // utility::vector1<core::conformation::ResidueCOP>
	residue_list_ = local_residue_list; // copy the non-const pointer(s) into the const pointer(s)
}

SAVE_AND_LOAD_SERIALIZABLE( core::conformation::parametric::Parameters );
CEREAL_REGISTER_TYPE( core::conformation::parametric::Parameters )

CEREAL_REGISTER_DYNAMIC_INIT( core_conformation_parametric_Parameters )
#endif // SERIALIZATION
