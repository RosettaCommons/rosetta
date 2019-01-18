// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file --path--/--class--.cc
/// @brief --brief--
/// @author --name-- (--email--)

// Unit headers
#include <--path--/--class--.hh>


// Basic/Utility headers

#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION


--namespace--
using namespace protocols::jd3;

//Constructor
--class--::--class--()
{}

//Destructor
--class--::~--class--() = default;

--end_namespace--

#ifdef    SERIALIZATION
template< class Archive >
void
--namespace_2colon--::--class--::load( Archive & arc ) {
	arc( cereal::base_class< protocols::jd3::JobSummary >( this ) );
	//arc( energy_ ); // core::Real
}

template< class Archive >
void
--namespace_2colon--::--class--::save( Archive & arc ) const {
	arc( cereal::base_class< protocols::jd3::JobSummary >( this ) );
	//arc( CEREAL_NVP( energy_ ) ); // core::Real
}

SAVE_AND_LOAD_SERIALIZABLE( --namespace_2colon--::--class-- );
CEREAL_REGISTER_TYPE( --namespace_2colon--::--class-- )

CEREAL_REGISTER_DYNAMIC_INIT(--namespace_underscore--_--class--)
#endif // SERIALIZATION
