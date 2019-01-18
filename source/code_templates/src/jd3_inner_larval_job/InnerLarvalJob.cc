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

#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION

--namespace--

//Constructor
--class--::--class--()
{}

--class--::--class--(core::Size nstruct, core::Size job_node):
	InnerLarvalJob(nstruct, job_node)
{}

//Destructor
--class--::~--class--() = default;

--class--::--class--( --class-- const & /*src*/ ) = default;


/// @brief Mutual comparison of this inner job to the other inner job
/// so that if either one thinks it's not the same as the other, then
/// it returns false.  Invokes the same() function on both this and other
bool
--class--::operator == ( InnerLarvalJob const & other ) const
{
	if ( InnerLarvalJob::operator == ( other ) ) {
		//auto const & other_std = static_cast< --class-- const & > ( other );
		//return some_var_ == other_std.some_var_;
		return true;
	}
	return false;
}

/// @details Since this is the base-class function, then by construction
/// other is equivalent to this.
/// @note classes derived from StandardInnerLarvalJob must perform dynamic casts
/// to ensure the other StandardInnerLarvalJob has the same type as them
bool
--class--::same_type( InnerLarvalJob const & other ) const
{
	return dynamic_cast< --class-- const * > ( &other );
}

void
--class--::show( std::ostream & out ) const
{
	out << "--class--::show stubbed out";
}

std::ostream &
operator<< ( std::ostream & out, const --class-- & inner_job )
{
	inner_job.show( out );
	return out;
}


--end_namespace--

#ifdef    SERIALIZATION

template< class Archive >
void
--namespace_2colon--::--class--::load( Archive & arc ) {
	arc( cereal::base_class< protocols::jd3::InnerLarvalJob >( this ) );
	//arc( energy_ ); // core::Real
}

template< class Archive >
void
--namespace_2colon--::--class--::save( Archive & arc ) const {
	arc( cereal::base_class< protocols::jd3::InnerLarvalJob >( this ) );
	//arc( CEREAL_NVP( energy_ ) ); // core::Real
}

SAVE_AND_LOAD_SERIALIZABLE( --namespace_2colon--::--class-- );
CEREAL_REGISTER_TYPE( --namespace_2colon--::--class-- )

CEREAL_REGISTER_DYNAMIC_INIT( --namespace_underscore--_--class-- )

#endif // SERIALIZATION
