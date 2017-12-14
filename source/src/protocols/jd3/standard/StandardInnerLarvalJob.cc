// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/jd3/standard/StandardInnerLarvalJob.cc
/// @brief  Method definitions for the InnerLarvalJob class
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Unit headers
#include <protocols/jd3/standard/StandardInnerLarvalJob.hh>

//C++ headers
#include <string>
#include <sstream>

// Utility headers
#include <utility/vector1.hh>


#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION

namespace protocols {
namespace jd3 {
namespace standard {

StandardInnerLarvalJob::StandardInnerLarvalJob() :
	InnerLarvalJob(),
	prelim_job_node_( 0 )
{}

StandardInnerLarvalJob::StandardInnerLarvalJob( core::Size nstruct ) :
	InnerLarvalJob( nstruct ),
	prelim_job_node_( 0 )
{}

StandardInnerLarvalJob::StandardInnerLarvalJob( core::Size nstruct, core::Size prelim_job_node ) :
	InnerLarvalJob( nstruct ),
	prelim_job_node_( prelim_job_node )
{}

StandardInnerLarvalJob::StandardInnerLarvalJob( StandardInnerLarvalJob const & /*src*/ ) = default;


StandardInnerLarvalJob::~StandardInnerLarvalJob() = default;

/// @brief Mutual comparison of this inner job to the other inner job
/// so that if either one thinks it's not the same as the other, then
/// it returns false.  Invokes the same() function on both this and other
bool
StandardInnerLarvalJob::operator == ( InnerLarvalJob const & other ) const
{
	if ( InnerLarvalJob::operator == ( other ) ) {
		auto const & other_std = static_cast< StandardInnerLarvalJob const & > ( other );
		return prelim_job_node_ == other_std.prelim_job_node_;
	}
	return false;
}

/// @details Since this is the base-class function, then by construction
/// other is equivalent to this.
/// @note classes derived from StandardInnerLarvalJob must perform dynamic casts
/// to ensure the other StandardInnerLarvalJob has the same type as them
bool
StandardInnerLarvalJob::same_type( InnerLarvalJob const & other ) const
{
	return dynamic_cast< StandardInnerLarvalJob const * > ( &other );
}

void
StandardInnerLarvalJob::show( std::ostream & out ) const
{
	out << "StandardInnerLarvalJob::show stubbed out";
}

std::ostream &
operator<< ( std::ostream & out, const StandardInnerLarvalJob & inner_job )
{
	inner_job.show( out );
	return out;
}

core::Size StandardInnerLarvalJob::prelim_job_node() const
{
	return prelim_job_node_;
}

void
StandardInnerLarvalJob::prelim_job_node( core::Size setting )
{
	prelim_job_node_ = setting;
}

} // namespace standard
} // namespace jd3
} // namespace protocols

#ifdef    SERIALIZATION

template< class Archive >
void
protocols::jd3::standard::StandardInnerLarvalJob::save( Archive & arc ) const {
	arc( cereal::base_class< protocols::jd3::InnerLarvalJob > ( this ) );
	arc( prelim_job_node_ );
}

template< class Archive >
void
protocols::jd3::standard::StandardInnerLarvalJob::load( Archive & arc ) {
	arc( cereal::base_class< protocols::jd3::InnerLarvalJob > ( this ) );
	arc( prelim_job_node_ );
}

SAVE_AND_LOAD_SERIALIZABLE( protocols::jd3::standard::StandardInnerLarvalJob );
CEREAL_REGISTER_TYPE( protocols::jd3::standard::StandardInnerLarvalJob )

CEREAL_REGISTER_DYNAMIC_INIT( protocols_jd3_standard_StandardInnerLarvalJob )
#endif // SERIALIZATION
