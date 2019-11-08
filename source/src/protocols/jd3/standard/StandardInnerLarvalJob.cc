// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/jd3/standard/StandardInnerLarvalJob.cc
/// @brief  Method definitions for the StandardInnerLarvalJob class
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Unit headers
#include <protocols/jd3/standard/StandardInnerLarvalJob.hh>

#include <utility/tag/Tag.hh>

// Basic headers
#include <basic/Tracer.hh>

#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION

namespace protocols {
namespace jd3 {
namespace standard {

static basic::Tracer TR( "protocols.jd3.standard.StandardInnerLarvalJob" );

StandardInnerLarvalJob::StandardInnerLarvalJob() :
	parent(),
	preliminary_job_node_( 0 )
{}

StandardInnerLarvalJob::StandardInnerLarvalJob( StandardInnerLarvalJob const & src ) :
	parent(src),
	preliminary_job_node_( src.preliminary_job_node_ )
{}


StandardInnerLarvalJob::StandardInnerLarvalJob(
	core::Size nstruct,
	JobDAGNodeID job_node,
	PrelimJobNodeID preliminary_job_node
) :
	parent( nstruct, job_node ),
	preliminary_job_node_( preliminary_job_node )
{}

StandardInnerLarvalJob::~StandardInnerLarvalJob() = default;

/// @brief Mutual comparison of this inner job to the other inner job
/// so that if either one thinks it's not the same as the other, then
/// it returns false.  Invokes the same() function on both this and other
bool
StandardInnerLarvalJob::operator == ( InnerLarvalJob const & other ) const
{
	if ( parent::operator == (other) ) {
		StandardInnerLarvalJob const & mrs_other = static_cast< StandardInnerLarvalJob const & > (other);
		return preliminary_job_node_ == mrs_other.preliminary_job_node_;
	}
	return false;
}

bool
StandardInnerLarvalJob::same_type( InnerLarvalJob const & other ) const
{
	return dynamic_cast< StandardInnerLarvalJob const * >( &other );
}

void
StandardInnerLarvalJob::show( std::ostream & out ) const
{
	out << "StandardInnerLarvalJob " << preliminary_job_node_ << " parent: ";
	parent::show(out);
}

std::ostream &
operator<< ( std::ostream & out, const StandardInnerLarvalJob & inner_job )
{
	inner_job.show( out );
	return out;
}

PrelimJobNodeID
StandardInnerLarvalJob::preliminary_job_node() const {
	return preliminary_job_node_;
}

} // namespace standard
} // namespace jd3
} // namespace protocols

#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
protocols::jd3::standard::StandardInnerLarvalJob::save( Archive & arc ) const {
	arc( cereal::base_class< jd3::InnerLarvalJob >( this ));
	arc( CEREAL_NVP( preliminary_job_node_ )); // core::Size
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
protocols::jd3::standard::StandardInnerLarvalJob::load( Archive & arc ) {
	arc( cereal::base_class< jd3::InnerLarvalJob >( this ));
	arc( preliminary_job_node_ ); // core::Size
}

SAVE_AND_LOAD_SERIALIZABLE( protocols::jd3::standard::StandardInnerLarvalJob );
CEREAL_REGISTER_TYPE( protocols::jd3::standard::StandardInnerLarvalJob )

CEREAL_REGISTER_DYNAMIC_INIT( protocols_jd3_standard_StandardInnerLarvalJob )
#endif // SERIALIZATION
