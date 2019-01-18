// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/jd3/chunk_library/ChunkLibraryInnerLarvalJob.cc
/// @brief  Method definitions for the InnerLarvalJob class
/// @author Andy Watkins (amw579@stanford.edu)

// Unit headers
#include <protocols/jd3/chunk_library/ChunkLibraryInnerLarvalJob.hh>

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
namespace chunk_library {

ChunkLibraryInnerLarvalJob::ChunkLibraryInnerLarvalJob() :
	InnerLarvalJob()
{}


ChunkLibraryInnerLarvalJob::ChunkLibraryInnerLarvalJob( core::Size nstruct, core::Size prelim_job_node ) :
	InnerLarvalJob( nstruct, prelim_job_node )
{}

ChunkLibraryInnerLarvalJob::ChunkLibraryInnerLarvalJob( ChunkLibraryInnerLarvalJob const & src ) :
	InnerLarvalJob( src )
{}


ChunkLibraryInnerLarvalJob::~ChunkLibraryInnerLarvalJob() = default;

/// @brief Mutual comparison of this inner job to the other inner job
/// so that if either one thinks it's not the same as the other, then
/// it returns false.  Invokes the same() function on both this and other
bool
ChunkLibraryInnerLarvalJob::operator == ( InnerLarvalJob const & other ) const
{
	if ( InnerLarvalJob::operator == ( other ) ) {
		//ChunkLibraryInnerLarvalJob const & other_std = static_cast< ChunkLibraryInnerLarvalJob const & > ( other );
		return true;
	}
	return false;
}

/// @details Since this is the base-class function, then by construction
/// other is equivalent to this.
/// @note classes derived from ChunkLibraryInnerLarvalJob must perform dynamic casts
/// to ensure the other ChunkLibraryInnerLarvalJob has the same type as them
bool
ChunkLibraryInnerLarvalJob::same_type( InnerLarvalJob const & other ) const
{
	return dynamic_cast< ChunkLibraryInnerLarvalJob const * > ( &other );
}

void
ChunkLibraryInnerLarvalJob::show( std::ostream & out ) const
{
	out << "ChunkLibraryInnerLarvalJob::show stubbed out";
}

std::ostream &
operator<< ( std::ostream & out, const ChunkLibraryInnerLarvalJob & inner_job )
{
	inner_job.show( out );
	return out;
}

} // namespace chunk_library
} // namespace jd3
} // namespace protocols

#ifdef    SERIALIZATION

template< class Archive >
void
protocols::jd3::chunk_library::ChunkLibraryInnerLarvalJob::save( Archive & arc ) const {
	arc( cereal::base_class< protocols::jd3::InnerLarvalJob > ( this ) );
}

template< class Archive >
void
protocols::jd3::chunk_library::ChunkLibraryInnerLarvalJob::load( Archive & arc ) {
	arc( cereal::base_class< protocols::jd3::InnerLarvalJob > ( this ) );
}

SAVE_AND_LOAD_SERIALIZABLE( protocols::jd3::chunk_library::ChunkLibraryInnerLarvalJob );
CEREAL_REGISTER_TYPE( protocols::jd3::chunk_library::ChunkLibraryInnerLarvalJob )

CEREAL_REGISTER_DYNAMIC_INIT( protocols_jd3_chunk_library_ChunkLibraryInnerLarvalJob )
#endif // SERIALIZATION
