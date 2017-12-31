// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/jd3/chunk_library/ChunkLibraryInnerLarvalJob.hh
/// @brief  Class definition for the ChunkLibraryInnerLarvalJob, which includes the index of
///         the preliminary-job node that this ILJ originated from, or "0" if none.
/// @author Andy Watkins (amw579@stanford.edu)

#ifndef INCLUDED_protocols_jd3_chunk_library_ChunkLibraryInnerLarvalJob_HH
#define INCLUDED_protocols_jd3_chunk_library_ChunkLibraryInnerLarvalJob_HH

// Unit headers
#include <protocols/jd3/InnerLarvalJob.hh>
#include <protocols/jd3/chunk_library/ChunkLibraryInnerLarvalJob.fwd.hh>

// Package headers
#include <protocols/jd3/chunk_library_inputters/ChunkLibraryInputSource.fwd.hh>

// Project headers
#include <core/types.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/tag/Tag.fwd.hh>

// Basic headers
#include <basic/datacache/ConstDataMap.fwd.hh>
#include <basic/resource_manager/JobOptions.fwd.hh>

//C++ headers
#include <string>

#include <utility/vector1.hh>


#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace protocols {
namespace jd3 {
namespace chunk_library {

class ChunkLibraryInnerLarvalJob : public InnerLarvalJob {
public:

	ChunkLibraryInnerLarvalJob();
	ChunkLibraryInnerLarvalJob( core::Size nstruct );
	ChunkLibraryInnerLarvalJob( core::Size nstruct, core::Size prelim_job_node );
	ChunkLibraryInnerLarvalJob( ChunkLibraryInnerLarvalJob const & src );

	virtual ~ChunkLibraryInnerLarvalJob();

	bool
	operator == ( InnerLarvalJob const & other ) const override;

	/// @brief returns true if this is the same type as other;
	/// does not call other.same_type()
	bool
	same_type( InnerLarvalJob const & other ) const override;

	void
	show( std::ostream & out ) const override;

	friend
	std::ostream &
	operator<< ( std::ostream & out, const InnerLarvalJob & inner_job );

	core::Size
	prelim_job_node() const;

	void prelim_job_node( core::Size setting );

	std::string job_tag() const override { return "S"; }

private:

	core::Size prelim_job_node_;

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

}; // InnerLarvalJob

} // namespace chunk_library
} // namespace jd3
} // namespace protocols

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( protocols_jd3_chunk_library_ChunkLibraryInnerLarvalJob )
#endif // SERIALIZATION


#endif //INCLUDED_protocols_jd3_chunk_library_ChunkLibraryInnerLarvalJob_HH
