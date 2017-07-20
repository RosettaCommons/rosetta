// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/jd3/full_model/FullModelInnerLarvalJob.hh
/// @brief  Class definition for the FullModelInnerLarvalJob, which includes the index of
///         the preliminary-job node that this ILJ originated from, or "0" if none.
/// @author Andy Watkins (amw579@stanford.edu)

#ifndef INCLUDED_protocols_jd3_full_model_FullModelInnerLarvalJob_HH
#define INCLUDED_protocols_jd3_full_model_FullModelInnerLarvalJob_HH

// Unit headers
#include <protocols/jd3/InnerLarvalJob.hh>
#include <protocols/jd3/full_model/FullModelInnerLarvalJob.fwd.hh>

// Package headers
#include <protocols/jd3/full_model_inputters/FullModelInputSource.fwd.hh>

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
namespace full_model {

class FullModelInnerLarvalJob : public InnerLarvalJob {
public:

	FullModelInnerLarvalJob();
	FullModelInnerLarvalJob( core::Size nstruct );
	FullModelInnerLarvalJob( core::Size nstruct, core::Size prelim_job_node );
	FullModelInnerLarvalJob( FullModelInnerLarvalJob const & src );

	virtual ~FullModelInnerLarvalJob();

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

private:

	core::Size prelim_job_node_;

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

}; // InnerLarvalJob

} // namespace full_model
} // namespace jd3
} // namespace protocols

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( protocols_jd3_full_model_FullModelInnerLarvalJob )
#endif // SERIALIZATION


#endif //INCLUDED_protocols_jd3_full_model_FullModelInnerLarvalJob_HH
