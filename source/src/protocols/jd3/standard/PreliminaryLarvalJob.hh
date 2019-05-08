// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/jd3/standard/PreliminaryLarvalJob.hh
/// @brief A simple class for that stores input information for each job defintion/input nstruct as a PreliminaryLarvalJob.
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)


#ifndef INCLUDED_protocols_jd3_standard_PreliminaryLarvalJob_hh
#define INCLUDED_protocols_jd3_standard_PreliminaryLarvalJob_hh

#include <protocols/jd3/standard/PreliminaryLarvalJob.fwd.hh>

//hh Needed by Binder
#include <protocols/jd3/InnerLarvalJob.hh>
#include <protocols/jd3/pose_inputters/PoseInputter.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/options/OptionCollection.hh>
#include <utility/tag/Tag.hh>

namespace protocols {
namespace jd3 {
namespace standard {

///@brief A PreliminaryLarvalJob is what gets created from the command-line options job_definition file for each input
///  and specified job in the Job Definition file.
///
/// The SJQ has a list of these after they are created in determine_preliminary_job_list()
///
struct PreliminaryLarvalJob
{
public:
	PreliminaryLarvalJob();
	~PreliminaryLarvalJob();
	PreliminaryLarvalJob( PreliminaryLarvalJob const & src );
	PreliminaryLarvalJob & operator = ( PreliminaryLarvalJob const & rhs );

	InnerLarvalJobOP inner_job;
	utility::tag::TagCOP job_tag;
	utility::options::OptionCollectionCOP job_options;
	pose_inputters::PoseInputterOP pose_inputter;
#ifdef    SERIALIZATION
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};


} //protocols
} //jd3
} //standard



#endif //INCLUDED_protocols_jd3_standard_PreliminaryLarvalJob_hh





