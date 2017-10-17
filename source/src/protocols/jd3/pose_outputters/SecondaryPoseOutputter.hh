// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/jd3/pose_outputters/SecondaryPoseOutputter.hh
/// @brief  Definition of the %SecondaryPoseOutputter class
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)


#ifndef INCLUDED_protocols_jd3_pose_outputters_SecondaryPoseOutputter_HH
#define INCLUDED_protocols_jd3_pose_outputters_SecondaryPoseOutputter_HH

//unit headers
#include <protocols/jd3/pose_outputters/SecondaryPoseOutputter.fwd.hh>

//package headers
#include <protocols/jd3/pose_outputters/PoseOutputter.hh>
#include <protocols/jd3/JobOutputIndex.fwd.hh>
#include <protocols/jd3/LarvalJob.fwd.hh>
#include <protocols/jd3/InnerLarvalJob.fwd.hh>

//project headers
#include <core/pose/Pose.fwd.hh>

// utility headers
#include <utility/options/OptionCollection.fwd.hh>
#include <utility/options/keys/OptionKey.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/tag/Tag.fwd.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>

namespace protocols {
namespace jd3 {
namespace pose_outputters {

/// @brief The %SecondaryPoseOutputter
class SecondaryPoseOutputter : public PoseOutputter
{
public:

	SecondaryPoseOutputter();
	virtual ~SecondaryPoseOutputter();

	/// @brief This responsibility does not belong in the SecondaryPoseOutputter, so it receives
	/// a no-op implementation in the base class -- it is primarily the responsibility of the
	/// PoseOutputter base class and will not / should not be asked of the SecondaryPoseOutputter.
	void determine_job_tag(
		utility::tag::TagCOP output_tag,
		utility::options::OptionCollection const & job_options,
		InnerLarvalJob & job
	) const override final;


	/// @brief This responsibility does not belong in the SecondaryPoseOutputter, so it receives
	/// a no-frills implementation in the base class -- it is primarily the responsibility of the
	/// PoseOutputter base class and will not / should not be asked of the SecondaryPoseOutputter.
	bool job_has_already_completed( LarvalJob const & job, utility::options::OptionCollection const & options ) const override final;

	/// @brief This responsibility does not belong in the SecondaryPoseOutputter, so it receives
	/// a no-op implementation in the base class -- it is primarily the responsibility of the
	/// PoseOutputter base class and will not / should not be asked of the SecondaryPoseOutputter.
	void mark_job_as_having_started( LarvalJob const & job, utility::options::OptionCollection const & options ) const override final;

};

} // namespace pose_outputters
} // namespace jd3
} // namespace protocols

#endif //INCLUDED_protocols_jd3_pose_outputters_SecondaryPoseOutputter_HH
