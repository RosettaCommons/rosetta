// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/jd2/InnerMultiThreadingJob.hh
/// @author James Thompson

#ifndef INCLUDED_protocols_jd2_InnerMultiThreadingJob_hh
#define INCLUDED_protocols_jd2_InnerMultiThreadingJob_hh

#include <protocols/jd2/InnerJob.hh>
#include <protocols/jd2/InnerMultiThreadingJob.fwd.hh>

#include <core/pose/Pose.fwd.hh>
#include <core/sequence/SequenceAlignment.fwd.hh>

#include <utility/pointer/ReferenceCount.hh>
#include <core/types.hh>
#include <utility/vector1.hh>

#include <string>

namespace protocols {
namespace jd2 {

class InnerMultiThreadingJob : public protocols::jd2::InnerJob {
public:
	/// @brief ctor.  Note that it takes only the input tag and max nstruct,
	/// pose instantiation is deferred until the pose is needed
	/* InnerMultiThreadingJob(
	std::string const & input_tag,
	core::Size nstruct_max,
	utility::vector1< core::sequence::SequenceAlignment > const & alns,
	utility::vector1< core::pose::Pose > const & templates
	); */

	/// @brief ctor.
	InnerMultiThreadingJob(
		core::pose::PoseCOP,
		std::string const & input_tag,
		core::Size nstruct_max,
		utility::vector1< core::sequence::SequenceAlignment > const & alns,
		utility::vector1< core::pose::Pose > const & templates
	);

	utility::vector1< core::sequence::SequenceAlignment > const &
	alignments() const;

	utility::vector1< core::pose::Pose > const &
	template_poses() const;

	virtual ~InnerMultiThreadingJob();

private:
	utility::vector1< core::sequence::SequenceAlignment > const & alns_;
	utility::vector1< core::pose::Pose > const & template_poses_;
};

} // namespace jd2
} // namespace protocols

#endif //INCLUDED_protocols_jd2_InnerMultiThreadingJob_HH
