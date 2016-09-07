// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/jd2/PoseInputStreamJobInputter.hh
/// @brief  header file for PoseInputStreamJobInputter class
/// @author James Thompson


#ifndef INCLUDED_protocols_jd2_PoseInputStreamJobInputter_hh
#define INCLUDED_protocols_jd2_PoseInputStreamJobInputter_hh

#include <protocols/jd2/JobInputter.hh>
#include <protocols/jd2/PoseInputStreamJobInputter.fwd.hh>
#include <protocols/jd2/Job.fwd.hh>
#include <protocols/jd2/JobsContainer.hh>
#include <core/import_pose/pose_stream/MetaPoseInputStream.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace jd2 {

/// @details This is the simplest implementation of JobInputter, which reads from -s/-l and SilentFile files.
class PoseInputStreamJobInputter : public protocols::jd2::JobInputter
{
public:

	PoseInputStreamJobInputter();

	~PoseInputStreamJobInputter() override;

	/// @brief This implementation simply calls fill_pose on the PoseInputStream
	/// object.
	void pose_from_job( core::pose::Pose & pose, JobOP job ) override;

	void fill_jobs( JobsContainer & jobs ) override;

	/// @brief Return the type of input source that the PoseInputStreamJobInputter is currently
	///  using.
	JobInputterInputSource::Enum input_source() const override;

private:
	core::chemical::ResidueTypeSetCOP rsd_set_;
	core::import_pose::pose_stream::MetaPoseInputStream input_;
}; // PoseInputStreamJobInputter

} // namespace jd2
} // namespace protocols

#endif //INCLUDED_protocols_jd2_PoseInputStreamJobInputter_HH
