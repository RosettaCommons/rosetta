// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

#ifndef INCLUDED_protocols_jd2_LazySilentFileJobInputter_HH
#define INCLUDED_protocols_jd2_LazySilentFileJobInputter_HH

#include <protocols/jd2/JobInputter.hh>
#include <protocols/jd2/Job.fwd.hh>
#include <protocols/jd2/JobsContainer.hh>
#include <core/io/silent/SilentStruct.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/io/silent/SilentFileData.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace jd2 {

class LazySilentFileJobInputter : public protocols::jd2::JobInputter
{
public:

	LazySilentFileJobInputter();

	virtual ~LazySilentFileJobInputter();

	virtual void pose_from_job( core::pose::Pose & pose, JobOP job );

	virtual core::io::silent::SilentStruct const& struct_from_job( JobOP job );

	virtual void fill_jobs( JobsContainer & jobs );

	virtual JobInputterInputSource::Enum input_source() const;

	core::io::silent::SilentFileData const& silent_file_data() const { return sfd_; };


private:
	core::io::silent::SilentFileData sfd_;
};

} //jd2
} //protocols


#endif //INCLUDED_protocols_jd2_LazySilentFileJobInputter_HH
