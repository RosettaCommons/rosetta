// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/jd2/MSDJobDistributor.cc
/// @brief  Job distributor subclass for running restrained multistate design
/// @brief  Takes in all input poses from the command line and passes them to any mover that
/// @brief  derives from VectorPoseMover, meaning that it is able to receive and operate on multiple poses
/// @author Alex Sevy (alex.sevy@gmail.com)

#ifndef INCLUDED_protocols_jd2_MSDJobDistributor_hh
#define INCLUDED_protocols_jd2_MSDJobDistributor_hh

#include <protocols/jd2/MSDJobDistributor.fwd.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/JobInputter.fwd.hh>
#include <protocols/jd2/JobOutputter.fwd.hh>
#include <protocols/jd2/Job.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <protocols/rosetta_scripts/ParsedProtocol.fwd.hh>

namespace protocols {
namespace jd2 {

class MSDJobDistributor : public JobDistributor {

public:


	~MSDJobDistributor() override;

	core::Size
	get_new_job_id() override;

	JobOP
	current_job() const override;

	void
	mark_current_job_id_for_repetition() override;

	void
	job_failed( core::pose::Pose & /*pose*/,
		bool /*will_retry*/ ) override;

	void handle_interrupt() override;

	void
	go( protocols::moves::MoverOP mover ) override;

	friend class JobDistributorFactory;

	bool apply_parsed_protocol( utility::vector1< core::pose::PoseOP > & working_poses,
		utility::vector1<protocols::rosetta_scripts::ParsedProtocolOP> & protocols,
		utility::vector1< core::Size > & pose_order );

private:
	std::map<core::Size, Jobs> job_map_;
	core::Size current_nstruct_;
	core::Size current_pose_;
	JobInputterOP job_inputter_;
	JobOutputterOP job_outputter_;
	bool randomize_input_;
};

} //jd2
} //protocols

#endif //INCLUDED_protocols_jd2_MSDJobDistributor_HH
