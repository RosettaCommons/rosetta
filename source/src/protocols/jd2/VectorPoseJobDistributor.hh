// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/jd2/VectorPoseJobDistributor.hh
/// @brief  Job distributor subclass for running RECON multistate design
/// @brief  Takes in all input poses from the command line and passes them to any mover or filter that
/// @brief  derives from VectorPoseMover or VectorPoseFilter, meaning that it is able to receive and
/// @brief  operate on multiple poses simultaneously. Only accessible through recon application.
/// @author Alex Sevy (alex.sevy@gmail.com)

#ifndef VectorPoseJOBDISTRIBUTOR_HH_
#define VectorPoseJOBDISTRIBUTOR_HH_

#ifdef USEMPI
#include <mpi.h>
#endif

// #include <devel/init.hh>
#include <protocols/jd2/VectorPoseJobDistributor.fwd.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <protocols/jd2/Job.hh>
// #include <protocols/minimization_packing/PackRotamersMover.fwd.hh>
#include <protocols/rosetta_scripts/ParsedProtocol.fwd.hh>

namespace protocols {
namespace jd2 {

/// @brief  Job distributor subclass for running RECON multistate design
/// Takes in all input poses from the command line and passes them to any mover or filter that
/// derives from VectorPoseMover or VectorPoseFilter, meaning that it is able to receive and
/// operate on multiple poses simultaneously. Only accessible through recon application.
class VectorPoseJobDistributor : public JobDistributor {

public:

	VectorPoseJobDistributor();

	virtual ~VectorPoseJobDistributor();

	core::Size get_new_job_id();

	JobOP current_job() const;

	void mark_current_job_id_for_repetition();

	void job_failed( core::pose::Pose &, bool );

	void handle_interrupt();

	bool apply_parsed_protocol_mpi( core::pose::PoseOP & pose,
		protocols::rosetta_scripts::ParsedProtocolOP & protocol );

	bool apply_parsed_protocol_serial( utility::vector1< core::pose::PoseOP > & working_poses,
		utility::vector1<protocols::rosetta_scripts::ParsedProtocolOP> & protocols,
		utility::vector1< core::Size > & pose_order );

	void go( protocols::moves::MoverOP mover );

	void go_serial( protocols::moves::MoverOP mover );

	void go_mpi( protocols::moves::MoverOP mover );

private:
	JobOP this_nodes_job_;
	core::Size number_jobs_;
	std::map<core::Size, Jobs> job_map_;
	core::Size current_nstruct_;
	// utility::vector1< core::Size > my_designable_residues_;
	// protocols::minimization_packing::PackRotamersMoverOP packer_;
	// core::scoring::ScoreFunctionOP sfxn_;
	core::Size current_pose_;
	JobInputterOP job_inputter_;
	JobOutputterOP job_outputter_;
	bool randomize_input_;
};

} // jd2
} // protocols



#endif /* VectorPoseJOBDISTRIBUTOR_HH_ */
