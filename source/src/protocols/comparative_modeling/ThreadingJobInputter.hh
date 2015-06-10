// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/comparative_modeling/ThreadingJobInputter.hh
/// @author Oliver Lange
/// @author James Thompson

#ifndef INCLUDED_protocols_comparative_modeling_ThreadingJobInputter_hh
#define INCLUDED_protocols_comparative_modeling_ThreadingJobInputter_hh

//unit headers
#include <protocols/jd2/JobInputter.hh>
#include <protocols/comparative_modeling/ThreadingJobInputter.fwd.hh>
#include <protocols/jd2/Job.fwd.hh>
#include <protocols/jd2/JobsContainer.hh>

//project headers
#include <core/pose/Pose.fwd.hh>
#ifdef WIN32
#include <core/sequence/SequenceAlignment.hh>
#endif

//utility headers
#include <map>

#include <core/sequence/SequenceAlignment.fwd.hh>
#include <utility/vector1.hh>

namespace protocols {
namespace comparative_modeling {

/// @details This is the simplest implementation of JobInputter
class ThreadingJobInputter : public protocols::jd2::JobInputter {
	typedef utility::vector1< core::sequence::SequenceAlignment > Alignments;
	typedef utility::vector1< core::pose::PoseOP > PoseOPs;
	typedef std::map< std::string, PoseOPs > PoseMap;
	typedef utility::vector1< std::string > FileList;
	typedef std::map< std::string, utility::vector1< core::Size > > ExtraResidues;
public:
	ThreadingJobInputter();
	/// @brief this function is responsible for filling the pose reference with
	/// the pose indicated by the job.  The Job object (within its InnerJob)
	/// contains a PoseCOP.  This function needs to either fill the pose
	/// reference from the InnerJob or, on first demand of a pose from that
	/// InnerJob, instantiate the pose, hand off a COP to the InnerJob, and fill
	/// the reference.
 	virtual void pose_from_job( core::pose::Pose & pose, protocols::jd2::JobOP job );

	/// @brief this function determines what jobs exist
	virtual void fill_jobs( protocols::jd2::JobsContainer & jobs );

	/// @brief Return the type of input source that the ThreadingJobInputter is
	/// currently using for template structures.
	virtual protocols::jd2::JobInputterInputSource::Enum input_source() const;

private:
	/// @brief Returns the number of templates in <template_poses_>
	size_t num_templates() const;

	Alignments alignments_;
	PoseMap template_poses_;
	ExtraResidues extra_residues_;
	protocols::jd2::JobInputterInputSource::Enum input_source_;
}; // ThreadingJobInputter

} // namespace comparative_modeling
} // namespace protocols

#endif //INCLUDED_protocols_comparative_modeling_ThreadingJobInputter_HH
