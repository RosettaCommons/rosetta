// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/comparative_modeling/ThreadingJob.hh
/// @brief  header file for ThreadingJob classes
/// @author

#ifndef INCLUDED_protocols_comparative_modeling_ThreadingJob_hh
#define INCLUDED_protocols_comparative_modeling_ThreadingJob_hh

//unit headers
#include <protocols/comparative_modeling/ThreadingJob.fwd.hh>
#include <protocols/jd2/InnerJob.hh>
#include <core/sequence/SequenceAlignment.hh>

//project headers

#include <protocols/jd2/Parser.fwd.hh> //for friendship

//utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <core/types.hh>

//C++ headers
#include <string>

#include <protocols/loops/Loops.fwd.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace comparative_modeling {

/// @details The InnerThreadingJob class is responsible for knowing input requirements for a given job - how many nstruct, and what the input is.  InnerThreadingJobs are relatively heavy; there is no need to duplicate a series of InnerThreadingJobs for each index into nstruct.  The companion ThreadingJob class handles the nstruct index and has a pointer to an InnerThreadingJob (which is shared across many ThreadingJobs).  InnerThreadingJob also holds a PoseOP to maintain the unmodified input pose for that job.
class ThreadingJob : public protocols::jd2::InnerJob {
public:
	/// @brief ctor.  Note that it takes only the input tag and max nstruct,
	/// pose instantiation is deferred until the pose is needed
	ThreadingJob(
		core::pose::PoseCOP template_pdb,
		core::sequence::SequenceAlignmentCOP alignment,
		std::string const & input_tag,
		core::Size nstruct_max
	);

	/// @brief the alignment for this Job
	core::sequence::SequenceAlignment const & alignment() const {
		return *alignment_;
	}

	/// @brief convenience : alignment id
	std::string alignment_id() const {
		return alignment().alignment_id();
	}

	/// @brief returns the "standard" loop definition (as conservative as possible)
	protocols::loops::Loops loops( core::Size nres ) const;

	/// @brief returns list of extra residues to steal
	utility::vector1< core::Size > const & extra_residues_to_steal() const;
	void extra_residues_to_steal( utility::vector1< core::Size > const & res );

	//NOTE:
	//get_pose() will give you the parent_pdb and not the starting model, the starting model is created from the ThreadingJobInputter
	// i have done it like this, since we can have multiple parent decoys for a given alignment read in via a silent file.
	// thus, to keep track which job should work on which parent-pdb + alignment I put them here.
	// alternative: have get_pose() refer to the starting_model
	// InnerJob --> contains a sufficiently uniq id ( or pointer ) to fish out the parent_decoy when starting the job.
private:
	core::sequence::SequenceAlignmentCOP alignment_; // alignment from input .ali file
	core::sequence::SequenceAlignmentCOP fasta2template_;
	core::Size nres_template_;
	// alignment from template seq (in .ali) to template pdb
	//    accounts for missing density in the input file
	utility::vector1< core::Size > extra_residues_to_steal_;
};

} // namespace comparative_modeling
} // namespace protocols

#endif //INCLUDED_protocols_comparative_modeling_ThreadingJob_HH
