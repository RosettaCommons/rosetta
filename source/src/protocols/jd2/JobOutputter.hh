// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/jd2/JobOutputter.hh
/// @brief  header file for JobOutputter class, part of August 2008 job distributor as planned at RosettaCon08
/// @author Steven Lewis smlewi@gmail.com


#ifndef INCLUDED_protocols_jd2_JobOutputter_hh
#define INCLUDED_protocols_jd2_JobOutputter_hh

//unit headers
#include <protocols/jd2/JobOutputter.fwd.hh>
#include <protocols/jd2/Job.fwd.hh>
#include <protocols/jd2/JobOutputterObserver.hh>
#ifdef WIN32
#include <protocols/jd2/Job.hh> // WIN32 INCLUDE
#endif

//project headers
#include <core/pose/Pose.fwd.hh>

#include <core/io/silent/SilentStruct.fwd.hh>

#include <protocols/evaluation/PoseEvaluator.fwd.hh>

//utility headers
#include <utility/pointer/ReferenceCount.hh>

//C++ headers
#include <string>

#include <utility/vector1.hh>
#include <set>

namespace protocols {
namespace jd2 {

/// @details the JobOutputter class is responsible for dealing with output, as well as determining what
/// jobs have already been output and what sort of name is associated with a job.  Derived classes will
/// be responsible for output such as PDBS, PDBS.tar.gz and silent files.
class JobOutputter : public utility::pointer::ReferenceCount
{
public:

	//constructor -- reads cmd-line to initialize evaluators
	JobOutputter();

	virtual ~JobOutputter();

	/// @brief this function is meant to be redefined in child classes to allow for flushing of memory buffers.
	/// Here's the long version: The SilentFileJobOutputter wanted to buffer output, but needed to guaruntee that
	/// the output would be flushed at end of runtime.  The original implementation was to A) bend over backward to ensure
	/// that the destructor was run (JobOutputter lives inside static JobDistributor, which was previously not destructed
	/// because it's static) and B) flush the buffers in the destructor.  This caused a problem because the buffer-flushing
	/// tried to use the Tracers, which had already been destructed...boom crash.
	///
	///  New solution: re-forbid use of destructors within the static JobDistributor system, and create a flush function to force this stuff out.  So here it is:
	virtual void flush();

	//////////////////////////////creating output functions/////////////////////////////////////////

	/// @brief this function takes a string and writes it to disk (separately from Tracer output).
	///use some sort of extention option system - default .dat?  .data?
	virtual
	void file( JobCOP job, std::string const & data ) = 0;

	/// @brief this function outputs the final result of a job.
	virtual
	void final_pose( JobOP job, core::pose::Pose const & pose, std::string const & tag = "" ) = 0;

	/// @brief this function is intended for saving mid-protocol poses; for example the final centroid structure in a combined centroid/fullatom protocol.
	virtual
	void other_pose( JobOP job, core::pose::Pose const & pose, std::string const & tag, int copy_count = -1, bool score_only = false ) = 0;

	/// @brief optionally pass a starting (reference) pose to a JobOutputter for later comparison purposes and/or as interface for initializing evaluators
	virtual
	void starting_pose( core::pose::Pose const & );

	/////////////////////////////////state of output functions/////////////////////////////////

	/// @brief this function is not used for output, but it belongs here since it needs to check the same output locations as the class normally writes to.  This class checks wherever output goes to see if the job's expected output already exists (on disk or whatever).  This is the most basic form of checkpointing.
	virtual
	bool job_has_completed( JobCOP job ) = 0;

	/// @brief this is the master function for determining the unique output identifier for a job - should this be defined in the base class?
	virtual
	std::string output_name( JobCOP job ) = 0;

	virtual
	std::string filename( JobCOP ) const {
		return "unknown_file";
	}

	///////////////////////////////// evaluator interface ////////////////////////////////////////////
public:
	void add_evaluation( evaluation::PoseEvaluatorOP );

	void set_evaluators( evaluation::PoseEvaluators const& );

	void clear_evaluators();

	evaluation::PoseEvaluatorsCOP evaluators() const;

	void evaluate( core::pose::Pose &pose, std::string tag, core::io::silent::SilentStruct &pss) const;

	/// @brief call all output_observers
	void call_output_observers( core::pose::Pose const& pose, JobOP job ) const;
	void set_defaults();

private:

	evaluation::PoseEvaluatorsOP evaluators_;

	//////////////////////////////// end evaluator interface /////////////////////////////////////////

protected:
	/// @brief this function provides the extended name, not just the output name.  e.g output_name returns 1UBQ_0001, this returns 1UBQ_0001.pdb
	//not necessary in the interface class - not all derived classes CAN implement this
	//virtual
	//std::string extended_name( JobCOP job ) = 0;

	/// @brief this function generates the affixed, numbered name for the job as prefix + input + suffix + number (affixes from options system).  This function is deliberately not virtual, this should be a common mechanism; your JobOutputter can not call it if it would like.
	std::string affixed_numbered_name( JobCOP job );

private:
	/// @brief operator= defined for the sake of the remove-headers-in-headers initiative.  As of this writing, there is no reason to actually call this function, so it is declared but UNIMPLEMENTED to force compiler errors if you try to do this (instead of allowing the compiler to autogenerate the code).  If you have a valid need for this function, feel free to implement it properly and make it a public function.
	JobOutputter & operator= (JobOutputter const & rhs);

	/// @brief copy ctor defined for the sake of the remove-headers-in-headers initiative.  As of this writing, there is no reason to actually call this function, so it is declared but UNIMPLEMENTED to force compiler errors if you try to do this (instead of allowing the compiler to autogenerate the code).  If you have a valid need for this function, feel free to implement it properly and make it a public function.
	JobOutputter(JobOutputter const & cpy);

	//options
	std::string prefix_;
	std::string suffix_;
	bool no_nstruct_label_;

}; // JobOutputter

} // namespace jd2
} // namespace protocols

#endif //INCLUDED_protocols_jd2_JobOutputter_HH
