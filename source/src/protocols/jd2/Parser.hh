// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/jd2/Parser.hh
/// @brief  header file for Parser class, part of August 2008 job distributor
/// @author Steven Lewis smlewi@gmail.com


#ifndef INCLUDED_protocols_jd2_Parser_hh
#define INCLUDED_protocols_jd2_Parser_hh

//unit headers
#include <protocols/moves/Mover.fwd.hh>
#include <protocols/jd2/Job.fwd.hh>
#include <protocols/jd2/Parser.fwd.hh>

//project headers
#include <core/pose/Pose.fwd.hh>

//utility headers
#include <utility/pointer/ReferenceCount.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace jd2 {

/// @details the Parser class can create a protocol (made of Movers) from an XML
/// file.  This interface class describes its functionality for the August '08
/// job distributor.

class Parser : public utility::pointer::ReferenceCount
{
public:
	typedef core::pose::Pose Pose;
	typedef protocols::moves::Mover Mover;
	typedef protocols::moves::MoverOP MoverOP;
public:

	~Parser() override;

	/// @brief generate_mover_from_job is the function called by the job
	/// distributor to request a mover.  It is defined in the base class (this
	/// class) and handles unpackaging the job and passes the work to
	/// generate_mover_from_pose.  The pose is repackaged into the JobOP so that
	/// jobs starting off that pose are properly modified.
	///
	/// @param [in]     job: this is the job we're working on (contains the input pose)
	/// @param [in]     pose: The pose that will be used in this job.
	/// @param [in/out] mover: this is a mover pointer reference; the function can choose
	///                 to overwrite the input with a new mover
	/// @param [in]     new_input: true if this is different input (a different pose)
	///                 from the last job distributor cycle
	/// @param [in]     allow_job_update: is true if the Parser is allowed to rewrite
	///                 the Pose inside the input Job after generate_mover_from_pose has
	///                 modified it.
	/// @param [in]     guarantee_new_mover: true if the input mover should be replaced
	///                 reguardless of whether or not new_input is true.
	void
	generate_mover_from_job(
		JobOP job,
		core::pose::Pose & pose,
		MoverOP & mover,
		bool new_input,
		bool allow_job_update = true,
		bool guarantee_new_mover = false
	);

	/// @brief generate_mover_from_pose is overloaded by derived classes to parse
	/// as desired.  The pose is passed by nonconst reference - the function is
	/// thus allowed to modify this pose (preferably not at all, but by adding
	/// constraints if necessary).  This pose will stay modified for jobs starting
	/// off that pose.  The function is expected to use its return bool to signal
	/// whether the pose has changed or not (and thus whether it needs to be
	/// repackaged back into the InnerJob).  This function should return
	/// IMMEDIATELY with false if it chooses not to modify the mover or pose.
	///
	/// @param [in]     JobOP job this is the job we're working on (contains the
	///                 input pose's name)
	/// @param [in/out] pose this is the "starting pose" for a series of jobs; the
	///                 function is allowed to modify it by adding constraints if
	///                 necessary
	/// @param [in/out]  mover this is a mover; the function can choose to overwrite
	///                 the input with a new mover
	/// @param [in]      bool new_input true if this is different input (a different
	///                 pose) from the last job distributor cycle
	/// @param [out]     bool the return value states whether or not the pose has
	///                 been changed - it will be repackaged into the InnerJob only
	///                 as necessary
	virtual
	bool
	generate_mover_from_pose(
		JobCOP job,
		Pose & pose,
		MoverOP & mover,
		bool new_input,
		std::string const & xml_file,
		bool guarantee_new_mover = false
	) = 0;

}; // Parser

} // namespace jd2
} // namespace protocols

#endif //INCLUDED_protocols_jd2_Parser_HH
