// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/enzdes/EnzdesJobInputter.hh
/// @author Florian Richter (floric@u.washington.edu ), Sagar Khare (khares@u.washington.edu)

#ifndef _INCLUDED_protocols_jd2_EnzdesJobInputter_hh_
#define _INCLUDED_protocols_jd2_EnzdesJobInputter_hh_

//unit headers
#include <protocols/jd2/PDBJobInputter.hh>
#include <protocols/jd2/EnzdesJobInputter.fwd.hh>

//package headers
#include <protocols/enzdes/EnzdesLoopsFile.fwd.hh>
#include <protocols/jd2/Job.fwd.hh>

//project headers
#include <core/pose/Pose.fwd.hh>

#include <utility/vector1.hh>


//utility headers

namespace protocols {
namespace jd2 {

class EnzdesJobInputter : public protocols::jd2::PDBJobInputter
{
public:

	EnzdesJobInputter();

	virtual ~EnzdesJobInputter();

	///@brief this function is responsible for filling the pose reference with the pose indicated by the job.  The Job object (within its InnerJob) contains a PoseCOP.  This function needs to either fill the pose reference from the InnerJob or, on first demand of a pose from that InnerJob, instantiate the pose, hand off a COP to the InnerJob, and fill the reference.  This implementation uses pose_from_pdb
 	virtual void pose_from_job( core::pose::Pose & pose, JobOP job );

private:

	//eventual loops file
	protocols::enzdes::EnzdesLoopsFileOP enz_loops_file_;

}; // EnzdesJobInputter

} // namespace jd2
} // namespace protocols

#endif  // _INCLUDED_protocols_jd2_EnzdesJobInputter_hh_
