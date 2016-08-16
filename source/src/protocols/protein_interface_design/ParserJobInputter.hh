// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/jd2/ParserJobInputter.hh
/// @brief  header file for PatchdockJobInputter. Feeds poses based on a starting structure and a patchdock output file describing the rigid-body transformations to apply to the pose
/// @author Sarel Fleishman sarelf@u.washington.edu


#ifndef INCLUDED_protocols_protein_interface_design_ParserJobInputter_hh
#define INCLUDED_protocols_protein_interface_design_ParserJobInputter_hh

//unit headers
#include <protocols/jd2/PDBJobInputter.hh>
#include <protocols/protein_interface_design/ParserJobInputter.fwd.hh>
#include <protocols/jd2/Job.fwd.hh>

//project headers
#include <core/pose/Pose.fwd.hh>

#include <utility/vector1.hh>


//utility headers

namespace protocols {
namespace protein_interface_design {

class ParserJobInputter : public protocols::jd2::PDBJobInputter
{
public:

	ParserJobInputter();

	virtual ~ParserJobInputter();

	/// @brief this function is responsible for filling the pose reference with the pose indicated by the job.  The Job object (within its InnerJob) contains a PoseCOP.  This function needs to either fill the pose reference from the InnerJob or, on first demand of a pose from that InnerJob, instantiate the pose, hand off a COP to the InnerJob, and fill the reference.  This implementation uses pose_from_pdb
	virtual void pose_from_job( core::pose::Pose & pose, protocols::jd2::JobOP job );

}; // ParserJobInputter

} // namespace protein_interface_design
} // namespace protocols

#endif //INCLUDED_protocols_protein_interface_design_ParserJobInputter_HH
