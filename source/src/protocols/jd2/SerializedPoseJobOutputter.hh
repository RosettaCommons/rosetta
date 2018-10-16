// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/jd2/SerializedPoseJobOutputter.hh
/// @brief  outputs serialized Poses that should only be read in by applications from the same git checkout on the same computer
/// @author Jack Maguire, jackmaguire1444@gmail.com

#ifndef INCLUDED_protocols_jd2_SerializedPoseJobOutputter_hh
#define INCLUDED_protocols_jd2_SerializedPoseJobOutputter_hh

//unit headers
#include <protocols/jd2/SerializedPoseJobOutputter.fwd.hh>
#include <protocols/jd2/wwPDBJobOutputter.hh>
#include <protocols/jd2/Job.fwd.hh>

//project headers
#include <core/pose/Pose.fwd.hh>

//C++ headers
#include <string>

#include <utility/vector1.hh>
#include <utility/io/ozstream.fwd.hh>
#include <iostream>


namespace protocols {
namespace jd2 {

class SerializedPoseJobOutputter : public protocols::jd2::wwPDBJobOutputter
{
public:

	SerializedPoseJobOutputter();

	~SerializedPoseJobOutputter() override;

protected:
	void dump_pose( JobCOP job, core::pose::Pose const & pose, utility::io::ozstream & out, std::string const &filename="" ) override;

}; // SerializedPoseJobOutputter

} // namespace jd2
} // namespace protocols

#endif //INCLUDED_protocols_jd2_SerializedPoseJobOutputter_HH
