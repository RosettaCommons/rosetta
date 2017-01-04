// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/jd3/pose_inputters/SilentFilePoseInputterCreator.hh
/// @brief  Creator class for the SilentFile Pose Inputter
/// @author Andy Watkins (amw579@stanford.edu)

#ifndef INCLUDED_protocols_jd3_pose_inputters_SilentFilePoseInputterCreator_hh
#define INCLUDED_protocols_jd3_pose_inputters_SilentFilePoseInputterCreator_hh

#include <protocols/jd3/pose_inputters/PoseInputterCreator.hh>

namespace protocols {
namespace jd3 {
namespace pose_inputters {

class SilentFilePoseInputterCreator : public PoseInputterCreator {
public:
	PoseInputterOP create_inputter() const override;
	std::string keyname() const override;
	void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const override;
	void list_options_read( utility::options::OptionKeyList & read_options ) const override;
};

} //pose_inputters
} //jd3
} //protocols

#endif //INCLUDED_protocols_jd3_pose_inputters_SilentFilePoseInputterCreator_hh

