// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/jd3/pose_outputters/SecondaryPoseOutputterCreator.hh
/// @brief  Declaration of the %SecondaryPoseOutputterCreator class
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)


#ifndef INCLUDED_protocols_jd3_pose_outputters_SecondaryPoseOutputterCreator_HH
#define INCLUDED_protocols_jd3_pose_outputters_SecondaryPoseOutputterCreator_HH

//unit headers
#include <protocols/jd3/pose_outputters/SecondaryPoseOutputterCreator.fwd.hh>

// Package headers
#include <protocols/jd3/pose_outputters/SecondaryPoseOutputter.fwd.hh>

// utility headers
#include <utility/options/keys/OptionKey.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>

// C++ headers
#include <string>

namespace protocols {
namespace jd3 {
namespace pose_outputters {

class SecondaryPoseOutputterCreator : public utility::pointer::ReferenceCount
{
public:
	virtual SecondaryPoseOutputterOP create_outputter() const = 0;
	virtual std::string keyname() const = 0;
	virtual void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const = 0;
	virtual void list_options_read( utility::options::OptionKeyList & read_options ) const = 0;
	virtual bool outputter_specified_by_command_line() const = 0;

};

} // namespace pose_outputters
} // namespace jd3
} // namespace protocols

#endif //INCLUDED_protocols_jd3_pose_outputters_SecondaryPoseOutputterCreator_HH
