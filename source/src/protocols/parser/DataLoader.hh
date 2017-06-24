// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/parser/DataLoader.hh
/// @brief  Declartion of the XML parser's DataLoader class
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_protocols_jd2_parser_DataLoader_hh
#define INCLUDED_protocols_jd2_parser_DataLoader_hh

// Package Headers
#include <basic/datacache/DataMap.fwd.hh>

// Project Headers
#include <core/pose/Pose.fwd.hh>

// Utility Headers
#include <utility/tag/Tag.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace parser {

/// @brief A class for loading arbitrary data into the XML parser's basic::datacache::DataMap.
class DataLoader : public utility::pointer::ReferenceCount
{
public:
	DataLoader();
	virtual ~DataLoader();

	/// @brief The DataLoader may load specific information into the basic::datacache::DataMap
	/// object using this interface.
	virtual
	void load_data(
		core::pose::Pose const & pose,
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data
	) const = 0;

};

} //namespace parser
} //namespace protocols

#endif
