// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief
/// @author Zhe Zhang (zhezhang1986@gmail.com)


#ifndef INCLUDED_devel_replica_docking_FnonnatFilterCreator_hh
#define INCLUDED_devel_replica_docking_FnonnatFilterCreator_hh

// Package Headers
#include <protocols/filters/FilterCreator.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>

// c++ headers
#include <string>

namespace devel {
namespace replica_docking {

class FnonnatFilterCreator : public protocols::filters::FilterCreator
{
public:
	virtual protocols::filters::FilterOP create_filter() const;
	virtual std::string keyname() const;
};

}
}

#endif

