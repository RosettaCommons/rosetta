// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/conformation/symmetry/SymmDataOptions.hh
/// @brief  load the SymmData data-structure, which is used to configure symmetric poses.
/// @author Matthew O'Meara (mattjomeara@gmail.com)

#ifndef INCLUDED_core_conformation_symmetry_SymmDataOptions_hh
#define INCLUDED_core_conformation_symmetry_SymmDataOptions_hh

//unit headers
#include <core/conformation/symmetry/SymmDataOptions.fwd.hh>
#include <basic/resource_manager/ResourceOptions.hh>

//package headers
#include <basic/resource_manager/types.hh>

//utility headers
#include <utility/pointer/ReferenceCount.hh>

#ifdef WIN32
#include <utility/tag/Tag.hh>
#else
#include <utility/tag/Tag.fwd.hh>
#endif

//C++ headers
#include <istream>

namespace core {
namespace conformation {
namespace symmetry {


class SymmDataOptions : public basic::resource_manager::ResourceOptions {
public:

	SymmDataOptions() {}

	SymmDataOptions(
		std::string const & name) :
		ResourceOptions(name)
	{}

	virtual ~SymmDataOptions() {}

	virtual
	void
	parse_my_tag(utility::tag::TagCOP) {}

	std::string
	type() const { return "SymmDataOptions"; }
};

} // namespace
} // namespace
} // namespace

#endif // include guard
