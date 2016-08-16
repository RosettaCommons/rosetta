// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/conformation/symmetry/SymmDataOptionsCreator.hh
/// @brief  load the SymmData data-structure, which is used to configure symmetric poses.
/// @author Matthew O'Meara (mattjomeara@gmail.com)

#ifndef INCLUDED_core_conformation_symmetry_SymmDataOptionsCreator_hh
#define INCLUDED_core_conformation_symmetry_SymmDataOptionsCreator_hh

//unit headers
#include <basic/resource_manager/ResourceOptionsCreator.hh>
#include <core/conformation/symmetry/SymmDataOptions.fwd.hh>

namespace core {
namespace conformation {
namespace symmetry {

class SymmDataOptionsCreator : public basic::resource_manager::ResourceOptionsCreator
{
public:
	SymmDataOptionsCreator();
	~SymmDataOptionsCreator();

	virtual
	basic::resource_manager::ResourceOptionsOP
	create_options() const;

	virtual
	std::string options_type() const;

};

} // namespace
} // namespace
} // namespace

#endif // include guard
