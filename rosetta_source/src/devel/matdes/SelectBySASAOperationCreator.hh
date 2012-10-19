// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   devel/matdes/SelectBySASAOperation.hh
/// @brief  Restrict design to residues matching user-specified SASA criteria in the monomeric, bound, or unbound state.
/// @author Jacob Bale (balej@uw.edu)

#ifndef INCLUDED_devel_matdes_SelectBySASAOperationCreator_hh
#define INCLUDED_devel_matdes_SelectBySASAOperationCreator_hh


#include <core/pack/task/operation/TaskOperationCreator.hh>

#include <string>

namespace devel {
namespace matdes {

class SelectBySASAOperationCreator : public core::pack::task::operation::TaskOperationCreator {
public:
	virtual core::pack::task::operation::TaskOperationOP create_task_operation() const;
	virtual std::string keyname() const { return "SelectBySASA"; }
};

} //namespace matdes
} //namespace devel
#endif

