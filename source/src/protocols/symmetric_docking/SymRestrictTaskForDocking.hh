// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file SymRestrictTaskForDocking.hh
/// @brief When passed to a PackerTask, pack/design is limited to the interface
/// @author Ingemar Andre

#ifndef INCLUDED_protocols_symmetric_docking_SymRestrictTaskForDocking_hh
#define INCLUDED_protocols_symmetric_docking_SymRestrictTaskForDocking_hh

#include <protocols/symmetric_docking/SymRestrictTaskForDocking.fwd.hh>
#include <core/pack/task/operation/TaskOperation.hh>

#include <core/types.hh>
#include <core/pack/task/PackerTask.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

#include <utility/vector1.hh>


// Symmetry

namespace protocols {
namespace symmetric_docking {

class SymRestrictTaskForDocking : public core::pack::task::operation::TaskOperation
{
public:
	SymRestrictTaskForDocking();

	SymRestrictTaskForDocking( core::scoring::ScoreFunctionCOP scorefxn, bool include_current, core::Real distance = 8  );

	~SymRestrictTaskForDocking() override;

	core::pack::task::operation::TaskOperationOP clone() const override;

	void apply( core::pose::Pose const & pose, core::pack::task::PackerTask & task ) const override;

private:
	core::scoring::ScoreFunctionCOP scorefxn_;
	bool include_current_;
	core::Real distance_;
};

} // namespace docking
} // namespace protocols

#endif
