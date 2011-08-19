// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is part of the Rosetta software suite and is made available under license.
// The Rosetta software is developed by the contributing members of the Rosetta Commons consortium.
// (C) 199x-2009 Rosetta Commons participating institutions and developers.
// For more information, see http://www.rosettacommons.org/.

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

// Symmetry
// AUTO-REMOVED #include <core/scoring/symmetry/SymmetricScoreFunction.fwd.hh>

namespace protocols {
namespace symmetric_docking {

class SymRestrictTaskForDocking : public core::pack::task::operation::TaskOperation
{
public:
	SymRestrictTaskForDocking();

	SymRestrictTaskForDocking( core::scoring::ScoreFunctionCOP scorefxn, bool include_current, core::Real distance = 8  );

	virtual ~SymRestrictTaskForDocking();

	virtual core::pack::task::operation::TaskOperationOP clone() const;

	virtual	void apply( core::pose::Pose const & pose, core::pack::task::PackerTask & task ) const;

private:
	core::scoring::ScoreFunctionCOP scorefxn_;
	bool include_current_;
	core::Real distance_;
};

} // namespace docking
} // namespace protocols

#endif
