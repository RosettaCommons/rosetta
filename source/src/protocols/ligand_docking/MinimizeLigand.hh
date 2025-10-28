// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/task/ResfileReader.hh
/// @brief  header of classes for resfile options
/// @author Gordon Lemmon

#ifndef INCLUDED_protocols_ligand_docking_MinimizeLigand_hh
#define INCLUDED_protocols_ligand_docking_MinimizeLigand_hh

// Unit Headers

// Package Headers
#include <protocols/ligand_docking/ResidueTorsionRestraints.fwd.hh>
#include <protocols/ligand_docking/MinimizeLigand.fwd.hh>
#include <utility/vector1.hh>
//
//// Project Headers
#include <protocols/moves/Mover.hh>

//// Utility Headers
//
//// STL Headers

///////////////////////////////////////////////////////////////////////

namespace protocols {
namespace ligand_docking {

/// @brief
class MinimizeLigand : public protocols::moves::Mover
{
public:
	MinimizeLigand();
	MinimizeLigand(std::string const & chain, core::Real degrees);
	~MinimizeLigand() override;
	MinimizeLigand(MinimizeLigand const & that);

	void apply( core::pose::Pose & pose ) override;

	/// @brief Remove the constraints added by this mover to the pose in apply, if any.
	void remove_constraints( core::pose::Pose & pose );

	std::string get_name() const override;

	utility::vector1<protocols::ligand_docking::ResidueTorsionRestraintsOP>::iterator begin();
	utility::vector1<protocols::ligand_docking::ResidueTorsionRestraintsOP>::iterator end();

private:
	std::string chain_;
	core::Real degrees_; // Size of one standard deviation
	utility::vector1< protocols::ligand_docking::ResidueTorsionRestraintsOP > ligand_torsion_restraints_;
};

} //namespace ligand_docking
} //namespace protocols

#endif
