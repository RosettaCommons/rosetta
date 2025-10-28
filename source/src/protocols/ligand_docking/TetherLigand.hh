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

#ifndef INCLUDED_protocols_ligand_docking_TetherLigand_hh
#define INCLUDED_protocols_ligand_docking_TetherLigand_hh

// Unit Headers

#include <core/scoring/constraints/Constraint.fwd.hh>

// Package Headers
#include <protocols/moves/Mover.hh>

//// Scripter Headers



///////////////////////////////////////////////////////////////////////

namespace protocols {
namespace ligand_docking {

/// @brief
class TetherLigand : public protocols::moves::Mover
{
public:
	TetherLigand();
	TetherLigand(std::string const & chain, const core::Real & angstroms);
	~TetherLigand() override;
	TetherLigand(TetherLigand const & that);

	void apply( core::pose::Pose & pose ) override;

	std::string get_name() const override;

	void release(core::pose::Pose & pose);

	core::scoring::constraints::ConstraintCOP const &
	get_ligand_tether() const;

private:
	std::string chain_;
	core::Real angstroms_; //size of one stdev for ligand restraint
	core::scoring::constraints::ConstraintCOP ligand_tether_;
};

core::scoring::constraints::ConstraintCOP
restrain_ligand_nbr_atom(
	core::Size const lig_id,
	core::Real const stddev_Angstroms,
	core::pose::Pose & pose
);

} //namespace ligand_docking
} //namespace protocols

#endif
