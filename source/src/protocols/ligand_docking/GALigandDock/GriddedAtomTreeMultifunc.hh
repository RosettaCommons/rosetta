// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/ligand_docking/GALigandDock/GridScorer.hh
///
/// @brief  Grid-based scoring for GA ligand docking
///         **NOTE** Currently only wraps scorefunction with some caching, not yet grid-based!
/// @author Hahnbeom Park and Frank DiMaio

#ifndef INCLUDED_protocols_ligand_docking_GALigandDock_GriddedGriddedAtomTreeMultifunc_hh
#define INCLUDED_protocols_ligand_docking_GALigandDock_GriddedGriddedAtomTreeMultifunc_hh

#include <protocols/ligand_docking/GALigandDock/GridScorer.hh>

#include <protocols/ligand_docking/GALigandDock/LigandConformer.hh>
#include <protocols/ligand_docking/GALigandDock/GridHash3D.hh>
#include <core/optimization/Multifunc.hh>
#include <core/optimization/MinimizerMap.hh>

namespace protocols {
namespace ligand_docking {
namespace ga_ligand_dock {

/// @brief Atom tree multifunction class
class GriddedAtomTreeMultifunc : public core::optimization::Multifunc {
public:

	GriddedAtomTreeMultifunc(
		LigandConformer & conf_in,
		core::pose::Pose & pose_in,
		GridScorer & scorefxn_in,
		core::optimization::MinimizerMap & min_map_in
	);

	~GriddedAtomTreeMultifunc() override;

public: // Methods
	core::Real
	operator ()( core::optimization::Multivec const & vars ) const override;

	void
	dfunc( core::optimization::Multivec const & vars, core::optimization::Multivec & dE_dvars ) const override;

	void
	dump( core::optimization::Multivec const & vars, core::optimization::Multivec const & vars2 ) const override;

private: // data

	LigandConformer & conf_;
	core::pose::Pose & pose_;
	GridScorer & sf_;
	core::optimization::MinimizerMap & min_map_;
}; // GriddedAtomTreeMultifunc


} // ga_ligand_dock
} // ligand_docking
} // protocols

#endif
