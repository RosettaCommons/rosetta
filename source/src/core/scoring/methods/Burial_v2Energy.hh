// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/Burial_v2Energy.hh
/// @brief  energy term use for score burial of a specific residue
/// @author TJ Brunette

#ifndef INCLUDED_core_scoring_methods_Burial_v2Energy_hh
#define INCLUDED_core_scoring_methods_Burial_v2Energy_hh

#include <core/scoring/methods/Burial_v2EnergyCreator.hh>
#include <core/scoring/methods/WholeStructureEnergy.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

#include <core/pose/Pose.fwd.hh>

#include <utility/vector1.hh>


namespace core {
namespace scoring {
namespace methods {


class Burial_v2Energy : public WholeStructureEnergy {
public:
	typedef WholeStructureEnergy parent;

	Burial_v2Energy() : WholeStructureEnergy( methods::EnergyMethodCreatorOP( new Burial_v2EnergyCreator ) ) {
		init_from_file();
	}

	/// clone
	virtual
	EnergyMethodOP
	clone() const;

	/////////////////////////////////////////////////////////////////////////////
	// scoring
	/////////////////////////////////////////////////////////////////////////////

	virtual void finalize_total_energy(
		pose::Pose & pose,
		ScoreFunction const &,
		EnergyMap & totals
	) const;

	virtual void indicate_required_context_graphs( utility::vector1< bool > & ) const {};

	virtual
	core::Size version() const;

private:
	core::Real using_atom_distance(core::pose::Pose const & pose) const;
	core::Real using_totalSasa(core::pose::Pose const & pose) const;
	void init_from_file();
	utility::vector1<core::Size> residue_ids_;
};

}
}
}

#endif // INCLUDED_core_scoring_methods_Burial_v2Energy_HH
