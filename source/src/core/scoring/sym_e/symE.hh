// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/symE/symE.hh
/// @brief  Header declarations for class symE
/// @author Sam DeLuca (samuel.l.deluca@vanderbilt.edu)

#ifndef INCLUDED_core_scoring_sym_e_symE_hh
#define INCLUDED_core_scoring_sym_e_symE_hh
#include <core/scoring/sym_e/symE.fwd.hh>
#include <core/scoring/methods/ContextDependentLRTwoBodyEnergy.hh>
#include <core/pose/Pose.fwd.hh>

#include <utility/vector1.hh>


namespace core {
namespace scoring {
namespace sym_e {

class symEnergy : public methods::ContextDependentLRTwoBodyEnergy {
public:
	typedef methods::ContextDependentLRTwoBodyEnergy parent;

public:
	//symE( methods::EnergyMethodOptions const & options );

	symEnergy();

	virtual methods::EnergyMethodOP clone() const;

	virtual void setup_for_scoring( pose::Pose & pose, ScoreFunction const & ) const;

	virtual void setup_for_packing(pose::Pose & pose, utility::vector1< bool > const &, utility::vector1< bool > const & ) const;

	virtual void indicate_required_context_graphs( utility::vector1< bool > & context_graphs_required ) const;

	virtual bool defines_intrares_energy( EnergyMap const & weights ) const;

	virtual void eval_intrares_energy(conformation::Residue const & rsd, pose::Pose const & pose, ScoreFunction const & sfxn, EnergyMap & emap) const;

	virtual methods::LongRangeEnergyType long_range_type() const;

	virtual bool defines_residue_pair_energy(const core::pose::Pose& pose, platform::Size res1, platform::Size res2) const;

	virtual void residue_pair_energy(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		pose::Pose const & pose,
		ScoreFunction const &,
		EnergyMap & emap
	) const;
	virtual
	core::Size version() const;
};


}  //symE
}  //scoring
}  //core

#endif /*symE_H_ */

