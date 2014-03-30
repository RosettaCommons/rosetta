//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

///@file core/scoring/Interface_/DDPscore.hh
///@brief Implementation of distance dependent interface score
///@author Hermann Zellner (hermann1.zellner@biologie.uni-regensburg.de)

#ifndef INCLUDED_core_scoring_interface_DDPscore_hh
#define INCLUDED_core_scoring_interface_DDPscore_hh

#include <core/scoring/methods/ContextDependentTwoBodyEnergy.hh>
#include <core/scoring/interface_/DDPscore.fwd.hh>
#include <core/scoring/interface_/DDPlookup.hh>
#include <core/pose/Pose.fwd.hh>

#include <core/scoring/ScoreFunction.fwd.hh>

#include <utility/vector1.hh>


namespace core{
namespace scoring{
namespace interface_{

class DDPscore : public methods::ContextDependentTwoBodyEnergy {

public:
	typedef methods::ContextDependentTwoBodyEnergy  parent;

public:
	DDPscore();

	virtual methods::EnergyMethodOP clone() const;

	virtual void setup_for_scoring( pose::Pose & pose, ScoreFunction const & ) const;

	virtual void setup_for_packing( pose::Pose & pose, utility::vector1< bool > const &, utility::vector1< bool > const & ) const;

	virtual void setup_for_derivatives( pose::Pose &pose, ScoreFunction const &  ) const;

	//virtual void setup_for_minimizing(pose::Pose & pose, ScoreFunction const & ,optimization::MinimizerMap const &) const;

	virtual void residue_pair_energy(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		const pose::Pose & pose,
		const ScoreFunction & sfxn,
		EnergyMap & emap
	) const;

	bool defines_intrares_energy(core::scoring::EnergyMap const &) const;

	virtual void eval_intrares_energy(
	                        const core::conformation::Residue &,
	                        const core::pose::Pose &,
	                        const core::scoring::ScoreFunction &,
	                        core::scoring::EnergyMap &
	                        ) const;

	virtual void indicate_required_context_graphs(
			utility::vector1< bool > & context_graphs_required
			) const;

	core::Distance atomic_interaction_cutoff() const;

private:
	DDPlookup const lookup_table_;
virtual
core::Size version() const;
}; // DDP
} // Interface_
} // scoring
} // core


#endif /* INCLUDED_core_scoring_Interface_DDPscore_HH_ */
