// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/sc/ElectrostaticSimilarityCalculator.hh
/// @brief  Headers for the ElectrostaticSimilarityCalculator
/// @details  The code closely follows Brian Coventry's implementation
///    of electrostatic complementarity with a added features
///    that in turn is based on the method:
///    McCoy, A. J., Epa, V. C., & Colman, P. M. (1997).
///    Electrostatic complementarity at protein/protein interfaces.
///    Journal of molecular biology, 268(2), 570-584.
/// @author   Andreas Scheck (andreas.scheck@epfl.ch)

#ifndef INCLUDED_core_scoring_sc_ElectrostaticSimilarityCalculator_hh
#define INCLUDED_core_scoring_sc_ElectrostaticSimilarityCalculator_hh

// Core Headers
#include <core/scoring/sc/ElectrostaticSimilarityCalculator.fwd.hh>
#include <core/scoring/sc/ShapeSimilarityCalculator.hh>
#include <core/id/AtomID_Map.hh>
#include <core/scoring/PoissonBoltzmannPotential.hh>
#include <core/select/residue_selector/ResidueSelector.fwd.hh>

#include <utility/vector0.hh>

namespace core {
namespace scoring {
namespace sc {

namespace ElectrostaticSimilarityDefaults {
static const bool PARTIALLY_SOLVATED = false;
static const bool IGNORE_RADIUS_RESOLUTION = 1;
}


struct ElectrostaticSimilarityResults {
	Real es_0_p;
	Real es_0_s;
	Real es_1_p;
	Real es_1_s;
	Real polarity_mag_avg0;
	Real polarity_avg0;
	Real polarity_mag_avg1;
	Real polarity_avg1;
};


class ElectrostaticSimilarityCalculator {

public:

	ElectrostaticSimilarityCalculator();
	~ElectrostaticSimilarityCalculator();

	int Init();
	void Reset();

	int Calc( core::pose::Pose const & pose );

	core::Size AddResidue( core::pose::Pose const & pose, int molecule, Size seqpos );

	ElectrostaticSimilarityResults const & GetResults() { return results_; }

public:
	void partially_solvated( bool partially_solvated ) { partially_solvated_ = partially_solvated; }
	void reference_pose( core::pose::PoseCOP reference_pose ) { reference_pose_ = reference_pose; }
	void selector1( select::residue_selector::ResidueSelectorCOP selector1 ) { selector1_ = selector1; }
	void selector2( select::residue_selector::ResidueSelectorCOP selector2 ) { selector2_ = selector2; }

private:
	id::AtomID_Map<bool>
	get_present_atoms(
		core::pose::Pose const & pose,
		id::AtomID_Map<bool> const & charged_side,
		id::AtomID_Map<bool> const & uncharged_side,
		id::AtomID_Map<bool> const & atoms_within_radius
	) const;


	void
	prepare_correlation_lists(
		std::vector<const DOT*> const & dots,
		PoissonBoltzmannPotential const & pb_a,
		PoissonBoltzmannPotential const & pb_b,
		utility::vector1<Real> & charges_from_a,
		utility::vector1<Real> & charges_from_b
	) const;


private:

	bool partially_solvated_;
	utility::vector0< id::AtomID_Map<bool> > molecule_atoms_;


	core::pose::PoseCOP reference_pose_;
	select::residue_selector::ResidueSelectorCOP selector1_;
	select::residue_selector::ResidueSelectorCOP selector2_;

	ShapeSimilarityCalculator sc_calc_;
	ElectrostaticSimilarityResults results_;
};


} // namespace sc
} // namespace scoring
} // namespace core


#endif
