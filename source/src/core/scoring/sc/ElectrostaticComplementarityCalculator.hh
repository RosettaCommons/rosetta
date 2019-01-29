// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/sc/ElectrostaticComplementarityCalculator.hh
/// @brief  Headers for the Electrostatic Complementarity Calculator
/// @author Brian Coventry (bcov@uw.edu)
/// @details Based on the method from
///          McCoy, A. J., Epa, V. C., & Colman, P. M. (1997).
///              Electrostatic complementarity at protein/protein interfaces1.
///              Journal of molecular biology, 268(2), 570-584.

#ifndef INCLUDED_core_scoring_sc_ElectrostaticComplementarityCalculator_hh
#define INCLUDED_core_scoring_sc_ElectrostaticComplementarityCalculator_hh


#include <core/scoring/sc/ElectrostaticComplementarityCalculator.fwd.hh>

#include <core/scoring/sc/ShapeComplementarityCalculator.hh>

#include <core/id/AtomID_Map.hh>
#include <core/scoring/PoissonBoltzmannPotential.hh>


namespace core {
namespace scoring {
namespace sc {

namespace ElectrostaticComplementarityDefaults {
static const Size IGNORE_RADIUS = 20;
static const Real INTERFACE_TRIM_RADIUS = 1.5;
static const bool PARTIALLY_SOLVATED = true;
static const bool IGNORE_RADIUS_RESOLUTION = 1;
}


struct ElectrostaticComplementarityResults {
	Real ec_0_p;
	Real ec_0_s;
	Real ec_1_p;
	Real ec_1_s;
	Real polarity_mag_avg0;
	Real polarity_avg0;
	Real polarity_mag_avg1;
	Real polarity_avg1;
};


class ElectrostaticComplementarityCalculator {

public:

	ElectrostaticComplementarityCalculator();
	~ElectrostaticComplementarityCalculator();

	int Init( core::pose::Pose const & pose );
	void Reset();

	int Calc( core::pose::Pose const & pose );
	int Calc( core::pose::Pose const & pose, core::Size jump_id );

	core::Size AddResidue( core::pose::Pose const & pose, int molecule, Size seqpos );

	ElectrostaticComplementarityResults const & GetResults() { return results_; }

public:
	void ignore_radius( Real ignore_radius ) { ignore_radius_ = ignore_radius; }
	void interface_trim_radius( Real interface_trim_radius ) { interface_trim_radius_ = interface_trim_radius; }
	void partially_solvated( bool partially_solvated ) { partially_solvated_ = partially_solvated; }

private:
	id::AtomID_Map<bool>
	get_present_atoms(
		core::pose::Pose const & pose,
		id::AtomID_Map<bool> const & charged_side,
		id::AtomID_Map<bool> const & uncharged_side,
		id::AtomID_Map<bool> const & atoms_within_radius
	) const;


	bool
	is_within_radius(
		numeric::xyzVector< Real > xyz,
		std::vector<const DOT*> const & dots
	) const;

	void
	get_atoms_within_radius(
		core::pose::Pose const & pose,
		id::AtomID_Map<bool> & atoms_within_radius,
		std::vector<const DOT*> const & dots,
		id::AtomID_Map<bool> const & molecule1,
		id::AtomID_Map<bool> const & molecule2
	) const;

	std::vector<const DOT*>
	trim_dots_to_resl( std::vector<const DOT*> const & dots ) const;

	void
	prepare_correlation_lists(
		std::vector<const DOT*> const & dots,
		PoissonBoltzmannPotential const & pb_a,
		PoissonBoltzmannPotential const & pb_b,
		utility::vector1<Real> & charges_from_a,
		utility::vector1<Real> & charges_from_b
	) const;


private:

	Real ignore_radius_;
	Real ignore_radius_resolution_;
	Real interface_trim_radius_;
	bool partially_solvated_;
	utility::vector0< id::AtomID_Map<bool> > molecule_atoms_;


	ShapeComplementarityCalculator sc_calc_;
	ElectrostaticComplementarityResults results_;
};


} // namespace sc
} // namespace scoring
} // namespace core


#endif
