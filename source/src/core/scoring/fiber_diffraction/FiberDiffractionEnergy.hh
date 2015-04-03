// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/fiber_diffraction/FiberDiffractionEnergy.hh
/// @brief  FiberDiffraction scoring and derivatives (all-atom) 
/// @author Wojciech Potrzebowski and Ingemar Andre


#ifndef INCLUDED_core_scoring_fiber_diffraction_FiberDiffractionEnergy_hh
#define INCLUDED_core_scoring_fiber_diffraction_FiberDiffractionEnergy_hh

#include <core/scoring/fiber_diffraction/FiberDiffractionEnergyCreator.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/methods/WholeStructureEnergy.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/AA.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/id/AtomID.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <numeric/xyzVector.hh>


#include <string>

namespace core {
namespace scoring {
namespace fiber_diffraction {

class FiberDiffractionEnergy : public methods::WholeStructureEnergy  {
public:
	typedef methods::WholeStructureEnergy  parent;

public:


	FiberDiffractionEnergy();


	FiberDiffractionEnergy( FiberDiffractionEnergy const & src ) :
		parent(src) {
		chi2_ = src.chi2_;
		scale_factor_ = src.scale_factor_;
		dchi2_d = src.dchi2_d;
		dchi2_d_cross_R = src.dchi2_d_cross_R;
	}

	virtual methods::EnergyMethodOP clone() const;

	virtual void finalize_total_energy(pose::Pose & pose,ScoreFunction const &,EnergyMap & totals) const;

	virtual void setup_for_scoring( pose::Pose & pose, ScoreFunction const & scorefxn ) const;

	virtual void setup_for_derivatives( pose::Pose & pose, ScoreFunction const & sf) const;

	virtual void
	eval_atom_derivative(
		id::AtomID const & id,
		pose::Pose const & pose,
		kinematics::DomainMap const &, // domain_map,
		ScoreFunction const & sfxn,
		EnergyMap const & weights,
		Vector & F1,
		Vector & F2
	) const;

	virtual void indicate_required_context_graphs(utility::vector1< bool > & /*context_graphs_required*/
    ) const {}

private:
	//////
	// PRIVATE DATA becuase the functions in this class are const we need to make variables mutable
	// store calculated intensities
	mutable utility::vector0 < utility::vector1 < core::Real > > I;
	// precomputed derivatives
	mutable utility::vector1< numeric::xyzVector< core::Real > > dchi2_d; 
	mutable utility::vector1< numeric::xyzVector< core::Real > > dchi2_d_cross_R;
	mutable std::map<  core::id::AtomID, core::Size > AtomID_to_atomnbr_;
	// saved chi2, c scaling factor
	mutable core::Real chi2_;
	mutable core::Real scale_factor_;
	mutable core::Real square_obs_, sum_obs_;
	mutable core::Real Rscale_factor;
	mutable core::Real Rsum_obs;
	mutable core::Real p_, c_;
	mutable core::Real a_, b_;
	virtual
	core::Size version() const;
};

}
}
}

#endif
