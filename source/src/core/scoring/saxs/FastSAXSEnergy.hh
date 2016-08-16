// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/saxs/SAXSEnergy.hh
/// @brief  FastSAX scoring of Stovgaard et al (BMC Bioinf. 2010)
/// @author Frank DiMaio


#ifndef INCLUDED_core_scoring_saxs_FastSAXSEnergy_hh
#define INCLUDED_core_scoring_saxs_FastSAXSEnergy_hh

// Package headers
#include <core/chemical/AA.hh>

#include <core/scoring/methods/WholeStructureEnergy.hh>

#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>

#include <string>

#include <utility/vector1.hh>
#include <numeric/xyzVector.hh>


namespace core {
namespace scoring {
namespace saxs {

class FastSAXSEnergy : public methods::WholeStructureEnergy  {
public:
	typedef methods::WholeStructureEnergy  parent;

public:


	FastSAXSEnergy();


	FastSAXSEnergy( FastSAXSEnergy const & src ) :
		parent(src),
		dchi2_dca(src.dchi2_dca),
		dchi2_dsc(src.dchi2_dsc),
		chi2(src.chi2),
		c(src.c),
		r_chi2(src.r_chi2)
	{}

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
	// PRIVATE DATA
	//////
	// precomputed bb/sc derivatives
	mutable utility::vector1< numeric::xyzVector< core::Real > > dchi2_dca;
	mutable utility::vector1< numeric::xyzVector< core::Real > > dchi2_dsc;

	// saved chi2, c scaling factor
	mutable core::Real chi2, c;

	// saved chi2 residuals
	mutable utility::vector1< core::Real > r_chi2;

	virtual
	core::Size version() const;
};


// a couple utility functions
// (0) remap core:chemical::AA enum to indices in formfactor table
core::Size aa2idx( core::chemical::AA aa );

// (1) load saxs spectrum from file (only do this once)
void load_fastsax_spectrum(
	core::Size &nq,
	utility::vector1< core::Real >::iterator &q,
	utility::vector1< core::Real >::iterator &i_obs,
	utility::vector1< core::Real >::iterator &i_sig);

// (2) load per-resiude form factors (resampled onto the same grid as the input spectrum)
void load_form_factors(
	core::Size nq,
	utility::vector1< core::Real >::iterator &q_in,
	utility::vector1< utility::vector1< core::Real > >::iterator &spectra);

}
}
}

#endif
