// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/core/energy_methods/ArgCationPiEnergy.hh
/// @brief  Cation pi term that specializes in bringing Arginine and rings together
/// @details Currently designed for canonical amino acids but easily extended beyond.
/// @author Brian Coventry (bcov@uw.edu)


#ifndef INCLUDED_core_energy_methods_ArgCationPiEnergy_hh
#define INCLUDED_core_energy_methods_ArgCationPiEnergy_hh

// Unit headers
#include <core/energy_methods/ArgCationPiEnergy.fwd.hh>

// Package headers
#include <core/scoring/methods/ContextIndependentTwoBodyEnergy.hh>
#include <core/scoring/methods/EnergyMethodCreator.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>

#include <core/chemical/ResidueType.hh>

// Utility headers
#include <utility/vector1.hh>
#include <utility/pointer/memory.hh>


namespace core {
namespace scoring {
namespace methods {


class ArgCationPiEnergy : public scoring::methods::ContextIndependentTwoBodyEnergy {

public:
	ArgCationPiEnergy( core::scoring::methods::EnergyMethodOptions const & options );

	scoring::methods::EnergyMethodOP
	clone() const override;

	void
	residue_pair_energy(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		pose::Pose const & pose,
		scoring::ScoreFunction const & sfxn,
		scoring::EnergyMap & emap
	) const override;

	Distance
	atomic_interaction_cutoff() const override;


	bool
	defines_score_for_residue_pair(
		conformation::Residue const & res1,
		conformation::Residue const & res2,
		bool res_moving_wrt_eachother
	) const override;


	void
	eval_residue_pair_derivatives(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		scoring::ResSingleMinimizationData const &,
		scoring::ResSingleMinimizationData const &,
		scoring::ResPairMinimizationData const & min_data,
		pose::Pose const & pose, // provides context
		scoring::EnergyMap const & weights,
		utility::vector1< scoring::DerivVectorPair > & r1_atom_derivs,
		utility::vector1< scoring::DerivVectorPair > & r2_atom_derivs
	) const override;


	bool
	defines_intrares_energy( scoring::EnergyMap const & weights ) const override;

	void
	eval_intrares_energy(
		conformation::Residue const & rsd,
		pose::Pose const & pose,
		scoring::ScoreFunction const & sfxn,
		scoring::EnergyMap & emap
	) const override;


	core::Size version() const override;

	// We don't want any of these
	void indicate_required_context_graphs(
		utility::vector1< bool > &
	) const override {}


	bool
	is_arg(
		core::chemical::ResidueType const & restype
	) const;

	bool
	is_phe(
		core::chemical::ResidueType const & restype
	) const;

	bool
	is_his(
		core::chemical::ResidueType const & restype
	) const;

	bool
	is_trp(
		core::chemical::ResidueType const & restype
	) const;

	bool
	is_tyr(
		core::chemical::ResidueType const & restype
	) const;

	bool
	is_pi(
		core::chemical::ResidueType const & restype
	) const;

	bool
	valid_res_pair(
		core::chemical::ResidueType const & restype1,
		core::chemical::ResidueType const & restype2
	) const;

	std::pair< std::pair< numeric::xyzVector<Real>, numeric::xyzVector<Real> >, utility::vector1<std::string> >
	get_ring(
		conformation::Residue const & pi,
		bool trp_little_ring
	) const;

	std::pair< std::pair< numeric::xyzVector<Real>, numeric::xyzVector<Real> >, utility::vector1<std::string> >
	get_ring_params(
		conformation::Residue const & pi,
		conformation::Residue const & arg
	) const;


private:

	bool his_can_be_pi_;


};

} // methods
} // scoring
} // core

#endif
