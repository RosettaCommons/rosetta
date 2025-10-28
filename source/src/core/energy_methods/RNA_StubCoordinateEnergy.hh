// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/rna/RNA_StubCoordinateEnergy.hh
/// @brief  Statistically derived rotamer pair potential class declaration
/// @author Rhiju Das

#ifndef INCLUDED_core_scoring_rna_RNA_StubCoordinateEnergy_HH
#define INCLUDED_core_scoring_rna_RNA_StubCoordinateEnergy_HH

// Unit Headers
#include <core/energy_methods/RNA_StubCoordinateEnergy.fwd.hh>

// Package headers
#include <core/scoring/methods/ContextIndependentLRTwoBodyEnergy.hh>
#include <core/scoring/func/Func.fwd.hh>

// Project headers
#include <core/kinematics/RT.fwd.hh>
#include <core/pose/rna/StubStubType.fwd.hh>

#include <numeric/xyzVector.hh>
#include <tuple>

namespace core {
namespace energy_methods {


class RNA_StubCoordinateEnergy : public core::scoring::methods::ContextIndependentLRTwoBodyEnergy {

public:
	typedef core::scoring::methods::ContextIndependentLRTwoBodyEnergy parent;
	typedef core::scoring::methods::EnergyMethodOP EnergyMethodOP;
	typedef numeric::xyzVector< Real > Vector;
	typedef numeric::xyzMatrix< Real > Matrix;

public:

	RNA_StubCoordinateEnergy();

	/// clone
	core::scoring::methods::EnergyMethodOP
	clone() const override;

	///
	void
	setup_for_scoring( pose::Pose & pose, core::scoring::ScoreFunction const & ) const override;

	void
	residue_pair_energy(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		pose::Pose const & pose,
		core::scoring::ScoreFunction const & sfxn,
		core::scoring::EnergyMap & emap
	) const override;

	void indicate_required_context_graphs( utility::vector1< bool > & ) const override {}

	bool
	defines_intrares_energy( core::scoring::EnergyMap const & /*weights*/ ) const override { return false; }

	void
	eval_intrares_energy(
		conformation::Residue const &,
		pose::Pose const &,
		core::scoring::ScoreFunction const &,
		core::scoring::EnergyMap & ) const override {};

	/// @brief Interface from the LongRangeTwoBodyEnergy base class; returns "true" if there's any non-zero
	/// or potentially non-zero interaction between a pair of residues in a pose.
	bool
	defines_residue_pair_energy(
		pose::Pose const & pose,
		Size res1,
		Size res2
	) const override;

	/// @brief Identification for this LR2B energy that links it with the
	/// long-range energy container that it stores in the Energies object
	core::scoring::methods::LongRangeEnergyType long_range_type() const override { return core::scoring::methods::rna_stub_coord_lr; }

	core::Size version() const override{ return 1; }

private:

	mutable Size res1_, res2_;
	std::tuple< utility::vector1<int>, utility::vector1<std::string>, utility::vector1< std::string > > jump_resnum_and_chain_;
	Vector target_xyz_in_takeoff_frame_;
	scoring::func::FuncOP func_;
	pose::rna::StubStubType stub_stub_type_;
	kinematics::RTOP reference_RT_;
};


} //scoring
} //core

#endif // INCLUDED_core_scoring_ScoreFunction_HH
