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
#include <core/scoring/rna/RNA_StubCoordinateEnergy.fwd.hh>

// Package headers
#include <core/scoring/methods/ContextIndependentLRTwoBodyEnergy.hh>
#include <core/scoring/func/Func.fwd.hh>

#include <numeric/xyzVector.hh>

namespace core {
namespace scoring {
namespace rna {


class RNA_StubCoordinateEnergy : public methods::ContextIndependentLRTwoBodyEnergy {

public:
	typedef methods::ContextIndependentLRTwoBodyEnergy parent;
	typedef methods::EnergyMethodOP EnergyMethodOP;
	typedef numeric::xyzVector< Real > Vector;
	typedef numeric::xyzMatrix< Real > Matrix;

public:

	RNA_StubCoordinateEnergy();

	/// clone
	virtual
	methods::EnergyMethodOP
	clone() const;

	///
	virtual
	void
	setup_for_scoring( pose::Pose & pose, ScoreFunction const & ) const;

	virtual
	void
	residue_pair_energy(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		pose::Pose const & pose,
		ScoreFunction const & sfxn,
		EnergyMap & emap
	) const;

	virtual
	void indicate_required_context_graphs( utility::vector1< bool > & ) const {}

	virtual
	bool
	defines_intrares_energy( EnergyMap const & /*weights*/ ) const { return false; }

	virtual
	void
	eval_intrares_energy(
		conformation::Residue const &,
		pose::Pose const &,
		ScoreFunction const &,
		EnergyMap & ) const {};

	/// @brief Interface from the LongRangeTwoBodyEnergy base class; returns "true" if there's any non-zero
	/// or potentially non-zero interaction between a pair of residues in a pose.
	virtual
	bool
	defines_residue_pair_energy(
		pose::Pose const & pose,
		Size res1,
		Size res2
	) const;

	/// @brief Identification for this LR2B energy that links it with the
	/// long-range energy container that it stores in the Energies object
	virtual methods::LongRangeEnergyType long_range_type() const { return methods::rna_stub_coord_lr; }

	virtual
	core::Size version() const{ return 1; }

private:

	Size takeoff_res_, landing_res_;
	Vector target_xyz_in_takeoff_frame_;
	scoring::func::FuncOP func_;

};


} //rna
} //scoring
} //core

#endif // INCLUDED_core_scoring_ScoreFunction_HH
