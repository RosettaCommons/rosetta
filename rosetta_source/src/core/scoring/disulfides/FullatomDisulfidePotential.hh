// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/disulfides/FullatomDisulfidePotential.hh
/// @brief  Fullatom Disulfide Potential class declaration
/// @author Andrew Leaver-Fay

#ifndef INCLUDED_core_scoring_disulfides_FullatomDisulfidePotential_hh
#define INCLUDED_core_scoring_disulfides_FullatomDisulfidePotential_hh

// Unit headers
#include <core/scoring/disulfides/FullatomDisulfidePotential.fwd.hh>

// Package headers
#include <core/scoring/disulfides/DisulfideAtomIndices.fwd.hh>
#include <core/scoring/EnergyMap.fwd.hh>
#include <core/scoring/constraints/Func.hh>
// AUTO-REMOVED #include <core/scoring/constraints/CircularHarmonicFunc.hh>
#include <core/scoring/constraints/CircularSigmoidalFunc.hh>

// Project headers
#include <core/types.hh>
#include <core/conformation/Residue.fwd.hh>
// AUTO-REMOVED #include <core/pose/Pose.fwd.hh>

// ObjexxFCL headers
// AUTO-REMOVED #include <ObjexxFCL/FArray1D.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <numeric/interpolation/Histogram.fwd.hh>

#include <utility/vector1.hh>


namespace core {
namespace scoring {
namespace disulfides {

class FullatomDisulfidePotential : public utility::pointer::ReferenceCount
{
public:
	typedef utility::pointer::ReferenceCount parent;

public:
	FullatomDisulfidePotential();
	virtual ~FullatomDisulfidePotential();

	void
	print_score_functions() const;

	void
	score_this_disulfide(
		conformation::Residue const & res1,
		conformation::Residue const & res2,
		DisulfideAtomIndices const & res1_atom_indices,
		DisulfideAtomIndices const & res2_atom_indices,
		Energy & distance_score_this_disulfide,
		Energy & csangles_score_this_disulfide,
		Energy & dihedral_score_this_disulfide,
		Energy & ca_dihedral_sc_this_disulf,
		bool & truefalse_fa_disulf
	) const;

	void
	get_disulfide_derivatives(
		conformation::Residue const & res1,
		conformation::Residue const & res2,
		DisulfideAtomIndices const & res1_atom_indices,
		DisulfideAtomIndices const & res2_atom_indices,
		Size const at1,
		EnergyMap const & weights,
		Vector & F1,
		Vector & F2
	) const;

private:

	void
	get_disulfide_params(
		conformation::Residue const & res1,
		conformation::Residue const & res2,
		DisulfideAtomIndices const & res1_atom_indices,
		DisulfideAtomIndices const & res2_atom_indices,
		Distance & dist_between_sulfurs, // dist between cys sulfurs
		Real & cs_bond_angle_1,
		Real & cs_bond_angle_2,
		Real & disulf_dihedral_angle, // dihedral (torsion) angle, cbeta-s-s-cbeta
		Real & disulf_ca_dihedral_angle_1,
		 // dihedral (torsion) angle, calpha1-cbeta1-s1-s2
		Real & disulf_ca_dihedral_angle_2 // dihedral (torsion) angle, calpha2-cbeta2-s2-s1
	) const;

private:
	// for interp300 this is 300, others it is 100.
	Real const disulf_ssdist_cutoff_;

	CBSG_Dihedral_FuncOP cbsg_dihedral_func_;
	SGSG_Dihedral_FuncOP sgsg_dihedral_func_;
	CB_Angle_FuncOP cb_angle_func_;
	SG_Dist_FuncOP sg_dist_func_;

};

class CBSG_Dihedral_Func : public constraints::Func
{
public:
	CBSG_Dihedral_Func();
	virtual ~CBSG_Dihedral_Func();

	virtual
	Real
	func( Real const ) const;

	constraints::FuncOP
	clone() const { return new CBSG_Dihedral_Func( *this ); };

	virtual
	Real
	dfunc( Real const ) const;
private:
	//Real const cbsg_pos_peak_;
	//Real const cbsg_pos_sd_;
	//Real const cbsg_neg_peak_;
	//Real const cbsg_neg_sd_;

	core::scoring::constraints::CircularSigmoidalFunc csf_cbang1_;
	core::scoring::constraints::CircularSigmoidalFunc csf_cbang2_;
	core::scoring::constraints::CircularSigmoidalFunc csf_cbang3_;

};

class SGSG_Dihedral_Func : public constraints::Func
{
public:
	SGSG_Dihedral_Func();

	~SGSG_Dihedral_Func();

	virtual
	Real
	func( Real const ) const;

	constraints::FuncOP
	clone() const { return new SGSG_Dihedral_Func( *this ); };

	virtual
	Real
	dfunc( Real const ) const;
private:
	/// Access the histogram for this Func
	//static
	//numeric::interpolation::HistogramCOP<core::Real,core::Real>::Type
	//fa_sgsg_dihedral_scores();

	core::scoring::constraints::CircularSigmoidalFunc csf_cbang1a_;
	core::scoring::constraints::CircularSigmoidalFunc csf_cbang2a_;
	core::scoring::constraints::CircularSigmoidalFunc csf_cbang1b_;
	core::scoring::constraints::CircularSigmoidalFunc csf_cbang2b_;
};

class CB_Angle_Func : public constraints::Func
{
public:
	CB_Angle_Func();

	~CB_Angle_Func();

	constraints::FuncOP
	clone() const { return new CB_Angle_Func( *this ); };

	virtual
	Real
	func( Real const ) const;


	virtual
	Real
	dfunc( Real const ) const;
private:
	/// Access the histogram for this Func
	//static
	//numeric::interpolation::HistogramCOP<core::Real,core::Real>::Type
	//fa_csang_scores();

	//core::scoring::constraints::CircularHarmonicFunc  chf_cbang_;
	//core::scoring::constraints::CircularSigmoidalFunc csf_cbang_;

	core::scoring::constraints::CircularSigmoidalFunc csf_cbang1_;
	core::scoring::constraints::CircularSigmoidalFunc csf_cbang2_;

	//core::scoring::constraints::CircularSigmoidalFunc csf_cbang2_;
	//core::scoring::constraints::CircularSigmoidalFunc csf_cbang3_;
};

class SG_Dist_Func : public constraints::Func
{
public:
	SG_Dist_Func();

	~SG_Dist_Func();

	constraints::FuncOP
	clone() const { return new SG_Dist_Func( *this ); };

	virtual
	Real
	func( Real const ) const;

	virtual
	Real
	dfunc( Real const ) const;

private:
	/// Access the histogram for this Func
	static
	numeric::interpolation::HistogramCOP<core::Real,core::Real>::Type
	fa_ssdist_scores();

	/// Access the histogram for this Func's derivative
	static
	numeric::interpolation::HistogramCOP<core::Real,core::Real>::Type
	fa_ssdist_scores_deriv();
};



} // namespace disulfides
} // namespace scoring
} // namespace core

#endif
