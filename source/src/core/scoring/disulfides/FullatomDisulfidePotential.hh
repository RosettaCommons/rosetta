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
#include <core/scoring/func/Func.hh>
#include <core/scoring/func/CircularSigmoidalFunc.hh>
#include <core/scoring/DerivVectorPair.hh>

// Project headers
#include <core/types.hh>
#include <core/conformation/Residue.fwd.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <numeric/interpolation/Histogram.fwd.hh>

#include <utility/vector1.hh>


namespace core {
namespace scoring {
namespace disulfides {

//fpd  May 6 2013 refitting
//fpd  Use the following formulations:
//fpd  For CBSG_Dihedral:
//   mixture of von mises with three components
//      score = sum( -A - kappa*cos(m-x) )
//      logAs = -15.8644 mean =  -72.2016 kappa = 13.3778
//      logAs = -16.9017 mean =   78.0303 kappa = 13.6370
//      logAs =  -7.0219 mean = -172.5505 kappa =  2.9327
//fpd  For SGSG_Dihedral:
//  mixture of von mises with two components
//     score = sum( -A - kappa*cos(m-x) )
//     logAs = -32.9599 mean = -86.0964 kappa = 30.9053
//     logAs = -23.3471 mean =  92.3915 kappa = 20.9805
//fpd  For CB angle:
//  single von mises
//     score = -A - kappa*cos(m-x)  [A handles normzalization of the distribution]
//     logAs = -419.8120 mean = 104.22 kappa = 419.7
//fpd  For SG-SG length:
// skewed normal distribution with:
//   location = 2.01
//      scale = 0.08
//      shape = 6.0 (transformed 'skewness')
//fpd  perhaps should be read from the DB?
struct FullatomDisulfideParams13 {
	FullatomDisulfideParams13() {
		d_location=2.01; d_scale=0.08; d_shape=6;

		a_logA=-419.8120; a_kappa=419.7; a_mu=104.22;

		dss_logA1=-32.9599; dss_kappa1=30.9053; dss_mu1=-86.0964;
		dss_logA2=-23.3471; dss_kappa2=20.9805; dss_mu2=92.3915;

		dcs_logA1=-15.8644; dcs_mu1=-72.2016;  dcs_kappa1=13.3778;
		dcs_logA2=-16.9017; dcs_mu2=78.0303;   dcs_kappa2=13.6370;
		dcs_logA3=-7.0219 ; dcs_mu3=-172.5505; dcs_kappa3=2.9327;
	}

	// distance
	Real d_location, d_scale, d_shape;
	// angle
	Real a_logA, a_kappa, a_mu;
	// SS dih
	Real dss_logA1, dss_kappa1, dss_mu1, dss_logA2, dss_kappa2, dss_mu2;
	// CS dih
	Real dcs_logA1, dcs_mu1, dcs_kappa1, dcs_logA2, dcs_mu2, dcs_kappa2, dcs_logA3, dcs_mu3, dcs_kappa3;
};


class FullatomDisulfidePotential : public utility::pointer::ReferenceCount
{
public:
	typedef utility::pointer::ReferenceCount parent;

public:
	FullatomDisulfidePotential();
	virtual ~FullatomDisulfidePotential();

	void
	print_score_functions() const;

	//fpd old version
	void
	score_this_disulfide_old(
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
	get_disulfide_derivatives_old(
		conformation::Residue const & res1,
		conformation::Residue const & res2,
		DisulfideAtomIndices const & res1_atom_indices,
		DisulfideAtomIndices const & res2_atom_indices,
		Size const at1,
		EnergyMap const & weights,
		Vector & F1,
		Vector & F2
	) const;

	//fpd new version
	void
	score_this_disulfide(
		conformation::Residue const & res1,
		conformation::Residue const & res2,
		DisulfideAtomIndices const & res1_atom_indices,
		DisulfideAtomIndices const & res2_atom_indices,
		Energy & score
	) const;

	void
	get_disulfide_derivatives(
		conformation::Residue const & res1,
		conformation::Residue const & res2,
		DisulfideAtomIndices const & res1_atom_indices,
		DisulfideAtomIndices const & res2_atom_indices,
		EnergyMap const & weights,
		utility::vector1< DerivVectorPair > & r1_atom_derivs,
		utility::vector1< DerivVectorPair > & r2_atom_derivs
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

	//fpd relative weighings and background probability on the new distribution
	Real wt_dihSS_, wt_dihCS_, wt_ang_, wt_len_;
	Real shift_;  // shift disulf down
	Real mest_;   // flatten
	FullatomDisulfideParams13 params_; // new parameters
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
	core::scoring::constraints::CircularSigmoidalFunc csf_cbang1_;
	core::scoring::constraints::CircularSigmoidalFunc csf_cbang2_;

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
