// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
//////////////////////////////////////////////////////////////////////
/// @begin Etable.hh
///
/// @brief
/// A class for generating the table for fa_atr/rep and fa_sol
///
/// @detailed
/// This class is called upon by the ScoringManager. Since actual calculating of the LJ potential
/// is time consuming if done multiple times, this class precomputes and discritizes the potential
/// (meaning that the potential is broken down into bins). Once the bins have been created, it will
/// smooth out the bins, for better interpolation.
///
///
/// @authors
/// Phil Bradley
/// Andrew Leaver-Fay
/// Steven Combs - comments and skipping of virtual atoms
///
/// @last_modified August 2013
/////////////////////////////////////////////////////////////////////////


#ifndef INCLUDED_core_scoring_etable_Etable_hh
#define INCLUDED_core_scoring_etable_Etable_hh

// Unit Headers
#include <core/scoring/etable/Etable.fwd.hh>

// Package Headers
#include <core/scoring/etable/EtableOptions.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/ScoreType.hh>

#include <core/scoring/etable/count_pair/CountPairFunction.fwd.hh>
#include <core/scoring/trie/RotamerTrieBase.fwd.hh>
#include <core/scoring/trie/TrieCountPairBase.fwd.hh>
#include <core/scoring/etable/etrie/EtableAtom.hh>

#include <core/chemical/AtomTypeSet.hh>
#include <core/conformation/Atom.hh>
#include <utility/pointer/access_ptr.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray1A.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/FArray3D.hh>

#include <utility/assert.hh>
#include <utility/pointer/ReferenceCount.hh>

#include <ObjexxFCL/FArray1.fwd.hh>
#include <ObjexxFCL/FArray2.fwd.hh>


// C++ Headers

namespace core {
namespace scoring {
namespace etable {

/// @brief %SplineParameters is a simple struct for holding the cubic spline polynomials used in
/// the etable to interpolate the lennard-jones attractive and
/// LK-solvation terms to zero smoothly.  These splines have
/// exactly two knots to represent them, and the same x values
/// are used for all the knots: thus the only parameters needed
/// are the y values at the knots, and the second-derivatives
/// for the polynomials at knots.
struct SplineParameters
{
	Real ylo,  yhi;
	Real y2lo, y2hi;
	SplineParameters() :
		ylo(0.0),
		yhi(0.0),
		y2lo(0.0),
		y2hi(0.0)
	{}
};

/// @brief the ExtraQuadraticRepulsion class is used to
/// add in extra repulsion for particular atom pairs, if needed,
/// (e.g. for OCbb/OCbb) where the functional form is:
/// fa_rep += (xhi - x)^2 * slope
/// for values of x between xlo and xhi, and
/// fa_rep += (x - xlo ) * extrapolated_slope + ylo
/// where extrapolated slope can be anything, but, to defined
/// a function with continuous derivatives, should be
/// extrapolated_slope = (xhi-xlo)^2*slope.  This is the
/// analytical implementation of the "modify_pot" behavior.
struct ExtraQuadraticRepulsion
{
	Real xlo, xhi;
	Real slope;
	Real extrapolated_slope, ylo;
	ExtraQuadraticRepulsion() :
		xlo( 0.0 ),
		xhi( 0.0 ),
		slope( 0.0 ),
		extrapolated_slope( 0.0 ),
		ylo( 0.0 )
	{}
};

/// @brief %EtableParamsOnePair describes all of the parameters for a particular
/// pair of atom types necessary to evaluate the Lennard-Jones and LK solvation
/// energies.
struct EtableParamsOnePair
{
	Real maxd2;
	bool hydrogen_interaction;
	Real ljrep_linear_ramp_d2_cutoff;
	Real lj_r6_coeff;
	Real lj_r12_coeff;
	Real lj_switch_intercept;
	Real lj_switch_slope;
	Real lj_minimum;
	Real lj_val_at_minimum;
	Real ljatr_spline_xlo;
	Real ljatr_spline_xhi;
	SplineParameters ljatr_spline_parameters;
	ExtraQuadraticRepulsion ljrep_extra_repulsion;
	bool ljrep_from_negcrossing;

	Real lk_coeff1;
	Real lk_coeff2;
	Real lk_min_dis2sigma_value;
	Real fasol_spline_close_start;
	Real fasol_spline_close_end;

	SplineParameters fasol_spline_close;  // mututal desolvation of atoms 1 and 2
	SplineParameters fasol_spline_far;    // mututal desolvation of atoms 1 and 2

	SplineParameters fasol_spline1_close; // desolvation of atom 1 by atom 2
	SplineParameters fasol_spline1_far;   // desolvation of atom 1 by atom 2

	SplineParameters fasol_spline2_close; // desolvation of atom 2 by atom 1
	SplineParameters fasol_spline2_far;   // desolvation of atom 2 by atom 1
	Real ljatr_final_weight;
	Real fasol_final_weight;

	EtableParamsOnePair() :
		maxd2( 0.0 ),
		hydrogen_interaction( false ),
		ljrep_linear_ramp_d2_cutoff(0.0),
		lj_r6_coeff(0.0),
		lj_r12_coeff(0.0),
		lj_switch_intercept(0.0),
		lj_switch_slope(0.0),
		lj_minimum(0.0),
		lj_val_at_minimum(0.0),
		ljatr_spline_xlo(0.0),
		ljatr_spline_xhi(0.0),
		ljrep_from_negcrossing(false),
		lk_coeff1(0.0),
		lk_coeff2(0.0),
		lk_min_dis2sigma_value(0.0),
		fasol_spline_close_start(0.0),
		fasol_spline_close_end(0.0),
		ljatr_final_weight(1.0),
		fasol_final_weight(1.0)
	{}


};


/// @brief Class definition for Etable
class Etable : public utility::pointer::ReferenceCount {

public:
	///@brief Automatically generated virtual destructor for class deriving directly from ReferenceCount
	virtual ~Etable();

	///  constructor
	Etable(
		chemical::AtomTypeSetCAP atom_set_in, // like etable namespace
		EtableOptions const & options,
		std::string const alternate_parameter_set = ""
	);

	Size
	n_atomtypes() const {
		return n_atomtypes_;
	}

	/// const access to the arrays
	ObjexxFCL::FArray3D< Real > const &
	ljatr() const
	{
		return ljatr_;
	}

	ObjexxFCL::FArray3D< Real > const &
	ljrep() const
	{
		return ljrep_;
	}

	ObjexxFCL::FArray3D< Real > const &
	solv1() const
	{
		return solv1_;
	}

	ObjexxFCL::FArray3D< Real > const &
	solv2() const
	{
		return solv2_;
	}

	/// const access to the deriv arrays
	ObjexxFCL::FArray3D< Real > const &
	dljatr() const
	{
		return dljatr_;
	}

	ObjexxFCL::FArray3D< Real > const &
	dljrep() const
	{
		return dljrep_;
	}

	/// @brief return the solvation derivative table for the desolvation of atom1 by atom2
	ObjexxFCL::FArray3D< Real > const &
	dsolv1() const
	{
		return dsolv1_;
	}

	/// @brief return the solvation derivative table that combines atom1 and atom2's desolvations
	ObjexxFCL::FArray3D< Real > const &
	dsolv() const
	{
		return dsolv_;
	}

	Real
	max_dis() const
	{
		return max_dis_;
	}

	Real
	get_safe_max_dis2() const
	{
		return safe_max_dis2;
	}

	int
	get_bins_per_A2() const
	{
		return bins_per_A2;
	}

	chemical::AtomTypeSetCAP
	atom_set() const
	{
		return atom_set_;
	}

	Real
	hydrogen_interaction_cutoff2() const
	{
		return hydrogen_interaction_cutoff2_;
	}

	Real
	max_heavy_heavy_cutoff() const {
		return max_dis_;
	}

	Real
	max_heavy_hydrogen_cutoff() const {
		return max_heavy_hydrogen_cutoff_;
	}

	Real
	max_hydrogen_hydrogen_cutoff() const {
		return max_hydrogen_hydrogen_cutoff_;
	}

	/// @brief The distance cutoff beyond which any pair of heavy-atoms is
	/// guaranteed to have an interaction energy of zero.  This function is
	/// used by the NeighborList
	Real
	nblist_dis2_cutoff_XX() const
	{
		return nblist_dis2_cutoff_XX_;
	}

	/// @brief The distance cutoff beyond which a hydrogen/heavy-atom pair is
	/// guaranteed to have an interaction energy of zero.  This function is used
	/// by the NeighborList
	Real
	nblist_dis2_cutoff_XH() const
	{
		return nblist_dis2_cutoff_XH_;
	}

	/// @brief The distance cutoff beyond which any hydrogen/hydrogen pair is guaranteed
	/// to have an interaction energy of zero.  This function is used by the NeighborList
	Real
	nblist_dis2_cutoff_HH() const
	{
		return nblist_dis2_cutoff_HH_;
	}

	/// @brief Returns the maximum lj radius for any non-hydrogen
	/// atom as defined by the atom-type-set used to create this Etable.
	Real
	max_non_hydrogen_lj_radius() const;

	/// @brief Returns the maximum lj radius for any hydrogen atom as
	/// defined by the input atom-type-set used to create this Etable.
	Real
	max_hydrogen_lj_radius() const;

	/// @brief Return the Lennard-Jones radius for an atom.
	inline
	Real
	lj_radius( int const i ) const
	{
		return lj_radius_[i];
	}

	/// @brief Return the Lennard-Jones well depth for an atom
	Real
	lj_wdepth( int const i ) const
	{
		return lj_wdepth_[i];
	}

	/// @brief Return the Lazardis Karplus DGFree value for an atom
	Real
	lk_dgfree( int const i ) const
	{
		return lk_dgfree_[i];
	}

	/// @brief Return the Lazaridis Karplus volume for an atom
	Real
	lk_volume( int const i ) const
	{
		return lk_volume_[i];
	}

	/// @brief Return the Lazaridis Karplus "lambda" value (correlation distance) for an atom
	Real
	lk_lambda( int const i ) const
	{
		return lk_lambda_[i];
	}

	/// @brief Use the analytic_etable_evaluation function to evaluate the energy
	/// of two atoms, but  evaluate the function at the old grid points and then
	/// interpolate between them the way the existing etable does (in square
	/// distance space). Useful for comparing the original etable evaluation with the
	/// analytic evaluation.
	void
	interpolated_analytic_etable_evaluation(
		conformation::Atom const & at1,
		conformation::Atom const & at2,
		Real & lj_atrE,
		Real & lj_repE,
		Real & fa_solE,
		Real & d2
	) const;

	/// @brief Use an analytic functional form of the etable to evaluate an atom-pair energy
	/// without reading from the enormous and uncachable tables.
	inline
	void
	analytic_etable_evaluation(
		conformation::Atom const & at1,
		conformation::Atom const & at2,
		Real & lj_atrE,
		Real & lj_repE,
		Real & fa_solE,
		Real & d2
	) const;

	/// @brief Use an analytic functional form of the etable to evaluate only the LK atom-pair energy
	/// computing the desolvation of atom1 by atom2 separately from the desolvation of atom2 by atom1.
	inline
	void
	analytic_lk_energy(
		conformation::Atom const & at1,
		conformation::Atom const & at2,
		Real & fa_solE1,
		Real & fa_solE2
	) const;

	/// @brief Analytically evaluate the energy derivatives for a pair of atoms
	inline
	void
	analytic_etable_derivatives(
		conformation::Atom const & at1,
		conformation::Atom const & at2,
		Real & dljatrE_ddis,
		Real & dljrepE_ddis,
		Real & dfasolE_ddis,
		Real & inv_d
	) const;

	/// @brief Analytically evaluate the LK solvation derivatives for a pair of atoms, separately
	/// computing the derivative for atom2's desolvation of atom1 (dfasolE1_ddis) and atom1's desolvation
	/// of atom2 (dfasolE2_ddis).
	inline
	void
	analytic_lk_derivatives(
		conformation::Atom const & at1,
		conformation::Atom const & at2,
		Real & dfasolE1_ddis,
		Real & dfasolE2_ddis,
		Real & inv_d
	) const;

	/// @brief Use an analytic functional form of the etable to evaluate an atom-pair energy
	/// without reading from the enormous and uncachable tables. This version is reserved for
	/// interactions with hydrogens (hydorgen/hydrogen and hydrogen heavyatom)
	//void
	//analytic_etable_evaluation_H(
	//	conformation::Atom const & at1,
	//	conformation::Atom const & at2,
	//	EtableParamsOnePair const & p,
	//	Real const dis2,
	//	Real & lj_atrE,
	//	Real & lj_repE,
	//	Real & fa_solE
	//) const;

private:
	// initialization routines

	void dimension_etable_arrays();
	void initialize_from_input_atomset( chemical::AtomTypeSetCAP atom_set_in );
	void calculate_nblist_distance_thresholds( EtableOptions const & options );
	void
	read_alternate_parameter_set(
		chemical::AtomTypeSetCAP atom_set_in,
		std::string const alternate_parameter_set
	);
	void calculate_hydrogen_atom_reach();
	void initialize_carbontypes_to_linearize_fasol();


private:
	inline
	Real
	analytic_ljrep_linearized(
		Real const dis,
		EtableParamsOnePair const & p
	) const;

	inline
	Real
	analytic_lj_generic_form(
		Real const dis2,
		Real const inv_dis2,
		EtableParamsOnePair const & p
	) const;

	inline
	Real
	analytic_ljatr_spline_ramp_to_zero(
		Real const dis,
		EtableParamsOnePair const & p
	) const;

	inline
	Real
	analytic_ljatr_spline_ramp_to_zero_deriv(
		Real const dis,
		EtableParamsOnePair const & p
	) const;


	/// @brief Evaluate the mututal desolvation energy as atom 1 and atom 2 approach.
	/// Combine the desolvation of atom 1 by atom 2 into the same value as the desolvation
	/// of atom 2 by atom 1.
	inline
	void
	analytic_lk_evaluation(
		int const atype1,
		int const atype2,
		EtableParamsOnePair const & p,
		Real const dis,
		Real const inv_dis2,
		Real & fa_solE
	) const;

	/// @brief Evaluate the LK solvation energy for a pair of atoms
	/// individually; fa_solE1 returns the desolvation of atom 1 by atom 2
	/// and fa_solE2 returns the desolvation of atom 2 by atom 1.
	inline
	void
	analytic_lk_evaluation_individual(
		int const atype1,
		int const atype2,
		EtableParamsOnePair const & p,
		Real const dis,
		Real const inv_dis2,
		Real & fa_solE1,
		Real & fa_solE2
	) const;


public:

	static
	inline
	Real
	eval_spline(
		Real const x,
		Real const xlo,
		Real const xhi,
		SplineParameters const & sp
	);

	static
	inline
	Real
	spline_deriv(
		Real const x,
		Real const xlo,
		Real const xhi,
		SplineParameters const & sp
	);

	static
	inline
	Real
	eval_spline(
		Real const x,
		Real const xlo,
		Real const xhi,
		Real const width,
		Real const invwidth,
		SplineParameters const & sp
	);

	static
	inline
	Real
	spline_deriv(
		Real const x,
		Real const xlo,
		Real const xhi,
		Real const width,
		Real const invwidth,
		SplineParameters const & sp
	);

	inline
	EtableParamsOnePair const &
	analytic_params_for_pair(
		Size atype1,
		Size atype2
	) const;

private:

	void
	output_etable(
		ObjexxFCL::FArray3D<Real> & etable,
		std::string label,
		std::ostream & out
	);

	void
	input_etable(
		ObjexxFCL::FArray3D<Real> & etable,
		const std::string label,
		std::istream & in
	);

	inline
	EtableParamsOnePair &
	analytic_params_for_pair(
		Size atype1,
		Size atype2
	);

public:
	Real get_lj_hbond_OH_donor_dis() const { return lj_hbond_OH_donor_dis; }
	Real get_lj_hbond_hdis() const { return lj_hbond_hdis; }

private:

	chemical::AtomTypeSetCAP atom_set_;

	// parameters:
	int const n_atomtypes_;

	// from options
	Real const max_dis_;
	int const bins_per_A2;
	Real const Wradius; // global mod to radii
	Real const lj_switch_dis2sigma; // actual value used for switch
	Real const max_dis2;
	int const etable_disbins;
	Real const lj_hbond_OH_donor_dis;
	Real const lj_hbond_dis;

	// hard-coded for now
	bool const lj_use_lj_deriv_slope;
	Real const lj_slope_intercept;
	bool const lj_use_hbond_radii;
	Real const lj_hbond_hdis;
	Real const lj_hbond_accOch_dis;
	Real const lj_hbond_accOch_hdis;
	bool const lj_use_water_radii;
	Real const lj_water_dis;
	Real const lj_water_hdis;
	Real const lk_min_dis2sigma;
	Real const min_dis;
	Real const min_dis2; // was double
	bool const add_long_range_damping;
	Real const long_range_damping_length;
	Real const epsilon;
	Real const safe_max_dis2;
	Real hydrogen_interaction_cutoff2_;
	Real max_heavy_hydrogen_cutoff_;
	Real max_hydrogen_hydrogen_cutoff_;
	Real nblist_dis2_cutoff_XX_; // for use by the old-style neighborlist
	Real nblist_dis2_cutoff_XH_; // for use by the old-style neighborlist
	Real nblist_dis2_cutoff_HH_; // for use by the old-style neighborlist
	Real max_non_hydrogen_lj_radius_;
	Real max_hydrogen_lj_radius_;

	// these three derived from other data
	Real lj_switch_sigma2dis;
	Real lj_switch_value2wdepth;
	Real lj_switch_slope_sigma2wdepth;

	ObjexxFCL::FArray2D< Real > lj_sigma_;
	ObjexxFCL::FArray2D< Real > lj_r6_coeff_;
	ObjexxFCL::FArray2D< Real > lj_r12_coeff_;
	ObjexxFCL::FArray2D< Real > lj_switch_intercept_;
	ObjexxFCL::FArray2D< Real > lj_switch_slope_;
	ObjexxFCL::FArray1D< Real > lk_inv_lambda2_;
	ObjexxFCL::FArray2D< Real > lk_coeff_;
	ObjexxFCL::FArray2D< Real > lk_min_dis2sigma_value_;

	/// OK: these values reflect the grid point with the lowest energy,
	/// which could be different from the sum of the lj radii or the lj well depths.
	/// These vals define the switch point from attractive to repulsive.  Repulsive
	/// interactions start counting at the lj_minima
	//ObjexxFCL::FArray2D< Real > lj_minima;
	//ObjexxFCL::FArray2D< Real > lj_vals_at_minima;

	/// Data needed to describe the splines for the ljatr term
	Real ljatr_spline_xlo;
	Real ljatr_spline_xhi;
	Real ljatr_spline_diff_xhi_xlo;
	Real ljatr_spline_diff_xhi_xlo_inv;
	//ObjexxFCL::FArray2D< std::pair< Real, Real > > ljatr_spline_xlo_xhi;
	//ObjexxFCL::FArray2D< SplineParameters > ljatr_spline_parameters;

	/// Add extra repulsion, if desired, for certain atom pair interactions
	//ObjexxFCL::FArray2D< ExtraQuadraticRepulsion > ljrep_extra_repulsion;

	/// Data needed to describe the splines for the fasol term:
	/// There are two splines needed: one to smooth the transition to where the fasol term goes flat
	/// as the distance becomes less than the sum of the van der Waals radii (the "close" spline),
	/// and a second to smooth the transition where the distance goes to the
	/// fa_max_dis and the energy goes to 0 (the far spline).  In between the close
	/// and far values, the exponential form of the energy function is used.
	//ObjexxFCL::FArray2D< std::pair< Real, Real > > fasol_spline_close_start_end; // starting and ending points for the lower spline
	//ObjexxFCL::FArray2D< SplineParameters > fasol_spline_close; // parameters for the lower spline
	Real fasol_spline_far_xlo;
	Real fasol_spline_far_xhi;
	Real fasol_spline_far_diff_xhi_xlo;
	Real fasol_spline_far_diff_xhi_xlo_inv;
	//ObjexxFCL::FArray2D< SplineParameters > fasol_spline_far;

	/// Should the repulsive component start at ljatr+ljrep = 0 (true) or
	/// when ljatr reaches its minimum (false )
	//ObjexxFCL::FArray2D< Size > ljrep_from_negcrossing;

	/// Data to turn off portions of the interactions between certain atom pairs
	/// e.g. Hydrogen atoms are repulsive, only, as are the REPLONLY atoms.
	//ObjexxFCL::FArray2D< Real > ljatr_final_weight;
	//ObjexxFCL::FArray2D< Real > fasol_final_weight;


	//
	utility::vector1< Real > lj_radius_;
	utility::vector1< Real > lj_wdepth_;
	utility::vector1< Real > lk_dgfree_;
	utility::vector1< Real > lk_volume_;
	utility::vector1< Real > lk_lambda_;

	// the etables themselves
	ObjexxFCL::FArray3D< Real > ljatr_;
	ObjexxFCL::FArray3D< Real > ljrep_;
	ObjexxFCL::FArray3D< Real > solv1_;
	ObjexxFCL::FArray3D< Real > solv2_;
	ObjexxFCL::FArray3D< Real > dljatr_;
	ObjexxFCL::FArray3D< Real > dljrep_;
	ObjexxFCL::FArray3D< Real > dsolv_;
	ObjexxFCL::FArray3D< Real > dsolv1_;

	utility::vector1< Size > carbon_types; // indices for the standard carbon types

	// upper-triangle of the n_atomtypes x n_atomtypes table.  N^2 - (N*(N-1))/2) entries.
	// indexed for at1, at2 w/ at1<at2:
	// ((at1-1)*n_atomtypes + (at2-1) + 1 - (at1*(at1-1)/2)
	utility::vector1< EtableParamsOnePair > analytic_parameters;

	bool slim_;

	// private methods

	// get the AtomType object corresponding to a give type index
	chemical::AtomType const &
	atom_type( int const type )
	{
		return (*atom_set_)[ type ];
	}

	//void smooth_etables();
	//void modify_pot();
	void make_pairenergy_table();

	int calculate_normal_disbins() const;

	void
	damp_long_range(
		int const normal_disbins,
		ObjexxFCL::FArray1A< Real > ljatr,
		ObjexxFCL::FArray1A< Real > dljatr,
		ObjexxFCL::FArray1A< Real > ljrep,
		ObjexxFCL::FArray1A< Real > dljrep,
		ObjexxFCL::FArray1A< Real > fasol1,
		ObjexxFCL::FArray1A< Real > fasol2,
		ObjexxFCL::FArray1A< Real > dfasol,
		ObjexxFCL::FArray1A< Real > dfasol1
	);

	void
	modify_pot_one_pair(
		Size atype1,
		Size atype2,
		ObjexxFCL::FArray1A< Real > ljatr,
		ObjexxFCL::FArray1A< Real > dljatr,
		ObjexxFCL::FArray1A< Real > ljrep,
		ObjexxFCL::FArray1A< Real > dljrep,
		ObjexxFCL::FArray1A< Real > fasol1,
		ObjexxFCL::FArray1A< Real > fasol2,
		ObjexxFCL::FArray1A< Real > dfasol,
		ObjexxFCL::FArray1A< Real > dfasol1
	);

	void
	assign_parameters_to_full_etables(
		Size atype1,
		Size atype2,
		ObjexxFCL::FArray1A< Real > ljatr,
		ObjexxFCL::FArray1A< Real > dljatr,
		ObjexxFCL::FArray1A< Real > ljrep,
		ObjexxFCL::FArray1A< Real > dljrep,
		ObjexxFCL::FArray1A< Real > fasol1,
		ObjexxFCL::FArray1A< Real > fasol2,
		ObjexxFCL::FArray1A< Real > dfasol,
		ObjexxFCL::FArray1A< Real > dfasol1
	);

	void
	smooth_etables_one_pair(
		Size atype1,
		Size atype2,
		ObjexxFCL::FArray1A< Real > ljatr,
		ObjexxFCL::FArray1A< Real > dljatr,
		ObjexxFCL::FArray1A< Real > ljrep,
		ObjexxFCL::FArray1A< Real > dljrep,
		ObjexxFCL::FArray1A< Real > fasol1,
		ObjexxFCL::FArray1A< Real > fasol2,
		ObjexxFCL::FArray1A< Real > dfasol,
		ObjexxFCL::FArray1A< Real > dfasol1
	);

	void
	zero_hydrogen_and_water_ljatr_one_pair(
		Size atype1,
		Size atype2,
		ObjexxFCL::FArray1A< Real > ljrep,
		ObjexxFCL::FArray1A< Real > ljatr,
		ObjexxFCL::FArray1A< Real > dljatr,
		ObjexxFCL::FArray1A< Real > fasol1,
		ObjexxFCL::FArray1A< Real > fasol2,
		ObjexxFCL::FArray1A< Real > dfasol,
		ObjexxFCL::FArray1A< Real > dfasol1
	);

	// helper functions
	void
	precalc_etable_coefficients(
		ObjexxFCL::FArray2< Real > & lj_sigma,
		ObjexxFCL::FArray2< Real > & lj_r6_coeff,
		ObjexxFCL::FArray2< Real > & lj_r12_coeff,
		ObjexxFCL::FArray2< Real > & lj_switch_intercept,
		ObjexxFCL::FArray2< Real > & lj_switch_slope,
		ObjexxFCL::FArray1< Real > & lk_inv_lambda2,
		ObjexxFCL::FArray2< Real > & lk_coeff,
		ObjexxFCL::FArray2< Real > & lk_min_dis2sigma_value
	);

	void
	calc_etable_value(
		Real dis2,
		int atype1,
		int atype2,
		Real & atrE,
		Real & d_atrE,
		Real & repE,
		Real & d_repE,
		Real & solvE1,
		Real & solvE2,
		Real & dsolvE1,
		Real & dsolvE2
	) const;

	/// @brief Analyticaly calculate the Lazaridis-Karplus solvation energy and derivatives
	/// for a pair of atoms at a given distance
	void
	lk_solv_energy_and_deriv(
		int atype1,
		int atype2,
		Real dis,
		Real & solvE1,
		Real & solvE2,
		Real & dsolvE1,
		Real & dsolvE2
	) const;

	//void
	//zero_hydrogen_and_water_ljatr();

};

//////////////////// Inline Etable Evaluators /////////////////
inline
void
Etable::analytic_etable_evaluation(
	conformation::Atom const & at1,
	conformation::Atom const & at2,
	Real & lj_atrE,
	Real & lj_repE,
	Real & fa_solE,
	Real & dis2
) const
{
	using namespace core;
	using namespace core::scoring::etable;
	dis2 = at1.xyz().distance_squared( at2.xyz() );

	int atype1 = at1.type() < at2.type() ? at1.type() : at2.type();
	int atype2 = at1.type() < at2.type() ? at2.type() : at1.type();
	lj_atrE = 0; lj_repE = 0; fa_solE = 0;

	// locals
	Real ljE;

	Real atrE = 0.;
	Real repE = 0.;

	//  ctsa - epsilon allows final bin value to be calculated
	if ( dis2 > max_dis2 + epsilon ) return;

	EtableParamsOnePair const & p = analytic_params_for_pair( atype1, atype2 );
	if ( dis2 < min_dis2 ) dis2 = min_dis2;

	if ( dis2 > p.maxd2 ) return;

	//if ( p.hydrogen_interaction ) {
	//	analytic_etable_evaluation_H(
	//		at1, at2, p, dis2, lj_atrE, lj_repE, fa_solE );
	//	return;
	//}

	Real const dis = std::sqrt(dis2);
	Real const inv_dis = 1.0/dis;
	Real const inv_dis2 = inv_dis * inv_dis;

	if ( dis2 < p.ljrep_linear_ramp_d2_cutoff ) { // dis * p.inv_lj_sigma < lj_switch_dis2sigma ) {
		//  ctsa - use linear ramp instead of lj when the dis/sigma  ratio drops below theshold
		ljE = analytic_ljrep_linearized( dis, p );
	} else if ( dis < p.ljatr_spline_xlo ) {
		//  ctsa - calc regular lennard-jones
		ljE = analytic_lj_generic_form( dis2, inv_dis2, p ); //  p.lj_r12_coeff * inv_dis12 + p.lj_r6_coeff * inv_dis6;
	} else if ( dis < p.ljatr_spline_xhi ) {
		ljE = analytic_ljatr_spline_ramp_to_zero( dis, p );
	} else {
		/// assuming ljatr_spline_xhi == LK distance cutoff, since this will skip lk evaluation
		return;
	}

	/// Divvy up the lennard-jones energies into attractive and repulsive components;
	/// for most atom pairs, the attractive component goes smoothly to a constant,
	/// as the atoms approach, and then the repulsive component takes over from there.
	if ( p.ljrep_from_negcrossing ) {
		// only for the REPLS and HREPS atom types: start repelling when the lennard-jones term
		// goes from being attractive to repulsive.
		if (ljE < 0 ) atrE = ljE;
		else repE = ljE;
	} else if ( dis < p.lj_minimum ) {
		atrE = p.lj_val_at_minimum;
		repE = ljE - p.lj_val_at_minimum;
	} else {
		atrE = ljE;
	}

	/// Some atom pairs include extra repulsion that's modeled as a quadratic function
	/// in some region and then linearized outside of that region.  Specifically, this code
	/// exists for the OCbb / Ocbb interactions.
	ExtraQuadraticRepulsion const & exrep = p.ljrep_extra_repulsion;
	if ( dis < exrep.xhi ) {
		if ( dis < exrep.xlo ) {
			repE += ( dis - exrep.xlo ) * exrep.extrapolated_slope + exrep.ylo;
		} else {
			repE += ( exrep.xhi - dis ) * ( exrep.xhi - dis ) * exrep.slope;
		}
	}

	/// Now zero out the attractive energy if necessry; this is done for hydrogen interactions
	lj_atrE = atrE * p.ljatr_final_weight;
	lj_repE = repE;

	/// Now handle solvation.
	analytic_lk_evaluation( atype1, atype2, p, dis, inv_dis2, fa_solE );
}

inline
void
Etable::analytic_lk_energy(
	conformation::Atom const & at1,
	conformation::Atom const & at2,
	Real & fa_solE1,
	Real & fa_solE2
) const
{
	using namespace core;
	using namespace core::scoring::etable;

	fa_solE1 = fa_solE2 = 0.0;

	Real dis2 = at1.xyz().distance_squared( at2.xyz() );

	//  ctsa - epsilon allows final bin value to be calculated
	if ( dis2 > max_dis2 + epsilon ) return;

	int atype1 = at1.type() < at2.type() ? at1.type() : at2.type();
	int atype2 = at1.type() < at2.type() ? at2.type() : at1.type();
	bool inorder = at1.type() < at2.type(); // otherwise, atype1 and atype2 are swapped

	EtableParamsOnePair const & p = analytic_params_for_pair( atype1, atype2 );
	if ( dis2 < min_dis2 ) dis2 = min_dis2;

	Real const dis = std::sqrt(dis2);
	Real const inv_dis = 1.0/dis;
	Real const inv_dis2 = inv_dis * inv_dis;

	analytic_lk_evaluation_individual( atype1, atype2, p, dis, inv_dis2, inorder ? fa_solE1 : fa_solE2, inorder ? fa_solE2 : fa_solE1 );
}

inline
void
Etable::analytic_etable_derivatives(
	conformation::Atom const & at1,
	conformation::Atom const & at2,
	Real & dljatrE_ddis,
	Real & dljrepE_ddis,
	Real & dfasolE_ddis,
	Real & inv_d
) const
{
	using namespace core;
	using namespace core::scoring::etable;
	Real dis2 = at1.xyz().distance_squared( at2.xyz() );

	int const atype1 = at1.type() < at2.type() ? at1.type() : at2.type();
	int const atype2 = at1.type() < at2.type() ? at2.type() : at1.type();
	dljatrE_ddis = dljrepE_ddis = dfasolE_ddis = 0.0; inv_d = 1;

	// locals
	Real dljE(1), inv_dis6(1);

	//  ctsa - epsilon allows final bin value to be calculated
	if ( dis2 > max_dis2 + epsilon ) return;

	EtableParamsOnePair const & p = analytic_params_for_pair( atype1, atype2 );
	if ( dis2 < min_dis2 ) dis2 = min_dis2;

	if ( dis2 > p.maxd2 ) return;

	//if ( p.hydrogen_interaction ) {
	//	analytic_etable_evaluation_H(
	//		at1, at2, p, dis2, lj_atrE, lj_repE, fa_solE );
	//	return;
	//}

	Real const dis = std::sqrt(dis2);
	Real const inv_dis = 1.0/dis; inv_d = inv_dis;
	Real const inv_dis2 = inv_dis * inv_dis;

	if ( dis2 < p.ljrep_linear_ramp_d2_cutoff ) { // dis * p.inv_lj_sigma < lj_switch_dis2sigma ) {
		//  ctsa - use linear ramp instead of lj when the dis/sigma  ratio drops below theshold
		dljE = p.lj_switch_slope;
	} else if ( dis < p.ljatr_spline_xlo ) {
		//  ctsa - calc regular lennard-jones
		inv_dis6 = inv_dis2 * inv_dis2 * inv_dis2;
		Real const inv_dis7 = inv_dis * inv_dis6;

		dljE = inv_dis7 * ( -12.*p.lj_r12_coeff * inv_dis6 - 6.* p.lj_r6_coeff );
	} else if ( dis < p.ljatr_spline_xhi ) {
		dljE = analytic_ljatr_spline_ramp_to_zero_deriv( dis, p );
	} else {
		/// assuming ljatr_spline_xhi == LK distance cutoff, since this will skip lk evaluation
		return;
	}

	/// Divvy up the lennard-jones energies into attractive and repulsive components;
	/// for most atom pairs, the attractive component goes smoothly to a constant,
	/// as the atoms approach, and then the repulsive component takes over from there.
	if ( p.ljrep_from_negcrossing ) {
		// only for the REPLS and HREPS atom types: start repelling when the lennard-jones term
		// goes from being attractive to repulsive.
		// so, calculate the energy to decide whether it's attractive or repulsive
		if ( dis2 < p.ljrep_linear_ramp_d2_cutoff ||  inv_dis6*( p.lj_r12_coeff*inv_dis6 + p.lj_r6_coeff ) > 0.0 ) {
			dljrepE_ddis = dljE;
		}
	} else if ( dis < p.lj_minimum ) {
		dljatrE_ddis = 0;
		dljrepE_ddis = dljE;
	} else {
		dljatrE_ddis = dljE;
	}

	/// Some atom pairs include extra repulsion that's modeled as a quadratic function
	/// in some region and then linearized outside of that region.  Specifically, this code
	/// exists for the OCbb / Ocbb interactions.
	ExtraQuadraticRepulsion const & exrep = p.ljrep_extra_repulsion;
	if ( dis < exrep.xhi ) {
		if ( dis < exrep.xlo ) {
			dljrepE_ddis += exrep.extrapolated_slope;
		} else {
			dljrepE_ddis += -2 * ( exrep.xhi - dis ) * exrep.slope;
		}
	}
	/// Finally, zero out the attractive energy derivatives if necessry; this is done for hydrogen interactions
	dljatrE_ddis *= p.ljatr_final_weight;

	/// Now handle solvation.
	/// a) At distances below p.fasol_spline_close_start, the value of fasol is held constant.
	/// b) Then there's a spline to smooth between this constant region and the exponential region.
	/// c) Then the exponential region.
	/// d) Then the spline to smooth between the exponential region and where the energy goes to zero.
	if ( dis < p.fasol_spline_close_start ) {
		dfasolE_ddis = 0;
	} else if ( dis < p.fasol_spline_close_end ) {
		dfasolE_ddis = spline_deriv( dis, p.fasol_spline_close_start, p.fasol_spline_close_end, p.fasol_spline_close );
		dfasolE_ddis *= p.fasol_final_weight;
		//Real const fasol_spline_knots_diff = p.fasol_spline_close_end - p.fasol_spline_close_start;
		//Real const fasol_spline_knots_diff_inv = 1/fasol_spline_knots_diff;
		//Real const a = ( p.fasol_spline_close_end - dis ) * fasol_spline_knots_diff_inv;
		//Real const b = ( dis - p.fasol_spline_close_start ) * fasol_spline_knots_diff_inv;
		//SplineParameters const & sp = p.fasol_spline_close;
		//fa_solE = a*sp.ylo + b*sp.yhi + ((a*a*a-a)*sp.y2lo + (b*b*b-b)*sp.y2hi)*fasol_spline_knots_diff*fasol_spline_knots_diff / 6;
	} else if ( dis < fasol_spline_far_xlo ) {
		/// exponential evaluation
		/// assert( atype1 <= atype2 ), which is accomplished at the top of this function.
		Real const dis_rad1 = dis - lj_radius(atype1);
		Real const x1 = ( dis_rad1 * dis_rad1 ) * lk_inv_lambda2_(atype1);
		Real const dis_rad2 = dis - lj_radius(atype2);
		Real const x2 = ( dis_rad2 * dis_rad2 ) * lk_inv_lambda2_(atype2);

		Real const solvE1 = std::exp(-x1) * p.lk_coeff1 * inv_dis2;
		Real const solvE2 = std::exp(-x2) * p.lk_coeff2 * inv_dis2;
		dfasolE_ddis = -2 * ( ( dis_rad1 * lk_inv_lambda2_(atype1) + inv_dis ) * solvE1 + ( dis_rad2 * lk_inv_lambda2_(atype2) + inv_dis ) * solvE2 );
		dfasolE_ddis *= p.fasol_final_weight;

	} else if ( dis < fasol_spline_far_xhi ) {
		dfasolE_ddis = spline_deriv( dis,
			fasol_spline_far_xlo, fasol_spline_far_xhi,
			fasol_spline_far_diff_xhi_xlo, fasol_spline_far_diff_xhi_xlo_inv,
			p.fasol_spline_far ) *
			p.fasol_final_weight;
	} else {
		dfasolE_ddis = 0;
	}
}

inline
void
Etable::analytic_lk_derivatives(
	conformation::Atom const & at1,
	conformation::Atom const & at2,
	Real & dfasolE1_ddis,
	Real & dfasolE2_ddis,
	Real & inv_d
) const
{
	using namespace core;
	using namespace core::scoring::etable;
	Real dis2 = at1.xyz().distance_squared( at2.xyz() );

	bool const inorder = at1.type() < at2.type(); // otherwise, atype1 and atype2 are swapped
	int const atype1 = inorder ? at1.type() : at2.type();
	int const atype2 = inorder ? at2.type() : at1.type();
	Real & dfasol1 = inorder ? dfasolE1_ddis : dfasolE2_ddis;
	Real & dfasol2 = inorder ? dfasolE2_ddis : dfasolE1_ddis;
	inv_d = 1;

	//  ctsa - epsilon allows final bin value to be calculated
	if ( dis2 > max_dis2 + epsilon ) return;

	EtableParamsOnePair const & p = analytic_params_for_pair( atype1, atype2 );
	if ( dis2 < min_dis2 ) dis2 = min_dis2;

	if ( dis2 > p.maxd2 ) return;
	Real const dis = std::sqrt(dis2);
	Real const inv_dis = 1.0/dis; inv_d = inv_dis;
	Real const inv_dis2 = inv_dis * inv_dis;

	if ( dis < p.fasol_spline_close_start ) {
		dfasol1 = 0;
		dfasol2 = 0;
	} else if ( dis < p.fasol_spline_close_end ) {
		dfasol1 = spline_deriv( dis, p.fasol_spline_close_start, p.fasol_spline_close_end, p.fasol_spline1_close );
		dfasol1 *= p.fasol_final_weight;
		dfasol2 = spline_deriv( dis, p.fasol_spline_close_start, p.fasol_spline_close_end, p.fasol_spline2_close );
		dfasol2 *= p.fasol_final_weight;
	} else if ( dis < fasol_spline_far_xlo ) {
		/// exponential evaluation
		Real const dis_rad1 = dis - lj_radius(atype1);
		Real const x1 = ( dis_rad1 * dis_rad1 ) * lk_inv_lambda2_(atype1);
		Real const dis_rad2 = dis - lj_radius(atype2);
		Real const x2 = ( dis_rad2 * dis_rad2 ) * lk_inv_lambda2_(atype2);

		Real const solvE1 = std::exp(-x1) * p.lk_coeff1 * inv_dis2;
		Real const solvE2 = std::exp(-x2) * p.lk_coeff2 * inv_dis2;
		dfasol1 = -2 * p.fasol_final_weight * ( dis_rad1 * lk_inv_lambda2_(atype1) + inv_dis ) * solvE1;
		dfasol2 = -2 * p.fasol_final_weight * ( dis_rad2 * lk_inv_lambda2_(atype2) + inv_dis ) * solvE2;

	} else if ( dis < fasol_spline_far_xhi ) {
		dfasol1 = spline_deriv( dis, fasol_spline_far_xlo, fasol_spline_far_xhi,
			fasol_spline_far_diff_xhi_xlo, fasol_spline_far_diff_xhi_xlo_inv,
			p.fasol_spline1_far ) * p.fasol_final_weight;
		dfasol2 = spline_deriv( dis, fasol_spline_far_xlo, fasol_spline_far_xhi,
			fasol_spline_far_diff_xhi_xlo, fasol_spline_far_diff_xhi_xlo_inv,
			p.fasol_spline2_far ) * p.fasol_final_weight;
	} else {
		dfasol1 = 0;
		dfasol2 = 0;
	}

}

//void
//Etable::analytic_etable_evaluation_H(
//	conformation::Atom const & at1,
//	conformation::Atom const & at2,
//	EtableParamsOnePair const & p,
//	Real const dis2,
//	Real & lj_atrE,
//	Real & lj_repE,
//	Real & fa_solE
//) const
//{
//	assert( dis2 <= p.maxd2 );
//  if ( dis2  < p.ljrep_linear_ramp_d2_cutoff ) {
//    //  ctsa - use linear ramp instead of lj when the dis/sigma
//    //    ratio drops below theshold
//		Real dis = std::sqrt( dis2 );
//		lj_repE = analytic_ljrep_linearized( dis, p ) - p.lj_val_at_minimum;
//  } else {
//    lj_repE = analytic_lj_generic_form( 1/dis2, p ) - p.lj_val_at_minimum;
//  }
//}

/// @details only to be called when the distance, dis, is less than the switch-point for
/// the lj_switch_dis2sigma value.
inline
Real
Etable::analytic_ljrep_linearized(
	Real const dis,
	EtableParamsOnePair const & p
) const
{
	assert( dis * dis < p.ljrep_linear_ramp_d2_cutoff );
	return  dis*p.lj_switch_slope + p.lj_switch_intercept;
}

/// @details: evaluate the standard Lennard-Jones 6-12 functional form.
/// Only call this function if the square distance is in the range
/// sqrt( p.ljrep_linear_ramp_d2_cutoff ) < dis < p.ljatr_spline_xhi
inline
Real
Etable::analytic_lj_generic_form(
	Real const ASSERT_ONLY( dis2 ),
	Real const inv_dis2,
	EtableParamsOnePair const & p
) const
{
	assert( dis2 >= p.ljrep_linear_ramp_d2_cutoff );
	assert( dis2 <= p.ljatr_spline_xhi * p.ljatr_spline_xhi );
	Real const inv_dis6  = inv_dis2 * inv_dis2 * inv_dis2;
	//Real const inv_dis12 = inv_dis6 * inv_dis6;

	return ( p.lj_r12_coeff * inv_dis6 + p.lj_r6_coeff ) * inv_dis6;
}

inline
Real
Etable::eval_spline(
	Real const x,
	Real const xlo,
	Real const xhi,
	SplineParameters const & sp
)
{
	assert( x >= xlo && x <= xhi );
	Real width = xhi - xlo;
	return eval_spline( x, xlo, xhi, width, 1/width, sp );
}

inline
Real
Etable::spline_deriv(
	Real const x,
	Real const xlo,
	Real const xhi,
	SplineParameters const & sp
)
{
	assert( x >= xlo && x <= xhi );
	Real width = xhi - xlo;
	return spline_deriv( x, xlo, xhi, width, 1/width, sp );
}

inline
Real
Etable::eval_spline(
	Real const x,
	Real const xlo,
	Real const xhi,
	Real const width,
	Real const invwidth,
	SplineParameters const & sp
)
{
	Real a = ( xhi - x ) * invwidth;
	Real b = ( x - xlo ) * invwidth;
	return a*sp.ylo + b*sp.yhi + ((a*a*a-a)*sp.y2lo + (b*b*b-b)*sp.y2hi)*width*width / 6;
}

inline
Real
Etable::spline_deriv(
	Real const x,
	Real const xlo,
	Real const xhi,
	Real const width,
	Real const invwidth,
	SplineParameters const & sp
)
{
	Real a = ( xhi - x ) * invwidth;
	Real b = ( x - xlo ) * invwidth;
	return -1*invwidth*sp.ylo + invwidth*sp.yhi + ((1 - 3*a*a )*sp.y2lo + (3*b*b - 1 )*sp.y2hi)*width / 6;
}

/// @details: evaluate the attractive component of the LJ term as it
/// ramps to zero.
/// Only call this function if the square distance is in the range
/// p.ljatr_spline_xlo < dis < p.ljatr_spline_xhi
inline
Real
Etable::analytic_ljatr_spline_ramp_to_zero(
	Real const dis,
	EtableParamsOnePair const & p
) const
{
	assert( dis >= p.ljatr_spline_xlo );
	assert( dis <= p.ljatr_spline_xhi );
	return eval_spline( dis, p.ljatr_spline_xlo, p.ljatr_spline_xhi, p.ljatr_spline_parameters );
}

/// @details: evaluate the derivative for the attractive component
/// of the LJ term as it ramps to zero.
/// Only call this function if the square distance is in the range
/// p.ljatr_spline_xlo < dis < p.ljatr_spline_xhi
inline
Real
Etable::analytic_ljatr_spline_ramp_to_zero_deriv(
	Real const dis,
	EtableParamsOnePair const & p
) const
{
	assert( dis >= p.ljatr_spline_xlo );
	assert( dis <= p.ljatr_spline_xhi );
	return spline_deriv( dis, p.ljatr_spline_xlo, p.ljatr_spline_xhi, p.ljatr_spline_parameters );
}

/// @details Note that atype1 must be less than or equal to atype2; the
/// EtableParamsOnePair are stored only once for atom types i and j.
void
Etable::analytic_lk_evaluation(
	int const atype1,
	int const atype2,
	EtableParamsOnePair const & p,
	Real const dis,
	Real const inv_dis2,
	Real & fa_solE
) const
{
	assert( atype1 <= atype2 );

	/// a) At distances below p.fasol_spline_close_start, the value of fasol is held constant.
	/// b) Then there's a spline to smooth between this constant region and the exponential region.
	/// c) Then the exponential region.
	/// d) Then the spline to smooth between the exponential region and where the energy goes to zero.
	if ( dis < p.fasol_spline_close_start ) {
		fa_solE = p.fasol_spline_close.ylo * p.fasol_final_weight;
	} else if ( dis < p.fasol_spline_close_end ) {
		fa_solE = eval_spline( dis, p.fasol_spline_close_start, p.fasol_spline_close_end, p.fasol_spline_close );
		fa_solE *= p.fasol_final_weight;
		//Real const fasol_spline_knots_diff = p.fasol_spline_close_end - p.fasol_spline_close_start;
		//Real const fasol_spline_knots_diff_inv = 1/fasol_spline_knots_diff;
		//Real const a = ( p.fasol_spline_close_end - dis ) * fasol_spline_knots_diff_inv;
		//Real const b = ( dis - p.fasol_spline_close_start ) * fasol_spline_knots_diff_inv;
		//SplineParameters const & sp = p.fasol_spline_close;
		//fa_solE = a*sp.ylo + b*sp.yhi + ((a*a*a-a)*sp.y2lo + (b*b*b-b)*sp.y2hi)*fasol_spline_knots_diff*fasol_spline_knots_diff / 6;
	} else if ( dis < fasol_spline_far_xlo ) {
		/// exponential evaluation
		/// assert( atype1 <= atype2 ), which is accomplished at the top of this function.
		Real const dis_rad1 = dis - lj_radius(atype1);
		Real const x1 = ( dis_rad1 * dis_rad1 ) * lk_inv_lambda2_(atype1);
		Real const dis_rad2 = dis - lj_radius(atype2);
		Real const x2 = ( dis_rad2 * dis_rad2 ) * lk_inv_lambda2_(atype2);

		fa_solE =  inv_dis2 * ( std::exp(-x1) * p.lk_coeff1 + std::exp(-x2) * p.lk_coeff2 );
		fa_solE *= p.fasol_final_weight;

	} else if ( dis < fasol_spline_far_xhi ) {
		fa_solE = eval_spline( dis,
			fasol_spline_far_xlo, fasol_spline_far_xhi,
			fasol_spline_far_diff_xhi_xlo, fasol_spline_far_diff_xhi_xlo_inv,
			p.fasol_spline_far );
		//Real const a = ( fasol_spline_far_xhi - dis ) * fasol_spline_far_diff_xhi_xlo_inv;
		//Real const b = ( dis - fasol_spline_far_xlo ) * fasol_spline_far_diff_xhi_xlo_inv;
		//SplineParameters const & sp = p.fasol_spline_far;
		//
		//fa_solE = a*sp.ylo + b*sp.yhi + ((a*a*a-a)*sp.y2lo + (b*b*b-b)*sp.y2hi)*fasol_spline_far_diff_xhi_xlo*fasol_spline_far_diff_xhi_xlo / 6;
		fa_solE *= p.fasol_final_weight;
	} else {
		fa_solE = 0;
	}
}

/// @details The etable parameters are stored such that the atom type of
/// atom 1 is less than or equal to the atom type of atom 2, so that
/// the caller of this function will need to keep track of which atom
/// is which.
void
Etable::analytic_lk_evaluation_individual(
	int const atype1,
	int const atype2,
	EtableParamsOnePair const & p,
	Real const dis,
	Real const inv_dis2,
	Real & fa_solE1,
	Real & fa_solE2
) const
{
	assert( atype1 <= atype2 );

	/// a) At distances below p.fasol_spline_close_start, the value of fasol is held constant.
	/// b) Then there's a spline to smooth between this constant region and the exponential region.
	/// c) Then the exponential region.
	/// d) Then the spline to smooth between the exponential region and where the energy goes to zero.
	if ( dis < p.fasol_spline_close_start ) {
		fa_solE1 = p.fasol_spline1_close.ylo * p.fasol_final_weight;
		fa_solE2 = p.fasol_spline2_close.ylo * p.fasol_final_weight;
	} else if ( dis < p.fasol_spline_close_end ) {
		fa_solE1 = p.fasol_final_weight * eval_spline( dis, p.fasol_spline_close_start, p.fasol_spline_close_end, p.fasol_spline1_close );
		fa_solE2 = p.fasol_final_weight * eval_spline( dis, p.fasol_spline_close_start, p.fasol_spline_close_end, p.fasol_spline2_close );
	} else if ( dis < fasol_spline_far_xlo ) {
		/// exponential evaluation
		Real const dis_rad1 = dis - lj_radius(atype1);
		Real const x1 = ( dis_rad1 * dis_rad1 ) * lk_inv_lambda2_(atype1);
		Real const dis_rad2 = dis - lj_radius(atype2);
		Real const x2 = ( dis_rad2 * dis_rad2 ) * lk_inv_lambda2_(atype2);

		fa_solE1 = p.fasol_final_weight * inv_dis2 * std::exp(-x1) * p.lk_coeff1;
		fa_solE2 = p.fasol_final_weight * inv_dis2 * std::exp(-x2) * p.lk_coeff2;

	} else if ( dis < fasol_spline_far_xhi ) {
		fa_solE1 = p.fasol_final_weight * eval_spline( dis,
			fasol_spline_far_xlo, fasol_spline_far_xhi,
			fasol_spline_far_diff_xhi_xlo, fasol_spline_far_diff_xhi_xlo_inv,
			p.fasol_spline1_far );

		fa_solE2 = p.fasol_final_weight * eval_spline( dis,
			fasol_spline_far_xlo, fasol_spline_far_xhi,
			fasol_spline_far_diff_xhi_xlo, fasol_spline_far_diff_xhi_xlo_inv,
			p.fasol_spline2_far );

	} else {
		fa_solE1 = fa_solE2 = 0;
	}
}

inline
EtableParamsOnePair const &
Etable::analytic_params_for_pair(
	Size atype1,
	Size atype2
) const
{
	Size i1 = atype1 < atype2 ? atype1 : atype2;
	Size i2 = atype1 < atype2 ? atype2 : atype1;
	Size index = (i1-1)*n_atomtypes_ + i2 - (i1*(i1-1)/2 );
	return analytic_parameters[ index ];
}

inline
EtableParamsOnePair &
Etable::analytic_params_for_pair(
	Size atype1,
	Size atype2
)
{
	Size i1 = atype1 < atype2 ? atype1 : atype2;
	Size i2 = atype1 < atype2 ? atype2 : atype1;
	Size index = (i1-1)*n_atomtypes_ + i2 - (i1*(i1-1)/2 );
	return analytic_parameters[ index ];
}



} // etable
} // scoring
} // core

#endif
