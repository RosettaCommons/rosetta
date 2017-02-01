// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.
//////////////////////////////////////////////////////////////////////
///
/// @brief
/// A class for generating the table for fa_atr/rep and fa_sol
///
/// @details
/// This class is called upon by the ScoringManager. Since actual calculating of the LJ potential
/// is time consuming if done multiple times, this class precomputes and discritizes the potential
/// (meaning that the potential is broken down into bins). Once the bins have been created, it will
/// smooth out the bins, for better interpolation.
///
///
/// @author
/// Phil Bradley
/// Andrew Leaver-Fay
/// Steven Combs - comments and skipping of virtual atoms
///
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

struct CubicPolynomial {
	Real c0, c1, c2, c3;
	CubicPolynomial() : c0(0), c1(0), c2(0), c3(0) {}
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
	Real ljatr_cubic_poly_xlo;
	Real ljatr_cubic_poly_xhi;
	CubicPolynomial ljatr_cubic_poly_parameters;
	ExtraQuadraticRepulsion ljrep_extra_repulsion;
	bool ljrep_from_negcrossing;

	Real lk_coeff1;
	Real lk_coeff2;
	Real lk_min_dis2sigma_value;
	Real fasol_cubic_poly_close_start;
	Real fasol_cubic_poly_close_end;

	CubicPolynomial fasol_cubic_poly_close;  // mututal desolvation of atoms 1 and 2
	Real fasol_cubic_poly_close_flat; // the fixed value used for distances beneath fasol_cubic_poly_close_start

	CubicPolynomial fasol_cubic_poly_far;    // mututal desolvation of atoms 1 and 2

	CubicPolynomial fasol_cubic_poly1_close; // desolvation of atom 1 by atom 2
	Real fasol_cubic_poly1_close_flat;       // the fixed value used for distances beneath fasol_cubic_poly_close_start
	CubicPolynomial fasol_cubic_poly1_far;   // desolvation of atom 1 by atom 2

	CubicPolynomial fasol_cubic_poly2_close; // desolvation of atom 2 by atom 1
	Real fasol_cubic_poly2_close_flat;       // the fixed value used for distances beneath fasol_cubic_poly_close_start
	CubicPolynomial fasol_cubic_poly2_far;   // desolvation of atom 2 by atom 1
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
		ljatr_cubic_poly_xlo(0.0),
		ljatr_cubic_poly_xhi(0.0),
		ljrep_from_negcrossing(false),
		lk_coeff1(0.0),
		lk_coeff2(0.0),
		lk_min_dis2sigma_value(0.0),
		fasol_cubic_poly_close_start(0.0),
		fasol_cubic_poly_close_end(0.0),
		fasol_cubic_poly_close_flat(0.0),
		fasol_cubic_poly1_close_flat(0.0),
		fasol_cubic_poly2_close_flat(0.0),
		ljatr_final_weight(1.0),
		fasol_final_weight(1.0)
	{}


};


/// @brief Class definition for Etable
class Etable : public utility::pointer::ReferenceCount {

public:
	/// @brief Automatically generated virtual destructor for class deriving directly from ReferenceCount
	virtual ~Etable();

	///  constructor
	Etable(
		chemical::AtomTypeSetCAP atom_set_in, // like etable namespace
		EtableOptions const & options,
		std::string const & alternate_parameter_set = ""
	);

	// Const Access methods for constants
	int
	n_atomtypes() const {
		return n_atomtypes_;
	}

	Real
	Wradius() const {
		return Wradius_;
	}

	Real
	lj_switch_dis2sigma() const {
		return lj_switch_dis2sigma_;
	}

	Real
	etable_disbins() const {
		return etable_disbins_;
	}


	bool
	lj_use_lj_deriv_slope() const {
		return lj_use_lj_deriv_slope_;
	}

	Real
	lj_slope_intercept() const {
		return lj_slope_intercept_;
	}

	bool
	lj_use_hbond_radii() const {
		return lj_use_hbond_radii_;
	}

	bool
	lj_use_water_radii() const {
		return lj_use_water_radii_;
	}

	Real
	lj_water_dis() const {
		return lj_water_dis_;
	}

	Real
	lj_water_hdis() const {
		return lj_water_hdis_;
	}

	Real
	lk_min_dis2sigma() const {
		return lk_min_dis2sigma_;
	}

	Real
	min_dis() const {
		return min_dis_;
	}

	Real
	min_dis2() const {
		return min_dis2_;
	}

	bool
	add_long_range_damping() const {
		return add_long_range_damping_;
	}

	Real
	long_range_damping_length() const {
		return long_range_damping_length_;
	}

	Real
	epsilon() const {
		return epsilon_;
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

	virtual
	ObjexxFCL::FArray3D< Real > const &
	solv1() const
	{
		return solv1_;
	}

	virtual
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
	virtual
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
	max_dis2() const
	{
		return max_dis2_;
	}

	Real
	get_safe_max_dis2() const
	{
		return safe_max_dis2_;
	}

	int
	get_bins_per_A2() const
	{
		return bins_per_A2_;
	}

	chemical::AtomTypeSetCAP
	atom_set() const
	{
		return atom_set_;
	}

	/// @brief Do hydrogens provide attractive forces or do they only repell?
	bool fa_hatr() const { return fa_hatr_; }

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

	/// @brief The square distance cutoff beyond which any pair of heavy-atoms is
	/// guaranteed to have an interaction energy of zero.  This function is
	/// used by the NeighborList
	virtual
	Real
	nblist_dis2_cutoff_XX() const
	{
		return nblist_dis2_cutoff_XX_;
	}

	/// @brief The square distance cutoff beyond which a hydrogen/heavy-atom pair is
	/// guaranteed to have an interaction energy of zero.  This function is used
	/// by the NeighborList
	virtual
	Real
	nblist_dis2_cutoff_XH() const
	{
		return nblist_dis2_cutoff_XH_;
	}

	/// @brief The square distance cutoff beyond which any hydrogen/hydrogen pair is guaranteed
	/// to have an interaction energy of zero.  This function is used by the NeighborList
	virtual
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
	virtual
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
	virtual
	Real
	lk_dgfree( int const i ) const
	{
		return lk_dgfree_[i];
	}

	/// @brief Return the Lazaridis Karplus volume for an atom
	virtual
	Real
	lk_volume( int const i ) const
	{
		return lk_volume_[i];
	}

	/// @brief Return the Lazaridis Karplus "lambda" value (correlation distance) for an atom
	virtual
	Real
	lk_lambda( int const i ) const
	{
		return lk_lambda_[i];
	}

	Real
	lk_inv_lambda2( int const i ) const {
		return lk_inv_lambda2_[i];
	}

	// get the AtomType object corresponding to a give type index
	chemical::AtomType const &
	atom_type( int const type )
	{
		chemical::AtomTypeSetCOP atom_set( atom_set_ );
		return (*atom_set)[ type ];
	}

	inline Real fasol_cubic_poly_far_xlo() const { return fasol_cubic_poly_far_xlo_; }
	inline Real fasol_cubic_poly_far_xhi() const { return fasol_cubic_poly_far_xhi_; }

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
	// conformation::Atom const & at1,
	// conformation::Atom const & at2,
	// EtableParamsOnePair const & p,
	// Real const dis2,
	// Real & lj_atrE,
	// Real & lj_repE,
	// Real & fa_solE
	//) const;

private:
	// initialization routines

	void dimension_etable_arrays();
	void initialize_from_input_atomset( chemical::AtomTypeSetCAP atom_set_in );
	void calculate_nblist_distance_thresholds( EtableOptions const & options );
	void
	read_alternate_parameter_set(
		chemical::AtomTypeSetCAP atom_set_in,
		std::string const & alternate_parameter_set
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
	analytic_ljatr_cubic_poly_ramp_to_zero(
		Real const dis,
		EtableParamsOnePair const & p
	) const;

	inline
	Real
	analytic_ljatr_cubic_poly_ramp_to_zero_deriv(
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
	CubicPolynomial
	cubic_polynomial_from_spline( Real xlo, Real xhi, SplineParameters const & sp );

	static
	inline
	Real
	eval_cubic_polynomial(
		Real const x,
		CubicPolynomial const & sp
	);

	static
	inline
	Real
	cubic_polynomial_deriv(
		Real const x,
		CubicPolynomial const & cp
	);

	inline
	EtableParamsOnePair const &
	analytic_params_for_pair(
		Size atype1,
		Size atype2
	) const;

public: // Interfaces for convenient IO

	virtual
	void
	output_etable(
		ObjexxFCL::FArray3D<Real> & etable,
		std::string label,
		std::ostream & out
	);

	virtual
	void
	input_etable(
		ObjexxFCL::FArray3D<Real> & etable,
		std::string const & label,
		std::istream & in
	);

private:

	inline
	EtableParamsOnePair &
	analytic_params_for_pair(
		Size atype1,
		Size atype2
	);

public:
	Real get_lj_hbond_OH_donor_dis() const { return lj_hbond_OH_donor_dis_; }
	Real get_lj_hbond_hdis() const { return lj_hbond_hdis_; }
	Real get_lj_hbond_dis() const { return lj_hbond_dis_; }

private:

	chemical::AtomTypeSetCAP atom_set_;

	// parameters:
	int const n_atomtypes_;

	// from options
	Real const max_dis_;
	int const bins_per_A2_;
	Real const Wradius_; // global mod to radii
	Real const lj_switch_dis2sigma_; // actual value used for switch
	Real const max_dis2_;
	int const etable_disbins_;
	Real const lj_hbond_OH_donor_dis_;
	Real const lj_hbond_dis_;

	// hard-coded for now
	bool const lj_use_lj_deriv_slope_;
	Real const lj_slope_intercept_;
	bool const lj_use_hbond_radii_;
	Real const lj_hbond_hdis_;
	bool const lj_use_water_radii_;
	Real const lj_water_dis_;
	Real const lj_water_hdis_;
	bool const enlarge_h_lj_wdepth_;
	Real const lk_min_dis2sigma_;
	bool const no_lk_polar_desolvation_;
	bool const proline_N_is_lk_nonpolar_;
	Real const min_dis_;
	Real const min_dis2_; // was double
	bool const add_long_range_damping_;
	Real const long_range_damping_length_;
	Real const epsilon_;
	Real const safe_max_dis2_;
	Real hydrogen_interaction_cutoff2_;
	Real max_heavy_hydrogen_cutoff_;
	Real max_hydrogen_hydrogen_cutoff_;
	Real nblist_dis2_cutoff_XX_; // for use by the old-style neighborlist
	Real nblist_dis2_cutoff_XH_; // for use by the old-style neighborlist
	Real nblist_dis2_cutoff_HH_; // for use by the old-style neighborlist
	Real max_non_hydrogen_lj_radius_;
	Real max_hydrogen_lj_radius_;
	bool fa_hatr_; // do hydrogens provide attractive forces, or are they for repulsion only?

	// these three derived from other data
	Real lj_switch_sigma2dis_;
	Real lj_switch_value2wdepth_;
	Real lj_switch_slope_sigma2wdepth_;

	ObjexxFCL::FArray2D< Real > lj_sigma_;
	ObjexxFCL::FArray2D< Real > lj_r6_coeff_;
	ObjexxFCL::FArray2D< Real > lj_r12_coeff_;
	ObjexxFCL::FArray2D< Real > lj_switch_intercept_;
	ObjexxFCL::FArray2D< Real > lj_switch_slope_;
	ObjexxFCL::FArray1D< Real > lk_inv_lambda2_;
	ObjexxFCL::FArray2D< Real > lk_coeff_;
	ObjexxFCL::FArray2D< Real > lk_min_dis2sigma_value_;

	/// Data needed to describe the splines for the fasol term:
	/// There are two cubic_polys needed: one to smooth the transition to where the fasol term goes flat
	/// as the distance becomes less than the sum of the van der Waals radii (the "close" cubic_poly),
	/// and a second to smooth the transition where the distance goes to the
	/// fa_max_dis and the energy goes to 0 (the far cubic_poly).  In between the close
	/// and far values, the exponential form of the energy function is used.
	Real fasol_cubic_poly_far_xlo_;
	Real fasol_cubic_poly_far_xhi_;

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

	utility::vector1< Size > carbon_types_; // indices for the standard carbon types

	// upper-triangle of the n_atomtypes x n_atomtypes table.  N^2 - (N*(N-1))/2) entries.
	// indexed for at1, at2 w/ at1<at2:
	// ((at1-1)*n_atomtypes + (at2-1) + 1 - (at1*(at1-1)/2)
	utility::vector1< EtableParamsOnePair > analytic_parameters_;

	bool slim_;

	// private methods

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
	if ( dis2 > max_dis2_ + epsilon_ ) return;

	EtableParamsOnePair const & p = analytic_params_for_pair( atype1, atype2 );
	if ( dis2 < min_dis2_ ) dis2 = min_dis2_;

	if ( dis2 > p.maxd2 ) return;

	Real const dis = std::sqrt(dis2);
	Real const inv_dis = 1.0/dis;
	Real const inv_dis2 = inv_dis * inv_dis;

	if ( dis2 < p.ljrep_linear_ramp_d2_cutoff ) { // dis * p.inv_lj_sigma < lj_switch_dis2sigma ) {
		//  ctsa - use linear ramp instead of lj when the dis/sigma  ratio drops below theshold
		ljE = analytic_ljrep_linearized( dis, p );
	} else if ( dis < p.ljatr_cubic_poly_xlo ) {
		//  ctsa - calc regular lennard-jones
		ljE = analytic_lj_generic_form( dis2, inv_dis2, p ); //  p.lj_r12_coeff * inv_dis12 + p.lj_r6_coeff * inv_dis6;
	} else if ( dis < p.ljatr_cubic_poly_xhi ) {
		ljE = analytic_ljatr_cubic_poly_ramp_to_zero( dis, p );
	} else {
		/// assuming ljatr_cubic_poly_xhi == LK distance cutoff, since this will skip lk evaluation
		return;
	}

	/// Divvy up the lennard-jones energies into attractive and repulsive components;
	/// for most atom pairs, the attractive component goes smoothly to a constant,
	/// as the atoms approach, and then the repulsive component takes over from there.
	if ( p.ljrep_from_negcrossing ) {
		// only for the REPLS and HREPS atom types: start repelling when the lennard-jones term
		// goes from being attractive to repulsive.
		if ( ljE < 0 ) atrE = ljE;
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
	if ( dis2 > max_dis2_ + epsilon_ ) return;

	int atype1 = at1.type() < at2.type() ? at1.type() : at2.type();
	int atype2 = at1.type() < at2.type() ? at2.type() : at1.type();
	bool inorder = at1.type() < at2.type(); // otherwise, atype1 and atype2 are swapped

	EtableParamsOnePair const & p = analytic_params_for_pair( atype1, atype2 );
	if ( dis2 < min_dis2_ ) dis2 = min_dis2_;

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

	//  ctsa - epsilon allows final_ bin value to be calculated
	if ( dis2 > max_dis2_ + epsilon_ ) return;

	EtableParamsOnePair const & p = analytic_params_for_pair( atype1, atype2 );
	if ( dis2 < min_dis2_ ) dis2 = min_dis2_;

	if ( dis2 > p.maxd2 ) return;

	Real const dis = std::sqrt(dis2);
	Real const inv_dis = 1.0/dis; inv_d = inv_dis;
	Real const inv_dis2 = inv_dis * inv_dis;

	if ( dis2 < p.ljrep_linear_ramp_d2_cutoff ) { // dis * p.inv_lj_sigma < lj_switch_dis2sigma ) {
		//  ctsa - use linear ramp instead of lj when the dis/sigma  ratio drops below theshold
		dljE = p.lj_switch_slope;
	} else if ( dis < p.ljatr_cubic_poly_xlo ) {
		//  ctsa - calc regular lennard-jones
		inv_dis6 = inv_dis2 * inv_dis2 * inv_dis2;
		Real const inv_dis7 = inv_dis * inv_dis6;

		dljE = inv_dis7 * ( -12.*p.lj_r12_coeff * inv_dis6 - 6.* p.lj_r6_coeff );
	} else if ( dis < p.ljatr_cubic_poly_xhi ) {
		dljE = analytic_ljatr_cubic_poly_ramp_to_zero_deriv( dis, p );
	} else {
		/// assuming ljatr_cubic_poly_xhi == LK distance cutoff, since this will skip lk evaluation
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
	/// a) At distances below p.fasol_cubic_poly_close_start, the value of fasol is held constant.
	/// b) Then there's a cubic_poly to smooth between this constant region and the exponential region.
	/// c) Then the exponential region.
	/// d) Then the cubic_poly to smooth between the exponential region and where the energy goes to zero.
	if ( dis < p.fasol_cubic_poly_close_start ) {
		dfasolE_ddis = 0;
	} else if ( dis < p.fasol_cubic_poly_close_end ) {
		dfasolE_ddis = cubic_polynomial_deriv( dis, p.fasol_cubic_poly_close );
		dfasolE_ddis *= p.fasol_final_weight;
	} else if ( dis < fasol_cubic_poly_far_xlo_ ) {
		/// exponential evaluation
		///debug_assert( atype1 <= atype2 ), which is accomplished at the top of this function.
		Real const dis_rad1 = dis - lj_radius(atype1);
		Real const x1 = ( dis_rad1 * dis_rad1 ) * lk_inv_lambda2_(atype1);
		Real const dis_rad2 = dis - lj_radius(atype2);
		Real const x2 = ( dis_rad2 * dis_rad2 ) * lk_inv_lambda2_(atype2);

		Real const solvE1 = std::exp(-x1) * p.lk_coeff1 * inv_dis2;
		Real const solvE2 = std::exp(-x2) * p.lk_coeff2 * inv_dis2;
		dfasolE_ddis = -2 * ( ( dis_rad1 * lk_inv_lambda2_(atype1) + inv_dis ) * solvE1 + ( dis_rad2 * lk_inv_lambda2_(atype2) + inv_dis ) * solvE2 );
		dfasolE_ddis *= p.fasol_final_weight;

	} else if ( dis < fasol_cubic_poly_far_xhi_ ) {
		dfasolE_ddis = cubic_polynomial_deriv( dis, p.fasol_cubic_poly_far ) * p.fasol_final_weight;
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
	if ( dis2 > max_dis2_ + epsilon_ ) return;

	EtableParamsOnePair const & p = analytic_params_for_pair( atype1, atype2 );
	if ( dis2 < min_dis2_ ) dis2 = min_dis2_;

	if ( dis2 > p.maxd2 ) return;
	Real const dis = std::sqrt(dis2);
	Real const inv_dis = 1.0/dis; inv_d = inv_dis;
	Real const inv_dis2 = inv_dis * inv_dis;

	if ( dis < p.fasol_cubic_poly_close_start ) {
		dfasol1 = 0;
		dfasol2 = 0;
	} else if ( dis < p.fasol_cubic_poly_close_end ) {
		dfasol1 = cubic_polynomial_deriv( dis, p.fasol_cubic_poly1_close ) * p.fasol_final_weight;
		dfasol2 = cubic_polynomial_deriv( dis, p.fasol_cubic_poly2_close ) * p.fasol_final_weight;
	} else if ( dis < fasol_cubic_poly_far_xlo_ ) {
		/// exponential evaluation
		Real const dis_rad1 = dis - lj_radius(atype1);
		Real const x1 = ( dis_rad1 * dis_rad1 ) * lk_inv_lambda2_(atype1);
		Real const dis_rad2 = dis - lj_radius(atype2);
		Real const x2 = ( dis_rad2 * dis_rad2 ) * lk_inv_lambda2_(atype2);

		Real const solvE1 = std::exp(-x1) * p.lk_coeff1 * inv_dis2;
		Real const solvE2 = std::exp(-x2) * p.lk_coeff2 * inv_dis2;
		dfasol1 = -2 * p.fasol_final_weight * ( dis_rad1 * lk_inv_lambda2_(atype1) + inv_dis ) * solvE1;
		dfasol2 = -2 * p.fasol_final_weight * ( dis_rad2 * lk_inv_lambda2_(atype2) + inv_dis ) * solvE2;

	} else if ( dis < fasol_cubic_poly_far_xhi_ ) {
		dfasol1 = cubic_polynomial_deriv( dis, p.fasol_cubic_poly1_far ) * p.fasol_final_weight;
		dfasol2 = cubic_polynomial_deriv( dis, p.fasol_cubic_poly2_far ) * p.fasol_final_weight;
	} else {
		dfasol1 = 0;
		dfasol2 = 0;
	}

}

/// @details only to be called when the distance, dis, is less than the switch-point for
/// the lj_switch_dis2sigma value.
inline
Real
Etable::analytic_ljrep_linearized(
	Real const dis,
	EtableParamsOnePair const & p
) const
{
	debug_assert( dis * dis < p.ljrep_linear_ramp_d2_cutoff );
	return  dis*p.lj_switch_slope + p.lj_switch_intercept;
}

/// @details: evaluate the standard Lennard-Jones 6-12 functional form.
/// Only call this function if the square distance is in the range
/// sqrt( p.ljrep_linear_ramp_d2_cutoff ) < dis < p.ljatr_cubic_poly_xhi
inline
Real
Etable::analytic_lj_generic_form(
	Real const ASSERT_ONLY( dis2 ),
	Real const inv_dis2,
	EtableParamsOnePair const & p
) const
{
	debug_assert( dis2 >= p.ljrep_linear_ramp_d2_cutoff );
	debug_assert( dis2 <= p.ljatr_cubic_poly_xhi * p.ljatr_cubic_poly_xhi );
	Real const inv_dis6  = inv_dis2 * inv_dis2 * inv_dis2;
	//Real const inv_dis12 = inv_dis6 * inv_dis6;

	return ( p.lj_r12_coeff * inv_dis6 + p.lj_r6_coeff ) * inv_dis6;
}

inline
Real
Etable::eval_cubic_polynomial(
	Real const x,
	CubicPolynomial const & cp
)
{
	return ((cp.c3*x+cp.c2)*x+cp.c1)*x + cp.c0;
}

inline
Real
Etable::cubic_polynomial_deriv(
	Real const x,
	CubicPolynomial const & cp
)
{
	return (3*cp.c3*x + 2*cp.c2)*x + cp.c1;
}

/// @details: evaluate the attractive component of the LJ term as it
/// ramps to zero.
/// Only call this function if the square distance is in the range
/// p.ljatr_cubic_poly_xlo < dis < p.ljatr_cubic_poly_xhi
inline
Real
Etable::analytic_ljatr_cubic_poly_ramp_to_zero(
	Real const dis,
	EtableParamsOnePair const & p
) const
{
	debug_assert( dis >= p.ljatr_cubic_poly_xlo );
	debug_assert( dis <= p.ljatr_cubic_poly_xhi );
	return eval_cubic_polynomial( dis, p.ljatr_cubic_poly_parameters );
}

/// @details: evaluate the derivative for the attractive component
/// of the LJ term as it ramps to zero.
/// Only call this function if the square distance is in the range
/// p.ljatr_cubic_poly_xlo < dis < p.ljatr_cubic_poly_xhi
inline
Real
Etable::analytic_ljatr_cubic_poly_ramp_to_zero_deriv(
	Real const dis,
	EtableParamsOnePair const & p
) const
{
	debug_assert( dis >= p.ljatr_cubic_poly_xlo );
	debug_assert( dis <= p.ljatr_cubic_poly_xhi );
	return cubic_polynomial_deriv( dis, p.ljatr_cubic_poly_parameters );
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
	debug_assert( atype1 <= atype2 );

	/// a) At distances below p.fasol_cubic_poly_close_start, the value of fasol is held constant.
	/// b) Then there's a cubic_poly to smooth between this constant region and the exponential region.
	/// c) Then the exponential region.
	/// d) Then the cubic_poly to smooth between the exponential region and where the energy goes to zero.
	if ( dis < p.fasol_cubic_poly_close_start ) {
		fa_solE = p.fasol_cubic_poly_close_flat * p.fasol_final_weight;
	} else if ( dis < p.fasol_cubic_poly_close_end ) {
		fa_solE = eval_cubic_polynomial( dis, p.fasol_cubic_poly_close );
		fa_solE *= p.fasol_final_weight;
	} else if ( dis < fasol_cubic_poly_far_xlo_ ) {
		/// exponential evaluation
		///debug_assert( atype1 <= atype2 ), which is accomplished at the top of this function.
		Real const dis_rad1 = dis - lj_radius(atype1);
		Real const x1 = ( dis_rad1 * dis_rad1 ) * lk_inv_lambda2_(atype1);
		Real const dis_rad2 = dis - lj_radius(atype2);
		Real const x2 = ( dis_rad2 * dis_rad2 ) * lk_inv_lambda2_(atype2);

		fa_solE =  inv_dis2 * ( std::exp(-x1) * p.lk_coeff1 + std::exp(-x2) * p.lk_coeff2 );
		fa_solE *= p.fasol_final_weight;

	} else if ( dis < fasol_cubic_poly_far_xhi_ ) {
		fa_solE = eval_cubic_polynomial( dis, p.fasol_cubic_poly_far );
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
	debug_assert( atype1 <= atype2 );

	/// a) At distances below p.fasol_cubic_poly_close_start, the value of fasol is held constant.
	/// b) Then there's a cubic_poly to smooth between this constant region and the exponential region.
	/// c) Then the exponential region.
	/// d) Then the cubic_poly to smooth between the exponential region and where the energy goes to zero.
	if ( dis < p.fasol_cubic_poly_close_start ) {
		fa_solE1 = p.fasol_cubic_poly1_close_flat * p.fasol_final_weight;
		fa_solE2 = p.fasol_cubic_poly2_close_flat * p.fasol_final_weight;
	} else if ( dis < p.fasol_cubic_poly_close_end ) {
		fa_solE1 = p.fasol_final_weight * eval_cubic_polynomial( dis, p.fasol_cubic_poly1_close );
		fa_solE2 = p.fasol_final_weight * eval_cubic_polynomial( dis, p.fasol_cubic_poly2_close );
	} else if ( dis < fasol_cubic_poly_far_xlo_ ) {
		/// exponential evaluation
		Real const dis_rad1 = dis - lj_radius(atype1);
		Real const x1 = ( dis_rad1 * dis_rad1 ) * lk_inv_lambda2_(atype1);
		Real const dis_rad2 = dis - lj_radius(atype2);
		Real const x2 = ( dis_rad2 * dis_rad2 ) * lk_inv_lambda2_(atype2);

		fa_solE1 = p.fasol_final_weight * inv_dis2 * std::exp(-x1) * p.lk_coeff1;
		fa_solE2 = p.fasol_final_weight * inv_dis2 * std::exp(-x2) * p.lk_coeff2;

	} else if ( dis < fasol_cubic_poly_far_xhi_ ) {
		fa_solE1 = p.fasol_final_weight * eval_cubic_polynomial( dis, p.fasol_cubic_poly1_far );
		fa_solE2 = p.fasol_final_weight * eval_cubic_polynomial( dis, p.fasol_cubic_poly2_far );
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
	return analytic_parameters_[ index ];
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
	return analytic_parameters_[ index ];
}


} // etable
} // scoring
} // core

#endif
