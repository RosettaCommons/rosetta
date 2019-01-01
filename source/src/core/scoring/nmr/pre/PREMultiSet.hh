// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    core/scoring/nmr/pre/PREMultiSet.hh
/// @brief   class that stores PRE data for one spinlabel site which can contain data for multiple
///          spinlabels and metals (i.e. multiple PRESingleSet objects)
/// @details last Modified: 10/12/16
/// @author  Georg Kuenze (georg.kuenze@vanderbilt.edu)

#ifndef INCLUDED_core_scoring_nmr_pre_PREMultiSet_HH
#define INCLUDED_core_scoring_nmr_pre_PREMultiSet_HH

// Unit headers
#include <core/scoring/nmr/pre/PREMultiSet.fwd.hh>

// Package headers
#include <core/scoring/nmr/pre/PRESingle.fwd.hh>
#include <core/scoring/nmr/pre/PRESingleSet.fwd.hh>
#include <core/scoring/nmr/NMRSpinlabel.fwd.hh>
#include <core/scoring/nmr/NMRGridSearch.fwd.hh>
#include <core/io/nmr/ParaIon.fwd.hh>
#include <core/scoring/nmr/util.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>

// Basic headers
#include <basic/Tracer.fwd.hh>

// Utility headers
#include <utility/vector1.hh>
#include <utility/fixedsizearray1.hh>

// Numeric headers
#include <numeric/xyzVector.fwd.hh>
#include <numeric/xyz.functions.fwd.hh>

// C++ headers
#include <iostream>
#include <string>

namespace core {
namespace scoring {
namespace nmr {
namespace pre {

class PREMultiSet {

public: // Typedefs

	typedef utility::vector1< Vector > CoordVector;
	// Vector of weight and radical atom xyz coordinates of N spinlabel conformers
	typedef utility::vector1< std::pair< Real, Vector > > WeightCoordVector;

	// Array with fixed size to hold predefined number of parameters
	// Can be rapidly created from the stack
	typedef utility::fixedsizearray1< Real, 6 > Vec6;

	typedef utility::vector1< utility::vector1< utility::vector1< Real > > > RealVector3D;
	typedef core::io::nmr::ParaIon    ParaIon;
	typedef core::io::nmr::ParaIonOP  ParaIonOP;
	typedef core::io::nmr::ParaIonCOP ParaIonCOP;

public: // Methods

	/// @brief construct PREMultiSet from vector of single datasets,
	///        the position of the spinlabel site, the type of paramagnetic
	///        ion and the dataset's weight.
	///        A default NMRSpinlabel object is created.
	PREMultiSet(
		utility::vector1<PRESingleSetOP> const & pre_singleset_vec,
		pose::Pose const & pose,
		Size const spinlabel_site,
		std::string const & iontype,
		Real const weight = 1.0
	);

	/// @brief construct PREMultiSet from vector of single datasets,
	///        the position of the spinlabel site, the type of paramagnetic
	///        ion, an NMRSpinlabel object and the dataset's weight.
	PREMultiSet(
		utility::vector1<PRESingleSetOP> const & pre_singleset_vec,
		pose::Pose const & pose,
		Size const spinlabel_site,
		std::string const & iontype,
		NMRSpinlabelOP spinlabel_ptr,
		Real const weight = 1.0
	);

	/// @brief construct PREMultiSet from vector of single datasets,
	///        the position of the spinlabel site, the type of paramagnetic
	///        ion, an NMRGridSearch object and the dataset's weight.
	PREMultiSet(
		utility::vector1<PRESingleSetOP> const & pre_singleset_vec,
		pose::Pose const & pose,
		Size const spinlabel_site,
		std::string const & iontype,
		NMRGridSearchOP gridsearch_ptr,
		Real const weight = 1.0
	);

	/// @brief copy constructor
	PREMultiSet(PREMultiSet const & other);

	/// @brief assignment operator
	PREMultiSet &
	operator=(PREMultiSet const & rhs);

	/// @brief destructor
	~PREMultiSet();

	/// @brief updates the spin coordinates every time the pose is changed
	///        make sure that this function is called before you call compute_pre_score()
	void update_spin_coordinates(pose::Pose const & pose);

	/// @brief calculates the PRE from a vector of spinlabel atom coordinates using
	///        the modified Solomon-Bloembergen (SBMF) equation and returns the
	///        weighted PRE score according to the single PRE value weighting scheme
	/// @details performs NLS fitting of the correlation times (tau_r, tau_e and tau_i)
	Real compute_pre_score_from_point_vector(WeightCoordVector const & points);

	/// @brief calculates the PRE from a single spinlabel atom position using the simplified
	///        Solomon-Bloembergen (SB) equation and returns the weighted PRE score
	///        according to the single PRE value weighting scheme
	/// @details performs NLS fitting of the correlation times (tau_r and tau_e) by default
	///          optionally, the spinlabel (sl) point position can be optimized too.
	Real compute_pre_score_from_single_point(
		Vector & point,
		bool opt_sl_pos=false
	);

	/// @brief determines the paraion/radical atom position of the spinlabel using a gridsearch or
	///        the dummy spinlabel rotamer ensemble and calculates the PRE score using
	///        the modified Solomon-Bloembergen (SBMF) or simplified Solomon-Bloembergen (SB) equation
	Real find_para_ion_position_and_compute_pre_score(
		pose::Pose const & pose,
		WeightCoordVector & spinlabel_counts_coords
	);

	/// @brief calculates and sets the xyz derivative of the PRE
	///        single PREs must be calculated in beforehand
	///        that's why make sure to call any of the three functions first:
	///        compute_pre_score_from_point_vector(), compute_pre_from_single_point() or
	///        find_para_ion_position_and_compute_pre_score()
	void
	set_atom_derivatives(
		pose::Pose & pose,
		WeightCoordVector const & spinlabel_counts_coords
	);

	// Getter and Setters
	utility::vector1<PRESingleSetOP> & get_pre_singleset_vec() { return pre_singleset_vec_; }
	utility::vector1<PRESingleSetOP> const & get_pre_singleset_vec() const { return pre_singleset_vec_; }
	Size get_total_number_pre() const { return total_number_pres_; }
	Size get_number_experiments() const { return number_experiments_; }
	Real get_weight() const { return weight_; }
	utility::vector1<Real> const & get_pre_values() const { return pre_values_; }
	utility::vector1<Real> const & get_pre_single_weights() const { return pre_single_weights_; }
	ParaIonCOP get_ion_type() const { return ion_type_; }
	Size get_spinlabel_site_rsd() const { return spinlabel_site_rsd_; }
	NMRSpinlabelOP get_spinlabel() { return spinlabel_; }
	NMRGridSearchOP get_gridsearch_iterator() { return gridsearch_iterator_; }
	utility::vector1< utility::vector1< CoordVector > > const & get_spin_coordinates() const { return spin_coordinates_; }
	NMR_VALUE_AVERAGING_TYPE get_averaging_type() const { return ave_type_; }
	Real get_protein_mass() const { return protein_mass_; }
	Real get_temperature() const { return temperature_; }
	Real get_tau_r() const { return tau_r_; }
	Real get_tau_c() const { return tau_c_; }
	Real get_tau_t() const { return tau_t_; }
	bool symmetric_pre_calc() const { return symmetric_pre_calc_; }
	Real get_tau_c_min() const { return tau_c_min_; }
	Real get_tau_c_max() const { return tau_c_max_; }
	bool optimize_paraion_position() const { return opt_para_ion_pos_; }

	void set_weight(Real weight) { weight_ = weight; }
	void set_averaging_type(std::string const & type);
	void set_spinlabel(NMRSpinlabelOP spinlabel_ptr) { spinlabel_ = spinlabel_ptr; }
	void set_gridsearch_iterator(NMRGridSearchOP gridsearch_ptr) { gridsearch_iterator_ = gridsearch_ptr; }
	void set_protein_mass(Real mass) {
		runtime_assert(mass > 0.0);
		protein_mass_ = mass;
		set_tau_r(calc_theoretical_tau_r());
	}
	void set_temperature(Real temperature) {
		runtime_assert(temperature > 0.0);
		temperature_ = temperature;
		set_tau_r(calc_theoretical_tau_r());
	}
	void set_tau_c(Real tau_c) { tau_c_ = tau_c; }
	void set_tau_t(Real tau_t) { tau_t_ = tau_t; }
	void set_tau_r(Real tau_r) {
		tau_r_ = tau_r;
		set_tau_c_limits(0.1*tau_r_, 10*tau_r_);
	}
	void set_tau_c_limits(Real tau_c_min, Real tau_c_max) {
		runtime_assert(tau_c_min <= tau_c_max);
		tau_c_min_ = tau_c_min;
		tau_c_max_ = tau_c_max;
	}

	void show(std::ostream & TR) const;

private: // Methods

	/// @brief default constructor
	PREMultiSet();

	void deep_copy_pre_singleset_vec(utility::vector1<PRESingleSetOP> const & other_vec);

	/// @brief utility function to initialize PREMultiSet
	void
	init_from_pre_singleset_vec(
		utility::vector1<PRESingleSetOP> const & pre_singleset_vec,
		pose::Pose const & pose
	);

	/// @brief set protein mass, temperature and correlation times to default values during initialization
	///        protein mass = 15kDa, temperature = 298 K
	void set_exp_default_conditions();

	/// @brief creates MTSL spinlabel as default
	void set_default_spinlabel();

	/// @brief register options
	void register_options();
	void init_from_cml();

	void create_para_ion_from_string(std::string const & iontype);

	/// @brief calculate rotational correlation time (in sec) of the protein
	///        using the Stokes-Einstein equation and the protein mass (in kDa == kg/mol)
	Real calc_theoretical_tau_r() const;

	/// @brief calculate sum of rotational correlation and electron relaxation time
	///        i.e. the total correlation time (in sec)
	Real calc_tau_c() const;

	/// @brief calculate electron spin frequency at given field strength in rad/s
	Real calc_omega_E(Real const & field_strength) const;

	/// @brief computes r6 and S2 vectors from given spinlabel atom coordinates and
	///        and current spin coordinates
	void update_r6_S2_values(WeightCoordVector const & spinlabel_counts_coords);

	/// @brief check if dimensionality of r6 and S2 vectors is consistent with averaging type
	///        and resize vectors if not.
	void resize_r6_S2_values();

	/// @brief computes the PRE score using the current r6 and S2 vectors and correlation times
	/// @params
	/// use_sb:   use only the simplified SB equation (e.g. if there is only one spinlabel)
	Real compute_pre_score_from_current_data(bool use_sb=false);

	/// @brief calculate the constant prefactor of the dipolar part of the PRE
	/// @params
	/// gamma_I:     gyromagnetic ratio of the nuclear spin (must be provided in rad/(s*T), dimension is 10^6)
	/// gJ:          electron Lande factor
	/// S:           total spin quantum number
	Real PRE_DD_prefactor(
		Real const & gamma_I,
		Real const & gJ,
		Real const & S,
		PRE_RATE_TYPE const & rate_type
	) const;

	/// @brief calculate the constant prefactor of the Curie part of the PRE
	/// @params
	/// gamma_I:     gyromagnetic ratio of the nuclear spin (must be provided in rad/(s*T), dimension is 10^6)
	/// gJ:          electron Lande factor
	/// S:           total spin quantum number
	/// B0:          magnetic field strength (in Tesla)
	/// T:           temperature (in Kelvin)
	Real PRE_Curie_prefactor(
		Real const & gamma_I,
		Real const & gJ,
		Real const & S,
		Real const & B0,
		Real const & T,
		PRE_RATE_TYPE const & rate_type
	) const;

	/// @brief calculate the ensemble averaged distance <r^-6> and generalized order parameter S2
	///        for a given spinlabel atom - nuclear spin connection vector as used in the SBMF equation
	/// @params
	/// spinlabel_counts_coords: vector of weights and para ion xyz coordinates of N distinct spinlabel conformers
	/// eq_spins_coords:         vector of xyz coordinates of nuclear spin(s) that give rise to the observed PRE
	/// one_over_r6:             averaged distance <r^-6> as return value
	/// S2:                      generalized order parameter S2 as return value
	void
	calc_r6_S2(
		WeightCoordVector const & spinlabel_counts_coords,
		CoordVector const & eq_spins_coords,
		Real & one_over_r6,
		Real & S2
	) const;

	/// @brief spectral density function for the simplified Solomon-Bloembergen (SB) equation
	/// @details A detailed description can be found in: Iwahara J, Schwieters CD & Clore GM, 2004, JACS 126,5879-5896
	/// @params
	/// omega_I:     nuclear spin resonance frequency (must be provided in rad/s, dimension is 10^6)
	/// tau_c:       correlation time (must be provided in s, typical dimension is 10^-9)
	/// one_over_r6: ensemble averaged distance <r^-6> (must be provided in Ang.^-6)
	Real
	J_SB(
		Real const & omega_I,
		Real const & tau_c,
		Real const & one_over_r6
	) const;

	/// @brief spectral density function for the modified Solomon-Bloembergen (SBMF) equation
	/// @details A detailed description can be found in: Iwahara J, Schwieters CD & Clore GM, 2004, JACS 126,5879-5896
	/// @params
	/// omega_I:     nuclear spin resonance frequency (must be provided in rad/s, dimension is 10^6)
	/// tau_c:       correlation time (must be provided in s, typical dimension is 10^-9)
	/// tau_t:       total correlation time (must be provided in s)
	/// one_over_r6: ensemble averaged distance <r^-6> (must be provided in Ang.^-6)
	/// S2:          generalized order parameter
	Real
	J_SBMF(
		Real const & omega_I,
		Real const & tau_c,
		Real const & tau_t,
		Real const & one_over_r6,
		Real const & S2
	) const;

	/// @brief calculate the R1 paramagnetic relaxation rate from the dipolar and Curie contribution
	/// @params:
	/// params[1] = prefac_dd:   prefactor of dipolar part of PRE
	/// params[2] = prefac_curie prefactor of Curie part of PRE
	/// params[3] = omega_I:     nuclear spin resonance frequency (must be provided in rad/s, dimension is 10^6)
	/// params[4] = tau_r:       rotational correlation time (must be provided in s, typical dimension is 10^-9)
	/// params[5] = tau_c:       correlation time (must be provided in s, typical dimension is 10^-9)
	/// params[6] = tau_t:       total correlation time (must be provided in s, typical dimension is 10^-10)
	/// one_over_r6:             ensemble averaged distance <r^-6> (must be provided in Ang.^-6)
	/// S2:                      generalized order parameter
	/// scaling:                 optional scaling factor
	/// use_sb:                  use only the simplified SB equation (e.g. if there is only one spinlabel)
	Real
	R1_Para(
		Vec6 const & params,
		Real const & one_over_r6,
		Real const & S2,
		Real const & scaling = 1.0,
		bool use_sb = false
	) const;

	/// @brief calculate the R2 paramagnetic relaxation rate from the dipolar and Curie contribution
	/// @params:
	/// params[1] = prefac_dd:   prefactor of dipolar part of PRE
	/// params[2] = prefac_curie prefactor of Curie part of PRE
	/// params[3] = omega_I:     nuclear spin resonance frequency (must be provided in rad/s, dimension is 10^6)
	/// params[4] = tau_r:       rotational correlation time (must be provided in s, typical dimension is 10^-9)
	/// params[5] = tau_c:       correlation time (must be provided in s, typical dimension is 10^-9)
	/// params[6] = tau_t:       total correlation time (must be provided in s, typical dimension is 10^-10)
	/// one_over_r6:             ensemble averaged distance <r^-6> (must be provided in Ang.^-6)
	/// S2:                      generalized order parameter
	/// scaling:                 optional scaling factor
	/// use_sb:                  use only the simplified SB equation (e.g. if there is only one spinlabel)
	Real
	R2_Para(
		Vec6 const & params,
		Real const & one_over_r6,
		Real const & S2,
		Real const & scaling = 1.0,
		bool use_sb = false
	) const;

	/// @brief calculate the xyz derivative of the ensemble averaged distance <r^-6>
	///        and the generalized order parameter S2 for a given spinlabel atom -
	///        nuclear spin connection vector
	/// @params
	/// spinlabel_counts_coords: vector of weights and para ion xyz coordinates of N distinct spinlabel conformers
	/// eq_spins_coords:         vector of xyz coordinates of nuclear spin(s) that give rise to the observed PRE
	/// dr6_dXYZ:                xyz derivative of <r^-6> as return value
	/// dS2_dXYZ:                xyz derivative of S2 as return value
	void
	calc_dr6_dS2_dXYZ(
		WeightCoordVector const & spinlabel_counts_coords,
		CoordVector const & eq_spins_coords,
		Vector & dr6_dXYZ,
		Vector & dS2_dXYZ
	) const;

	/// @brief xyz derivative of the spectral density function for the simplified Solomon-Bloembergen (SB) equation
	/// @params
	/// omega_I:     nuclear spin resonance frequency (must be provided in rad/s, dimension is 10^6)
	/// tau_c:       correlation time (must be provided in s, typical dimension is 10^-9)
	/// dr6_dXYZ:    xyz derivative of <r^-6>
	Vector
	dJ_dXYZ_SB(
		Real const & omega_I,
		Real const & tau_c,
		Vector const & dr6_dXYZ
	) const;

	/// @brief xyz derivative of the spectral density function for the modified Solomon-Bloembergen (SBMF) equation
	/// @params
	/// omega_I:     nuclear spin resonance frequency (must be provided in rad/s, dimension is 10^6)
	/// tau_c:       correlation time (must be provided in s, typical dimension is 10^-9)
	/// tau_t:       total correlation time (must be provided in s, typical dimension is 10^-10)
	/// one_over_r6: ensemble averaged distance <r^-6> (must be provided in Ang.^-6)
	/// S2:          generalized order parameter
	/// dr6_dXYZ:    xyz derivative of <r^-6>
	/// dS2_dXYZ:    xyz derivative of S2
	Vector
	dJ_dXYZ_SBMF(
		Real const & omega_I,
		Real const & tau_c,
		Real const & tau_t,
		Real const & one_over_r6,
		Real const & S2,
		Vector const & dr6_dXYZ,
		Vector const & dS2_dXYZ
	) const;

	/// @brief calculate the xyz derivative of the R1 paramagnetic relaxation rate
	/// @params
	/// params[1] = prefac_dd:   prefactor of dipolar part of PRE
	/// params[2] = prefac_curie prefactor of Curie part of PRE
	/// params[3] = omega_I:     nuclear spin resonance frequency (must be provided in rad/s, dimension is 10^6)
	/// params[4] = tau_r:       rotational correlation time (must be provided in s, typical dimension is 10^-9)
	/// params[5] = tau_c:       correlation time (must be provided in s, typical dimension is 10^-9)
	/// params[6] = tau_t:       total correlation time (must be provided in s, typical dimension is 10^-10)
	/// one_over_r6:             ensemble averaged distance <r^-6> (must be provided in Ang.^-6)
	/// S2:                      generalized order parameter
	/// dr6_dXYZ:                xyz derivative of <r^-6>
	/// dS2_dXYZ:                xyz derivative of S2
	/// scaling:                 optional scaling factor
	/// use_sb:                  use only the simplified SB equation (e.g. if there is only one spinlabel)
	Vector
	dR1_dXYZ(
		Vec6 const & params,
		Real const & one_over_r6,
		Real const & S2,
		Vector const & dr6_dXYZ,
		Vector const & dS2_dXYZ,
		Real const & scaling = 1.0,
		bool use_sb = false
	) const;

	/// @brief calculate the xyz derivative of the R2 paramagnetic relaxation rate
	/// @params
	/// params[1] = prefac_dd:   prefactor of dipolar part of PRE
	/// params[2] = prefac_curie prefactor of Curie part of PRE
	/// params[3] = omega_I:     nuclear spin resonance frequency (must be provided in rad/s, dimension is 10^6)
	/// params[4] = tau_r:       rotational correlation time (must be provided in s, typical dimension is 10^-9)
	/// params[5] = tau_c:       correlation time (must be provided in s, typical dimension is 10^-9)
	/// params[6] = tau_t:       total correlation time (must be provided in s, typical dimension is 10^-10)
	/// one_over_r6:             ensemble averaged distance <r^-6> (must be provided in Ang.^-6)
	/// S2:                      generalized order parameter
	/// dr6_dXYZ:                xyz derivative of <r^-6>
	/// dS2_dXYZ:                xyz derivative of S2
	/// scaling:                 optional scaling factor
	/// use_sb:                  use only the simplified SB equation (e.g. if there is only one spinlabel)
	Vector
	dR2_dXYZ(
		Vec6 const & params,
		Real const & one_over_r6,
		Real const & S2,
		Vector const & dr6_dXYZ,
		Vector const & dS2_dXYZ,
		Real const & scaling = 1.0,
		bool use_sb = false
	) const;

	/// @brief pre error function used in the lmmin function to optimize tau
	///        * par is an array of fit parameters [tau_c, tau_t]
	///        * data is a pointer to the PREMultiSet object i.e. to all data needed
	///        for PRE calculation and NLS fitting
	///        * fvc is an array holding the residuals of the fit calculation
	friend
	void
	pre_erf_opt_tau(
		Real const *par,
		int m_dat,
		void const *data,
		Real *fvec,
		int */*info*/
	);

	/// @brief pre error function used in the lmmin function to optimize tau and the para ion position
	///        * par is an array of fit parameters [xM, yM, zM, tau_c]
	///        * data is a pointer to the PREMultiSet object i.e. to all data needed
	///        for PRE calculation and NLS fitting
	///        * fvc is an array holding the residuals of the fit calculation
	friend
	void
	pre_erf_opt_tau_xyz(
		Real const *par,
		int m_dat,
		void const *data,
		Real *fvec,
		int */*info*/
	);

private: // Data

	utility::vector1< PRESingleSetOP > pre_singleset_vec_;
	Size total_number_pres_;
	Size number_experiments_;
	Real weight_;
	utility::vector1<Real> pre_values_;
	utility::vector1<Real> pre_single_weights_;

	// To avoid calculating the average SL-spin distance and
	// the order parameter multiple times during scoring (e.g
	// during the multiple NLS fits) we want to calculate them
	// only once after the pose has changed and store them as
	// private data members here.
	// In case we perform mean averaging of the PRE, the inner
	// vector has just dimension one. If we perform sum averaging
	// the inner vector's size equals the number of components
	// that we want to sum up.
	RealVector3D one_over_r6_values_;
	RealVector3D s2_values_;

	// Type of paramagnetic ion or radical atom that gives rise to the PRE
	ParaIonOP ion_type_;

	// Residue number of the spinlabel site
	Size spinlabel_site_rsd_;

	// The paraion/radical position can be inferred from the spinlabel structure
	// or be found by performing a grid search
	NMRSpinlabelOP spinlabel_;
	NMRGridSearchOP gridsearch_iterator_;

	// To avoid retrieving the spin coordinates every time we compute the PRE and fit the
	// correlation times, we store the spin coordinates here in a lookup table.
	// The outer vector runs over the number of PREs, the middle vector over symmetrical spins
	// (in case the pose is not symmetric this vector has only size 1) and the inner vector runs
	// over equivalent spins (e.g. methyl protons)
	utility::vector1< utility::vector1< CoordVector > > spin_coordinates_;

	// Experimental conditions in addition to those in PRESingleSet
	Real protein_mass_;    // in kDa
	Real temperature_;     // in Kelvin

	Real tau_r_;            // rotational correlation time; approximated from the protein molecular mass
	Real tau_c_;            // sum of rotational correlation time and electron relaxation time 1/tau_c = 1/tau_r + 1/tau_e
	Real tau_t_;            // total correlation time
	// sum of tau_c and correlation of internal spinlabel motion 1/tau_t = 1/tau_c + 1/tau_i
	// to be optimized by non-linear least squares fitting
	// Values the user might optionally provide
	// Otherwise default values are used
	bool symmetric_pre_calc_;
	NMR_VALUE_AVERAGING_TYPE ave_type_;
	Size nls_repeats_;
	Real tau_c_min_;  // lower and upper limit of the correlation time tau_c used to constrain the NLS fit
	Real tau_c_max_;
	bool opt_para_ion_pos_;

};

} // namespace pre
} // namespace nmr
} // namespace scoring
} // namespace core

#endif // INCLUDED_core_scoring_nmr_pre_PREMultiSet_HH
