// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    core/scoring/nmr/pre/PREMultiSet.cc
/// @brief   implementation of class PREMultiSet
/// @details last Modified: 10/12/16
/// @author  Georg Kuenze (georg.kuenze@vanderbilt.edu)

// Unit headers
#include <core/scoring/nmr/pre/PREMultiSet.hh>

// Package headers
#include <core/scoring/nmr/pre/PRESingleSet.hh>
#include <core/scoring/nmr/pre/PRESingle.hh>
#include <core/scoring/nmr/NMRSpinlabel.hh>
#include <core/scoring/nmr/NMRDummySpinlabelEnsemble.hh>
#include <core/scoring/nmr/NMRDummySpinlabelVoxelGrid.hh>
#include <core/scoring/nmr/NMRGridSearch.hh>
#include <core/scoring/nmr/util.hh>
#include <core/io/nmr/ParaIon.hh>
#include <core/io/nmr/ParamagneticDatabaseHandler.hh>
#include <core/io/nmr/util.hh>

// Project headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/id/AtomID.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/pose/symmetry/util.hh>

// Basic headers
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/nmr.OptionKeys.gen.hh>

// Utility headers
#include <utility/exit.hh>
#include <utility/vector1.hh>

// Numeric headers
#include <numeric/xyzVector.hh>
#include <numeric/constants.hh>
#include <numeric/nls/lmmin.hh>
#include <numeric/random/random.hh>
#include <numeric/xyz.functions.hh>

// C++ headers
#include <cmath>
#include <limits>
#include <iostream>
#include <iomanip>
#include <string>
#include <algorithm>
#include <numeric>

// Boost headers
#include <boost/algorithm/string.hpp>

namespace core {
namespace scoring {
namespace nmr {
namespace pre {

static basic::Tracer TR( "core.scoring.nmr.pre.PREMultiSet" );

/// @brief construct PREMultiSet from vector of single datasets,
///        the position of the spinlabel site, the type of paramagnetic
///        ion and the dataset's weight.
///        A default NMRSpinlabel object is created.
PREMultiSet::PREMultiSet(
	utility::vector1<PRESingleSetOP> const & pre_singleset_vec,
	pose::Pose const & pose,
	Size const spinlabel_site,
	std::string const & iontype,
	Real const weight
) :
	pre_singleset_vec_(pre_singleset_vec),
	number_experiments_(pre_singleset_vec.size()),
	weight_(weight),
	spinlabel_site_rsd_(spinlabel_site),
	gridsearch_iterator_(nullptr),
	symmetric_pre_calc_(false),
	ave_type_(MEAN)
{
	runtime_assert_msg(pre_singleset_vec.size() > 0, "ERROR while trying to construct PREMultiSet. Vector of single PRE experiments has zero size.");
	register_options();

	// initialize members from command line
	init_from_cml();

	// resize and fill pre_values and pre_weights from vector of PRESingleSets, resize spin_coordinates
	init_from_pre_singleset_vec(pre_singleset_vec, pose);

	// Create ParaIon object
	create_para_ion_from_string(iontype);

	// Set MTSL spinlabel as default
	set_default_spinlabel();

	// Set default experimental conditions and correlation times
	set_exp_default_conditions();
}

/// @brief construct PREMultiSet from vector of single datasets,
///        the position of the spinlabel site, the type of paramagnetic
///        ion, an NMRSpinlabel object and the dataset's weight.
PREMultiSet::PREMultiSet(
	utility::vector1<PRESingleSetOP> const & pre_singleset_vec,
	pose::Pose const & pose,
	Size const spinlabel_site,
	std::string const & iontype,
	NMRSpinlabelOP spinlabel_ptr,
	Real const weight
) :
	pre_singleset_vec_(pre_singleset_vec),
	number_experiments_(pre_singleset_vec.size()),
	weight_(weight),
	spinlabel_site_rsd_(spinlabel_site),
	spinlabel_(new NMRSpinlabel(*spinlabel_ptr)),
	gridsearch_iterator_(nullptr),
	symmetric_pre_calc_(false),
	ave_type_(MEAN)
{
	runtime_assert_msg(pre_singleset_vec.size() > 0, "ERROR while trying to construct PREMultiSet. Vector of single PRE experiments has zero size.");
	register_options();

	// initialize members from command line
	init_from_cml();

	// resize and fill pre_values and pre_weights from vector of PRESingleSets, resize spin_coordinates
	init_from_pre_singleset_vec(pre_singleset_vec, pose);

	// Create ParaIon object
	create_para_ion_from_string(iontype);

	// Set default experimental conditions and correlation times
	set_exp_default_conditions();
}

/// @brief construct PREMultiSet from vector of single datasets,
///        the position of the spinlabel site, the type of paramagnetic
///        ion, an NMRGridSearch object and the dataset's weight.
PREMultiSet::PREMultiSet(
	utility::vector1<PRESingleSetOP> const & pre_singleset_vec,
	pose::Pose const & pose,
	Size const spinlabel_site,
	std::string const & iontype,
	NMRGridSearchOP gridsearch_ptr,
	Real const weight
) :
	pre_singleset_vec_(pre_singleset_vec),
	number_experiments_(pre_singleset_vec.size()),
	weight_(weight),
	spinlabel_site_rsd_(spinlabel_site),
	spinlabel_(nullptr),
	gridsearch_iterator_(new NMRGridSearch(*gridsearch_ptr)),
	symmetric_pre_calc_(false),
	ave_type_(MEAN)
{
	runtime_assert_msg(pre_singleset_vec.size() > 0, "ERROR while trying to construct PREMultiSet. Vector of single PRE experiments has zero size.");
	register_options();

	// initialize members from command line
	init_from_cml();

	// resize and fill pre_values and pre_weights from vector of PRESingleSets, resize spin_coordinates
	init_from_pre_singleset_vec(pre_singleset_vec, pose);

	// Create ParaIon object
	create_para_ion_from_string(iontype);

	// Set default experimental conditions and correlation times
	set_exp_default_conditions();
}

/// @brief copy constructor
PREMultiSet::PREMultiSet(PREMultiSet const & other) :
	total_number_pres_(other.total_number_pres_),
	number_experiments_(other.number_experiments_),
	weight_(other.weight_),
	pre_values_(other.pre_values_),
	pre_single_weights_(other.pre_single_weights_),
	one_over_r6_values_(other.one_over_r6_values_),
	s2_values_(other.s2_values_),
	ion_type_( other.ion_type_ ? new ParaIon( *(other.ion_type_) ) : nullptr ),
	spinlabel_site_rsd_(other.spinlabel_site_rsd_),
	spinlabel_( other.spinlabel_ ? new NMRSpinlabel( *(other.spinlabel_) ) : nullptr ),
	gridsearch_iterator_( other.gridsearch_iterator_ ? new NMRGridSearch( *(other.gridsearch_iterator_) ) : nullptr ),
	spin_coordinates_(other.spin_coordinates_),
	protein_mass_(other.protein_mass_),
	temperature_(other.temperature_),
	tau_r_(other.tau_r_),
	tau_c_(other.tau_c_),
	tau_t_(other.tau_t_),
	symmetric_pre_calc_(other.symmetric_pre_calc_),
	ave_type_(other.ave_type_),
	nls_repeats_(other.nls_repeats_),
	tau_c_min_(other.tau_c_min_),
	tau_c_max_(other.tau_c_max_),
	opt_para_ion_pos_(other.opt_para_ion_pos_)
{
	// Make deep copy of PRESingleSet vector
	deep_copy_pre_singleset_vec(other.pre_singleset_vec_);
}

/// @brief assignment operator
PREMultiSet &
PREMultiSet::operator=(PREMultiSet const & rhs) {
	if ( this != &rhs ) {
		total_number_pres_ = rhs.total_number_pres_;
		number_experiments_ = rhs.number_experiments_;
		weight_ = rhs.weight_;
		pre_values_ = rhs.pre_values_;
		pre_single_weights_ = rhs.pre_single_weights_;
		one_over_r6_values_ = rhs.one_over_r6_values_;
		s2_values_ = rhs.s2_values_;
		ion_type_ = rhs.ion_type_ ? ParaIonOP( new ParaIon( *(rhs.ion_type_) ) ) : nullptr;
		spinlabel_site_rsd_ = rhs.spinlabel_site_rsd_;
		spinlabel_ = rhs.spinlabel_ ? NMRSpinlabelOP( new NMRSpinlabel( *(rhs.spinlabel_) ) ) : nullptr;
		gridsearch_iterator_ = rhs.gridsearch_iterator_ ? NMRGridSearchOP( new NMRGridSearch( *(rhs.gridsearch_iterator_) ) ) : nullptr;
		spin_coordinates_ = rhs.spin_coordinates_;
		protein_mass_ = rhs.protein_mass_;
		temperature_ = rhs.temperature_;
		tau_r_ = rhs.tau_r_;
		tau_c_ = rhs.tau_c_;
		tau_t_ = rhs.tau_t_;
		symmetric_pre_calc_ = rhs.symmetric_pre_calc_;
		ave_type_ = rhs.ave_type_;
		nls_repeats_ = rhs.nls_repeats_;
		tau_c_min_ = rhs.tau_c_min_;
		tau_c_max_ = rhs.tau_c_max_;
		opt_para_ion_pos_ = rhs.opt_para_ion_pos_;
		// Make deep copy of PRESingleSet vector
		pre_singleset_vec_.clear();
		deep_copy_pre_singleset_vec(rhs.pre_singleset_vec_);
	}
	return *this;
}

/// @brief destructor
PREMultiSet::~PREMultiSet() {}

/// @brief utility function to initialize PREMultiSet
void
PREMultiSet::init_from_pre_singleset_vec(
	utility::vector1<PRESingleSetOP> const & pre_singleset_vec,
	pose::Pose const & pose
)
{
	using pose::symmetry::is_symmetric;
	using pose::symmetry::symmetry_info;

	pre_singleset_vec_ = pre_singleset_vec;
	// If we perform PRE calculation with automatic deduction of protein symmetry,
	// we resize the middle vector of the spin_coordinates_ array to the number of symmetric subunits
	Size num_subunits(1);
	if ( symmetric_pre_calc_ && is_symmetric( pose ) ) {
		conformation::symmetry::SymmetryInfoCOP syminfo_ptr = symmetry_info( pose );
		num_subunits = syminfo_ptr->subunits();
	}

	// Resize vectors for PRE values, single weights, <r-6>, S2 and spin coordinates
	total_number_pres_ = 0;
	for ( auto& p : pre_singleset_vec_ ) {
		total_number_pres_ += p->get_number_pre();
	}
	pre_values_.resize(total_number_pres_);
	pre_single_weights_.resize(total_number_pres_);
	one_over_r6_values_.resize(total_number_pres_);
	s2_values_.resize(total_number_pres_);
	spin_coordinates_.resize(total_number_pres_);

	// Fill vector for PRE values and single weights
	Size index_offset(0);
	for ( Size i = 1; i <= number_experiments_; ++i ) {
		utility::vector1< PRESingle > const & pre_single_vec = pre_singleset_vec_[i]->get_pre_single_vec();

		Real max_pre(0);          // largest PRE of current experiment i
		Size no_pres = pre_singleset_vec_[i]->get_number_pre(); // number of PREss of current experiment i; increments the index_offset

		for ( Size j = 1; j <= no_pres; ++j ) {
			// We don't have to do additional scaling of the PRE data here, since this is already done during instantiation of the PRESingleSet.
			// Thus, we simply copy the PRE values from the PRESingleSet.
			pre_values_[index_offset + j] = pre_single_vec[j].get_pre_exp();
			if ( std::abs(pre_values_[index_offset + j]) > max_pre ) {
				max_pre = std::abs(pre_values_[index_offset + j]);
			}

			// Resize spin_coordinates array, <r-6> and S2 vectors
			spin_coordinates_[index_offset + j].resize(num_subunits);
			one_over_r6_values_[index_offset + j].resize(num_subunits);
			s2_values_[index_offset + j].resize(num_subunits);

			for ( Size su = 1; su <= num_subunits; ++su ) {
				spin_coordinates_[index_offset + j][su].resize(pre_single_vec[j].get_protein_spins().size());
				// If we perform simple mean averaging the most inner vector has size 1
				// If we perform sum averaging (e.g. we store explicitly the coords of
				// symmetric spins in the most inner vector of spin_coordinates_ the
				// third dimension of <r-6> and S2 vectors should have the same dimension
				if ( ave_type_ == MEAN ) {
					one_over_r6_values_[index_offset + j][su].resize(1);
					s2_values_[index_offset + j][su].resize(1);
				} else if ( ave_type_ == SUM ) {
					one_over_r6_values_[index_offset + j][su].resize(pre_single_vec[j].get_protein_spins().size());
					s2_values_[index_offset + j][su].resize(pre_single_vec[j].get_protein_spins().size());
				}
			}
		} // number PREs per experiment

		Real pre_err(1.0);
		Real pre_weight(1.0);
		for ( Size j = 1; j <= no_pres; ++j ) {
			if ( pre_singleset_vec_[i]->get_single_pre_weighting_scheme() == CONST ) {
				pre_weight = 1.0;
			} else {
				pre_err = pre_single_vec[j].get_pre_err();
				runtime_assert_msg(pre_err > 1e-6, "ERROR in calculating single PRE weights. Experimental PRE error is unreasonably small (< 1e-6) and would produce a large weight (> 10e+12) for the chosen weighting scheme. Check PRE datafile.");
				if ( pre_singleset_vec_[i]->get_single_pre_weighting_scheme() == SIGMA ) {
					pre_weight = (1.0/(pre_err * pre_err));
				} else if ( pre_singleset_vec_[i]->get_single_pre_weighting_scheme() == OBSIG ) {
					pre_weight = (std::abs(pre_values_[index_offset + j])/(pre_err * pre_err * max_pre));
				}
			}
			pre_single_weights_[index_offset + j] = pre_weight;
			pre_singleset_vec_[i]->pre_single_vec_[j].set_weight(pre_weight);
		} // number PREs per experiment

		// increment the index offset
		index_offset += no_pres;

	} // number of PRE experiments
}

/// @brief updates the spin coordinates every time the pose is changed
///        make sure that this function is called before you call compute_score()
void
PREMultiSet::update_spin_coordinates(pose::Pose const & pose) {
	using pose::symmetry::is_symmetric;
	using pose::symmetry::symmetry_info;

	// Here we split behavior. If summetric_pre_calc_ is true and the pose is indeed symmetric
	// we deduce the symmetric spins from the pose and fill their coordinates in the middle vector of spin_coordinates_
	Size num_subunits(1);
	conformation::symmetry::SymmetryInfoCOP syminfo_ptr;
	if ( symmetric_pre_calc_ && is_symmetric( pose ) ) {
		syminfo_ptr = symmetry_info( pose );
		num_subunits = syminfo_ptr->subunits();
	}

	Size index_offset(0);
	// The outer vector runs over the number of experiments * No. PREs
	for ( Size i = 1; i <= number_experiments_; ++i ) {
		Size no_pres = pre_singleset_vec_[i]->get_number_pre();

		for ( Size j = 1; j <= no_pres; ++j ) {

			// Number of equivalent spins over which we perform averaging of the PRE
			Size num_eq_spins(pre_singleset_vec_[i]->get_pre_single_vec()[j].get_protein_spins().size());

			for ( Size k = 1; k <= num_eq_spins; ++k ) {

				Size spin_rsd_asu(pre_singleset_vec_[i]->get_pre_single_vec()[j].get_protein_spins()[k].rsd());
				utility::vector1< Size > rsds_for_pre_all_subunits( 1, spin_rsd_asu );

				if ( symmetric_pre_calc_ && is_symmetric( pose ) && syminfo_ptr ) {
					runtime_assert_msg(syminfo_ptr->is_asymmetric_seqpos(spin_rsd_asu),
						"ERROR: For PRE calculation with automatic symmetry deduction, the residue number in PRE datafile must refer only to asymmetric subunit.");
					utility::vector1< Size > symm_spin_rsds = syminfo_ptr->bb_clones(spin_rsd_asu);
					rsds_for_pre_all_subunits.resize(num_subunits);
					std::copy(symm_spin_rsds.begin(), symm_spin_rsds.end(), rsds_for_pre_all_subunits.begin()+1);
				}
				for ( Size su = 1; su <= num_subunits; ++su ) {
					spin_coordinates_[index_offset + j][su][k] = pose.residue(rsds_for_pre_all_subunits[su]).atom(pre_singleset_vec_[i]->get_pre_single_vec()[j].get_protein_spins()[k].atomno()).xyz();

				} // number subunits

			} //number equivalent spins

		} // number PREs per experiment

		// increment index offset
		index_offset += no_pres;

	} // number of PRE experiments
}

/// @brief pre error function used in the lmmin function top optimize
///        * par is an array of fit parameters [tau_c, tau_t]
///        * data is a pointer to the PREMultiSet object i.e. to all data needed
///        for PRE calculation and NLS fitting
///        * fvc is an array holding the residuals of the fit calculation
void
pre_erf_opt_tau(
	Real const *par,
	int m_dat,
	void const *data,
	Real *fvec,
	int */*info*/
)
{
	PREMultiSet * pre_multiset_ptr = static_cast< PREMultiSet * >(const_cast< void* >(data));
	Real * nonconst_par = const_cast< Real* >(par);

	// Make sure that ion type is set
	if ( !pre_multiset_ptr->ion_type_ ) {
		utility_exit_with_message("ERROR while trying to calculate the PRE in PREMultiSet. Ion type is not set.");
	}

	// Constrain correlation times to reasonable ranges
	// tau_c (in range tc_min - tc_max)
	if ( !( pre_multiset_ptr->tau_c_min_ <= par[0] ) || !( par[0] <= pre_multiset_ptr->tau_c_max_ ) ) {
		nonconst_par[0] = 0.5*(pre_multiset_ptr->tau_c_max_ - pre_multiset_ptr->tau_c_min_)*std::tanh(par[0])
			+(pre_multiset_ptr->tau_c_max_ + pre_multiset_ptr->tau_c_min_)/2.0;
	}
	// tau_t (in range 0 - tau_c)
	if ( !( 0.0 <= par[1]) || !( par[1] <= nonconst_par[0] ) ) {
		nonconst_par[1] = 0.5*nonconst_par[0]*std::tanh(par[1])+(0.5*nonconst_par[0]);
	}

	PREMultiSet::Vec6 pre_params;
	// params[1] = prefac_dd
	// params[2] = prefac_curie
	// params[3] = omega_I
	// params[4] = tau_r
	// params[5] = tau_c
	// params[6] = tau_t

	Real scal(1.0); // scaling factor per PRE experiment (normalization by experiment StdDev)

	// We need to keep track of with which index a new experiment in the PRE value and weights vector starts
	// so that we can update the parameter set for the PRE equation. The start index is 1 plus
	// the number of PREs of the previous experiments ( e.g. 1, 101, 201, ... ).
	utility::vector1< Size > start_index(pre_multiset_ptr->number_experiments_);
	for ( Size i = 1; i <= pre_multiset_ptr->number_experiments_; ++i ) {
		if ( i == 1 ) { start_index[i] = 1; }
		else { start_index[i] = start_index[i-1] + pre_multiset_ptr->pre_singleset_vec_[i-1]->get_number_pre(); }
	}
	runtime_assert(start_index.back() == ( pre_multiset_ptr->total_number_pres_ - pre_multiset_ptr->pre_singleset_vec_.back()->get_number_pre() + 1));

	// i is the index for the ith PRE in the pre_values_ vector
	// j is the index for the jth experiment in the PRESingleSet vector
	// k is the index for the start_index vector (to figure out when to change the Params for the PRE eq.)
	for ( Size i = 1, j = 0, k = 1, i_end = static_cast< Size >(m_dat); i <= i_end; ++i ) { // loop over the total number of PREs
		fvec[i-1] = pre_multiset_ptr->pre_values_[i];
		// update parameters for next experiment and increase index j
		if ( i == start_index[k] ) {
			++j; // increment index for PRESingleSet vector
			k < pre_multiset_ptr->number_experiments_ ? ++k : k = pre_multiset_ptr->number_experiments_;

			Real gamma_I(pre_multiset_ptr->pre_singleset_vec_[j]->calc_gamma_I());
			Real gJ(pre_multiset_ptr->ion_type_->get_gJ());
			Real S(pre_multiset_ptr->ion_type_->get_S());
			Real B0(pre_multiset_ptr->pre_singleset_vec_[j]->get_field_strength() / 42.576); // B0 in Tesla

			pre_params[1] = pre_multiset_ptr->PRE_DD_prefactor(gamma_I, gJ, S, pre_multiset_ptr->pre_singleset_vec_[j]->get_pre_rate_type());
			pre_params[2] = pre_multiset_ptr->PRE_Curie_prefactor(gamma_I, gJ, S, B0, pre_multiset_ptr->get_temperature(), pre_multiset_ptr->pre_singleset_vec_[j]->get_pre_rate_type());
			pre_params[3] = pre_multiset_ptr->pre_singleset_vec_[j]->calc_omega_I();
			pre_params[4] = pre_multiset_ptr->tau_r_;
			pre_params[5] = nonconst_par[0];
			pre_params[6] = nonconst_par[1];

			scal = 1.0/pre_multiset_ptr->pre_singleset_vec_[j]->get_scaling_factor();
		}

		for ( Size su = 1, su_end = pre_multiset_ptr->spin_coordinates_[i].size(); su <= su_end; ++su ) {  // loop over spins in symmetric subunits (if symmetric_pre_calc_ is true)
			// otherwise middle vector in spin_coordinates_ has dimension 1
			if ( pre_multiset_ptr->ave_type_ == MEAN ) {
				// calculate average PRE for a group of degenerate spins (e.g. Methyls) and calculate the residual
				if ( pre_multiset_ptr->pre_singleset_vec_[j]->get_pre_rate_type() == R1_PARA ) {
					fvec[i-1] -= pre_multiset_ptr->R1_Para(pre_params, pre_multiset_ptr->one_over_r6_values_[i][su][1], pre_multiset_ptr->s2_values_[i][su][1], scal);
				} else if ( pre_multiset_ptr->pre_singleset_vec_[j]->get_pre_rate_type() == R2_PARA ) {
					fvec[i-1] -= pre_multiset_ptr->R2_Para(pre_params, pre_multiset_ptr->one_over_r6_values_[i][su][1], pre_multiset_ptr->s2_values_[i][su][1], scal);
				}
			} else if ( pre_multiset_ptr->ave_type_ == SUM ) {
				// if spins in the third dimension are not equivalent but e.g. symmetric spins sum up their PRE and calculate the residual
				for ( Size kk = 1, kk_end = pre_multiset_ptr->spin_coordinates_[i][su].size(); kk <= kk_end; ++kk ) {
					if ( pre_multiset_ptr->pre_singleset_vec_[j]->get_pre_rate_type() == R1_PARA ) {
						fvec[i-1] -= pre_multiset_ptr->R1_Para(pre_params, pre_multiset_ptr->one_over_r6_values_[i][su][kk], pre_multiset_ptr->s2_values_[i][su][kk], scal);
					} else if ( pre_multiset_ptr->pre_singleset_vec_[j]->get_pre_rate_type() == R2_PARA ) {
						fvec[i-1] -= pre_multiset_ptr->R2_Para(pre_params, pre_multiset_ptr->one_over_r6_values_[i][su][kk], pre_multiset_ptr->s2_values_[i][su][kk], scal);
					}
				}

			} // number equivalent spins

		} // number subunits
		// Weight by single PRE value weight
		//fvec[i-1] *= std::sqrt(pre_multiset_ptr->pre_single_weights_[i]);
		fvec[i-1] *= pre_multiset_ptr->pre_single_weights_[i];

	} // total number number PREs
}

/// @brief pre error function used in the lmmin function to optimize tau and the para ion position
///        * par is an array of fit parameters [xM, yM, zM, tau_c]
///        * data is a pointer to the PREMultiSet object i.e. to all data needed
///        for PRE calculation and NLS fitting
///        * fvc is an array holding the residuals of the fit calculation
void
pre_erf_opt_tau_xyz(
	Real const *par,
	int m_dat,
	void const *data,
	Real *fvec,
	int */*info*/
)
{
	PREMultiSet * pre_multiset_ptr = static_cast< PREMultiSet * >(const_cast< void* >(data));
	Real * nonconst_par = const_cast< Real* >(par);

	// Make sure that ion type and grid search iterator are set
	if ( !pre_multiset_ptr->ion_type_ ) {
		utility_exit_with_message("ERROR while trying to calculate the PRE in PREMultiSet. Ion type is not set.");
	}
	if ( !pre_multiset_ptr->gridsearch_iterator_ ) {
		utility_exit_with_message("ERROR while trying to calculate the PRE in PREMultiSet. GridSearch iterator is not set.");
	}

	// Constrain the para ion coordinates to the boundaries of the grid search box
	Real xMin(pre_multiset_ptr->gridsearch_iterator_->get_grid_search_center().x() - pre_multiset_ptr->gridsearch_iterator_->get_grid_max_radius());
	Real xMax(pre_multiset_ptr->gridsearch_iterator_->get_grid_search_center().x() + pre_multiset_ptr->gridsearch_iterator_->get_grid_max_radius());
	Real yMin(pre_multiset_ptr->gridsearch_iterator_->get_grid_search_center().y() - pre_multiset_ptr->gridsearch_iterator_->get_grid_max_radius());
	Real yMax(pre_multiset_ptr->gridsearch_iterator_->get_grid_search_center().y() + pre_multiset_ptr->gridsearch_iterator_->get_grid_max_radius());
	Real zMin(pre_multiset_ptr->gridsearch_iterator_->get_grid_search_center().z() - pre_multiset_ptr->gridsearch_iterator_->get_grid_max_radius());
	Real zMax(pre_multiset_ptr->gridsearch_iterator_->get_grid_search_center().z() + pre_multiset_ptr->gridsearch_iterator_->get_grid_max_radius());

	if ( !(xMin <= par[0]) || !(par[0] <= xMax) ) {
		nonconst_par[0] = 0.5*(xMax-xMin)*std::tanh(par[0])+(xMax+xMin)/2.0;
	}
	if ( !(yMin <= par[1]) || !(par[1] <= yMax) ) {
		nonconst_par[1] = 0.5*(yMax-yMin)*std::tanh(par[1])+(yMax+yMin)/2.0;
	}
	if ( !(zMin <= par[2]) || !(par[2] <= zMax) ) {
		nonconst_par[2] = 0.5*(zMax-zMin)*std::tanh(par[2])+(zMax+zMin)/2.0;
	}

	// Constrain correlation times to reasonable ranges
	// tau_c (in range tc_min - tc_max)
	if ( !( pre_multiset_ptr->tau_c_min_ <= par[3] ) || !( par[3] <= pre_multiset_ptr->tau_c_max_ ) ) {
		nonconst_par[3] = 0.5*(pre_multiset_ptr->tau_c_max_ - pre_multiset_ptr->tau_c_min_)*std::tanh(par[3])
			+(pre_multiset_ptr->tau_c_max_ + pre_multiset_ptr->tau_c_min_)/2.0;
	}

	PREMultiSet::Vec6 pre_params;
	// params[1] = prefac_dd
	// params[2] = prefac_curie
	// params[3] = omega_I
	// params[4] = tau_r
	// params[5] = tau_c
	// params[6] = tau_t

	PREMultiSet::WeightCoordVector paraion_position(1, std::pair<Real,Vector>(1.0, Vector(nonconst_par[0], nonconst_par[1], nonconst_par[2])));

	Real scal(1.0); // scaling factor per PRE experiment (normalization by experiment StdDev)
	bool use_only_sb(true); // Since we have only one spinlabel atom (para ion), we use the simplified SB equation.

	// We need to keep track of with which index a new experiment in the PRE value and weights vector starts
	// so that we can update the parameter set for the PRE equation. The start index is 1 plus
	// the number of PREs of the previous experiments ( e.g. 1, 101, 201, ... ).
	utility::vector1< Size > start_index(pre_multiset_ptr->number_experiments_);
	for ( Size i = 1; i <= pre_multiset_ptr->number_experiments_; ++i ) {
		if ( i == 1 ) { start_index[i] = 1; }
		else { start_index[i] = start_index[i-1] + pre_multiset_ptr->pre_singleset_vec_[i-1]->get_number_pre(); }
	}
	runtime_assert(start_index.back() == ( pre_multiset_ptr->total_number_pres_ - pre_multiset_ptr->pre_singleset_vec_.back()->get_number_pre() + 1));

	// i is the index for the ith PRE in the pre_values_ vector
	// j is the index for the jth experiment in the PRESingleSet vector
	// k is the index for the start_index vector (to figure out when to change the Params for the PRE eq.)
	for ( Size i = 1, j = 0, k = 1, i_end = static_cast< Size >(m_dat); i <= i_end; ++i ) { // loop over the total number of PREs
		fvec[i-1] = pre_multiset_ptr->pre_values_[i];
		// update parameters for next experiment and increase index j
		if ( i == start_index[k] ) {
			++j; // increment index for PRESingleSet vector
			k < pre_multiset_ptr->number_experiments_ ? ++k : k = pre_multiset_ptr->number_experiments_;

			Real gamma_I(pre_multiset_ptr->pre_singleset_vec_[j]->calc_gamma_I());
			Real gJ(pre_multiset_ptr->ion_type_->get_gJ());
			Real S(pre_multiset_ptr->ion_type_->get_S());
			Real B0(pre_multiset_ptr->pre_singleset_vec_[j]->get_field_strength() / 42.576); // B0 in Tesla

			pre_params[1] = pre_multiset_ptr->PRE_DD_prefactor(gamma_I, gJ, S, pre_multiset_ptr->pre_singleset_vec_[j]->get_pre_rate_type());
			pre_params[2] = pre_multiset_ptr->PRE_Curie_prefactor(gamma_I, gJ, S, B0, pre_multiset_ptr->get_temperature(), pre_multiset_ptr->pre_singleset_vec_[j]->get_pre_rate_type());
			pre_params[3] = pre_multiset_ptr->pre_singleset_vec_[j]->calc_omega_I();
			pre_params[4] = pre_multiset_ptr->tau_r_;
			pre_params[5] = nonconst_par[3];
			pre_params[6] = nonconst_par[3]; // Note that tau_t isn't actually used and we set it to tau_c

			scal = 1.0/pre_multiset_ptr->pre_singleset_vec_[j]->get_scaling_factor();
		}

		for ( Size su = 1, su_end = pre_multiset_ptr->spin_coordinates_[i].size(); su <= su_end; ++su ) {  // loop over spins in symmetric subunits (if symmetric_pre_calc_ is true)
			// otherwise middle vector in spin_coordinates_ has dimension 1
			if ( pre_multiset_ptr->ave_type_ == MEAN ) {

				// First, calculate <r-6> and S2 from new para ion xyz coordinates
				pre_multiset_ptr->calc_r6_S2(paraion_position, pre_multiset_ptr->spin_coordinates_[i][su],
					pre_multiset_ptr->one_over_r6_values_[i][su][1], pre_multiset_ptr->s2_values_[i][su][1]);

				// Second, calculate the residual
				// Here, we calculate the average PRE for a group of equivalent spins (e.g. Methyls)
				if ( pre_multiset_ptr->pre_singleset_vec_[j]->get_pre_rate_type() == R1_PARA ) {
					fvec[i-1] -= pre_multiset_ptr->R1_Para(pre_params, pre_multiset_ptr->one_over_r6_values_[i][su][1],
						pre_multiset_ptr->s2_values_[i][su][1], scal, use_only_sb);

				} else if ( pre_multiset_ptr->pre_singleset_vec_[j]->get_pre_rate_type() == R2_PARA ) {
					fvec[i-1] -= pre_multiset_ptr->R2_Para(pre_params, pre_multiset_ptr->one_over_r6_values_[i][su][1],
						pre_multiset_ptr->s2_values_[i][su][1], scal, use_only_sb);
				}
			} else if ( pre_multiset_ptr->ave_type_ == SUM ) {
				for ( Size kk = 1, kk_end = pre_multiset_ptr->spin_coordinates_[i][su].size(); kk <= kk_end; ++kk ) {

					// First, calculate <r-6> and S2 from new para ion xyz coordinates
					pre_multiset_ptr->calc_r6_S2(paraion_position, PREMultiSet::CoordVector(1, pre_multiset_ptr->spin_coordinates_[i][su][kk]),
						pre_multiset_ptr->one_over_r6_values_[i][su][kk], pre_multiset_ptr->s2_values_[i][su][kk]);

					// Second, calculate the residual
					// Here, we sum up the PRE of the spins in the third dimension of spin_coordinates (if they are not equivalent but e.g. hold symmetric spins)
					if ( pre_multiset_ptr->pre_singleset_vec_[j]->get_pre_rate_type() == R1_PARA ) {
						fvec[i-1] -= pre_multiset_ptr->R1_Para(pre_params, pre_multiset_ptr->one_over_r6_values_[i][su][kk],
							pre_multiset_ptr->s2_values_[i][su][kk], scal, use_only_sb);
					} else if ( pre_multiset_ptr->pre_singleset_vec_[j]->get_pre_rate_type() == R2_PARA ) {
						fvec[i-1] -= pre_multiset_ptr->R2_Para(pre_params, pre_multiset_ptr->one_over_r6_values_[i][su][kk],
							pre_multiset_ptr->s2_values_[i][su][kk], scal, use_only_sb);
					}
				}

			} // number equivalent spins

		} // number subunits
		// Weight by single PRE value weight
		//fvec[i-1] *= std::sqrt(pre_multiset_ptr->pre_single_weights_[i]);
		fvec[i-1] *= pre_multiset_ptr->pre_single_weights_[i];

	} // total number number PREs
}

/// @brief calculates the PRE from a vector of spinlabel atom coordinates using
///        the modified Solomon-Bloembergen (SBMF) equation and returns the
///        weighted PRE score according to the single PRE value weighting scheme
/// @details performs NLS fitting of the correlation times (tau_c and tau_t)
Real
PREMultiSet::compute_pre_score_from_point_vector(WeightCoordVector const & points) {

	// Make sure that ion type is set
	if ( !ion_type_ ) {
		utility_exit_with_message("ERROR while trying to compute PREs in PREMultiSet. Ion type is not set.");
	}

	// First, update <r-6> and S2 vectors
	update_r6_S2_values(points);

	// Set initial values for fit parameters (tau_c, tau_t)
	Size number_fit_params(2);
	utility::fixedsizearray1<Real,2> fit_params;
	utility::fixedsizearray1<Real,2> fit_params_best;
	// tau_c_min <= tau_c <= tau_c_max
	fit_params[1] = tau_c_min_ + (tau_c_max_ - tau_c_min_) * numeric::random::uniform();
	// 0.0 <= tau_t <= tau_c
	fit_params[2] = fit_params[1] * numeric::random::uniform();
	// Set best fit params to initial guess so that they are not uninitialized in case lmmin fails
	fit_params_best=fit_params;

	// definition of auxiliary parameters for lmmin
	numeric::nls::lm_status_struct status;
	Real bestnorm = std::numeric_limits< Real >::infinity();

	// Perform nonlinear least squares fitting
	for ( Size i = 1; i <= nls_repeats_; ++i ) {
		TR.Trace << "PRE NLS fitting" << std::endl;
		// Random starting values
		fit_params[1] = tau_c_min_ + (tau_c_max_ - tau_c_min_) * numeric::random::uniform();
		fit_params[2] = fit_params[1] * numeric::random::uniform();
		TR.Trace << "Initialize tau_c and tau_t to random values: " << fit_params[1] << " " << fit_params[2] << std::endl;
		numeric::nls::lmmin( number_fit_params, &fit_params[1], pre_values_.size(), (const void*) this, pre_erf_opt_tau, &status, numeric::nls::lm_printout_std);
		TR.Trace << "Iteration: " << i << " status.fnorm: " << status.fnorm << " bestnorm: "<< bestnorm << std::endl;
		//save to best fitting parameter
		if ( status.fnorm < bestnorm ) {
			TR.Trace << "status.fnorm: " << status.fnorm << " replaced bestnorm: " << bestnorm << std::endl;
			bestnorm=status.fnorm;
			fit_params_best = fit_params;
		}
	}
	bool use_sb(points.size() > 1 ? false : true);
	// Set correlation times as obtained by fit
	set_tau_c(fit_params_best[1]);
	if ( !use_sb ) {
		set_tau_t(fit_params_best[2]);
	} else {
		set_tau_t(fit_params_best[1]);
	}

	// Calculate PRE and score from current r6, S2 vectors and fitted correlation times
	Real total_score=compute_pre_score_from_current_data(use_sb);
	if ( TR.Trace.visible() ) {
		TR.Trace << "PRE score for " << number_experiments_ << " experiment(s) with " << (spinlabel_ ? "spinlabel "+spinlabel_->get_code() : "gridsearch")
			<< " at position " << spinlabel_site_rsd_ << " after NLS fit: " << total_score << std::endl;
		if ( spinlabel_ ) {
			spinlabel_->show(TR);
		}
		TR.Trace << "tau_r     (nsec) = " << std::scientific << tau_r_*1.0e+9 << std::endl;
		TR.Trace << "tau_c     (nsec) = " << std::scientific << tau_c_*1.0e+9 << std::endl;
		TR.Trace << "tau_t     (nsec) = " << std::scientific << tau_t_*1.0e+9 << std::endl;
	}
	return total_score;
}

/// @brief calculates the PRE from a single spinlabel atom position using the simplified
///        Solomon-Bloembergen (SB) equation and returns the weighted PRE score
///        according to the single PRE value weighting scheme
/// @details performs NLS fitting of the correlation time (tau_c) by default
///          optionally, the point position can be optimized too.
Real
PREMultiSet::compute_pre_score_from_single_point(
	Vector & point,
	bool opt_sl_pos)
{
	Real total_score(0);
	if ( opt_sl_pos ) {
		if ( !gridsearch_iterator_ ) {
			TR.Warning << "Optimization of single spinlabel position only possible in combination with grid search, but no NMRGridSeach object is set." << std::endl;
			TR.Warning << "Compute PRE from single spinlabel position without further optimization." << std::endl;

			total_score += compute_pre_score_from_point_vector(WeightCoordVector(1, std::pair<Real,Vector>(1.0,point)));

		} else {
			// Make sure that ion type is set
			if ( !ion_type_ ) {
				utility_exit_with_message("ERROR while trying to compute PREs in PREMultiSet. Ion type is not set.");
			}

			// Set initial values for fit parameters (xM, yM, zM, tau_c)
			Size number_fit_params(4);
			utility::fixedsizearray1<Real,4> fit_params;
			utility::fixedsizearray1<Real,4> fit_params_best;

			fit_params[1] = point.x();
			fit_params[2] = point.y();
			fit_params[3] = point.z();
			// tau_c_min <= tau_c <= tau_c_max
			fit_params[4] = tau_c_min_ + (tau_c_max_ - tau_c_min_) * numeric::random::uniform();
			// Set best fit params to initial guess so that they are not uninitialized in case lmmin fails
			fit_params_best=fit_params;

			// definition of auxiliary parameters for lmmin
			numeric::nls::lm_status_struct status;
			Real bestnorm = std::numeric_limits< Real >::infinity();

			// Perform nonlinear least squares fitting
			for ( Size i = 1; i <= nls_repeats_; ++i ) {
				TR.Trace << "PRE NLS fitting" << std::endl;
				// Random starting values
				fit_params[4] = tau_c_min_ + (tau_c_max_ - tau_c_min_) * numeric::random::uniform();
				TR.Trace << "Initialize tau_c to random value: " << fit_params[4] << std::endl;
				numeric::nls::lmmin( number_fit_params, &fit_params[1], pre_values_.size(), (const void*) this, pre_erf_opt_tau_xyz, &status, numeric::nls::lm_printout_std);
				TR.Trace << "Iteration: " << i << " status.fnorm: " << status.fnorm << " bestnorm: "<< bestnorm << std::endl;
				//save to best fitting parameter
				if ( status.fnorm < bestnorm ) {
					TR.Trace << "status.fnorm: " << status.fnorm << " replaced bestnorm: " << bestnorm << std::endl;
					bestnorm=status.fnorm;
					fit_params_best = fit_params;
				}
			}
			// Set spinlabel atom position as determined by the fit
			point.x() = fit_params_best[1];
			point.y() = fit_params_best[2];
			point.z() = fit_params_best[3];
			gridsearch_iterator_->set_best_grid_point(point);

			// Set correlation times as obtained by fit
			set_tau_c(fit_params_best[4]);
			set_tau_t(fit_params_best[4]);

			// Calculate PRE and score from current r6, S2 vectors and fitted correlation times
			update_r6_S2_values(WeightCoordVector(1,std::pair<Real,Vector>(1,point)));
			bool use_sb(true);
			total_score += compute_pre_score_from_current_data(use_sb);

			if ( TR.Trace.visible() ) {
				TR.Trace << "PRE score for " << number_experiments_ << " experiment(s) with " << (spinlabel_ ? "spinlabel "+spinlabel_->get_code() : "gridsearch")
					<< " at position " << spinlabel_site_rsd_ << " after NLS fit: " << total_score << std::endl;
				if ( spinlabel_ ) {
					spinlabel_->show(TR);
				}
				TR.Trace << "tau_r     (nsec) = " << std::scientific << tau_r_*1.0e+9 << std::endl;
				TR.Trace << "tau_c     (nsec) = " << std::scientific << tau_c_*1.0e+9 << std::endl;
				TR.Trace << "tau_t     (nsec) = " << std::scientific << tau_t_*1.0e+9 << std::endl;
			}
		}
	} else {
		total_score += compute_pre_score_from_point_vector(WeightCoordVector(1, std::pair<Real,Vector>(1,point)));
	}
	return total_score;
}

/// @brief determines the para ion/radical atom position of the spinlabel using a gridsearch or
///        the dummy spinlabel conformer ensemble and calculates the PRE score using
///        the modified Solomon-Bloembergen (SBMF) equation
Real
PREMultiSet::find_para_ion_position_and_compute_pre_score(
	pose::Pose const & pose,
	WeightCoordVector & spinlabel_counts_coords)
{
	Real total_score(0);

	update_spin_coordinates(pose);

	if ( spinlabel_ ) {
		// Note that we are using a simple distance calculation to the NBR_ATOM for spinlabel conformer filtering here.
		if ( spinlabel_->get_highres_conformer_filter_type() != NMRSpinlabel::DISTANCE ) {
			TR.Debug << "Will compute PRE score in PREMultiSet by filtering NMRDummySpinlabelConformers based on distance to neighborhood residues. ";
			TR.Debug << "Setting the conformer filter type to anything else than DISTANCE will be ignored." << std::endl;
		}

		// Filter dummy ensemble by neighbor count calculation and get the spinlabel coordinates
		// We can simply call the member method of the NMRSpinlabel. Internally, it performs the
		// clash filter, looks up the coordinates of the radical atom and transforms their coordinates
		// into the coordinate frame of the spinlabel site. It also performs clustering if the number of
		// spinlabel conformers exceeds the maximum ensemble size
		runtime_assert_msg( spinlabel_->get_dummy_ensemble(), "ERROR while trying to calculate PRE score. NMRDummySpinlabelEnsemble not set. Check if this spinlabel type has database file in Rosetta database." );
		WeightCoordVector spinlabel_atom_coords = spinlabel_->filter_spinlabel_ensemble_by_distance_check(pose, spinlabel_site_rsd_);

		// Finally calculate the score
		total_score = compute_pre_score_from_point_vector(spinlabel_atom_coords);
		spinlabel_counts_coords = spinlabel_atom_coords;

	} else if ( gridsearch_iterator_ ) {
		Vector current_paraion_coords;
		Real current_score(0);
		Real best_score(std::numeric_limits<Real>::max());
		gridsearch_iterator_->set_grid_search_center( pose );
		Vector best_paraion_coords(gridsearch_iterator_->get_grid_search_center());

		// Perform NLS fitting of para ion position and correlation times
		if ( opt_para_ion_pos_ ) {
			total_score = compute_pre_score_from_single_point(best_paraion_coords, true);
			gridsearch_iterator_->set_best_grid_point(best_paraion_coords);
		} else {
			// Run the gridsearch of para ion position and fit only the correlation times
			while ( gridsearch_iterator_->valid_next_grid_point(current_paraion_coords) ) {
				current_score = compute_pre_score_from_single_point(current_paraion_coords, false);
				if ( current_score < best_score ) {
					best_score = current_score;
					best_paraion_coords = current_paraion_coords;
					gridsearch_iterator_->set_best_grid_point(best_paraion_coords);
				}
			}
			total_score = compute_pre_score_from_single_point(best_paraion_coords, false);
		}
		// Vector of counts and radical atom xyz coordinates. Note that with the gridsearch method we consider only one point instead of a vector of coordinates.
		spinlabel_counts_coords = WeightCoordVector(1, std::pair<Real,Vector>(1,best_paraion_coords));
	} else {
		utility_exit_with_message("ERROR while trying to calculate PREMultiSet score. Spinlabel and gridsearch iterator are not set.");
	}
	return total_score;
}

/// @brief calculates and sets the xyz derivative of the PRE
///        single PREs must be calculated in beforehand
///        that's why make sure to call any of the three functions first:
///        compute_pre_score_from_point_vector(), compute_pre_from_single_point() or
///        find_para_ion_position_and_compute_pre_score()
void
PREMultiSet::set_atom_derivatives(
	pose::Pose & pose,
	WeightCoordVector const & spinlabel_counts_coords
)
{
	using pose::symmetry::is_symmetric;
	using pose::symmetry::symmetry_info;

	// Make sure that ion type is set
	if ( !ion_type_ ) {
		utility_exit_with_message("ERROR while trying to compute PRE derivatives in PREMultiSet. Ion type is not set.");
	}

	// Is pose symmetric
	conformation::symmetry::SymmetryInfoCOP syminfo_ptr;
	Size num_subunits(1);
	if ( is_symmetric( pose ) ) {
		syminfo_ptr = symmetry_info( pose );
		debug_assert(syminfo_ptr);
		num_subunits = syminfo_ptr->subunits();
	}
	Size index_offset(0);
	Vec6 pre_params;
	Real scal(1.0);
	Vector dr6(0.0), dS2(0.0);
	Real one_over_r6(0.0), S2(0.0);
	bool use_only_sb(spinlabel_ ? false : true);

	// Calculate and set xyz derivatives
	for ( Size i = 1; i <= number_experiments_; ++i ) {
		Size no_pres(pre_singleset_vec_[i]->get_number_pre());

		Real gamma_I(pre_singleset_vec_[i]->calc_gamma_I());
		Real gJ(ion_type_->get_gJ());
		Real S(ion_type_->get_S());
		Real B0(pre_singleset_vec_[i]->get_field_strength() / 42.576);

		pre_params[1] = PRE_DD_prefactor(gamma_I, gJ, S, pre_singleset_vec_[i]->get_pre_rate_type());
		pre_params[2] = PRE_Curie_prefactor(gamma_I, gJ, S, B0, temperature_, pre_singleset_vec_[i]->get_pre_rate_type());
		pre_params[3] = pre_singleset_vec_[i]->calc_omega_I();
		pre_params[4] = tau_r_;
		pre_params[5] = tau_c_;
		pre_params[6] = tau_t_;

		scal = 1.0 / pre_singleset_vec_[i]->get_scaling_factor();

		for ( Size j = 1; j <= no_pres; ++j ) {
			Real calc_pre = pre_singleset_vec_[i]->get_pre_single_vec()[j].get_pre_calc();
			Real diff = calc_pre - pre_values_[index_offset + j];

			// Here we iterate over the number of equivalent (degenerate) spins that produce the same PRE
			for ( Size k = 1, k_end = pre_singleset_vec_[i]->get_pre_single_vec()[j].get_protein_spins().size(); k <= k_end; ++k ) {
				Vector dPRE_dXYZ(0.0, 0.0, 0.0);

				// If symmetric_pre_calc_ is true i.e. we deduce the symmetric spins automatically
				// we iterate over the number of symmetric subunits and transform the spinlabel coordinates
				// into the frame of the other subunits to calculate the derivative with
				// respect to all spinlabels in the whole structure that give rise to one PRE value
				Size residue_from = pre_singleset_vec_[i]->get_pre_single_vec()[j].get_protein_spins()[k].rsd();
				utility::vector1< Size > residues_to(1, residue_from);

				if ( is_symmetric( pose ) && syminfo_ptr ) {
					residues_to.resize(num_subunits);
					if ( !syminfo_ptr->is_asymmetric_seqpos(residue_from) ) {
						residue_from = syminfo_ptr->get_asymmetric_seqpos(residue_from);
						residues_to[1] = residue_from;
					}
					utility::vector1< Size > symm_rsds = syminfo_ptr->bb_clones(residue_from);
					std::copy(symm_rsds.begin(), symm_rsds.end(), residues_to.begin()+1);
				}
				for ( Size su = 1; su <= num_subunits; ++su ) {
					// Get symmetric spinlabel coordinates
					WeightCoordVector spinlabel_counts_coords_symm;
					for ( WeightCoordVector::const_iterator iter = spinlabel_counts_coords.begin(),
							last  = spinlabel_counts_coords.end(); iter != last; ++iter ) {
						spinlabel_counts_coords_symm.push_back( std::make_pair( (*iter).first, apply_vector_transformation( pose, (*iter).second, residue_from, residues_to[su] ) ) );
					}

					// Here we have to switch behavior again, depending on if the coordinates in the inner vector of spin_coordinates_
					// refer to the degenerate spins (for which we perform mean averaging of their PRE) or symmetric spins (for which
					// we perform summation of their PRE)
					if ( ave_type_ == MEAN ) {
						// Calculate <r-6>, S2 and their derivatives.
						// We cannot reuse the private data because this time
						// we calculate the PRE (and its derivative) relative
						// to the spinlabel located on the symmetric subunits.
						calc_r6_S2(spinlabel_counts_coords_symm, spin_coordinates_[index_offset + j][1], one_over_r6, S2);
						calc_dr6_dS2_dXYZ(spinlabel_counts_coords_symm, spin_coordinates_[index_offset + j][1], dr6, dS2);
						if ( pre_singleset_vec_[i]->get_pre_rate_type() == R1_PARA ) {
							dPRE_dXYZ += dR1_dXYZ(pre_params, one_over_r6, S2, dr6, dS2, scal, use_only_sb);

						} else if ( pre_singleset_vec_[i]->get_pre_rate_type() == R2_PARA ) {
							dPRE_dXYZ += dR2_dXYZ(pre_params, one_over_r6, S2, dr6, dS2, scal, use_only_sb);

						}
					} else if ( ave_type_ == SUM ) {
						calc_r6_S2(spinlabel_counts_coords_symm, PREMultiSet::CoordVector(1, spin_coordinates_[index_offset + j][1][k]), one_over_r6, S2);
						calc_dr6_dS2_dXYZ(spinlabel_counts_coords_symm, PREMultiSet::CoordVector(1, spin_coordinates_[index_offset + j][1][k]), dr6, dS2);
						if ( pre_singleset_vec_[i]->get_pre_rate_type() == R1_PARA ) {
							dPRE_dXYZ += dR1_dXYZ(pre_params, one_over_r6, S2, dr6, dS2, scal, use_only_sb);

						} else if ( pre_singleset_vec_[i]->get_pre_rate_type() == R2_PARA ) {
							dPRE_dXYZ += dR2_dXYZ(pre_params, one_over_r6, S2, dr6, dS2, scal, use_only_sb);

						}
					}
					dPRE_dXYZ *= (diff * pre_single_weights_[index_offset + j]);
				} // number of subunits

				pre_singleset_vec_[i]->pre_single_vec_[j].set_atom_derivatives(k, dPRE_dXYZ.x(), dPRE_dXYZ.y(), dPRE_dXYZ.z());
			} // equivalent spins

		} // number of PREs per experiment
		// increment index offset
		index_offset += no_pres;

	} // number of experiments
}

void
PREMultiSet::show(std::ostream & tracer) const {
	// Store old iostream manipulator flags
	std::ios oldState(nullptr);
	oldState.copyfmt(tracer);

	tracer << "   * * * PREMultiSet Summary Report * * *   " << std::endl;
	tracer << "Spinlabel site: " << spinlabel_site_rsd_ << std::endl;
	tracer << "No experiments: " << number_experiments_ << std::endl;
	tracer << "No PREs:        " << total_number_pres_ << std::endl;
	tracer << "PRE Spinlabel:  ";
	if ( spinlabel_ ) {
		tracer << std::endl;
		tracer << " - Name           = " << spinlabel_->get_name() << std::endl;
		tracer << " - 3-Letter code  = " << spinlabel_->get_code() << std::endl;
		tracer << " - Radical ion    = " << spinlabel_->get_radical_atom() << std::endl;
		tracer << " - Ensemble size  = " << spinlabel_->get_current_ensemble_size() << std::endl;
	} else {
		tracer << "None" << std::endl;
	}
	tracer << "Experimental conditions: " << std::endl;
	for ( Size i = 1; i <= number_experiments_; ++i ) {
		tracer << " - " << pre_singleset_vec_[i]->get_dataset_name() << ": [ " << pre_singleset_vec_[i]->get_number_pre() << " PREs, "
			<< pre_singleset_vec_[i]->pre_rate_type_to_string() << ", " << pre_singleset_vec_[i]->get_field_strength() << " MHz ]" << std::endl;
	}
	tracer << std::setprecision(2) << std::fixed << std::right;
	tracer << " - MW     (kDa)   = " << std::setw(9) << protein_mass_ << std::endl;
	tracer << " - T      (K)     = " << std::setw(9) << temperature_ << std::endl;
	tracer << "Optimized correlation time: " << std::endl;
	tracer << " - tau_r  (nsec)  = " << std::setw(9) << std::scientific << tau_r_*1.0e+9 << std::endl;
	tracer << " - tau_c  (nsec)  = " << std::setw(9) << std::scientific << tau_c_*1.0e+9 << std::endl;
	tracer << " - tau_t  (nsec)  = " << std::setw(9) << std::scientific << tau_t_*1.0e+9 << std::endl;
	tracer.copyfmt(oldState);
}

void
PREMultiSet::set_averaging_type(std::string const & type) {
	ave_type_ = convert_string_to_averaging_type(type);
	resize_r6_S2_values();
}

/// @brief register options
void
PREMultiSet::register_options() {
	using namespace basic::options;
	option.add_relevant(OptionKeys::nmr::pre::use_symmetry_calc);
	option.add_relevant(OptionKeys::nmr::pre::nls_repeats);
	option.add_relevant(OptionKeys::nmr::pre::optimize_para_ion_position);
}

void
PREMultiSet::init_from_cml() {
	using namespace basic::options;
	symmetric_pre_calc_ = option[ basic::options::OptionKeys::nmr::pre::use_symmetry_calc ]();
	nls_repeats_ = option[ basic::options::OptionKeys::nmr::pre::nls_repeats ]();
	opt_para_ion_pos_ = option[ basic::options::OptionKeys::nmr::pre::optimize_para_ion_position ]();
}

void
PREMultiSet::deep_copy_pre_singleset_vec(utility::vector1<PRESingleSetOP> const & other_vec) {
	pre_singleset_vec_.resize(other_vec.size());
	for ( Size i = 1; i <= pre_singleset_vec_.size(); ++i ) {
		pre_singleset_vec_[i] = PRESingleSetOP( new PRESingleSet( *(other_vec[i]) ) );
	}
}

void
PREMultiSet::create_para_ion_from_string(std::string const & iontype) {
	using namespace core::io::nmr;
	std::map< std::string, ParaIon > const & ion_data_table = ParamagneticDatabaseHandler::get_instance()->get_ion_data_table();
	std::map< std::string, ParaIon >::const_iterator iter = ion_data_table.find( iontype );
	if ( iter == ion_data_table.end() ) {
		utility_exit_with_message("ERROR while trying to create paramagnetic ion of PREMultiSet. " + iontype + " not found in paramagnetic ion database.");
	}
	ion_type_ = ParaIonOP( new ParaIon( (*iter).second));
}

/// @brief set default values for tau_r_, tau_c_, tau_c_min_ and tau_c_max_
void
PREMultiSet::set_exp_default_conditions() {
	protein_mass_ = 15.0;
	temperature_ = 298.0;
	tau_r_ = calc_theoretical_tau_r();
	tau_c_ = calc_tau_c();
	tau_c_min_ = 0.1 * tau_r_;
	tau_c_max_ = 10.0* tau_r_;
	tau_t_ = tau_c_;
}

/// @brief creates MTSL spinlabel as default
void
PREMultiSet::set_default_spinlabel() {
	spinlabel_ = NMRSpinlabelOP( new NMRSpinlabel("fa_standard", "R1A") );
}

/// @brief calculate rotational correlation time (in sec) of the protein
///        using the Stokes-Einstein equation and the protein mass (in kDa == kg/mol)
Real
PREMultiSet::calc_theoretical_tau_r() const {
	Real const NA(6.022140857);   // 10^+23 (1/mol)
	Real const kb(1.38064852);   // 10^-23 (m^2*kg)/(s^2 * K)
	Real const rho(1300.0);    // (kg/m^3)
	Real const viscosity(0.0035);      // (kg/s*m)
	runtime_assert_msg(temperature_ > 0.0 , "ERROR in calculating tau r. Temperature must be a positive value.");
	runtime_assert_msg(protein_mass_ > 0.0, "ERROR in calculating tau r. Protein mass must be a positive value.");
	return (viscosity * protein_mass_)/(rho * NA * kb * temperature_);
}

/// @brief calculate sum of rotational correlation and electron relaxation time
///        i.e. the total correlation time (in sec)
Real
PREMultiSet::calc_tau_c() const {
	if ( !ion_type_ ) {
		utility_exit_with_message("ERROR while trying to calculate tau_c in PREMultiSet. Ion type is not set.");
	}
	Real tau_e = ion_type_->get_tau_e()*1.0e-012;
	runtime_assert_msg(tau_e > 0.0, "ERROR in calculating tau c. Electron relaxation time tau e must be a positive value.");
	return (1.0/(1.0/tau_r_ + 1.0/tau_e));
}

/// @brief calculate electron spin frequency at given field strength in rad/s
Real
PREMultiSet::calc_omega_E(Real const & field_strength) const {
	//Real const mu_b(5.788381e-5);  // Bohr magneton in eV/T
	//Real const hbar(6.582119e-16); // Planck constant in eV*s
	if ( !ion_type_ ) {
		utility_exit_with_message("ERROR while trying to calculate electron resonance frequency in PREMultiSet. Ion type is not set.");
	}
	Real const mu_b_over_h(8.804e+010); // mu_b / hbar
	Real gJ(ion_type_->get_gJ());
	runtime_assert_msg(field_strength > 0.0, "ERROR in calculating nuclear spin frequency. Magnetic field strength must be a positive value.");
	return gJ * mu_b_over_h * field_strength;
}

/// @brief computes r6 and S2 vectors from given spinlabel atom coordinates and
///        and current spin coordinates
void
PREMultiSet::update_r6_S2_values(WeightCoordVector const & spinlabel_counts_coords) {
	Size index_offset(0);
	for ( Size i = 1; i <= number_experiments_; ++i ) {
		Size no_pres(pre_singleset_vec_[i]->get_number_pre());
		for ( Size j = 1; j <= no_pres; ++j ) {
			for ( Size su = 1, su_end = spin_coordinates_[index_offset + j].size(); su <= su_end; ++su ) {
				if ( ave_type_ == MEAN ) {
					runtime_assert_msg(one_over_r6_values_[index_offset + j][su].size() == 1,
						"ERROR: Length of inner <r-6> value vector != 1. Incompatible with averaging type MEAN.");
					runtime_assert_msg(s2_values_[index_offset + j][su].size() == 1,
						"ERROR: Length of inner S2 value vector != 1. Incompatible with averaging type MEAN.");
					calc_r6_S2(spinlabel_counts_coords, spin_coordinates_[index_offset + j][su],
						one_over_r6_values_[index_offset + j][su][1], s2_values_[index_offset + j][su][1]);
				} else if ( ave_type_ == SUM ) {
					for ( Size k = 1, k_end = spin_coordinates_[index_offset + j][su].size(); k <= k_end; ++k ) {
						calc_r6_S2(spinlabel_counts_coords, PREMultiSet::CoordVector(1, spin_coordinates_[index_offset + j][su][k]),
							one_over_r6_values_[index_offset + j][su][k], s2_values_[index_offset + j][su][k]);
					}
				}
			}
		} // increment index offset
		index_offset += no_pres;
	}
}

/// @brief check if dimensionality of r6 and S2 vectors is consistent with averaging type
///        and resize vectors if not.
void
PREMultiSet::resize_r6_S2_values() {
	Size index_offset(0);
	for ( Size i = 1; i <= number_experiments_; ++i ) {
		Size no_pres = pre_singleset_vec_[i]->get_number_pre();
		for ( Size j = 1; j <= no_pres; ++j ) {
			for ( Size su = 1, su_end = spin_coordinates_[index_offset + j].size(); su <= su_end; ++su ) {
				if ( ave_type_ == MEAN ) {
					one_over_r6_values_[index_offset + j][su].resize(1);
					s2_values_[index_offset + j][su].resize(1);
				} else if ( ave_type_ == SUM ) {
					one_over_r6_values_[index_offset + j][su].resize(pre_singleset_vec_[i]->pre_single_vec_[j].get_protein_spins().size());
					s2_values_[index_offset + j][su].resize(pre_singleset_vec_[i]->pre_single_vec_[j].get_protein_spins().size());
				}
			}
		}
	}
}

/// @brief computes the PRE score using the current r6 and S2 vectors and correlation times
/// @params
/// use_sb:   use only the simplified SB equation (e.g. if there is only one spinlabel)
Real
PREMultiSet::compute_pre_score_from_current_data(bool use_sb) {
	Real total_score(0);
	Size index_offset(0);
	Vec6 pre_params;
	// params[1] = prefac_dd
	// params[2] = prefac_curie
	// params[3] = omega_I
	// params[4] = tau_r
	// params[5] = tau_c
	// params[6] = tau_t
	Real scal(1.0);

	for ( Size i = 1; i <= number_experiments_; ++i ) {
		Size no_pres(pre_singleset_vec_[i]->get_number_pre());
		Real singleset_score(0);

		Real gamma_I(pre_singleset_vec_[i]->calc_gamma_I());
		Real gJ(ion_type_->get_gJ());
		Real S(ion_type_->get_S());
		Real B0(pre_singleset_vec_[i]->get_field_strength() / 42.576);

		pre_params[1] = PRE_DD_prefactor(gamma_I, gJ, S, pre_singleset_vec_[i]->get_pre_rate_type());
		pre_params[2] = PRE_Curie_prefactor(gamma_I, gJ, S, B0, temperature_, pre_singleset_vec_[i]->get_pre_rate_type());
		pre_params[3] = pre_singleset_vec_[i]->calc_omega_I();
		pre_params[4] = tau_r_;
		pre_params[5] = tau_c_;
		pre_params[6] = tau_t_;

		scal = 1.0 / pre_singleset_vec_[i]->get_scaling_factor();

		// loop over PREs in experiment i
		for ( Size j = 1; j <= no_pres; ++j ) {
			Real calc_pre(0);
			// loop over spins in symmetric subunits (if symmetric_pre_calc_ is true)
			// otherwise middle vector in spin_coordinates_ has dimension 1
			for ( Size su = 1, su_end = spin_coordinates_[index_offset + j].size(); su <= su_end; ++su ) {

				if ( ave_type_ == MEAN ) {
					// calculate average PRE for a group of degenerate spins (e.g. Methyls)
					if ( pre_singleset_vec_[i]->get_pre_rate_type() == R1_PARA ) {
						calc_pre += R1_Para(pre_params, one_over_r6_values_[index_offset + j][su][1], s2_values_[index_offset + j][su][1], scal, use_sb);

					} else if ( pre_singleset_vec_[i]->get_pre_rate_type() == R2_PARA ) {
						calc_pre += R2_Para(pre_params, one_over_r6_values_[index_offset + j][su][1], s2_values_[index_offset + j][su][1], scal, use_sb);

					}
				} else if ( ave_type_ == SUM ) {
					// if spins in the third dimension are not equivalent but e.g. symmetric spins, sum up their PRE
					for ( Size k = 1, k_end = spin_coordinates_[index_offset + j][su].size(); k <= k_end; ++k ) {
						if ( pre_singleset_vec_[i]->get_pre_rate_type() == R1_PARA ) {
							calc_pre += R1_Para(pre_params, one_over_r6_values_[index_offset + j][su][k], s2_values_[index_offset + j][su][k], scal, use_sb);

						} else if ( pre_singleset_vec_[i]->get_pre_rate_type() == R2_PARA ) {
							calc_pre += R2_Para(pre_params, one_over_r6_values_[index_offset + j][su][k], s2_values_[index_offset + j][su][k], scal, use_sb);

						}
					} // number equivalent spins
				}
			} // number subunits
			pre_singleset_vec_[i]->pre_single_vec_[j].set_pre_calc(calc_pre);
			Real diff = calc_pre - pre_values_[index_offset + j];
			singleset_score += diff * diff * pre_single_weights_[index_offset + j];
		} // number PREs per experiment
		//total_score += std::sqrt(singleset_score) * pre_singleset_vec_[i]->get_weight();
		total_score += singleset_score * pre_singleset_vec_[i]->get_weight();

		// increment index offset
		index_offset += no_pres;

	} // number of PRE experiments
	return total_score;
}

/// @brief calculate the constant prefactor of the dipolar part of the PRE
/// @params
/// gamma_I:     gyromagnetic ratio of the nuclear spin (must be provided in rad/(s*T), dimension is 10^6)
/// gJ:          electron Lande factor
/// S:           total spin quantum number
Real
PREMultiSet::PRE_DD_prefactor(
	Real const & gamma_I,
	Real const & gJ,
	Real const & S,
	PRE_RATE_TYPE const & rate_type
) const
{
	// Some physical constants
	Real const mu_0 = 1.2566370614; // vacuum permeability (in N/A^2 * 10^-6)
	Real const mu_b = 9.274009994;  // Bohr magneton (in J/T * 10^-24)

	// Precompute some terms
	Real mu_0_2 = mu_0 * mu_0;
	Real mu_b_2 = mu_b * mu_b;
	Real pi_2 = numeric::constants::r::pi * numeric::constants::r::pi;
	Real gamma_I_2 = gamma_I * gamma_I;
	Real gJ_2 = gJ * gJ;
	Real C(mu_0_2/(16.0*pi_2) * gamma_I_2 * gJ_2 * mu_b_2 * S*(S+1));

	return rate_type == R1_PARA ? 2.0/15.0 * C : 1.0/15.0 * C;
}

/// @brief calculate the constant prefactor of the Curie part of the PRE
/// @params
/// gamma_I:     gyromagnetic ratio of the nuclear spin (must be provided in rad/(s*T), dimension is 10^6)
/// gJ:          electron Lande factor
/// S:           total spin quantum number
/// B0:          magnetic field strength (in Tesla)
/// T:           temperature (in Kelvin)
Real
PREMultiSet::PRE_Curie_prefactor(
	Real const & gamma_I,
	Real const & gJ,
	Real const & S,
	Real const & B0,
	Real const & T,
	PRE_RATE_TYPE const & rate_type
) const
{
	// Some physical constants
	Real const mu_0 = 1.2566370614; // vacuum permeability (in N/A^2 * 10^-6)
	Real const mu_b = 9.274009994;  // Bohr magneton (in J/T * 10^-24)
	Real const kb = 13.8064852;     // Boltzmann constant ((m^2 * kg)/(s^2 * K) * 10^-24)

	// Precompute some terms
	Real mu_0_2 = mu_0 * mu_0;
	Real mu_b_2 = mu_b * mu_b;
	Real mu_b_4 = mu_b_2 * mu_b_2;
	Real kb_2 = kb * kb;
	Real pi_2 = numeric::constants::r::pi * numeric::constants::r::pi;
	Real gamma_I_2 = gamma_I * gamma_I;
	Real gJ_2 = gJ * gJ;
	Real gJ_4 = gJ_2 * gJ_2;
	Real B0_2 = B0 * B0;
	Real T_2 = T * T;
	Real C(mu_0_2/(16.0*pi_2) * (gamma_I_2 * B0_2 * gJ_4 * mu_b_4 * S*S*(S+1)*(S+1)) / (9.0 * kb_2 * T_2));

	return rate_type == R1_PARA ? 2.0/5.0 * C : 1.0/5.0 * C;
}

/// @brief calculate the ensemble averaged distance <r^-6> and generalized order parameter S2
///        for a given spinlabel atom - nuclear spin connection vector as used in the SBMF equation
/// @params
/// spinlabel_counts_coords: vector of weights and para ion xyz coordinates of N distinct spinlabel conformers
/// eq_spins_coords:         vector of xyz coordinates of nuclear spin(s) that give rise to the observed PRE
/// one_over_r6:             averaged distance <r^-6> as return value
/// S2:                      generalized order parameter S2 as return value
void
PREMultiSet::calc_r6_S2(
	WeightCoordVector const & spinlabel_counts_coords,
	CoordVector const & eq_spins_coords,
	Real & one_over_r6,
	Real & S2
) const
{
	S2 = 0;
	Real S2_angular(0.0);
	Real S2_radial_numerator(0.0);
	Real S2_radial_denominator(0.0);
	Real S2_angular_per_spin(0.0);
	Real S2_radial_numerator_per_spin(0.0);
	Real S2_radial_denominator_per_spin(0.0);

	one_over_r6 = 0;
	Real one_over_r6_per_spin(0.0);

	// Sum of weights of spinlabel conformers
	Real total_sum_weights(0.0);
	// Center of spinlabel radical atom coordinates weighted by the number of occurrences
	Vector spinlabel_coords_center(0.0, 0.0, 0.0);

	for ( std::pair< Real, Vector > const & sl : spinlabel_counts_coords ) {
		total_sum_weights += sl.first;
		spinlabel_coords_center += (sl.first * sl.second);
	}
	spinlabel_coords_center /= total_sum_weights;

	// Center of equivalent spins
	Vector eq_spins_coords_center(numeric::center_of_mass(eq_spins_coords));
	// Vector and euclidean distance between the two centers
	Real Xc(spinlabel_coords_center.x() - eq_spins_coords_center.x());
	Real Yc(spinlabel_coords_center.y() - eq_spins_coords_center.y());
	Real Zc(spinlabel_coords_center.z() - eq_spins_coords_center.z());
	Real rc(std::sqrt(Xc*Xc + Yc*Yc + Zc*Zc));

	// Iterate over the number of spinlabel conformers (weighted by the number of their occurrences)
	for ( std::pair< Real, Vector > const & sl : spinlabel_counts_coords ) {
		Real xSL(sl.second.x());
		Real ySL(sl.second.y());
		Real zSL(sl.second.z());
		S2_angular_per_spin = 0.0;
		S2_radial_numerator_per_spin = 0.0;
		S2_radial_denominator_per_spin = 0.0;
		one_over_r6_per_spin = 0.0;

		// Iterate over equivalent spins
		for ( Size i = 1, i_end = eq_spins_coords.size(); i <= i_end; ++i ) {
			Real xS(eq_spins_coords[i].x());
			Real yS(eq_spins_coords[i].y());
			Real zS(eq_spins_coords[i].z());

			Real Xi(xSL - xS);
			Real Yi(ySL - yS);
			Real Zi(zSL - zS);
			Real ri(std::sqrt(Xi*Xi + Yi*Yi + Zi*Zi));
			Real ri3(ri * ri * ri);
			Real ri6(ri3 * ri3);
			Real cosine_ci = (Xc*Xi + Yc*Yi + Zc*Zi)/(rc * ri);

			S2_angular_per_spin += (1.5*(cosine_ci * cosine_ci) - 0.5);
			S2_radial_numerator_per_spin += (1.0 / ri3);
			S2_radial_denominator_per_spin += (1.0 / ri6);
			one_over_r6_per_spin += (1.0 / ri6);
		}
		S2_angular += (S2_angular_per_spin * sl.first);
		S2_radial_numerator += (S2_radial_numerator_per_spin * sl.first);
		S2_radial_denominator += (S2_radial_denominator_per_spin * sl.first);
		one_over_r6 += (one_over_r6_per_spin * sl.first);
	}
	S2_angular /= (eq_spins_coords.size() * total_sum_weights);
	Real S2_radial = (S2_radial_numerator * S2_radial_numerator)/(eq_spins_coords.size() * total_sum_weights * S2_radial_denominator);

	one_over_r6 /= (total_sum_weights * eq_spins_coords.size());
	S2 = S2_angular * S2_radial;
}

/// @brief spectral density function for the simplified Solomon-Bloembergen (SB) equation
/// @details A detailed description can be found in: Iwahara J, Schwieters CD & Clore GM, 2004, JACS 126,5879-5896
/// @params
/// omega_I:     nuclear spin resonance frequency (must be provided in rad/s, dimension is 10^6)
/// tau_c:       correlation time (must be provided in s, typical dimension is 10^-9)
/// one_over_r6: ensemble averaged distance <r^-6> (must be provided in Ang.^-6)
Real
PREMultiSet::J_SB(
	Real const & omega_I,
	Real const & tau_c,
	Real const & one_over_r6
) const
{
	return one_over_r6 * (tau_c / (1.0 + (omega_I*tau_c)*(omega_I*tau_c)));
}

/// @brief spectral density function for the modified Solomon-Bloembergen (SBMF) equation
/// @details A detailed description can be found in: Iwahara J, Schwieters CD & Clore GM, 2004, JACS 126,5879-5896
/// @params
/// omega_I:     nuclear spin resonance frequency (must be provided in rad/s, dimension is 10^6)
/// tau_c:       correlation time (must be provided in s, typical dimension is 10^-9)
/// tau_t:       total correlation time (must be provided in s)
/// one_over_r6: ensemble averaged distance <r^-6> (must be provided in Ang.^-6)
/// S2:          generalized order parameter
Real
PREMultiSet::J_SBMF(
	Real const & omega_I,
	Real const & tau_c,
	Real const & tau_t,
	Real const & one_over_r6,
	Real const & S2
) const
{
	return one_over_r6 * ( (S2 * tau_c) / (1.0 + (omega_I*tau_c)*(omega_I*tau_c)) + ((1.0-S2)*tau_t) / (1.0 + (omega_I*tau_t)*(omega_I*tau_t)));
}

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
PREMultiSet::R1_Para(
	Vec6 const & params,
	Real const & one_over_r6,
	Real const & S2,
	Real const & scaling,
	bool use_sb
) const
{
	// Dipolar part
	Real J_omega_tau_c_t = use_sb ? J_SB(params[3], params[5], one_over_r6) : J_SBMF(params[3], params[5], params[6], one_over_r6, S2);
	Real pre_dd = scaling * params[1] * 3.0 * J_omega_tau_c_t;

	// Curie part
	Real J_omega_tau_r = J_SB(params[3], params[4], one_over_r6);
	Real J_omega_tau_c = J_SB(params[3], params[5], one_over_r6);
	Real pre_curie = scaling * params[2] * ( 3.0 * J_omega_tau_r - 3.0 * J_omega_tau_c);

	return pre_dd + pre_curie;
}

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
PREMultiSet::R2_Para(
	Vec6 const & params,
	Real const & one_over_r6,
	Real const & S2,
	Real const & scaling,
	bool use_sb
) const
{
	// Dipolar part
	Real J_0_tau_c_t     = use_sb ? J_SB(0.0, params[5], one_over_r6) : J_SBMF(0.0, params[5], params[6], one_over_r6, S2);
	Real J_omega_tau_c_t = use_sb ? J_SB(params[3], params[5], one_over_r6) : J_SBMF(params[3], params[5], params[6], one_over_r6, S2);
	Real pre_dd          = scaling * params[1] * (4.0 * J_0_tau_c_t + 3.0 * J_omega_tau_c_t);

	// Curie part
	Real J_0_tau_r      = J_SB(0.0, params[4], one_over_r6);
	Real J_omega_tau_r  = J_SB(params[3], params[4], one_over_r6);
	Real J_0_tau_c      = J_SB(0.0, params[5], one_over_r6);
	Real J_omega_tau_c  = J_SB(params[3], params[5], one_over_r6);
	Real pre_curie      = scaling * params[2] * (4.0 * J_0_tau_r + 3.0 * J_omega_tau_r - 4.0 * J_0_tau_c - 3.0 * J_omega_tau_c);

	return pre_dd + pre_curie;
}

/// @brief calculate the xyz derivative of the ensemble averaged distance <r^-6>
///        and the generalized order parameter S2 for a given spinlabel atom -
///        nuclear spin connection vector
/// @params
/// spinlabel_counts_coords: vector of weights and para ion xyz coordinates of N distinct spinlabel conformers
/// eq_spins_coords:         vector of xyz coordinates of nuclear spin(s) that give rise to the observed PRE
/// dr6_dXYZ:                xyz derivative of <r^-6> as return value
/// dS2_dXYZ:                xyz derivative of S2 as return value
void
PREMultiSet::calc_dr6_dS2_dXYZ(
	WeightCoordVector const & spinlabel_counts_coords,
	CoordVector const & eq_spins_coords,
	Vector & dr6_dXYZ,
	Vector & dS2_dXYZ
) const
{
	// For calculating dS2_dXYZ and dr6_dXYZ, we need to calculate the derivatives of different sum terms
	// which we initialize to 0.0 here
	Real S2_angular(0.0);
	Real S2_angular_per_spin(0.0);
	Vector dS2_angular_dXYZ(0.0, 0.0, 0.0);
	Vector dS2_angular_dXYZ_per_spin(0.0, 0.0, 0.0);

	Real S2_radial_numerator(0.0);
	Real S2_radial_numerator_per_spin(0.0);
	Vector dS2_radial_numerator_dXYZ(0.0, 0.0, 0.0);
	Vector dS2_radial_numerator_dXYZ_per_spin(0.0, 0.0, 0.0);

	Real S2_radial_denominator(0.0);
	Real S2_radial_denominator_per_spin(0.0);
	Vector dS2_radial_denominator_dXYZ(0.0, 0.0, 0.0);
	Vector dS2_radial_denominator_dXYZ_per_spin(0.0, 0.0, 0.0);

	dr6_dXYZ = Vector(0.0, 0.0, 0.0);
	Vector dr6_dXYZ_per_spin(0.0, 0.0, 0.0);

	// Sum of weights of spinlabel conformers
	Real total_sum_weights(0.0);
	// Center of spinlabel radical atom coordinates weighted by the number of occurrences
	Vector spinlabel_coords_center(0.0, 0.0, 0.0);
	for ( std::pair< Real, Vector > const & sl : spinlabel_counts_coords ) {
		total_sum_weights += sl.first;
		spinlabel_coords_center += (sl.first * sl.second);
	}
	spinlabel_coords_center /= total_sum_weights;

	// Center of equivalent spins
	Vector eq_spins_coords_center(numeric::center_of_mass(eq_spins_coords));

	// vector and euclidean distance between the SL center and spin
	Real Xc(spinlabel_coords_center.x() - eq_spins_coords_center.x());
	Real Yc(spinlabel_coords_center.y() - eq_spins_coords_center.y());
	Real Zc(spinlabel_coords_center.z() - eq_spins_coords_center.z());
	Real Xc2(Xc * Xc);
	Real Yc2(Yc * Yc);
	Real Zc2(Zc * Zc);
	Real rc(std::sqrt(Xc2 + Yc2 + Zc2));
	Real rc2(rc * rc);

	// Iterate over the number of spinlabel conformers (weighted by the number of their occurrences)
	for ( std::pair< Real, Vector > const & sl : spinlabel_counts_coords ) {
		Real xSL(sl.second.x());
		Real ySL(sl.second.y());
		Real zSL(sl.second.z());

		dS2_angular_dXYZ_per_spin = 0.0;
		dS2_radial_numerator_dXYZ_per_spin = 0.0;
		dS2_radial_denominator_dXYZ_per_spin = 0.0;
		S2_angular_per_spin = 0.0;
		S2_radial_numerator_per_spin = 0.0;
		S2_radial_denominator_per_spin = 0.0;
		dr6_dXYZ_per_spin = 0.0;

		for ( Size i = 1, i_end = eq_spins_coords.size(); i <= i_end; ++i ) {
			Real xS(eq_spins_coords[i].x());
			Real yS(eq_spins_coords[i].y());
			Real zS(eq_spins_coords[i].z());

			Real Xi(xSL - xS);
			Real Yi(ySL - yS);
			Real Zi(zSL - zS);
			Real Xi2(Xi * Xi);
			Real Yi2(Yi * Yi);
			Real Zi2(Zi * Zi);

			Real ri(std::sqrt(Xi2 + Yi2 + Zi2));
			Real ri2(ri  * ri );
			Real ri3(ri  * ri2);
			Real ri5(ri2 * ri3);
			Real ri6(ri3 * ri3);
			Real ri8(ri3 * ri5);
			Real cosine_ci = (Xc*Xi + Yc*Yi + Zc*Zi)/(rc * ri);

			dS2_angular_dXYZ_per_spin.x() += ( (3.0 * (Xi*Xc + Yi*Yc + Zi*Zc)
				* ( (Xi+Xc) * (ri2*rc2) - (Xi*Xc + Yi*Yc + Zi*Zc) * (Xi*(Xc2 + Yc2 + Zc2) + Xc*(Xi2 + Yi2 + Zi2)))) / (ri2*ri2*rc2*rc2) );

			dS2_angular_dXYZ_per_spin.y() += ( (3.0 * (Xi*Xc + Yi*Yc + Zi*Zc)
				* ( (Yi+Yc) * (ri2*rc2) - (Xi*Xc + Yi*Yc + Zi*Zc) * (Yi*(Xc2 + Yc2 + Zc2) + Yc*(Xi2 + Yi2 + Zi2)))) / (ri2*ri2*rc2*rc2) );

			dS2_angular_dXYZ_per_spin.z() += ( (3.0 * (Xi*Xc + Yi*Yc + Zi*Zc)
				* ( (Zi+Zc) * (ri2*rc2) - (Xi*Xc + Yi*Yc + Zi*Zc) * (Zi*(Xc2 + Yc2 + Zc2) + Zc*(Xi2 + Yi2 + Zi2)))) / (ri2*ri2*rc2*rc2) );

			dS2_radial_numerator_dXYZ_per_spin.x() += (-3.0*Xi)/(ri5);
			dS2_radial_numerator_dXYZ_per_spin.y() += (-3.0*Yi)/(ri5);
			dS2_radial_numerator_dXYZ_per_spin.z() += (-3.0*Zi)/(ri5);

			dS2_radial_denominator_dXYZ_per_spin.x() += (-6.0*Xi)/(ri8);
			dS2_radial_denominator_dXYZ_per_spin.y() += (-6.0*Yi)/(ri8);
			dS2_radial_denominator_dXYZ_per_spin.z() += (-6.0*Zi)/(ri8);

			S2_angular_per_spin += (1.5*(cosine_ci * cosine_ci) - 0.5);
			S2_radial_numerator_per_spin += (1.0 / ri3);
			S2_radial_denominator_per_spin += (1.0 / ri6);

			dr6_dXYZ_per_spin.x() += (-6.0*Xi / ri8);
			dr6_dXYZ_per_spin.y() += (-6.0*Yi / ri8);
			dr6_dXYZ_per_spin.z() += (-6.0*Zi / ri8);
		}

		dS2_angular_dXYZ += (dS2_angular_dXYZ_per_spin * sl.first);
		dS2_radial_numerator_dXYZ += (dS2_radial_numerator_dXYZ_per_spin * sl.first);
		dS2_radial_denominator_dXYZ += (dS2_radial_denominator_dXYZ_per_spin * sl.first);

		S2_angular += (S2_angular_per_spin * sl.first);
		S2_radial_numerator += (S2_radial_numerator_per_spin * sl.first);
		S2_radial_denominator += (S2_radial_denominator_per_spin * sl.first);

		dr6_dXYZ += (dr6_dXYZ_per_spin * sl.first);
	}

	dS2_angular_dXYZ /= (total_sum_weights * eq_spins_coords.size());
	Vector dS2_radial_dXYZ = ( (2.0 * dS2_radial_numerator_dXYZ * S2_radial_numerator * S2_radial_denominator
		- S2_radial_numerator * S2_radial_numerator * dS2_radial_denominator_dXYZ)
		/ (total_sum_weights * eq_spins_coords.size() * S2_radial_denominator * S2_radial_denominator) );

	S2_angular /= (eq_spins_coords.size() * total_sum_weights);
	Real S2_radial = (S2_radial_numerator * S2_radial_numerator)/(total_sum_weights * eq_spins_coords.size() * S2_radial_denominator);

	dS2_dXYZ = dS2_angular_dXYZ * S2_radial + S2_angular * dS2_radial_dXYZ;
	dr6_dXYZ /= (total_sum_weights * eq_spins_coords.size());
}

/// @brief xyz derivative of the spectral density function for the simplified Solomon-Bloembergen (SB) equation
/// @params
/// omega_I:     nuclear spin resonance frequency (must be provided in rad/s, dimension is 10^6)
/// tau_c:       correlation time (must be provided in s, typical dimension is 10^-9)
/// dr6_dXYZ:    xyz derivative of <r^-6>
Vector
PREMultiSet::dJ_dXYZ_SB(
	Real const & omega_I,
	Real const & tau_c,
	Vector const & dr6_dXYZ
) const
{
	Real C1(tau_c / (1.0 + (omega_I*omega_I*tau_c*tau_c)));
	Vector dJ_dXYZ = dr6_dXYZ * C1;
	return dJ_dXYZ;
}

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
PREMultiSet::dJ_dXYZ_SBMF(
	Real const & omega_I,
	Real const & tau_c,
	Real const & tau_t,
	Real const & one_over_r6,
	Real const & S2,
	Vector const & dr6_dXYZ,
	Vector const & dS2_dXYZ
) const
{
	Real C1(tau_c / (1.0 + (omega_I*omega_I*tau_c*tau_c)));
	Real C2(tau_t / (1.0 + (omega_I*omega_I*tau_t*tau_t)));
	Vector dJ_dXYZ = dr6_dXYZ * (C1*S2 + C2*(1.0 - S2)) + one_over_r6 * (C1*dS2_dXYZ - C2*dS2_dXYZ);
	return dJ_dXYZ;
}

/// @brief calculate the xyz derivative of the R1 paramagnetic relaxation
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
PREMultiSet::dR1_dXYZ(
	Vec6 const & params,
	Real const & one_over_r6,
	Real const & S2,
	Vector const & dr6_dXYZ,
	Vector const & dS2_dXYZ,
	Real const & scaling,
	bool use_sb
) const
{
	// Dipolar part
	Vector dJ_omega_tau_c_t = use_sb ? dJ_dXYZ_SB(params[3], params[5], dr6_dXYZ) : dJ_dXYZ_SBMF(params[3], params[5], params[6], one_over_r6, S2, dr6_dXYZ, dS2_dXYZ);
	//Vector dpre_dXYZ_dd = -scaling * params[1] * 3.0 * dJ_omega_tau_c_t;
	Vector dpre_dXYZ_dd = scaling * params[1] * 3.0 * dJ_omega_tau_c_t;

	// Curie part
	Vector dJ_omega_tau_r = dJ_dXYZ_SB(params[3], params[4], dr6_dXYZ);
	Vector dJ_omega_tau_c = dJ_dXYZ_SB(params[3], params[5], dr6_dXYZ);
	//Vector dpre_dXYZ_curie = -scaling * params[2] * ( 3.0 * dJ_omega_tau_r - 3.0 * dJ_omega_tau_c);
	Vector dpre_dXYZ_curie = scaling * params[2] * ( 3.0 * dJ_omega_tau_r - 3.0 * dJ_omega_tau_c);

	return dpre_dXYZ_dd + dpre_dXYZ_curie;
}

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
PREMultiSet::dR2_dXYZ(
	Vec6 const & params,
	Real const & one_over_r6,
	Real const & S2,
	Vector const & dr6_dXYZ,
	Vector const & dS2_dXYZ,
	Real const & scaling,
	bool use_sb
) const
{
	// Dipolar part
	Vector dJ_0_tau_c_t = use_sb ? dJ_dXYZ_SB(0.0, params[5], dr6_dXYZ) : dJ_dXYZ_SBMF(0.0, params[5], params[6], one_over_r6, S2, dr6_dXYZ, dS2_dXYZ);
	Vector dJ_omega_tau_c_t = use_sb ? dJ_dXYZ_SB(params[3], params[5], dr6_dXYZ) : dJ_dXYZ_SBMF(params[3], params[5], params[6], one_over_r6, S2, dr6_dXYZ, dS2_dXYZ);
	//Vector dpre_dXYZ_dd = -scaling * params[1] * (4.0 * dJ_0_tau_c_t + 3.0 * dJ_omega_tau_c_t);
	Vector dpre_dXYZ_dd = scaling * params[1] * (4.0 * dJ_0_tau_c_t + 3.0 * dJ_omega_tau_c_t);

	// Curie part
	Vector dJ_0_tau_r = dJ_dXYZ_SB(0.0, params[4], dr6_dXYZ);
	Vector dJ_omega_tau_r = dJ_dXYZ_SB(params[3], params[4], dr6_dXYZ);
	Vector dJ_0_tau_c = dJ_dXYZ_SB(0.0, params[5], dr6_dXYZ);
	Vector dJ_omega_tau_c = dJ_dXYZ_SB(params[3], params[5], dr6_dXYZ);
	//Vector dpre_dXYZ_curie = -scaling * params[2] * (4.0 * dJ_0_tau_r + 3.0 * dJ_omega_tau_r - 4.0 * dJ_0_tau_c - 3.0 * dJ_omega_tau_c);
	Vector dpre_dXYZ_curie = scaling * params[2] * (4.0 * dJ_0_tau_r + 3.0 * dJ_omega_tau_r - 4.0 * dJ_0_tau_c - 3.0 * dJ_omega_tau_c);

	return dpre_dXYZ_dd + dpre_dXYZ_curie;
}

} // namespace pre
} // namespace nmr
} // namespace scoring
} // namespace core
