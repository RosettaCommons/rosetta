// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    protocols/nmr/pre/PREMover.cc
/// @brief   Mover that converts NMR PRE rates into CB-CB atom pair distances and
///          assigns them as atom pair constraints to the pose. In case of degenerate spins
///          (e.g. degenerate protons or symmetric spins) an AmbiguousNMRDistanceConstraint
///          is created.
///          As constraints function a SplineFunc is fitted to a PRE distance histogram.
///          By default, the histogram file is read from the spinlabel database and the generated
///          spline potential converts the HN PRE distances into CB-CB distance constraints.
///          Alternatively, the user can provide a different histogram file to generate CB-CB distance
///          constraints from PRE rates collected for a different type of PRE nucleus e.g. 15N or 1HA.
/// @details last Modified: 12/16/16
///          Note that the calculation of PRE distances is simplified and takes into account only
///          the dipolar and in part Curie relaxation so far. Cross-correlated relaxation is neglected.
///          This approach is reasonable for radicals or paramagnetic ions with a nearly isotropic g-tensor
///          (e.g. nitroxide, Mn2+, Cu2+) which are commonly used in NMR PRE experiments. However, for ions
///          with an anisotropic g-tensor (e.g. lanthanides) Curie and cross-correlated relaxation become
///          more prominent and a direct conversion of PRE rates to distances is not easily possible any longer.
/// @author  Georg Kuenze (georg.kuenze@vanderbilt.edu)

// Unit headers
#include <protocols/nmr/pre/PREMover.hh>
#include <protocols/nmr/pre/PREMoverCreator.hh>

// Package headers
#include <core/scoring/nmr/pre/PREData.hh>
#include <core/scoring/nmr/pre/PREMultiSet.hh>
#include <core/scoring/nmr/pre/PRESingleSet.hh>
#include <core/scoring/nmr/pre/PRESingle.hh>
#include <core/scoring/nmr/NMRSpinlabel.hh>
#include <core/io/nmr/ParaIon.hh>

// Project headers
#include <basic/datacache/DataMap.hh>
#include <core/pose/util.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/conformation/Residue.hh>
#include <core/import_pose/import_pose.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/types.hh>
#include <core/id/AtomID.hh>
#include <core/chemical/AA.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/Constraint.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/AmbiguousNMRDistanceConstraint.hh>
#include <core/scoring/constraints/ConstraintIO.hh>
#include <core/scoring/func/Func.hh>
#include <core/scoring/func/SplineFunc.hh>
#include <protocols/filters/Filter.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <protocols/minimization_packing/MinMover.hh>
#include <protocols/moves/mover_schemas.hh>

// Utility headers
#include <utility/excn/Exceptions.hh>
#include <utility/vector1.hh>
#include <utility/string_util.hh>
#include <utility/tag/Tag.hh>
#include <utility/file/file_sys_util.hh>
#include <utility/tag/XMLSchemaGeneration.hh>

// Numeric headers
#include <numeric/constants.hh>

// Basic headers
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/keys/nmr.OptionKeys.gen.hh>

// ObjexxFCL headers
#include <ObjexxFCL/format.hh>

// C++ headers
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <map>
#include <set>
#include <cmath>
#include <numeric>
#include <algorithm>

namespace protocols {
namespace nmr {
namespace pre {

static basic::Tracer TR( "protocols.nmr.pre.PREMover" );

/// @brief Default constructor
PREMover::PREMover() :
	protocols::moves::Mover( PREMover::mover_name() ),
	histogram_files_(),
	pre_data_(),
	weighted_average_(false),
	minimize_(false)
{
	sfxn_ = core::scoring::get_score_function();
}

/// @brief Construct PREMover from PRE data input file
PREMover::PREMover(
	std::string const & pre_data_file,
	Pose const & pose
) :
	protocols::moves::Mover( PREMover::mover_name() ),
	pre_data_( PREDataOP( new PREData( pre_data_file, pose) ) ),
	weighted_average_(false),
	minimize_(false)
{
	using namespace core::scoring::nmr::pre;

	utility::vector1< PREMultiSetOP > & multiset_vec = pre_data_->get_pre_multiset_vec();
	Real bin_size(0.5);
	for ( auto & p : multiset_vec ) {
		if ( !p->get_spinlabel() ) {
			utility_exit_with_message("ERROR while trying to create PREMover from PREData. No NMRSpinlabel exists in PREData. Set \"spinlabel_type\" option in PRE data input file.");
		} else {
			std::string sl_name = p->get_spinlabel()->get_code();
			std::string sl_histogram_file = p->get_spinlabel()->get_distance_potential_histogram_file();

			// It is possible that the spinlabel has no histogram file stored in the database.
			// In this case, the spinlabel in PREData (PREMultiSet) can still be created,
			// but it cannot provide a PRE distance potential and cannot be used to calculate
			// atom pair constraints. Therefore, we check if the file exists and print a warning
			// if not. Alternatively, we request the user to parse his/her own histogram file.
			// If the histogram file is not set before calling the apply() method, the program is exited.

			if ( !utility::file::file_exists(sl_histogram_file) ) {
				TR.Warning << "Cannot not find histogram file " << sl_histogram_file << " for spinlabel " << sl_name << "." << std::endl;
				TR.Warning << "Histogram is needed for derivation of a spline distance potential and transformation of PRE data into atom pair constraints." << std::endl;
				TR.Warning << "Set histogram file for spinlabel " << sl_name << " through XML script or set_histogram_file() method prior to calling PREMover.apply()." << std::endl;
			}

			if ( histogram_files_.find(sl_name) == histogram_files_.end() ) {
				histogram_files_[ sl_name ] = std::make_pair(sl_histogram_file, bin_size);
			}
		}
	}
}

/// @brief Copy constructor
PREMover::PREMover(PREMover const & other) :
	protocols::moves::Mover( other ),
	histogram_files_(other.histogram_files_),
	pre_data_( new PREData( *(other.pre_data_) ) ),
	sfxn_( other.sfxn_->clone() ),
	weighted_average_(other.weighted_average_),
	minimize_(other.minimize_)
{ }

/// @brief Copy assignment
PREMover &
PREMover::operator=(PREMover const & rhs) {
	if ( this == &rhs ) {
		return *this;
	}
	return *( new PREMover( *this ) );
}

/// @brief destructor
PREMover::~PREMover() {}

std::string
PREMoverCreator::keyname() const {
	return PREMover::mover_name();
}

/// @brief Get the name of this mover
std::string
PREMover::get_name() const {
	return mover_name();
}

std::string
PREMover::mover_name() {
	return "PREMover";
}

protocols::moves::MoverOP
PREMoverCreator::create_mover() const {
	return protocols::moves::MoverOP(new PREMover);
}

/// @brief Make a deep copy of this mover
protocols::moves::MoverOP
PREMover::clone() const {
	return protocols::moves::MoverOP(new PREMover(*this));
}

/// @brief Create a fresh instance of this mover
protocols::moves::MoverOP
PREMover::fresh_instance() const {
	return protocols::moves::MoverOP(new PREMover);
}

void
PREMover::add_histogram_file(
	std::string const & spinlabel_name,
	std::string const & histogram_file,
	Real bin_size
)
{
	if ( !utility::file::file_exists(histogram_file) ) {
		utility_exit_with_message("ERROR: Histogram file " + histogram_file + " for spinlabel " + spinlabel_name + " does not exist.");
	} else {
		histogram_files_[ spinlabel_name ] = std::make_pair(histogram_file, bin_size);
	}
}

/// @brief Calculate distance from R2 relaxation rate
/// @details Considers dipolar and Curie relaxation
/// @params
/// params[1] = gamma_I:     gyromagnetic ratio of the nuclear spin (must be provided in rad/(s*T), dimension is 10^6)
/// params[2] = gJ:          electron Lande factor
/// params[3] = S:           total spin quantum number
/// params[4] = omega_I:     nuclear spin resonance frequency (must be provided in rad/s, dimension is 10^6)
/// params[5] = tau_c:       total correlation time (must be provided in s, typical dimension is 10^-9)
/// params[6] = tau_r:       rotational correlation time (must be provided in s, typical dimension is 10^-9)
/// params[7] = B0:          magnetic field strength (in Tesla)
/// params[8] = T:           temperature (in K)
/// R2:                      R2 relaxation rate (in Hz)
core::Real
PREMover::R2_to_dist_dd_curie(
	Vec8 const & params,
	Real const R2
)
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
	Real gamma_I_2 = params[1] * params[1];
	Real gJ_2 = params[2] * params[2];
	Real gJ_4 = gJ_2 * gJ_2;
	Real B0_2 = params[7] * params[7];
	Real T_2 = params[8] * params[8];

	Real prefac_dd = 1.0/15.0 * mu_0_2/(16.0*pi_2) * gamma_I_2 * gJ_2 * mu_b_2 * params[3]*(params[3]+1);
	Real prefac_curie = 1.0/5.0 * mu_0_2/(16.0*pi_2) * (gamma_I_2 * B0_2 * gJ_4 * mu_b_4 * params[3]*params[3]*(params[3]+1)*(params[3]+1)) / (9.0*kb_2*T_2);

	Real r6_dd = 1.0 / R2 * prefac_dd * (4.0 * params[5] + 3.0 * params[5] / (1.0 + (params[4]*params[5])*(params[4]*params[5]) ) );
	Real r6_curie = 1.0 / R2 * prefac_curie * (4.0 * params[6] + 3.0 * params[6] / (1.0 + (params[4]*params[6])*(params[4]*params[6]))
		- 4.0 * params[5] - 3.0 * params[5] / (1.0 + (params[4]*params[5])*(params[4]*params[5])) );

	runtime_assert(r6_dd > 0.0);
	runtime_assert(r6_curie > 0.0);
	Real r = exp(log(r6_dd + r6_curie)/6.0);
	return r;
}

/// @brief Calculate distance from R1 relaxation rate
/// @details Considers dipolar and Curie relaxation
/// @params
/// params[1] = gamma_I:     gyromagnetic ratio of the nuclear spin (must be provided in rad/(s*T), dimension is 10^6)
/// params[2] = gJ:          electron Lande factor
/// params[3] = S:           total spin quantum number
/// params[4] = omega_I:     nuclear spin resonance frequency (must be provided in rad/s, dimension is 10^6)
/// params[5] = tau_c:       total correlation time (must be provided in s, typical dimension is 10^-9)
/// params[6] = tau_r:       rotational correlation time (must be provided in s, typical dimension is 10^-9)
/// params[7] = B0:          magnetic field strength (in Tesla)
/// params[8] = T:           temperature (in K)
/// R1:                      R1 relaxation rate (in Hz)
core::Real
PREMover::R1_to_dist_dd_curie(
	Vec8 const & params,
	Real const R1
)
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
	Real gamma_I_2 = params[1] * params[1];
	Real gJ_2 = params[2] * params[2];
	Real gJ_4 = gJ_2 * gJ_2;
	Real B0_2 = params[7] * params[7];
	Real T_2 = params[8] * params[8];

	Real prefac_dd = 2.0/15.0 * mu_0_2/(16.0*pi_2) * gamma_I_2 * gJ_2 * mu_b_2 * params[3]*(params[3]+1);
	Real prefac_curie = 2.0/5.0 * mu_0_2/(16.0*pi_2) * (gamma_I_2 * B0_2 * gJ_4 * mu_b_4 * params[3]*params[3]*(params[3]+1)*(params[3]+1)) / (9.0*kb_2*T_2);

	Real r6_dd = 1.0 / R1 * prefac_dd * (3.0 * params[5] / (1.0 + (params[4]*params[5])*(params[4]*params[5]) ) );
	Real r6_curie = 1.0 / R1 * prefac_curie * (3.0 * params[6] / (1.0 + (params[4]*params[6])*(params[4]*params[6]))
		- 3.0 * params[5] / (1.0 + (params[4]*params[5])*(params[4]*params[5])) );

	runtime_assert(r6_dd > 0.0);
	runtime_assert(r6_curie > 0.0);
	Real r = exp(log(r6_dd + r6_curie)/6.0);
	return r;
}

/// @brief Calculate distances from relaxation rates and map them
///        to their respective spinlabel and protein residue(s)
void
PREMover::pre_data_to_distances(
	core::scoring::nmr::pre::PREData & pre_data,
	SpinlabelToPREDistances & all_sl_distances
)
{
	using namespace core::scoring::nmr;
	using namespace core::scoring::nmr::pre;

	utility::vector1< PREMultiSetOP > & multiset_vec = pre_data.get_pre_multiset_vec();
	Size num_sl_sites = pre_data.get_number_spinlabel_sites();
	Vec8 exp_params(8);

	for ( Size i = 1; i <= num_sl_sites; ++i ) {

		// Some basic setup
		utility::vector1< PRESingleSetOP > & singleset_vec = multiset_vec[i]->get_pre_singleset_vec();
		Size sl_rsd_num = multiset_vec[i]->get_spinlabel_site_rsd();
		Size num_experiments = multiset_vec[i]->get_number_experiments();

		// Map of average PRE distances and weights for this spinlabel site
		// Note that we calculate the average value on-the-fly
		PREDistances sl_av_distances;
		PREDistances::iterator sl_av_distances_iter;

		// total sum of all experiment weights; for calculating weighted average PRE distances
		Real sum_expt_wts(0);

		for ( Size j = 1; j <= num_experiments; ++j ) {

			utility::vector1< PRESingle > const & single_pres_vec = singleset_vec[j]->get_pre_single_vec();
			Size num_pres = singleset_vec[j]->get_number_pre();
			sum_expt_wts += singleset_vec[j]->get_weight();

			// Experimental conditions
			exp_params[1] = singleset_vec[j]->calc_gamma_I();
			exp_params[2] = multiset_vec[i]->get_ion_type()->get_gJ();
			exp_params[3] = multiset_vec[i]->get_ion_type()->get_S();
			exp_params[4] = singleset_vec[j]->calc_omega_I();
			exp_params[5] = multiset_vec[i]->get_tau_c();
			exp_params[6] = multiset_vec[i]->get_tau_r();
			exp_params[7] = singleset_vec[j]->get_field_strength() / 42.576; // field strength in Tesla
			exp_params[8] = multiset_vec[i]->get_temperature();

			for ( Size k = 1; k <= num_pres; ++k ) {
				// Residue IDs for this PRE
				std::set< Size > rsds;
				for ( auto const & id : single_pres_vec[k].get_protein_spins() ) {
					rsds.insert(id.rsd());
				}
				Real pre = single_pres_vec[k].get_pre_exp();
				Real err = single_pres_vec[k].get_pre_err();
				if ( err > pre ) {
					std::ostringstream ss;
					for ( auto const & rsd : rsds ) { ss << rsd << " "; }
					TR.Warning << "PRE error = " << err << " for residue(s) " << ss.str() << "is bigger than PRE value = " << pre << ".";
					TR.Warning << "Set PRE error to PRE value / 2.0." << std::endl;
					err = pre / 2.0;
				}
				// Calculate PRE distance
				Real upper_dist(0.0);
				Real lower_dist(0.0);
				if ( singleset_vec[j]->get_pre_rate_type() == R1_PARA ) {
					upper_dist = R1_to_dist_dd_curie(exp_params, pre-err);
					lower_dist = R1_to_dist_dd_curie(exp_params, pre+err);
				} else if ( singleset_vec[j]->get_pre_rate_type() == R2_PARA ) {
					upper_dist = R2_to_dist_dd_curie(exp_params, pre-err);
					lower_dist = R2_to_dist_dd_curie(exp_params, pre+err);
				}
				Real pre_dist = 0.5 * (upper_dist + lower_dist);
				Real pre_dist_tol = 0.5 * (upper_dist - lower_dist);

				if ( weighted_average_ ) {
					pre_dist *= singleset_vec[j]->get_weight();
					pre_dist_tol *= singleset_vec[j]->get_weight();
				}

				if ( TR.Debug.visible() ) {
					std::ostringstream ss;
					for ( auto const & rsd : rsds ) { ss << rsd << " "; }
					TR.Debug << "PRE distance between SL residue " << sl_rsd_num << " and residue(s) " << ss.str() << ": " << pre_dist << " +/- " << pre_dist_tol << std::endl;
				}

				if ( j == 1 ) {
					sl_av_distances.push_back(PREDistanceRecord(rsds, pre_dist, pre_dist_tol));
				} else {
					// Is Residue set with same ID in PREDistance vector?
					for ( sl_av_distances_iter = sl_av_distances.begin();
							sl_av_distances_iter != sl_av_distances.end(); ++sl_av_distances_iter ) {
						if ( rsds == (*sl_av_distances_iter).rsds() ) {
							// Found same residue set
							// Add PRE distance and calculate average on-the-fly
							(*sl_av_distances_iter).add_and_average(pre_dist, pre_dist_tol);
							break;
						}
					}
					// Did not find same residue set
					// Create a new PREDistanceRecord object and
					// add it to the PREDistances vector
					if ( sl_av_distances_iter == sl_av_distances.end() ) {
						sl_av_distances.push_back(PREDistanceRecord(rsds, pre_dist, pre_dist_tol));
					}
				}
			} // number of PREs

		} // number of experiments

		if ( weighted_average_ ) {
			for ( sl_av_distances_iter = sl_av_distances.begin();
					sl_av_distances_iter != sl_av_distances.end(); ++sl_av_distances_iter ) {
				(*sl_av_distances_iter).set_dist( (*sl_av_distances_iter).get_dist() / (sum_expt_wts / (*sl_av_distances_iter).count() ) );
				(*sl_av_distances_iter).set_tol( (*sl_av_distances_iter).get_tol() / (sum_expt_wts / (*sl_av_distances_iter).count() ) );
			}
		}

		// Finally set PREDistances vector entry in SpinlabelToPREDistances map
		all_sl_distances[ sl_rsd_num ] = sl_av_distances;

	} // number of spinlabel sites
}

/// @brief Calculate CB-CB distances from PRE rates and append
///        them as atom pair distance constraints to the pose
void
PREMover::apply(Pose & pose) {
	using namespace core::id;
	using namespace core::scoring::nmr::pre;
	using namespace core::scoring::func;
	using namespace core::scoring::constraints;
	using namespace core::scoring;

	// PRE distances for all spinlabel sites
	SpinlabelToPREDistances sl_all_distances;
	ConstraintSetOP csts( new ConstraintSet);

	if ( !pre_data_ ) {
		utility_exit_with_message("PREData member of PREMover is not set. Must be set prior to \"apply()\" function.");
	}
	pre_data_to_distances( *pre_data_, sl_all_distances );

	utility::vector1< PREMultiSetOP > & multiset_vec = pre_data_->get_pre_multiset_vec();
	Size num_sl_sites = pre_data_->get_number_spinlabel_sites();

	for ( Size i = 1; i <= num_sl_sites; ++i ) {
		if ( !multiset_vec[i]->get_spinlabel() ) {
			utility_exit_with_message("ERROR while calculating PRE distances constraints with PREMover. No NMRSpinlabel exists in PREData. Set \"spinlabel_type\" option in PRE data input file.");
		}

		std::string sl_name = multiset_vec[i]->get_spinlabel()->get_code();
		if ( histogram_files_.count(sl_name) == 0 ) {
			utility_exit_with_message("ERROR while calculating PRE distances constraints with PREMover. No histogram file for spinlabel residue " + sl_name + " found.");
		}

		// Get AtomID for this spinlabel site
		Size sl_site = multiset_vec[i]->get_spinlabel_site_rsd();
		AtomID sl_site_atomid;
		if ( !core::chemical::is_canonical_L_aa_or_gly( pose.residue(sl_site).aa() ) ) {
			utility_exit_with_message("ERROR while calculating PRE distance constraints with PREMover. Residue "
				+ utility::to_string(sl_site) + " is not a canonical L-amino acid.");
		}
		if ( pose.residue(sl_site).has("CB") ) {
			sl_site_atomid = AtomID(pose.residue(sl_site).type().atom_index("CB"), sl_site);
		} else if ( pose.residue(sl_site).has("2HA") ) {
			sl_site_atomid = AtomID(pose.residue(sl_site).type().atom_index("2HA"), sl_site);
		} else {
			utility_exit_with_message("ERROR while calculating PRE distance constraints with PREMover. Residue "
				+ utility::to_string(sl_site) + " has no CB or 2HA atom.");
		}

		// PRE distances for this spinlabel site
		PREDistances const & sl_one_distances = sl_all_distances[sl_site];

		// Determine the maximal error and give distances with a lower error a higher weight
		Real max_tol(0);
		for ( auto const & d : sl_one_distances ) {
			if ( d.get_tol() > max_tol ) { max_tol = d.get_tol(); }
		}

		for ( Size j = 1, j_end = sl_one_distances.size(); j <= j_end; ++j ) {
			utility::vector1< Size > nuc_spins_rsds( sl_one_distances[j].rsds().begin(), sl_one_distances[j].rsds().end() );
			utility::vector1< AtomID > nuc_spins_atomids;
			for ( auto const & rsd : nuc_spins_rsds ) {
				if ( !core::chemical::is_canonical_L_aa_or_gly( pose.residue(rsd).aa() ) ) {
					utility_exit_with_message("ERROR while calculating PRE distance constraints with PREMover. Residue " + utility::to_string(rsd) + " is not a canonical L-amino acid.");
				}
				if ( pose.residue(rsd).has("CB") ) {
					nuc_spins_atomids.push_back( AtomID(pose.residue(rsd).type().atom_index("CB"), rsd) );
				} else if ( pose.residue(rsd).has("2HA") ) {
					nuc_spins_atomids.push_back( AtomID(pose.residue(rsd).type().atom_index("2HA"), rsd) );
				} else {
					utility_exit_with_message("ERROR while calculating PRE distance constraints with PREMover. Residue " + utility::to_string(rsd) + " has no CB or 2HA atom.");
				}
			}

			Real wt = (max_tol * max_tol) / (sl_one_distances[j].get_tol() * sl_one_distances[j].get_tol());

			// Set up csts function from string stream object
			// KB_potential_description << histogram_file << exp_val << weight << bin_size
			SplineFuncOP spline( new SplineFunc);
			std::stringstream ss;
			ss << "difference " << histogram_files_[sl_name].first << " " <<  sl_one_distances[j].get_dist() << " "
				<< wt << " " << histogram_files_[sl_name].second;
			TR.Debug << ss.str() << std::endl;
			spline->read_data(ss);

			// Single nuclear spin -> AtomPair constraint
			if ( nuc_spins_rsds.size() == 1 ) {
				AtomPairConstraintOP pre_cst( new AtomPairConstraint(sl_site_atomid, nuc_spins_atomids[1], spline) );
				csts->add_constraint(pre_cst);
			} else {
				// Multiple nuclear spins -> AmbiguousNMRDistance constraint
				AmbiguousNMRDistanceConstraintOP pre_cst( new AmbiguousNMRDistanceConstraint(utility::vector1<AtomID>(1, sl_site_atomid), nuc_spins_atomids, spline) );
				csts->add_constraint(pre_cst);
			}

		} // PRE csts for this spinlabel

	} // all spinlabel

	pose.constraint_set(csts);

	// Show PRE constraints and write them to file
	if ( TR.Trace.visible() ) {
		TR.Trace << "Write PRE distance constraints to file \"pre_distances.csts\"." << std::endl;
		csts->show_definition(TR.Trace, pose);
		std::string pre_csts_file("pre_distances.csts");
		ConstraintIO::get_instance()->write_constraints(pre_csts_file, *csts, pose);
	}

	if ( minimize_ ) {
		ScoreFunctionOP sfxn_copy = sfxn_->clone();
		if ( sfxn_copy->has_zero_weight(atom_pair_constraint) ) {
			TR.Warning << "PREMover scorefunction not set up for minimization with PRE constraints. Atom pair constraint weight is zero." << std::endl;
			TR.Warning << "Making copy of scorefunction and setting atom_pair_constraint weight to 1.0." << std::endl;
			sfxn_copy->set_weight(atom_pair_constraint,1.0);
		}
		// Fast minimization after setting the constraints
		core::kinematics::MoveMapOP mm( new core::kinematics::MoveMap() );
		for ( Size i = 1; i <= pose.total_residue(); ++i ) {
			mm->set_bb(i, true);
			mm->set_chi(i, true);
		}
		mm->set_jump(1, true);
		protocols::minimization_packing::MinMoverOP min_mover(new protocols::minimization_packing::MinMover(mm, sfxn_copy,"lbfgs_armijo_nonmonotone", 0.01, true));
		min_mover->apply(pose);
	}
}

void
PREMover::show(std::ostream& tracer) const {
	tracer << " * * * PREMover Data * * * " << std::endl;
	if ( pre_data_ ) {
		pre_data_->show(tracer);
	}
	tracer << "Spinlabel histograms:" << std::endl;
	for ( SpinlabelHistogramMap::const_iterator iter = histogram_files_.begin(), last = histogram_files_.end(); iter != last; ++iter ) {
		tracer << (*iter).first << ": " << (*iter).second.first << std::endl;
	}
}

/// @brief Parse tags of XML script
void
PREMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & datamap,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	core::pose::Pose const & pose
)
{
	try {
		// PREData
		if ( tag->hasOption("pre_input_file") ) {
			std::string filename;
			filename = tag->getOption< std::string >( "pre_input_file", "" );
			if ( filename == "" ) {
				utility_exit_with_message("ERROR: No PRE data input filename provided for PREMover.");
			} else {
				pre_data_ = PREDataOP(new PREData(filename, pose));
			}
		}
		// ScoreFunction for optional minimization
		if ( tag->hasOption("scorefxn") ) {
			sfxn_ = protocols::rosetta_scripts::parse_score_function(tag, "scorefxn", datamap);
		}
		if ( tag->hasOption("weighted_average") ) {
			weighted_average_ = tag->getOption< bool >("weighted_average", false);
		}
		if ( tag->hasOption("minimize_w_csts") ) {
			minimize_ = tag->getOption< bool >("minimize_w_csts", false);
		}
		//Spinlabel histogram files
		utility::vector1< utility::tag::TagCOP > const branch_tags( tag->getTags() );
		utility::vector1< utility::tag::TagCOP >::const_iterator tag_it;
		for ( tag_it = branch_tags.begin(); tag_it != branch_tags.end(); ++tag_it ) {
			if ( (*tag_it)->getName() == "Histograms" ) {
				std::string sl_name = (*tag_it)->getOption<std::string>( "spinlabel" );
				std::string histogram_file = (*tag_it)->getOption<std::string>( "histogram_file" );
				Real bin_size = (*tag_it)->getOption<Real>( "bin_size" );
				add_histogram_file(sl_name, histogram_file, bin_size);
			}
		}
	} catch ( utility::excn::RosettaScriptsOptionError & excn ) {
		TR << "caught exception " << excn.msg() << std::endl;
	}
}

/// @brief Create XML schema definition for PREMover
void
PREMover::provide_xml_schema(utility::tag::XMLSchemaDefinition & xsd) {
	using namespace utility::tag;

	// Basic attributes
	AttributeList attlist;
	attlist
		+ XMLSchemaAttribute::required_attribute( "pre_input_file", xs_string, "PRE data input file" )
		+ XMLSchemaAttribute::attribute_w_default( "minimize_w_csts", xsct_rosetta_bool, "Do one round of minimization with PRE distance constraints", "false" )
		+ XMLSchemaAttribute( "scorefxn", xs_string, "Use this scorefunction when minimizing with PRE distance constraints" )
		+ XMLSchemaAttribute::attribute_w_default( "weighted_average", xsct_rosetta_bool, "Calculate average distance from multiple PRE experiments for the same spinlabel site using single experiment weights", "false" );

	// attributes for Histograms subelement
	AttributeList histogams_subelement_attributes;
	histogams_subelement_attributes
		+ XMLSchemaAttribute::required_attribute( "spinlabel", xs_string, "Three-letter spinlabel code")
		+ XMLSchemaAttribute::required_attribute( "histogram_file", xs_string, "Path to spinlabel histogram file used for transformation of PREs into CB-CB atom pair constraints" )
		+ XMLSchemaAttribute::attribute_w_default( "bin_size", xsct_real, "Histogram bin size", "0.5" );

	XMLSchemaSimpleSubelementList subelements;
	subelements.add_simple_subelement( "Histograms", histogams_subelement_attributes, "Provide spinlabel histograms for transformation of PREs into CB-CB atom pair constraints if not available in Rosetta database" );
	protocols::moves::xsd_type_definition_w_attributes_and_repeatable_subelements( xsd, mover_name(),
		"This Mover converts PRE rates into CB-CB atom pair constraints and attaches them to the pose", attlist, subelements );
}

void
PREMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const {
	PREMover::provide_xml_schema( xsd );
}

} // namespace pre
} // namespace nmr
} // namespace protocols
