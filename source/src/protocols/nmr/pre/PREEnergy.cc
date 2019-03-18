// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    protocols/nmr/pre/PREEnergy.cc
/// @brief   Implementation of class PREEnergy
/// @details last Modified: 10/23/16
/// @author  Georg Kuenze (georg.kuenze@vanderbilt.edu)

// Unit headers
#include <protocols/nmr/pre/PREEnergy.hh>
#include <protocols/nmr/pre/PREEnergyCreator.hh>

// Project headers
#include <core/scoring/nmr/NMRDataFactory.hh>
#include <core/scoring/nmr/pre/PRESingle.hh>
#include <core/scoring/nmr/pre/PRESingleSet.hh>
#include <core/scoring/nmr/pre/PREMultiSet.hh>
#include <core/scoring/nmr/pre/PREData.hh>
#include <core/scoring/nmr/NMRSpinlabel.hh>
#include <protocols/nmr/nmrspinlabel_util.hh>
#include <core/scoring/nmr/NMRDummySpinlabelEnsemble.hh>
#include <core/scoring/nmr/NMRGridSearch.hh>
#include <core/scoring/nmr/util.hh>

// Package headers
#include <core/pose/Pose.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/pose/symmetry/util.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/ResidueType.hh>
#include <core/scoring/methods/WholeStructureEnergy.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <basic/datacache/BasicDataCache.hh>
#include <core/id/AtomID_Map.hh>
#include <core/id/AtomID.hh>

// Basic headers
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/nmr.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>

// Utility headers
#include <utility/exit.hh>
#include <utility/vector1.hh>
#include <utility/graph/Graph.hh>
#include <utility/string_util.hh>

// Numeric headers
#include <numeric/xyzVector.hh>
#include <numeric/xyzMatrix.hh>
#include <numeric/HomogeneousTransform.hh>
#include <numeric/constants.hh>

// boost headers
#include <boost/unordered/unordered_map.hpp>
#include <boost/functional/hash.hpp>

// C++ headers
#include <string>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <map>
#include <cmath>
#include <limits>
#include <algorithm>

namespace protocols {
namespace nmr {
namespace pre {

static basic::Tracer TR( "protocols.nmr.pre.PREEnergy" );

/// @brief Instantiate a new PREEnergy
core::scoring::methods::EnergyMethodOP
PREEnergyCreator::create_energy_method( core::scoring::methods::EnergyMethodOptions const & ) const {
	return core::scoring::methods::EnergyMethodOP( new PREEnergy );
}

// @brief Return the set of score types claimed by the EnergyMethod
// this EnergyMethodCreator creates in its create_energy_method() function
core::scoring::ScoreTypes
PREEnergyCreator::score_types_for_method() const {
	core::scoring::ScoreTypes sts;
	sts.push_back( core::scoring::nmr_pre );
	return sts;
}

/// @brief default constructor
PREEnergy::PREEnergy() :
	core::scoring::methods::WholeStructureEnergy( core::scoring::methods::EnergyMethodCreatorOP( new PREEnergyCreator ) )
{ }

// @brief copy constructor
PREEnergy::PREEnergy(PREEnergy const & other) :
	core::scoring::methods::WholeStructureEnergy( other ),
	atom_id_to_pre_xyz_deriv_map_(other.atom_id_to_pre_xyz_deriv_map_)
{ }

/// @brief destructor
PREEnergy::~PREEnergy() { }

core::scoring::methods::EnergyMethodOP
PREEnergy::clone() const {
	return core::scoring::methods::EnergyMethodOP( new PREEnergy );
}

/// @brief Return PREData from pose. Create PREData if not present and attach them to the pose.
core::scoring::nmr::pre::PREData &
PREEnergy::get_pre_data_from_pose(Pose & pose) const {
	using namespace core::pose::datacache;
	using namespace core::scoring::nmr;
	using namespace core::scoring::nmr::pre;

	if ( pose.data().has( CacheableDataType::NMR_PRE_DATA ) ) {
		return *( utility::pointer::static_pointer_cast< PREData > ( pose.data().get_ptr( CacheableDataType::NMR_PRE_DATA ) ) );
	} else {
		std::string inputfile;
		if ( basic::options::option[ basic::options::OptionKeys::nmr::pre::input_file ].active() ) {
			inputfile = basic::options::option[ basic::options::OptionKeys::nmr::pre::input_file ]();
		} else {
			utility_exit_with_message( "No PRE input file given. You must provide the input file on the command line \"-nmr::pre::input_file filename\"." );
		}
		PREDataOP pre_data_all_ptr(utility::pointer::static_pointer_cast< PREData >( NMRDataFactory::get_instance()->get_nmr_data("PRE", inputfile, pose) ) );
		pose.data().set( CacheableDataType::NMR_PRE_DATA, pre_data_all_ptr );
		return *pre_data_all_ptr;
	}
}

/// @brief Calculate the PRE score and write it into the pose's total energy EnergyMap
void
PREEnergy::finalize_total_energy(
	Pose & pose,
	ScoreFunction const & /*sxfn*/,
	EnergyMap & totals
) const
{
	totals[ core::scoring::nmr_pre ] = calculate_total_score( pose );
}

/// @brief Calculate the total PRE score from PREData retrieved from the pose
core::Real
PREEnergy::calculate_total_score(Pose & pose) const {
	using namespace core::scoring::nmr;
	using namespace core::scoring::nmr::pre;
	using WeightCoordVector = core::scoring::nmr::NMRSpinlabel::WeightCoordVector;

	PREData & pre_data_all = get_pre_data_from_pose(pose);

	// Some basic setup
	Real total_score(0);
	Size number_spinlabel_sites = pre_data_all.get_number_spinlabel_sites();
	utility::vector1< PREMultiSetOP > multiset_vec = pre_data_all.get_pre_multiset_vec();
	runtime_assert(number_spinlabel_sites == multiset_vec.size());

	// Loop over the spinlabel sites (i.e. the PREMultiSets)
	for ( Size i = 1; i<= number_spinlabel_sites; ++i ) {

		// Prepare PRE data
		Real score_multiset(0);
		multiset_vec[i]->update_spin_coordinates(pose);
		Size spinlabel_position = multiset_vec[i]->get_spinlabel_site_rsd();
		// Vector of weights and radical atom xyz coordinates for N distinct spinlabel conformers
		WeightCoordVector spinlabel_wghts_coords;

		if ( multiset_vec[i]->get_spinlabel() ) {
			NMRSpinlabelOP spinlabel = multiset_vec[i]->get_spinlabel();
			// Here we split behavior depending on the representation of the pose and the computation type
			if ( pose.is_fullatom() && spinlabel->get_highres_conformer_filter_type() == NMRSpinlabel::BUMP_ENERGY ) {

				spinlabel_wghts_coords = filter_spinlabel_ensemble_by_packerenergy(pose, *spinlabel, spinlabel_position);

				// PRE calculation with implicit dummy spinlabel ensemble
			} else {
				// We need to create and filter the dummy spinlabel ensemble only for the ASU.
				// The PREMultiSet is taking care of to also calculate and sum up the PRE for
				// the other subunits in case of a symmetric pose.

				// Filter dummy ensemble by clash score calculation and get the spinlabel coordinates
				// We can simply call the member method of the NMRSpinlabel. Internally, it performs the
				// clash filter, looks up the coordinates of the radical atom, transforms their coordinates
				// into the coordinate frame of the spinlabel site and performs additional binning if the
				// vector size exceeds the maximal number of spinlabel conformers
				runtime_assert_msg( spinlabel->get_dummy_ensemble(),
					"ERROR while trying to calculate PRE score. Dummy ensemble of NMRSpinlabel not set. Check if database file for this NMRDummySpinlabel exists in Rosetta database." );
				spinlabel_wghts_coords = spinlabel->filter_spinlabel_ensemble_by_distance_check(pose, spinlabel_position);
			}
			score_multiset = multiset_vec[i]->compute_pre_score_from_point_vector(spinlabel_wghts_coords);

		} else if ( multiset_vec[i]->get_gridsearch_iterator() ) {

			Vector current_metal_coords;
			Real current_score(0);
			Real best_score(std::numeric_limits<Real>::max());
			multiset_vec[i]->get_gridsearch_iterator()->set_grid_search_center( pose );
			Vector best_metal_coords(multiset_vec[i]->get_gridsearch_iterator()->get_grid_search_center());

			// Perform NLS fitting of para ion position and correlation times
			if ( multiset_vec[i]->optimize_paraion_position() ) {
				score_multiset = multiset_vec[i]->compute_pre_score_from_single_point(best_metal_coords, true);
			} else {
				// Run the gridsearch of metal ion position and fit only the correlation times
				while ( multiset_vec[i]->get_gridsearch_iterator()->valid_next_grid_point(current_metal_coords) ) {
					current_score = multiset_vec[i]->compute_pre_score_from_single_point(current_metal_coords, false);
					if ( current_score < best_score ) {
						best_score = current_score;
						best_metal_coords = current_metal_coords;
						multiset_vec[i]->get_gridsearch_iterator()->set_best_grid_point(best_metal_coords);
					}
				}
				score_multiset = multiset_vec[i]->compute_pre_score_from_single_point(best_metal_coords, false);
			}
			spinlabel_wghts_coords = WeightCoordVector(1, std::pair<Real,Vector>(1, best_metal_coords));

		} else {
			utility_exit_with_message("ERROR during PRE score calculation. No NMRSpinlabel or gridsearch iterator are set for PRE dataset at spinlabel site " + utility::to_string(spinlabel_position));
		}

		// Finally set the atom derivatives
		multiset_vec[i]->set_atom_derivatives(pose, spinlabel_wghts_coords);

		TR.Trace << "PRE Score for spinlabel position " << spinlabel_position << " " << score_multiset << " (unweighted) or "
			<< score_multiset * multiset_vec[i]->get_weight() << " (weighted)." << std::endl;
		total_score += score_multiset * multiset_vec[i]->get_weight();
	}
	TR.Trace << "Total PRE score " << total_score << std::endl;
	return total_score;
}

/// @brief Called at the beginning of atom tree minimization, this method
///        allows the derived class the opportunity to initialize pertinent data
///        that will be used during minimization.
///        Here, the function creates and updates the atom_id_to_pre_xyz_deriv_map_ which
///        is needed by the eval_atom_derivative() function.
void
PREEnergy::setup_for_minimizing(
	Pose & pose,
	ScoreFunction const & /*sxfn*/,
	core::kinematics::MinimizerMapBase const & /*minmap*/
) const
{
	using namespace core::scoring::nmr::pre;

	PREData & pre_data_all = get_pre_data_from_pose(pose);

	// Set to zero because we are using operator+=
	Vector fij(0 ,0 ,0);
	atom_id_to_pre_xyz_deriv_map_.default_value(fij);
	atom_id_to_pre_xyz_deriv_map_.fill_with(fij);

	utility::vector1< PREMultiSetOP > & multiset_vec = pre_data_all.get_pre_multiset_vec();
	Size number_spinlabel_sites(pre_data_all.get_number_spinlabel_sites());
	utility::vector1<PRESingle>::const_iterator iter;

	for ( Size i = 1; i <= number_spinlabel_sites; ++i ) {
		utility::vector1< PRESingleSetOP > const & singleset_vec = multiset_vec[i]->get_pre_singleset_vec();
		Size number_experiments(multiset_vec[i]->get_number_experiments());
		core::conformation::symmetry::SymmetryInfoCOP syminfo_ptr;
		if ( multiset_vec[i]->symmetric_pre_calc() && core::pose::symmetry::is_symmetric(pose) ) {
			syminfo_ptr = core::pose::symmetry::symmetry_info( pose );
		}

		for ( Size j = 1; j <= number_experiments; ++j ) {
			utility::vector1<PRESingle> const & single_pre_vec = singleset_vec[j]->get_pre_single_vec();

			for ( iter = single_pre_vec.begin(); iter != single_pre_vec.end(); ++iter ) {
				utility::vector1< core::id::AtomID > const & protein_spins = iter->get_protein_spins();
				utility::vector1< Vector > const & atom_derivatives = iter->get_atom_derivatives();
				runtime_assert_msg(protein_spins.size() == atom_derivatives.size(), "ERROR in setup of atom tree minimization from PRE derivatives. Vector of AtomIDs and of derivatives have unequal length.");
				for ( Size k = 1, k_end = protein_spins.size(); k <= k_end; ++k ) {
					fij = atom_id_to_pre_xyz_deriv_map_.get( protein_spins[k] );  // returns (0,0,0) if atom is not present in map, otherwise the corresponding derivative
					fij += atom_derivatives[k];                         // accumulates the derivative for that particular atom (for all pre data)
					atom_id_to_pre_xyz_deriv_map_.set( protein_spins[k], fij);   // sets the accumulated derivative and resizes the map if necessary

				}
				// Do another pass to set PRE derivatives for symmetric residues
				if ( syminfo_ptr ) {
					for ( Size k = 1, k_end = protein_spins.size(); k <= k_end; ++k ) {
						utility::vector1< Size >protein_spins_symm_rsd = syminfo_ptr->bb_clones(protein_spins[k].rsd());

						for ( Size l = 1, l_end = protein_spins_symm_rsd.size(); l <= l_end; ++l ) {
							core::id::AtomID symm_spin(protein_spins[k].atomno(), protein_spins_symm_rsd[l]);

							// return the derivative of that atom, if not present in map, return the default value which is (0,0,0)
							Vector fij_symm_spin = atom_id_to_pre_xyz_deriv_map_.get(symm_spin);

							// rotate derivative vector into the frame of the symmetric subunit
							// and accumulate the contribution from all pre experiments
							fij_symm_spin += core::scoring::nmr::apply_vector_rotation(pose, atom_derivatives[k], protein_spins[k].rsd(), protein_spins_symm_rsd[l]);

							// set the accumulated derivative and resize the map if necessary
							atom_id_to_pre_xyz_deriv_map_.set(symm_spin, fij_symm_spin);
						}
					}
				}
			}
		}
	}
}


/// @brief Evaluate the xyz derivative of the PRE for an atom in the pose.
///        Called during the atomtree derivative calculation, atom_tree_minimize.cc,
///        through the ScoreFunction::eval_atom_derivative intermediary.
///        F1 and F2 should not be zeroed, rather, setup_for_minimizing() accumulates the PRE
///        contribution from the xyz derivatives of atom id
void
PREEnergy::eval_atom_derivative(
	core::id::AtomID const & id,
	Pose const & pose,
	core::kinematics::DomainMap const & /*domain_map*/,
	ScoreFunction const & /*sfxn*/,
	EnergyMap const & weights,
	Vector & F1,
	Vector & F2
) const
{
	if ( !atom_id_to_pre_xyz_deriv_map_.has(id) ) {
		// only checks if id residue number is in range of AtomID_Map
		return;
	}
	Vector fij = atom_id_to_pre_xyz_deriv_map_.get(id); //  returns (0,0,0) if id is not present in map, otherwise the accumulated derivative of the PRE for this atom

	Vector atom_x = pose.xyz(id);
	Vector const f2( -fij );
	Vector const atom_y = atom_x - f2;   // a "fake" atom in the direction of the gradient
	Vector const f1( atom_x.cross( atom_y ) );

	F1 += weights[ core::scoring::nmr_pre ] * f1;
	F2 += weights[ core::scoring::nmr_pre ] * f2;
}

/// @brief register options
void
PREEnergy::register_options() {
	using namespace basic::options;
	option.add_relevant(OptionKeys::nmr::pre::show_info);
}

/// @brief Indicate in the context-graphs-required list which
///        context-graphs this energy method requires that the Pose
///        maintain when doing neighbor evaluation. Context graphs are allowed.
void
PREEnergy::indicate_required_context_graphs(utility::vector1< bool > &) const { }

/// @brief Return the version of the energy method
core::Size
PREEnergy::version() const { return 1; }

/// @brief show additional information of the energy method
void
PREEnergy::show_additional_info(
	std::ostream & tracer,
	Pose & pose,
	bool verbose
) const
{
	using namespace core::scoring::nmr;
	using namespace core::scoring::nmr::pre;

	if ( basic::options::option[ basic::options::OptionKeys::nmr::pre::show_info ].user() ) {
		verbose = basic::options::option[ basic::options::OptionKeys::nmr::pre::show_info ]();
	}

	// Some basic setup
	PREData & pre_data_all = get_pre_data_from_pose(pose);
	Real pre_score = calculate_total_score(pose);
	Size number_spinlabel_sites = pre_data_all.get_number_spinlabel_sites();
	utility::vector1< Real > scores_all_spinlabels(number_spinlabel_sites, 0.0);
	utility::vector1< PREMultiSetOP > & multiset_vec = pre_data_all.get_pre_multiset_vec();
	runtime_assert(number_spinlabel_sites == multiset_vec.size());

	tracer << " * * * * * * *       PREEnergy Info       * * * * * * * " << std::endl;
	tracer << "Total Score: " << pre_score << std::endl;
	tracer << "Is norm of " << number_spinlabel_sites << " spinlabel sites." << std::endl;

	Real sum_all_dev_square(0);
	Real sum_all_exp_square(0);
	Size total_number_pre(0);
	Size total_number_experiments(0);

	// Loop over PREMultiSets
	for ( Size i = 1; i <= number_spinlabel_sites; ++i ) {
		Real sum_dev_square(0);
		Real sum_exp_square(0);

		utility::vector1< PRESingleSetOP > & singleset_vec = multiset_vec[i]->get_pre_singleset_vec();
		Size number_experiments(multiset_vec[i]->get_number_experiments());
		runtime_assert(number_experiments == singleset_vec.size());
		Size number_pres_per_spinlabel(multiset_vec[i]->get_total_number_pre());
		utility::vector1<Real> const & pre_values = multiset_vec[i]->get_pre_values();
		utility::vector1<Real> const & pre_single_weights = multiset_vec[i]->get_pre_single_weights();
		runtime_assert_msg(number_pres_per_spinlabel == pre_values.size(), "ERROR in PREEnergy's show_additional_info() function. Number of PREs and length of PRE value vector are not the same");
		runtime_assert_msg(number_pres_per_spinlabel == pre_single_weights.size(), "ERROR in PREEnergy's show_additional_info() function. Number of PREs and length of PRE weights vector are not the same");

		Size index_offset(0);

		// Loop over PRESingleSets
		for ( Size j = 1; j <= number_experiments; ++j ) {

			Real single_score(0);
			Size number_pres_per_experiment(singleset_vec[j]->get_number_pre());
			utility::vector1< PRESingle > const & pre_single_vec = singleset_vec[j]->get_pre_single_vec();
			runtime_assert(number_pres_per_experiment == pre_single_vec.size());
			Real scal(1.0 / singleset_vec[j]->get_scaling_factor());

			tracer << "PRE dataset: " << singleset_vec[j]->get_dataset_name() << std::endl;
			if ( verbose ) {
				tracer << " * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * " << std::endl;
				tracer << " * * * * * * * * * * * * * * * *              exp vs. calc PRE             * * * * * * * * * * * * * * * * * " << std::endl;
				tracer << " * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * " << std::endl;
				tracer << " Spin_Resid_Atom     Obs_PRE     Scal_Obs_PRE     Calc_PRE     Scal_Calc_PRE     Deviation     Abs_Deviation " << std::endl;
				tracer << " ----------------------------------------------------------------------------------------------------------- " << std::endl;
			}

			// Loop over PRE values per experiment
			for ( Size k = 1; k <= number_pres_per_experiment; ++k ) {
				Real calc_pre(pre_single_vec[k].get_pre_calc());
				Real pre_dev(calc_pre - pre_values[ index_offset + k ]);
				Real pre_abs_dev(std::abs(pre_dev));
				sum_dev_square += pre_dev * pre_dev;
				sum_exp_square += pre_values[ index_offset + k ] * pre_values[ index_offset + k ];
				single_score += pre_dev * pre_dev * pre_single_weights[ index_offset + k ];

				if ( verbose ) {
					std::ostringstream converter;
					for ( core::id::AtomID const & id : pre_single_vec[k].get_protein_spins() ) {
						converter << std::setw(4) << id.rsd() << " " << pose.residue(id.rsd()).atom_name(id.atomno()) << " ";
					}
					tracer << " ( " << std::right << converter.str() << ")";
					tracer << std::setw(12)<< std::fixed << std::setprecision(3) << pre_values[ index_offset + k ] / scal; // Revert normalization of PRE data that was applied during creation of PRESingleSet
					tracer << std::setw(17)<< std::fixed << std::setprecision(3) << pre_values[ index_offset + k ];
					tracer << std::setw(13)<< std::fixed << std::setprecision(3) << calc_pre / scal;       // Revert normalization of PRE data that was applied during creation of PRESingleSet
					tracer << std::setw(18)<< std::fixed << std::setprecision(3) << calc_pre;
					tracer << std::setw(14)<< std::fixed << std::setprecision(3) << pre_dev;
					tracer << std::setw(18)<< std::fixed << std::setprecision(3) << pre_abs_dev << std::endl;
				}
			} // PREs per experiment
			if ( verbose ) {
				tracer << " ----------------------------------------------------------------------------------------------------------- " << std::endl;
			}

			//scores_all_spinlabels[i] += std::sqrt(single_score) * singleset_vec[j]->get_weight();
			scores_all_spinlabels[i] += single_score * singleset_vec[j]->get_weight();
			tracer << "Weighted score for PRE dataset " << singleset_vec[j]->get_dataset_name() << " " << single_score << std::endl;

			// Increment index offset
			index_offset += number_pres_per_experiment;

		} // PREs per spinlabel position

		tracer << "Overall weighted PRE score for spinlabel position " << multiset_vec[i]->get_spinlabel_site_rsd() << ": " << scores_all_spinlabels[i] << std::endl;
		Real Q_fac_per_sl_site(std::sqrt(sum_dev_square/sum_exp_square));
		Real Rmsd_per_sl_site(std::sqrt(sum_dev_square/number_pres_per_spinlabel));
		tracer << "PRE Q-factor and Rmsd for spinlabel position " << multiset_vec[i]->get_spinlabel_site_rsd() << ": " << Q_fac_per_sl_site << ", " << Rmsd_per_sl_site << std::endl;
		sum_all_dev_square += sum_dev_square;
		sum_all_exp_square += sum_exp_square;
		total_number_pre += number_pres_per_spinlabel;
		total_number_experiments += number_experiments;

	} // all PRE data

	Real total_Q_fac(std::sqrt(sum_all_dev_square/sum_all_exp_square));
	Real total_Rmsd(std::sqrt(sum_all_dev_square/total_number_pre));
	tracer << "Total PRE Q-factor and Rmsd for " << number_spinlabel_sites << " spinlabel positions and " << total_number_experiments
		<< " experiments: " << total_Q_fac << ", " << total_Rmsd << std::endl;
}

} // namespace pcs
} // namespace nmr
} // namespace protocols
