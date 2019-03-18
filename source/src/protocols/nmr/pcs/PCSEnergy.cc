// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    protocols/nmr/pcs/PCSEnergy.cc
/// @brief   Implementation of class PCSEnergy
/// @details last Modified: 07/19/16
/// @author  Georg Kuenze (georg.kuenze@vanderbilt.edu)

// Unit headers
#include <protocols/nmr/pcs/PCSEnergy.hh>
#include <protocols/nmr/pcs/PCSEnergyCreator.hh>

// Project headers
#include <core/scoring/nmr/NMRDataFactory.hh>
#include <core/scoring/nmr/NMRGridSearch.hh>
#include <core/scoring/nmr/NMRSpinlabel.hh>
#include <protocols/nmr/nmrspinlabel_util.hh>
#include <core/scoring/nmr/NMRDummySpinlabelEnsemble.hh>
#include <core/scoring/nmr/pcs/PCSSingle.hh>
#include <core/scoring/nmr/pcs/PCSSingleSet.hh>
#include <core/scoring/nmr/pcs/PCSMultiSet.hh>
#include <core/scoring/nmr/pcs/PCSData.hh>
#include <core/scoring/nmr/pcs/PCSTensor.hh>
#include <core/scoring/nmr/util.hh>
#include <protocols/nmr/pcs/PCSTensorOptimizer.hh>

// Package headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/pose/symmetry/util.hh>
#include <core/conformation/Residue.hh>
#include <core/scoring/methods/WholeStructureEnergy.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <basic/datacache/BasicDataCache.hh>
#include <core/optimization/Minimizer.hh>
#include <core/optimization/MinimizerOptions.hh>
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
#include <utility/string_util.hh>

// Numeric headers
#include <numeric/xyzVector.hh>
#include <numeric/xyzMatrix.hh>
#include <numeric/constants.hh>
#include <numeric/HomogeneousTransform.hh>
#include <numeric/geometry/BoundingBox.hh>

// ObjexxFCL headers
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/FArray1D.hh>

// C++ headers
#include <string>
#include <limits>
#include <iostream>
#include <sstream>
#include <iomanip>

namespace protocols {
namespace nmr {
namespace pcs {

static basic::Tracer TR( "protocols.nmr.pcs.PCSEnergy" );

/// @brief Instantiate a new PCSEnergy
core::scoring::methods::EnergyMethodOP
PCSEnergyCreator::create_energy_method(core::scoring::methods::EnergyMethodOptions const &) const {
	return core::scoring::methods::EnergyMethodOP( new PCSEnergy );
}

/// @brief Return the set of score types claimed by the EnergyMethod
/// this EnergyMethodCreator creates in its create_energy_method() function
core::scoring::ScoreTypes
PCSEnergyCreator::score_types_for_method() const {
	core::scoring::ScoreTypes sts;
	sts.push_back( core::scoring::nmr_pcs );
	return sts;
}

/// @brief default constructor
PCSEnergy::PCSEnergy() :
	core::scoring::methods::WholeStructureEnergy( core::scoring::methods::EnergyMethodCreatorOP( new PCSEnergyCreator ) )
{ }

/// @brief copy constructor
PCSEnergy::PCSEnergy(PCSEnergy const & other) :
	core::scoring::methods::WholeStructureEnergy( other ),
	atom_id_to_pcs_xyz_deriv_map_(other.atom_id_to_pcs_xyz_deriv_map_)
{ }

/// @brief destructor
PCSEnergy::~PCSEnergy() { }

core::scoring::methods::EnergyMethodOP
PCSEnergy::clone() const {
	return core::scoring::methods::EnergyMethodOP( new PCSEnergy );
}

/// @brief Return PCSData from pose. Create PCSData if not present and attach them to the pose.
core::scoring::nmr::pcs::PCSData &
PCSEnergy::get_pcs_data_from_pose(Pose & pose) const {
	using namespace core::pose::datacache;
	using namespace core::scoring::nmr;
	using namespace core::scoring::nmr::pcs;

	if ( pose.data().has( CacheableDataType::NMR_PCS_DATA ) ) {
		return *( utility::pointer::static_pointer_cast< PCSData > ( pose.data().get_ptr( CacheableDataType::NMR_PCS_DATA ) ) );
	} else {
		std::string inputfile;
		if ( basic::options::option[ basic::options::OptionKeys::nmr::pcs::input_file ].active() ) {
			inputfile = basic::options::option[ basic::options::OptionKeys::nmr::pcs::input_file ]();
		} else {
			utility_exit_with_message( "No PCS input file given. You must provide the input file on the command line \"-nmr::pcs::input_file filename\"." );
		}
		PCSDataOP pcs_data_all_ptr(utility::pointer::static_pointer_cast< PCSData >( NMRDataFactory::get_instance()->get_nmr_data("PCS", inputfile, pose) ) );
		pose.data().set( CacheableDataType::NMR_PCS_DATA, pcs_data_all_ptr );
		return *pcs_data_all_ptr;
	}
}

/// @brief Calculate the PCS score and write it into the pose's total energy EnergyMap
void
PCSEnergy::finalize_total_energy(
	Pose & pose,
	ScoreFunction const & /*sxfn*/,
	EnergyMap & totals
) const
{
	totals[ core::scoring::nmr_pcs ] = calcualate_total_score_and_tensors( pose );
}

/// @brief Calculate the total PCS score from PCSData retrieved from the pose
core::Real
PCSEnergy::calcualate_total_score_and_tensors(Pose & pose) const {
	using namespace core::scoring::nmr;
	using namespace core::scoring::nmr::pcs;

	PCSData & pcs_data_all = get_pcs_data_from_pose(pose);

	// Some basic setup
	Real total_score(0);
	utility::vector1< PCSMultiSetOP > & multiset_vec = pcs_data_all.get_pcs_multiset_vec();
	Size number_tags(pcs_data_all.get_number_tags());

	// Loop over all tagging sites
	for ( Size i = 1; i <= number_tags; ++i ) {

		// Split behavior depending on if we are going to solve the tensor or
		// calculate the score from fixed tensor values.
		if ( multiset_vec[i]->tensors_fixed() ) {
			Size number_metal_ions(multiset_vec[i]->get_number_metal_ions());
			Size tag_residue_number(multiset_vec[i]->get_tag_residue_number());
			utility::vector1< PCSSingleSetOP > & singleset_vec = multiset_vec[i]->get_pcs_singleset_vec();
			Real score_per_tag(0);

			for ( Size j = 1; j <= number_metal_ions; ++j ) {
				singleset_vec[j]->update_spin_coordinates(pose);
				score_per_tag += singleset_vec[j]->compute_pcs_values_and_score_from_tensor() * singleset_vec[j]->get_weight();
				singleset_vec[j]->set_atom_derivatives( pose );
			}
			if ( TR.Trace.visible() ) {
				TR.Trace << "PCS Score for tagging site at residue " << tag_residue_number << " " << score_per_tag << " (unweighted) "
					<< "or " << score_per_tag*multiset_vec[i]->get_weight() << " (weighted)." << std::endl;
				for ( Size j = 1; j <= number_metal_ions; ++j ) {
					TR.Trace << "Dataset: " << singleset_vec[j]->get_dataset_name() << std::endl;
					singleset_vec[j]->get_tensor()->show_tensor_stats(TR.Trace, true);
				}
			}
			total_score += score_per_tag * multiset_vec[i]->get_weight();
			// Solve tensor and calculate score using either the NMRSpinlabel or GridSearch method
		} else {
			// Calculate with spinlabel
			if ( multiset_vec[i]->get_spinlabel() ) {
				total_score += calculate_score_and_tensors_with_spinlabel(pose, *(multiset_vec[i]), pcs_data_all.optimize_tensors());
				// Calculate with grid search
			} else if ( multiset_vec[i]->get_gridsearch_iterator() ) {
				total_score += calculate_score_and_tensors_with_gridsearch(pose, *(multiset_vec[i]), pcs_data_all.optimize_tensors());
				// Raise error when no spinlabel or grid search
			} else {
				utility_exit_with_message("ERROR during PCS score calculation. No NMR spinlabel or gridsearch iterator are set for PCS dataset at spinlabel site "
					+ utility::to_string(multiset_vec[i]->get_tag_residue_number()));
			}
		}
	}
	TR.Trace << "Total PCS score " << total_score << std::endl;
	return total_score;
}

/// @bried Calculate PCS score for input PCSMultiSet with spinlabel
core::Real
PCSEnergy::calculate_score_and_tensors_with_spinlabel(
	Pose & pose,
	PCSMultiSet & pcs_data,
	bool optimize_tensors
) const
{
	using namespace core::scoring::nmr;
	using namespace core::scoring::nmr::pcs;
	using namespace core::pose::symmetry;
	using WeightCoordVector = NMRSpinlabel::WeightCoordVector;

	// Some basic setup
	Size number_metal_ions(pcs_data.get_number_metal_ions());
	Size tag_residue_number(pcs_data.get_tag_residue_number());
	Real score_this_tag(0);
	utility::vector1< PCSSingleSetOP > & singleset_vec = pcs_data.get_pcs_singleset_vec();
	runtime_assert(number_metal_ions==singleset_vec.size());
	utility::vector1< PCSSingleSetOP > singlesets_for_svd;
	utility::vector1< PCSSingleSetOP > singlesets_for_nls;

	// Update spin coordinates and group PCSSingleSets according to their computation type
	for ( Size i(1); i <= number_metal_ions; ++i ) {
		singleset_vec[ i ]->update_spin_coordinates(pose);
		if ( singleset_vec[ i ] ->get_computation_type() == PCSSingleSet::SVD ) {
			singlesets_for_svd.push_back( singleset_vec[ i ] );
		} else {
			singlesets_for_nls.push_back( singleset_vec[ i ] );
		}
	}

	NMRSpinlabelOP nmrspinlabel(pcs_data.get_spinlabel());
	WeightCoordVector spinlabel_wghts_coords;

	// Pose is in full atom mode and spinlabel clash score is approximated by calculating bump energy
	if ( pose.is_fullatom() && nmrspinlabel->get_highres_conformer_filter_type() == NMRSpinlabel::BUMP_ENERGY ) {

		spinlabel_wghts_coords = filter_spinlabel_ensemble_by_packerenergy(pose, *nmrspinlabel, tag_residue_number);

	} else {
		// Filter dummy ensemble by clash score calculation and get the spinlabel coordinates
		// We can simply call the member method of the NMRSpinlabel. Internally, it performs the
		// clash filter, looks up the coordinates of the radical atom, transforms their coordinates
		// into the coordinate frame of the spinlabel site and performs additional binning if the
		// vector size exceeds the maximal number of spinlabel conformers
		runtime_assert_msg( nmrspinlabel->get_dummy_ensemble(),
			"ERROR while trying to calculate PCS score with implicit spinlabel model. Dummy ensemble of NMRSpinlabel not set. Check if NMRSpinlabel database file exists in Rosetta database." );
		spinlabel_wghts_coords = nmrspinlabel->filter_spinlabel_ensemble_by_distance_check(pose, tag_residue_number);
	}

	// Average metal ion coordinates and bounding box of the complete xyz range
	Vector sl_av_metal_coords(0.0);
	Real totwght(0.0);
	numeric::geometry::BoundingBox<Vector> bbox;
	bbox.set_lower(Vector(std::numeric_limits< Real >::max()));
	bbox.set_upper(Vector(std::numeric_limits< Real >::min()));
	for ( Size k = 1, k_end = spinlabel_wghts_coords.size(); k <= k_end; ++k ) {
		totwght += spinlabel_wghts_coords[k].first;
		sl_av_metal_coords += spinlabel_wghts_coords[k].first * spinlabel_wghts_coords[k].second;
		bbox.add(spinlabel_wghts_coords[k].second);
	}
	sl_av_metal_coords /= totwght;

	// Solve by SVD
	if ( singlesets_for_svd.size() ) {
		Real score_svd = calculate_score_and_tensors_by_svd(pose, singlesets_for_svd, sl_av_metal_coords, tag_residue_number, optimize_tensors);
		score_this_tag += score_svd;
	}
	// Solve by NLS
	if ( singlesets_for_nls.size() ) {
		utility::fixedsizearray1<Real,6> metal_coord_ranges_for_nls;
		metal_coord_ranges_for_nls[1] = bbox.lower().x();
		metal_coord_ranges_for_nls[2] = bbox.lower().y();
		metal_coord_ranges_for_nls[3] = bbox.lower().z();
		metal_coord_ranges_for_nls[4] = bbox.upper().x();
		metal_coord_ranges_for_nls[5] = bbox.upper().y();
		metal_coord_ranges_for_nls[6] = bbox.upper().z();
		Real score_nls = calculate_score_and_tensors_by_nls(pose, singlesets_for_nls, sl_av_metal_coords, metal_coord_ranges_for_nls, tag_residue_number);
		score_this_tag += score_nls;
	}
	TR.Trace << "PCS Score for spinlabel site at residue " << tag_residue_number << " " << score_this_tag << " (unweighted) "
		<< "or " << score_this_tag*pcs_data.get_weight() << " (weighted)." << std::endl;

	return score_this_tag*pcs_data.get_weight();
}

/// @bried Calculate PCS score for input PCSMultiSet with grid search
core::Real
PCSEnergy::calculate_score_and_tensors_with_gridsearch(
	Pose & pose,
	PCSMultiSet & pcs_data,
	bool optimize_tensors
) const
{
	using namespace core::scoring::nmr;
	using namespace core::scoring::nmr::pcs;

	// Some basic setup
	Size number_metal_ions(pcs_data.get_number_metal_ions());
	Size tag_residue_number(pcs_data.get_tag_residue_number());
	Real score_this_tag(0);
	utility::vector1< PCSSingleSetOP > & singleset_vec = pcs_data.get_pcs_singleset_vec();
	runtime_assert(number_metal_ions==singleset_vec.size());
	utility::vector1< PCSSingleSetOP > singlesets_for_svd;
	utility::vector1< PCSSingleSetOP > singlesets_for_nls;

	// Update spin coordinates and group PCSSingleSets according to their computation type
	for ( Size i(1); i <= number_metal_ions; ++i ) {
		singleset_vec[ i ]->update_spin_coordinates(pose);
		if ( singleset_vec[ i ] ->get_computation_type() == PCSSingleSet::SVD ) {
			singlesets_for_svd.push_back( singleset_vec[ i ] );
		} else {
			singlesets_for_nls.push_back( singleset_vec[ i ] );
		}
	}

	NMRGridSearchOP gridsearch_iterator(pcs_data.get_gridsearch_iterator());
	gridsearch_iterator->set_grid_search_center(pose);
	Vector grid_search_center = gridsearch_iterator->get_grid_search_center();

	// Solve by SVD
	if ( singlesets_for_svd.size() ) {
		Vector current_metal_coords;
		Vector best_metal_coords;
		Real score_svd(0);
		Real best_score_svd(std::numeric_limits<Real>::max());
		// Run grid search
		while ( gridsearch_iterator->valid_next_grid_point(current_metal_coords) ) {
			score_svd = 0;
			for ( Size j(1); j <= singlesets_for_svd.size() ; ++j ) {
				score_svd += singlesets_for_svd[ j ]->solve_tensor_and_compute_score_by_svd(current_metal_coords) * singlesets_for_svd[ j ]->get_weight();
				if ( score_svd > best_score_svd ) { // if one dataset gives already a score higher than the best score, there is no need to calculate the scores of the remaining datasets
					break;
				}
			}
			if ( score_svd < best_score_svd ) {
				best_score_svd = score_svd;
				best_metal_coords = current_metal_coords;
				gridsearch_iterator->set_best_grid_point(best_metal_coords);
			}
		}
		// Do score calculation by SVD one more time to update calculated PCS values and optionally do optimization
		best_score_svd = calculate_score_and_tensors_by_svd(pose, singlesets_for_svd, best_metal_coords, tag_residue_number, optimize_tensors);
		score_this_tag += best_score_svd;
	}
	// Solve by NLS
	if ( singlesets_for_nls.size() ) {
		utility::fixedsizearray1<Real,6> metal_coord_ranges_for_nls;
		metal_coord_ranges_for_nls[1] = grid_search_center.x() - gridsearch_iterator->get_grid_max_radius();
		metal_coord_ranges_for_nls[2] = grid_search_center.y() - gridsearch_iterator->get_grid_max_radius();
		metal_coord_ranges_for_nls[3] = grid_search_center.z() - gridsearch_iterator->get_grid_max_radius();
		metal_coord_ranges_for_nls[4] = grid_search_center.x() + gridsearch_iterator->get_grid_max_radius();
		metal_coord_ranges_for_nls[5] = grid_search_center.y() + gridsearch_iterator->get_grid_max_radius();
		metal_coord_ranges_for_nls[6] = grid_search_center.z() + gridsearch_iterator->get_grid_max_radius();
		Real score_nls = calculate_score_and_tensors_by_nls(pose, singlesets_for_nls, grid_search_center, metal_coord_ranges_for_nls, tag_residue_number);
		score_this_tag += score_nls;
	}
	TR.Trace << "PCS Score for spinlabel site at residue " << tag_residue_number << " " << score_this_tag << " (unweighted) "
		<< "or " << score_this_tag*pcs_data.get_weight() << " (weighted)." << std::endl;

	return score_this_tag*pcs_data.get_weight();
}

/// @brief Minimize input PCS tensors with line search
core::Real
PCSEnergy::minimize_tensors(
	utility::vector1< Real > & tensor_params_for_optimization,
	utility::vector1< core::scoring::nmr::pcs::PCSSingleSetOP > & singlesets_for_optimization
) const
{
	PCSTensorOptimizer optimizer(singlesets_for_optimization);
	core::optimization::MinimizerOptions options( "lbfgs_armijo_nonmonotone", 0.00001, true, false, false );
	core::optimization::Minimizer minimizer(optimizer, options );
	Real optimized_score = minimizer.run( tensor_params_for_optimization );
	return optimized_score;
}

/// @brief Calculate PCS score for input PCSSingleSets at input metal position by SVD
core::Real
PCSEnergy::calculate_score_and_tensors_by_svd(
	Pose & pose,
	utility::vector1< core::scoring::nmr::pcs::PCSSingleSetOP > & pcs_singlesets,
	Vector & metal_coords,
	Size tag_residue_number,
	bool optimize_tensors
) const
{
	Real score_svd(0);
	for ( Size i(1); i <= pcs_singlesets.size(); ++i ) {
		score_svd += pcs_singlesets[ i ]->solve_tensor_and_compute_score_by_svd(metal_coords) * pcs_singlesets[ i ]->get_weight();
		pcs_singlesets[ i ]->set_atom_derivatives( pose );
	}

	// Do tensor optimization
	if ( optimize_tensors ) {
		utility::vector1< Real > tensor_params_for_optimization;
		tensor_params_for_optimization.push_back(metal_coords.x());
		tensor_params_for_optimization.push_back(metal_coords.y());
		tensor_params_for_optimization.push_back(metal_coords.z());

		TR.Trace << "SVD Score for tagging site at residue " << tag_residue_number << " before optimization = " << score_svd << std::endl;
		TR.Trace << "Tensors before optimization" << std::endl;
		for ( Size i(1); i <= pcs_singlesets.size(); ++i ) {
			tensor_params_for_optimization.push_back( pcs_singlesets[ i ]->get_tensor_const()->get_T_xx() );
			tensor_params_for_optimization.push_back( pcs_singlesets[ i ]->get_tensor_const()->get_T_xy() );
			tensor_params_for_optimization.push_back( pcs_singlesets[ i ]->get_tensor_const()->get_T_xz() );
			tensor_params_for_optimization.push_back( pcs_singlesets[ i ]->get_tensor_const()->get_T_yy() );
			tensor_params_for_optimization.push_back( pcs_singlesets[ i ]->get_tensor_const()->get_T_yz() );

			TR.Trace << "Dataset: " << pcs_singlesets[ i ]->get_dataset_name() << std::endl;
			pcs_singlesets[ i ]->get_tensor_const()->show_tensor_stats(TR.Trace, false);
		}

		// Calculate optimized tensors and score. Note that the score value already includes the singleset weights.
		Real optimized_score = minimize_tensors(tensor_params_for_optimization, pcs_singlesets);

		TR.Trace << "SVD Score for tagging site at residue " << tag_residue_number << " after optimization = " << optimized_score << std::endl;
		if ( optimized_score < score_svd ) { // If score improved after optimization set new tensor parameters
			score_svd = optimized_score;
			TR.Trace << "Tensors after optimization" << std::endl;
			Size k(1);
			for ( Size i(1); i <= pcs_singlesets.size(); ++i ) {
				core::scoring::nmr::pcs::PCSTensorOP tensor = pcs_singlesets[ i ]->get_tensor();
				tensor->set_metal_center(tensor_params_for_optimization[1], tensor_params_for_optimization[2], tensor_params_for_optimization[3]);
				tensor->set_T_xx(tensor_params_for_optimization[3 + 5*(k-1) + 1]);
				tensor->set_T_xy(tensor_params_for_optimization[3 + 5*(k-1) + 2]);
				tensor->set_T_xz(tensor_params_for_optimization[3 + 5*(k-1) + 3]);
				tensor->set_T_yy(tensor_params_for_optimization[3 + 5*(k-1) + 4]);
				tensor->set_T_yz(tensor_params_for_optimization[3 + 5*(k-1) + 5]);
				++k;
				pcs_singlesets[ i ]->set_atom_derivatives( pose ); // Set atom derivatives
				pcs_singlesets[ i ]->compute_pcs_values_and_score_from_tensor( *tensor ); // Update calculated PCS values with new tensor
				TR.Trace << "Dataset: " << pcs_singlesets[ i ]->get_dataset_name() << std::endl;
				tensor->show_tensor_stats(TR.Trace, false);
			}
		} else {
			TR.Warning << "Warning: Score after optimization higher = " << std::scientific << std::setprecision(12) << optimized_score << ". Keep the old score = " << score_svd << " and tensor parameter." << std::endl;
		}
	} else {
		if ( TR.Trace.visible() ) {
			TR.Trace << "SVD Score for tagging site at residue " << tag_residue_number << " " << score_svd << std::endl;
			for ( Size i(1); i <= pcs_singlesets.size(); ++i ) {
				pcs_singlesets[ i ]->get_tensor_const()->show_tensor_stats(TR.Trace, false);
			}
		}
	}
	return score_svd;
}

/// @brief Calculate PCS score for input PCSSingleSets at input metal position by NLS
core::Real
PCSEnergy::calculate_score_and_tensors_by_nls(
	Pose & pose,
	utility::vector1< core::scoring::nmr::pcs::PCSSingleSetOP > & pcs_singlesets,
	Vector & metal_coords,
	utility::fixedsizearray1<Real,6> const & metal_coord_range,
	Size tag_residue_number
) const
{
	Real score_nls(0);
	for ( Size i(1); i <= pcs_singlesets.size(); ++i ) {
		pcs_singlesets[ i ]->set_metal_coord_bounds(metal_coord_range);
		score_nls += pcs_singlesets[ i ]->solve_tensor_and_compute_score_by_nls(metal_coords) * pcs_singlesets[ i ]->get_weight();
		pcs_singlesets[ i ]->set_atom_derivatives( pose ); // Set atom derivatives
	}
	if ( TR.Trace.visible() ) {
		TR.Trace << "NLS Score for tagging site at residue " << tag_residue_number << " " << score_nls << std::endl;
		for ( Size i(1); i <= pcs_singlesets.size(); ++i ) {
			core::scoring::nmr::pcs::PCSTensorCOP tensor = pcs_singlesets[ i ]->get_tensor_const();
			TR.Trace << "Dataset: " << pcs_singlesets[ i ]->get_dataset_name() << std::endl;
			tensor->show_tensor_stats(TR.Trace, true);
		}
	}
	return score_nls;
}

/// @brief Called at the beginning of atom tree minimization, this method
///        allows the derived class the opportunity to initialize pertinent data
///        that will be used during minimization.
///        Here, the function creates and updates the atom_id_to_pcs_xyz_deriv_map_ which
///        is needed by the eval_atom_derivative() function.
void
PCSEnergy::setup_for_minimizing(
	Pose & pose ,
	ScoreFunction const & /*sxfn*/,
	core::kinematics::MinimizerMapBase const & /*minmap*/
) const
{
	using namespace core::scoring::nmr::pcs;

	PCSData & pcs_data_all = get_pcs_data_from_pose(pose);

	// Set to zero because we are using operator+=
	Vector fij(0 ,0 ,0);
	atom_id_to_pcs_xyz_deriv_map_.default_value(fij);
	atom_id_to_pcs_xyz_deriv_map_.fill_with(fij);

	utility::vector1< PCSMultiSetOP > & multiset_vec = pcs_data_all.get_pcs_multiset_vec();
	Size number_tags(pcs_data_all.get_number_tags());
	utility::vector1<PCSSingle>::const_iterator iter;

	for ( Size i = 1; i <= number_tags; ++i ) {
		utility::vector1< PCSSingleSetOP > & singleset_vec = multiset_vec[i]->get_pcs_singleset_vec();
		Size number_metals = multiset_vec[i]->get_number_metal_ions();

		for ( Size j = 1; j <= number_metals; ++j ) {
			utility::vector1<PCSSingle> const & single_pcs_vec = singleset_vec[j]->get_single_pcs_vec();
			core::conformation::symmetry::SymmetryInfoCOP syminfo_ptr;
			if ( singleset_vec[j]->symmetric_pcs_calc() && core::pose::symmetry::is_symmetric(pose) ) {
				syminfo_ptr = core::pose::symmetry::symmetry_info( pose );
			}

			for ( iter = single_pcs_vec.begin(); iter != single_pcs_vec.end(); ++iter ) {
				utility::vector1< core::id::AtomID > const & protein_spins = iter->get_protein_spins();
				utility::vector1< Vector > const & atom_derivatives = iter->get_atom_derivatives();
				runtime_assert_msg(protein_spins.size() == atom_derivatives.size(), "ERROR in setup of atom tree minimization from PCS derivatives. Vector of AtomIDs and of derivatives have unequal length.");
				for ( Size k = 1, k_end = protein_spins.size(); k <= k_end; ++k ) {
					fij = atom_id_to_pcs_xyz_deriv_map_.get( protein_spins[k] );  // returns (0,0,0) if atom is not present in map, otherwise the corresponding derivative
					fij += atom_derivatives[k];                         // accumulates the derivative for that particular atom (for all pcs data)
					atom_id_to_pcs_xyz_deriv_map_.set( protein_spins[k], fij);   // sets the accumulated derivative and resizes the map if necessary
				}

				// Do another pass to set PCS derivatives for symmetric residues
				if ( syminfo_ptr ) {
					for ( Size k = 1, k_end = protein_spins.size(); k <= k_end; ++k ) {
						utility::vector1< Size >protein_spins_symm_rsd = syminfo_ptr->bb_clones(protein_spins[k].rsd());
						for ( Size l = 1, l_end = protein_spins_symm_rsd.size(); l <= l_end; ++l ) {
							core::id::AtomID symm_spin(protein_spins[k].atomno(), protein_spins_symm_rsd[l]);

							// return the derivative of that atom, if not present in map, return the default value which is (0,0,0)
							Vector fij_symm_spin = atom_id_to_pcs_xyz_deriv_map_.get(symm_spin);

							// rotate derivative vector into the frame of the symmetric subunit
							// and accumulate the contribution from all pcs experiments
							fij_symm_spin += core::scoring::nmr::apply_vector_rotation(pose, atom_derivatives[k], protein_spins[k].rsd(), protein_spins_symm_rsd[l]);

							// set the accumulated derivative and resize the map if necessary
							atom_id_to_pcs_xyz_deriv_map_.set(symm_spin, fij_symm_spin);
						}
					}
				}
			}
		}
	}
}

/// @brief Evaluate the xyz derivative of the PCS for an atom in the pose.
///        Called during the atomtree derivative calculation, atom_tree_minimize.cc,
///        through the ScoreFunction::eval_atom_derivative intermediary.
///        F1 and F2 should not be zeroed, rather, setup_for_minimizing() accumulates the PCS
///        contribution from the xyz derivatives of atom id
void
PCSEnergy::eval_atom_derivative(
	core::id::AtomID const & id,
	Pose const & pose,
	core::kinematics::DomainMap const & /*domain_map*/,
	ScoreFunction const & /*sfxn*/,
	EnergyMap const & weights,
	Vector & F1,
	Vector & F2
) const
{
	if ( !atom_id_to_pcs_xyz_deriv_map_.has(id) ) {
		// only checks if id residue number is in range of AtomID_Map
		return;
	}
	Vector fij = atom_id_to_pcs_xyz_deriv_map_.get(id); //  returns (0,0,0) if id is not present in map, otherwise the accumulated derivative of the PCS for this atom
	Vector atom_x = pose.xyz(id);
	Vector const f2( -fij );
	Vector const atom_y = atom_x - f2;   // a "fake" atom in the direction of the gradient
	Vector const f1( atom_x.cross( atom_y ) );

	F1 += weights[ core::scoring::nmr_pcs ] * f1;
	F2 += weights[ core::scoring::nmr_pcs ] * f2;
}

/// @brief register options
void
PCSEnergy::register_options() {
	using namespace basic::options;
	option.add_relevant(OptionKeys::nmr::pcs::show_info);
}

/// @brief Indicate in the context-graphs-required list which
///        context-graphs this energy method requires that the Pose
///        maintain when doing neighbor evaluation. Context graphs are allowed.
void
PCSEnergy::indicate_required_context_graphs(utility::vector1< bool > &) const { }

/// @brief Return the version of the energy method
core::Size
PCSEnergy::version() const { return 1; }

/// @brief show additional information of the energy method
void
PCSEnergy::show_additional_info(
	std::ostream & tracer,
	Pose & pose,
	bool verbose
) const
{
	using namespace core::scoring::nmr;
	using namespace core::scoring::nmr::pcs;

	if ( basic::options::option[ basic::options::OptionKeys::nmr::pcs::show_info ].user() ) {
		verbose = basic::options::option[ basic::options::OptionKeys::nmr::pcs::show_info ]();
	}

	// Some basic setup
	PCSData & pcs_data_all = get_pcs_data_from_pose(pose);
	Real pcs_score = calcualate_total_score_and_tensors(pose);
	Size number_tags = pcs_data_all.get_number_tags();
	utility::vector1< utility::vector1< PCSTensor > > tensors_all_tags_and_metal_ions(number_tags);
	utility::vector1< Real > scores_all_tags(number_tags, 0.0); // set all score to zero because we are going to use +=
	utility::vector1< PCSMultiSetOP > & multiset_vec = pcs_data_all.get_pcs_multiset_vec();

	tracer << " * * * * * * *       PCSEnergy Info       * * * * * * * " << std::endl;
	tracer << "Total Score: " << pcs_score << std::endl;
	tracer << "Is norm of " << number_tags << " tagging sites." << std::endl;

	// Back-calculate PCSs from determined tensor params
	// Show calc vs exp PCS and tensor params
	// Loop over PCSMultiSets
	Real sum_all_dev_square(0);
	Real sum_all_exp_square(0);
	Size number_all_pcs(0);
	Size total_number_experiments(0);
	for ( Size i = 1; i <= number_tags; ++i ) {
		Size number_metal_ions(multiset_vec[i]->get_number_metal_ions());
		tensors_all_tags_and_metal_ions[i].resize(number_metal_ions);
		utility::vector1< PCSSingleSetOP > & singleset_vec = multiset_vec[i]->get_pcs_singleset_vec();

		Real sum_dev_square(0);
		Real sum_exp_square(0);
		Size number_pcs_per_tag(0);
		// Loop over PCSSingleSets
		for ( Size j = 1; j <= number_metal_ions; ++j ) {

			// Update xyz coordinates of active spins and do some basic checks
			singleset_vec[j]->update_spin_coordinates(pose);
			utility::vector1< utility::vector1< utility::vector1< Vector > > > const & spin_coordinates = singleset_vec[j]->get_spin_coordinates();
			ObjexxFCL::FArray1D<Real> const & exp_pcs_values = singleset_vec[j]->get_pcs_values();
			ObjexxFCL::FArray1D<Real> const & pcs_single_weights = singleset_vec[j]->get_pcs_single_weights();
			Size npcs(singleset_vec[j]->get_number_pcs());
			number_pcs_per_tag += npcs;
			runtime_assert_msg(npcs == spin_coordinates.size(), "ERROR in PCSEnergy's show_additional_info() function. Number of PCS and length of spin coordinate vector are not the same");
			runtime_assert_msg(npcs == exp_pcs_values.size(), "ERROR in PCSEnergy's show_additional_info() function. Number of PCS and length of PCS value vector are not the same");
			runtime_assert_msg(npcs == pcs_single_weights.size(), "ERROR in PCSEnergy's show_additional_info() function. Number of PCS and length of PCS weights vector are not the same");

			// Bring all tensors in unified tensor frame representation
			PCSTensor tensor( *( singleset_vec[j]->get_tensor() ) );
			tensor.diagonalize_tensor();
			tensor.reorder_tensor();
			tensors_all_tags_and_metal_ions[i][j] = tensor;
			Real xM(tensor.get_metal_center().x());
			Real yM(tensor.get_metal_center().y());
			Real zM(tensor.get_metal_center().z());
			Real Xax(tensor.get_ax());
			Real Xrh(tensor.get_rh());
			Real alpha(tensor.get_alpha());
			Real beta(tensor.get_beta());
			Real gamma(tensor.get_gamma());

			Real scal(1.0 / singleset_vec[j]->get_scaling_factor());
			Vector euler_angles(alpha, beta, gamma);
			numeric::xyzMatrix< Real > rotM(rotation_matrix_from_euler_angles(euler_angles, tensor.get_euler_convention()));
			Real single_score(0);

			tracer << "PCS dataset: " << singleset_vec[j]->get_dataset_name() << std::endl;
			if ( verbose ) {
				tracer << " * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * " << std::endl;
				tracer << " * * * * * * * * * * * * * * * *              exp vs. calc PCS             * * * * * * * * * * * * * * * * * " << std::endl;
				tracer << " * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * " << std::endl;
				tracer << " Spin_Resid_Atom     Obs_PCS     Scal_Obs_PCS     Calc_PCS     Scal_Calc_PCS     Deviation     Abs_Deviation " << std::endl;
				tracer << " ----------------------------------------------------------------------------------------------------------- " << std::endl;
			}

			// Calculate all PCSs for this PCSSingleSet
			for ( Size k = 1; k <= npcs; ++k ) {
				Real calc_pcs(0);
				Size n_subunits(spin_coordinates[k].size());
				for ( Size l = 1; l <= n_subunits; ++l ) {
					Real calc_pcs_per_subunit(0);
					Size n_eq_spins(spin_coordinates[k][l].size());
					for ( Size m = 1; m <= n_eq_spins; ++m ) {
						Real x(spin_coordinates[k][l][m].x() - xM);
						Real y(spin_coordinates[k][l][m].y() - yM);
						Real z(spin_coordinates[k][l][m].z() - zM);

						Real x_t(rotM(1,1)*x + rotM(1,2)*y + rotM(1,3)*z);
						Real y_t(rotM(2,1)*x + rotM(2,2)*y + rotM(2,3)*z);
						Real z_t(rotM(3,1)*x + rotM(3,2)*y + rotM(3,3)*z);

						Real r2(x_t*x_t + y_t*y_t + z_t*z_t);
						Real r5(r2 * r2 * std::sqrt(r2));
						Real value_1_12_PI_r5(10000.0 / (12.0 * numeric::constants::d::pi * r5));

						calc_pcs_per_subunit += scal * value_1_12_PI_r5 * (Xax * (3.0 * z_t*z_t - r2) + Xrh * 1.5 * (x_t*x_t - y_t*y_t));
					} // no equivalent spins (e.g. CH3 protons)
					if ( singleset_vec[j]->get_averaging_type() == core::scoring::nmr::MEAN ) { calc_pcs_per_subunit /= n_eq_spins; }
					calc_pcs += calc_pcs_per_subunit;
				} // no protein subunits
				Real pcs_dev(calc_pcs - exp_pcs_values(k));
				sum_dev_square += pcs_dev*pcs_dev;
				sum_exp_square += exp_pcs_values(k)*exp_pcs_values(k);
				Real pcs_abs_dev(std::abs(pcs_dev));
				single_score += pcs_dev * pcs_dev * pcs_single_weights(k);

				if ( verbose ) {
					std::ostringstream converter;
					for ( core::id::AtomID const & id : singleset_vec[j]->get_single_pcs_vec()[k].get_protein_spins() ) {
						converter << std::setw(4) << id.rsd() << " " << pose.residue(id.rsd()).atom_name(id.atomno()) << " ";
					}
					tracer << " ( " << std::right << converter.str() << ")";
					tracer << std::setw(12)<< std::fixed << std::setprecision(3) << exp_pcs_values(k) / scal;  // Revert normalization of PCS data that was applied during creation of PCSSingleSet
					tracer << std::setw(17)<< std::fixed << std::setprecision(3) << exp_pcs_values(k);
					tracer << std::setw(13)<< std::fixed << std::setprecision(3) << calc_pcs / scal;    // Revert normalization of PCS data that was applied during creation of PCSSingleSet
					tracer << std::setw(18)<< std::fixed << std::setprecision(3) << calc_pcs;
					tracer << std::setw(14)<< std::fixed << std::setprecision(3) << pcs_dev;
					tracer << std::setw(18)<< std::fixed << std::setprecision(3) << pcs_abs_dev << std::endl;
				}
			}
			if ( verbose ) {
				tracer << " ----------------------------------------------------------------------------------------------------------- " << std::endl;
			}
			//scores_all_tags[i] += std::sqrt(single_score) * singleset_vec[j]->get_weight();
			scores_all_tags[i] += single_score * singleset_vec[j]->get_weight();
			tracer << "Weighted score for dataset " << singleset_vec[j]->get_dataset_name() << " " << single_score << std::endl;
			tensor.show_tensor_stats(tracer, true);
		}
		tracer << "Overall weighted score for tagging site " << multiset_vec[i]->get_tag_residue_number() << ": " << scores_all_tags[i] << std::endl;
		Real Q_fac_per_tag(std::sqrt(sum_dev_square/sum_exp_square));
		Real Rmsd_per_tag(std::sqrt(sum_dev_square/number_pcs_per_tag));
		tracer << "PCS Q-factor and Rmsd for tagging site " << multiset_vec[i]->get_tag_residue_number() << ": " << Q_fac_per_tag << ", " << Rmsd_per_tag << std::endl;
		sum_all_dev_square += sum_dev_square;
		sum_all_exp_square += sum_exp_square;
		number_all_pcs += number_pcs_per_tag;
		total_number_experiments += number_metal_ions;
	}
	Real total_Q_fac(std::sqrt(sum_all_dev_square/sum_all_exp_square));
	Real total_Rmsd(std::sqrt(sum_all_dev_square/number_all_pcs));
	tracer << "Total PCS Q-factor and Rmsd for " << number_tags << " tagging sites and " << total_number_experiments
		<< " experiments: " << total_Q_fac << ", " << total_Rmsd << std::endl;
}


} // namespace pcs
} // namespace nmr
} // namespace protocols

