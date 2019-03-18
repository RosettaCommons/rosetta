// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    core/scoring/nmr/pcs/PCSMultiSet.cc
/// @brief   Implementation of class PCSMultiSet
/// @details last Modified: 06/30/16
/// @author  Georg Kuenze (georg.kuenze@vanderbilt.edu)

// Unit headers
#include <core/scoring/nmr/pcs/PCSMultiSet.hh>

// Package headers
#include <core/scoring/nmr/pcs/PCSSingleSet.hh>
#include <core/scoring/nmr/pcs/PCSTensor.hh>
#include <core/scoring/nmr/NMRGridSearch.hh>
#include <core/scoring/nmr/NMRSpinlabel.hh>
#include <core/scoring/nmr/NMRDummySpinlabelEnsemble.hh>
#include <core/scoring/nmr/NMRDummySpinlabelVoxelGrid.hh>

// Project headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/id/AtomID.hh>

// Basic headers
#include <basic/Tracer.hh>

// Utility headers
#include <utility/exit.hh>
#include <utility/vector1.hh>
#include <utility/fixedsizearray1.hh>

// Numeric headers
#include <numeric/xyzVector.hh>
#include <numeric/xyzMatrix.hh>
#include <numeric/HomogeneousTransform.hh>
#include <numeric/geometry/BoundingBox.hh>

// C++ headers
#include <cmath>
#include <limits>
#include <string>
#include <iostream>
#include <algorithm>
#include <numeric>

namespace core {
namespace scoring {
namespace nmr {
namespace pcs {

static basic::Tracer TR( "core.scoring.nmr.pcs.PCSMultiSet" );

/// @brief construct from a vector of PCSSingleSets and the
///        position of the tagging site. Set a default NMRSpinlabel.
PCSMultiSet::PCSMultiSet(
	utility::vector1<PCSSingleSetOP> const & pcs_singleset_vec,
	Size const spinlabel_site,
	Real const weight
) :
	pcs_singleset_vec_(pcs_singleset_vec),
	number_metal_ions_(pcs_singleset_vec.size()),
	tag_residue_number_(spinlabel_site),
	weight_(weight),
	gridsearch_iterator_(nullptr),
	spinlabel_(nullptr),
	fixed_tensors_(false)
{
	set_default_spinlabel();
}

/// @brief construct from a vector of PCSSingleSets, the position
///        of the tagging site and an NMRSpinlabel object
PCSMultiSet::PCSMultiSet(
	utility::vector1<PCSSingleSetOP> const & pcs_singleset_vec,
	Size const spinlabel_site,
	NMRSpinlabelOP spinlabel_ptr,
	Real const weight
) :
	pcs_singleset_vec_(pcs_singleset_vec),
	number_metal_ions_(pcs_singleset_vec.size()),
	tag_residue_number_(spinlabel_site),
	weight_(weight),
	gridsearch_iterator_(nullptr),
	spinlabel_(new NMRSpinlabel(*spinlabel_ptr)),
	fixed_tensors_(false)
{ }

/// @brief construct from a vector of PCSSingleSets, the position
///        of the tagging site and an NMRGridSearch object
PCSMultiSet::PCSMultiSet(
	utility::vector1<PCSSingleSetOP> const & pcs_singleset_vec,
	Size const spinlabel_site,
	NMRGridSearchOP gridsearch_ptr,
	Real const weight
) :
	pcs_singleset_vec_(pcs_singleset_vec),
	number_metal_ions_(pcs_singleset_vec.size()),
	tag_residue_number_(spinlabel_site),
	weight_(weight),
	gridsearch_iterator_(new NMRGridSearch(*gridsearch_ptr)),
	spinlabel_(nullptr),
	fixed_tensors_(false)
{ }

/// @brief copy constructor
PCSMultiSet::PCSMultiSet(PCSMultiSet const & other) :
	number_metal_ions_(other.number_metal_ions_),
	tag_residue_number_(other.tag_residue_number_),
	weight_(other.weight_),
	gridsearch_iterator_( other.gridsearch_iterator_ ? new NMRGridSearch( *(other.gridsearch_iterator_) ) : nullptr ),
	spinlabel_( other.spinlabel_ ? new NMRSpinlabel( *(other.spinlabel_) ) : nullptr ),
	fixed_tensors_(other.fixed_tensors_)
{
	// Make deep copy of PCSSingleSet vector
	deep_copy_pcs_singleset_vec(other.pcs_singleset_vec_);
}

/// @brief assignment operator
PCSMultiSet&
PCSMultiSet::operator=(PCSMultiSet const & rhs) {
	if ( this != &rhs ) {
		number_metal_ions_ = rhs.number_metal_ions_;
		tag_residue_number_ = rhs.tag_residue_number_;
		weight_ = rhs.weight_;
		gridsearch_iterator_ = rhs.gridsearch_iterator_ ? NMRGridSearchOP( new NMRGridSearch( *(rhs.gridsearch_iterator_) ) ) : nullptr;
		spinlabel_ = rhs.spinlabel_ ? NMRSpinlabelOP( new NMRSpinlabel( *(rhs.spinlabel_) ) ) : nullptr;
		fixed_tensors_ = rhs.fixed_tensors_;

		// Make deep copy of PCSSingleSet vector
		pcs_singleset_vec_.clear();
		deep_copy_pcs_singleset_vec(rhs.pcs_singleset_vec_);
	}
	return *this;
}

/// @brief destructor
PCSMultiSet::~PCSMultiSet() { }

/// @brief compute the score of all PCS datasets given the input pose.
///        the PCS tensor is solved by SVD or NLS. In the first case,
///        the metal ion position can be determined by a grid search.
///        Alternatively, if the spinlabel data member has been set,
///        the metal ion position can be inferred from the conformational
///        ensemble of the spinlabel. Here, the metal coordinates are
///        calculated as the weighted average across the ensemble.
Real
PCSMultiSet::compute_score(
	pose::Pose & pose,
	utility::vector1< PCSTensorCOP > & singleset_tensors
)
{
	using WeightCoordVector = NMRSpinlabel::WeightCoordVector;

	Real total_score(0);
	utility::vector1< PCSTensorCOP > singleset_tensors_temp(number_metal_ions_);

	// Split behavior depending on if we are going to solve the tensor or
	// calculate the score from fixed tensor values.
	if ( fixed_tensors_ ) {
		for ( Size i = 1; i <= number_metal_ions_; ++i ) {
			pcs_singleset_vec_[i]->update_spin_coordinates( pose );
			total_score += pcs_singleset_vec_[i]->compute_pcs_values_and_score_from_tensor() * pcs_singleset_vec_[i]->get_weight();
			pcs_singleset_vec_[i]->set_atom_derivatives( pose );
			singleset_tensors_temp[i] = pcs_singleset_vec_[i]->get_tensor_const();
		}
	} else {
		// Split behavior depending on if we are going to use the spinlabel ensemble or the gridsearch

		utility::fixedsizearray1<Real,6> metal_coord_ranges_for_nls; // metal coordinate ranges to constrain the NLS fitting

		if ( spinlabel_ ) {

			WeightCoordVector spinlabel_wghts_coords;

			// Note that we are using a simple distance calculation to the NBR_ATOM for spinlabel conformer filtering here.
			if ( spinlabel_->get_highres_conformer_filter_type() != NMRSpinlabel::DISTANCE ) {
				TR.Warning << "Computation of PCS score in PCSMultiSet with NMRDummySpinlabelEnsemble uses DISTANCE-based conformer filter. ";
				TR.Warning << "Setting the conformer filter type to anything else than DISTANCE will be ignored." << std::endl;
			}

			// Perform spinlabel clash filter
			runtime_assert_msg( spinlabel_->get_dummy_ensemble(), "ERROR while trying to calculate PCS score. NMRDummySpinlabelEnsemble not set. Check if this spinlabel type has database file in Rosetta database." );
			spinlabel_wghts_coords = spinlabel_->filter_spinlabel_ensemble_by_distance_check(pose, tag_residue_number_);

			// Average metal ion coordinates and bounding box to keep
			// track of the xyz range that the spinlabel samples
			Vector metal_coo_mean(0.0);
			Real totwght(0.0);
			numeric::geometry::BoundingBox<Vector> bbox;
			bbox.set_lower(Vector(std::numeric_limits< Real >::max()));
			bbox.set_upper(Vector(std::numeric_limits< Real >::min()));

			for ( Size k = 1, k_end = spinlabel_wghts_coords.size(); k <= k_end; ++k ) {
				totwght += spinlabel_wghts_coords[k].first;
				metal_coo_mean += spinlabel_wghts_coords[k].first * spinlabel_wghts_coords[k].second;
				bbox.add(spinlabel_wghts_coords[k].second);
			}
			metal_coo_mean /= totwght;
			metal_coord_ranges_for_nls[1] = bbox.lower().x();
			metal_coord_ranges_for_nls[2] = bbox.lower().y();
			metal_coord_ranges_for_nls[3] = bbox.lower().z();
			metal_coord_ranges_for_nls[4] = bbox.upper().x();
			metal_coord_ranges_for_nls[5] = bbox.upper().y();
			metal_coord_ranges_for_nls[6] = bbox.upper().z();

			for ( Size i = 1; i <= number_metal_ions_; ++i ) {
				pcs_singleset_vec_[i]->update_spin_coordinates( pose );
				pcs_singleset_vec_[i]->set_metal_coord_bounds(metal_coord_ranges_for_nls);
				if ( pcs_singleset_vec_[i]->get_computation_type() == PCSSingleSet::SVD ) {
					total_score += pcs_singleset_vec_[i]->solve_tensor_and_compute_score_by_svd(metal_coo_mean) * pcs_singleset_vec_[i]->get_weight();
				} else if ( pcs_singleset_vec_[i]->get_computation_type() == PCSSingleSet::NLS ) {
					total_score += pcs_singleset_vec_[i]->solve_tensor_and_compute_score_by_nls(metal_coo_mean) * pcs_singleset_vec_[i]->get_weight();
				}
				pcs_singleset_vec_[i]->set_atom_derivatives( pose );
				singleset_tensors_temp[i] = pcs_singleset_vec_[i]->get_tensor_const();
			}
		} else if ( gridsearch_iterator_ ) {

			// Some basic setup
			Real score_svd(0);
			Real best_score_svd(std::numeric_limits<Real>::max());

			Vector current_metal_coords;             // coordinates that get updated during the grid search
			Vector best_metal_coords;              // best coordinates obtained by grid search and SVD
			gridsearch_iterator_->set_grid_search_center( pose );       // update coordinates of the grid search center
			Vector grid_search_center(gridsearch_iterator_->get_grid_search_center());  // grid search center as initial guess for NLS fitting

			metal_coord_ranges_for_nls[1] = grid_search_center.x() - gridsearch_iterator_->get_grid_max_radius();
			metal_coord_ranges_for_nls[2] = grid_search_center.y() - gridsearch_iterator_->get_grid_max_radius();
			metal_coord_ranges_for_nls[3] = grid_search_center.z() - gridsearch_iterator_->get_grid_max_radius();
			metal_coord_ranges_for_nls[4] = grid_search_center.x() + gridsearch_iterator_->get_grid_max_radius();
			metal_coord_ranges_for_nls[5] = grid_search_center.y() + gridsearch_iterator_->get_grid_max_radius();
			metal_coord_ranges_for_nls[6] = grid_search_center.z() + gridsearch_iterator_->get_grid_max_radius();

			// Which datasets in pcs_singleset_vec_ for SVD and NLS calculation?
			// Find the corresponding indices and score calculation
			utility::vector1<Size> idx_pcs_sets_with_svd_calc;
			utility::vector1<Size> idx_pcs_sets_with_nls_calc;
			utility::vector1<Size>::const_iterator idx;

			for ( Size i = 1; i <= number_metal_ions_; ++i ) {
				// retrieve the spin coordinates from the pose to build matrix A
				pcs_singleset_vec_[i]->update_spin_coordinates( pose );
				if ( pcs_singleset_vec_[i]->get_computation_type() == PCSSingleSet::SVD ) {
					idx_pcs_sets_with_svd_calc.push_back(i);
				} else {
					idx_pcs_sets_with_nls_calc.push_back(i);
				}
			}

			// Grid search and SVD
			while ( gridsearch_iterator_->valid_next_grid_point(current_metal_coords) ) {
				score_svd = 0;
				for ( idx = idx_pcs_sets_with_svd_calc.begin(); idx != idx_pcs_sets_with_svd_calc.end(); ++idx ) {
					score_svd += pcs_singleset_vec_[*idx]->solve_tensor_and_compute_score_by_svd(current_metal_coords) * pcs_singleset_vec_[*idx]->get_weight();
					if ( score_svd > best_score_svd ) { // if one dataset gives already a score higher than the best score, there is no need to calculate the scores of the remaining datasets
						break;
					}
				}
				if ( score_svd < best_score_svd ) {
					best_score_svd = score_svd;
					best_metal_coords = current_metal_coords;
					gridsearch_iterator_->set_best_grid_point(best_metal_coords);
				}
			}

			// Set PCSTensors to best values as calculated by SVD and set atom derivatives
			for ( idx = idx_pcs_sets_with_svd_calc.begin(); idx != idx_pcs_sets_with_svd_calc.end(); ++idx ) {
				best_score_svd = pcs_singleset_vec_[*idx]->solve_tensor_and_compute_score_by_svd(best_metal_coords) * pcs_singleset_vec_[*idx]->get_weight();
				pcs_singleset_vec_[*idx]->set_atom_derivatives( pose );
				singleset_tensors_temp[*idx] = pcs_singleset_vec_[*idx]->get_tensor_const();
			}

			total_score += best_score_svd;

			// NLS fitting
			for ( idx = idx_pcs_sets_with_nls_calc.begin(); idx != idx_pcs_sets_with_nls_calc.end(); ++idx ) {
				pcs_singleset_vec_[*idx]->set_metal_coord_bounds(metal_coord_ranges_for_nls);
				total_score += pcs_singleset_vec_[*idx]->solve_tensor_and_compute_score_by_nls(grid_search_center) * pcs_singleset_vec_[*idx]->get_weight();
				pcs_singleset_vec_[*idx]->set_atom_derivatives( pose );
				singleset_tensors_temp[*idx] = pcs_singleset_vec_[*idx]->get_tensor_const();
			}

		} else {
			utility_exit_with_message("ERROR while trying to calculate PCSMultiSet score. Spinlabel and gridsearch iterator are not set.");
		}
	}
	if ( TR.Trace.visible() ) {
		TR.Trace << "PCS score for " << number_metal_ions_ << " experiment(s) with " << (spinlabel_ ? "spinlabel "+spinlabel_->get_code(): "gridsearch")
			<< " at position " << tag_residue_number_ << ": " << total_score << std::endl;
		if ( spinlabel_ ) {
			spinlabel_->show(TR);
		}
		for (  Size i = 1; i <= number_metal_ions_; ++i ) {
			TR.Trace << "Dataset: " << pcs_singleset_vec_[i]->get_dataset_name() << std::endl;
			if ( pcs_singleset_vec_[i]->get_computation_type() == PCSSingleSet::SVD && !fixed_tensors_ ) {
				singleset_tensors_temp[i]->show_tensor_stats(TR.Trace, false);
			} else {
				singleset_tensors_temp[i]->show_tensor_stats(TR.Trace, true);
			}
		}
	}
	singleset_tensors = singleset_tensors_temp;
	return total_score;
}

void
PCSMultiSet::show(std::ostream & tracer) const {
	auto sum_calc = [](Size const & a, PCSSingleSetOP const & b) { return a + b->get_number_pcs(); };
	Size total_pcs = std::accumulate(pcs_singleset_vec_.begin(), pcs_singleset_vec_.end(), 0, sum_calc);
	tracer << "   * * * PCSMultiSet Summary Report * * *   " << std::endl;
	tracer << "Spinlabel site: " << tag_residue_number_ << std::endl;
	tracer << "No experiments: " << number_metal_ions_ << std::endl;
	tracer << "No PCSs:        " << total_pcs << std::endl;
	tracer << "PCS Spinlabel:  ";
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
	for ( Size i = 1; i <= number_metal_ions_; ++i ) {
		tracer << " - " << pcs_singleset_vec_[i]->get_dataset_name() << ": [ " << pcs_singleset_vec_[i]->get_number_pcs() << " PCSs, "
			<< pcs_singleset_vec_[i]->get_metal_ion_label() << " ]" << std::endl;
	}
}

void
PCSMultiSet::deep_copy_pcs_singleset_vec(utility::vector1<PCSSingleSetOP> const & other_vec) {
	pcs_singleset_vec_.resize(other_vec.size());
	for ( Size i = 1; i <= pcs_singleset_vec_.size(); ++i ) {
		pcs_singleset_vec_[i] = PCSSingleSetOP( new PCSSingleSet( *(other_vec[i]) ) );
	}
}

/// @brief creates MTSL spinlabel as default
void
PCSMultiSet::set_default_spinlabel() {
	spinlabel_ = NMRSpinlabelOP( new NMRSpinlabel("fa_standard", "R1A") );
}

Size
PCSMultiSet::get_total_number_pcs() const {
	auto sum = [](Size const & a, PCSSingleSetOP const & b) { return a + b->get_number_pcs(); };
	return std::accumulate(pcs_singleset_vec_.begin(), pcs_singleset_vec_.end(), 0, sum);
}

} // namespace pcs
} // namespace nmr
} // namespace scoring
} // namespace core
