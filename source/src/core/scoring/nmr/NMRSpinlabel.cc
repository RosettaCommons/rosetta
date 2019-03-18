// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    core/scoring/nmr/NMRSpinlabel.cc
/// @brief   Implementation of class NMRSpinlabel
/// @details last Modified: 08/31/16
/// @author  Georg Kuenze (georg.kuenze@vanderbilt.edu)

// Unit headers
#include <core/scoring/nmr/NMRSpinlabel.hh>

// Package headers
#include <core/scoring/nmr/NMRDummySpinlabelEnsemble.hh>
#include <core/io/nmr/SpinlabelDatabaseHandler.hh>
#include <core/scoring/nmr/util.hh>

// Project headers
#include <core/types.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueType.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/rms_util.hh>
#include <numeric/agglomerative_hierarchical_clustering.hh>
#include <numeric/ClusteringTreeNode.hh>

// Basic headers
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/nmr.OptionKeys.gen.hh>

// Utility headers
#include <utility/vector1.hh>
#include <utility/exit.hh>
#include <utility/string_util.hh>

// C++ headers
#include <iostream>
#include <algorithm>
#include <string>
#include <numeric>

// Boost headers
#include <boost/algorithm/string.hpp>

namespace core {
namespace scoring {
namespace nmr {

static basic::Tracer TR( "core.scoring.nmr.NMRSpinlabel" );

/// @brief construct from strings of residue type set and residue type;
///        the other properties (e.g. radical atom) are looked up in the database
NMRSpinlabel::NMRSpinlabel(
	std::string const & residue_type_set,
	std::string const & residue_type
) :
	utility::pointer::ReferenceCount(),
	residue_type_( new chemical::ResidueType( *((chemical::ChemicalManager::get_instance()->residue_type_set(residue_type_set))->get_representative_type_name3(residue_type)) ) ),
	radical_atom_(""),
	distance_potential_histogram_file_(""),
	dummy_ensemble_(nullptr),
	weights_coordinates_table_(),
	highres_conformer_filter_type_(DISTANCE),
	boltzmann_factor_(2.0)
{
	using namespace core::io::nmr;

	if ( !residue_type_ ) {
		utility_exit_with_message("ERROR while creating NMRSpinlabel. Residue type " + residue_type + " not found.");
	}

	// Create dummy ensemble and set remaining data members
	SpinlabelDatabaseEntry const & spinlabel_data = SpinlabelDatabaseHandler::get_instance()->get_spinlabel_data( residue_type );
	dummy_ensemble_ = NMRDummySpinlabelEnsembleOP( new NMRDummySpinlabelEnsemble(spinlabel_data.ensemble_conformers(), *residue_type_) );
	radical_atom_ = spinlabel_data.radical_atom();
	distance_potential_histogram_file_ = spinlabel_data.distance_potential_histogram();
	init_radical_atom_weights_and_coordinates();
	init_spinlabel_name();
	register_options();
	init_from_cml();
}

NMRSpinlabel::NMRSpinlabel(
	chemical::ResidueTypeCOP residue_type,
	std::string const & radical_atom,
	NMRDummySpinlabelEnsembleCOP dummy_ensemble
) :
	utility::pointer::ReferenceCount(),
	residue_type_( new chemical::ResidueType( *residue_type ) ),
	radical_atom_(radical_atom),
	distance_potential_histogram_file_(""),
	dummy_ensemble_( new NMRDummySpinlabelEnsemble( *dummy_ensemble ) ),
	weights_coordinates_table_(),
	highres_conformer_filter_type_(DISTANCE),
	boltzmann_factor_(2.0)
{
	init_radical_atom_weights_and_coordinates();
	init_spinlabel_name();
	register_options();
	init_from_cml();
}

/// @brief copy constructor
NMRSpinlabel::NMRSpinlabel(NMRSpinlabel const & other) :
	utility::pointer::ReferenceCount(other),
	residue_type_( new chemical::ResidueType( *(other.residue_type_) )),
	radical_atom_(other.radical_atom_),
	distance_potential_histogram_file_(other.distance_potential_histogram_file_),
	dummy_ensemble_( other.dummy_ensemble_ ? new NMRDummySpinlabelEnsemble( *(other.dummy_ensemble_) ) : nullptr ),
	weights_coordinates_table_(other.weights_coordinates_table_),
	max_ensemble_size_(other.max_ensemble_size_),
	highres_conformer_filter_type_(other.highres_conformer_filter_type_),
	boltzmann_factor_(other.boltzmann_factor_)
{
	init_spinlabel_name();
}

/// @brief assignment operator
NMRSpinlabel & NMRSpinlabel::operator=(NMRSpinlabel const & rhs) {
	if ( this != &rhs ) {
		residue_type_ = chemical::ResidueTypeCOP( new chemical::ResidueType( *(rhs.residue_type_) ) );
		radical_atom_ = rhs.radical_atom_;
		distance_potential_histogram_file_ = rhs.distance_potential_histogram_file_;
		dummy_ensemble_ = rhs.dummy_ensemble_ ? NMRDummySpinlabelEnsembleOP( new NMRDummySpinlabelEnsemble( *(rhs.dummy_ensemble_) ) ) : nullptr;
		weights_coordinates_table_ = rhs.weights_coordinates_table_;
		max_ensemble_size_ = rhs.max_ensemble_size_;
		highres_conformer_filter_type_ = rhs.highres_conformer_filter_type_;
		boltzmann_factor_ = rhs.boltzmann_factor_;
		init_spinlabel_name();
	}
	return *this;
}

/// @brief destructor
NMRSpinlabel::~NMRSpinlabel() { }

void
NMRSpinlabel::init_spinlabel_name() {
	name_ = residue_type_->name();
	three_letter_code_ = residue_type_->name3();
}

void
NMRSpinlabel::init_radical_atom_weights_and_coordinates() {
	// Set initial weights and coordinates of radical atom
	utility::vector1< NMRDummySpinlabelConformerOP > const & ndsl_all_conformers = dummy_ensemble_->get_conformer_table();
	for ( Size i(1); i <= ndsl_all_conformers.size(); ++i ) {
		NMRDummySpinlabelAtomTable const & atom_table = ndsl_all_conformers[i]->get_atom_table();
		if ( atom_table.find(radical_atom_) != atom_table.end() ) {
			weights_coordinates_table_.push_back( std::make_pair( static_cast<Real>(ndsl_all_conformers[i]->get_nobs()),
				atom_table.at(radical_atom_).get_coordinates() ) );
		} else {
			utility_exit_with_message("Unable to find radical atom " + radical_atom_ + " in NMRDummySpinlabel representation");
		}
	}
}

void
NMRSpinlabel::show( std::ostream & tracer ) {
	tracer << " * * * NMRSpinlabel data * * * " << std::endl;
	tracer << "Full name        : " << name_ << std::endl;
	tracer << "Three-letter code: " << three_letter_code_ << std::endl;
	tracer << "Radical ion/atom : " << radical_atom_ << std::endl;
	tracer << "Ensemble size    : " << weights_coordinates_table_.size() << std::endl;
}

void
NMRSpinlabel::set_weights_and_coordinates(WeightCoordVector const & weights_coords) {
	weights_coordinates_table_.clear();
	weights_coordinates_table_ = weights_coords;
}

/// @brief filter dummy spinlabel ensemble given the neighborhood
///        of a particular target residue in the pose. Return a vector of
///        each spinlabel conformer's weight and radical atom coordinates.
///        Performs also clustering of coordinates internally such that the
///        vector size does not exceed the maximal number of spinlabel conformers.
NMRSpinlabel::WeightCoordVector
NMRSpinlabel::filter_spinlabel_ensemble_by_distance_check(
	pose::Pose const & pose,
	Size const & target_resid
)
{
	auto sl_sum = []( Size const previous, NMRDummySpinlabelConformerOP current)
		{ return previous + Size(!current->has_clash()); };

	// Number of NMRDummySpinlabelConformers before filter
	TR.Debug << "Filter NMRDummySpinlabelConformers on target site " << target_resid << " by their neighbor count" << std::endl;
	if ( TR.Debug.visible() ) {
		Size nonclash = std::accumulate(dummy_ensemble_->get_conformer_table().begin(),
			dummy_ensemble_->get_conformer_table().end(), 0, sl_sum);
		TR.Debug << "# Non-clashing conformers before filter: " << nonclash << " / " << dummy_ensemble_->get_ensemble_size() << std::endl;
	}

	// Filter dummy spinlabel with voxelgrid
	Real radius(residue_type_->nbr_radius() > 12.0 ? residue_type_->nbr_radius() : 12.0);
	runtime_assert_msg(dummy_ensemble_, "ERROR: NMRSpinlabel dummy ensemble not set.");
	dummy_ensemble_->clash_check(pose, target_resid, radius);

	// Number of NMRDummySpinlabelConformers after filter
	if ( TR.Debug.visible() ) {
		Size nonclash = std::accumulate(dummy_ensemble_->get_conformer_table().begin(),
			dummy_ensemble_->get_conformer_table().end(), 0, sl_sum);
		TR.Debug << "# Non-clashing conformers after filter: " << nonclash << " / " << dummy_ensemble_->get_ensemble_size() << std::endl;
	}

	WeightCoordVector radical_atom_wghts_coords;
	numeric::HomogeneousTransform<Real> Q(dummy_ensemble_->coordinate_transform_onto_target_site(pose, target_resid));

	// Get non-clashing dummy spinlabel conformers and pull out their elements from the RMSD matrix
	utility::vector1< NMRDummySpinlabelConformerOP > & ndsl_all_conformers = dummy_ensemble_->get_conformer_table();
	utility::vector1<utility::vector1<Real>> const & rmsd_mat = dummy_ensemble_->get_rmsd_mat();
	runtime_assert_msg(ndsl_all_conformers.size()==rmsd_mat.size(), "Vector of NMRDummySpinlabelConformers and RMSD matrix have different length.");
	utility::vector1<utility::vector1<Real>> reduced_rmsd_mat;
	reduced_rmsd_mat.reserve(rmsd_mat.size());

	// Add RMSD values for pairs of non-clashing spinlabel conformers into reduced RMSD matrix
	for ( Size i(1); i <= ndsl_all_conformers.size(); ++i ) {
		if ( !(ndsl_all_conformers[i]->has_clash()) ) {
			NMRDummySpinlabelAtomTable const & atom_table = ndsl_all_conformers[i]->get_atom_table();
			if ( atom_table.find(radical_atom_) != atom_table.end() ) {
				radical_atom_wghts_coords.push_back( std::make_pair( static_cast<Real>(ndsl_all_conformers[i]->get_nobs()), // cast number of obs to weight
					Q * atom_table.at(radical_atom_).get_coordinates() ) );
			} else {
				utility_exit_with_message("Unable to find radical atom " + radical_atom_ + " in NMRDummySpinlabel representation");
			}
			utility::vector1<Real> rmsd_row;
			rmsd_row.reserve(rmsd_mat[i].size());
			for ( Size j(1); j <= ndsl_all_conformers.size(); ++j ) {
				if ( !(ndsl_all_conformers[j]->has_clash()) ) {
					rmsd_row.push_back(rmsd_mat[i][j]);
				}
			}
			reduced_rmsd_mat.push_back(rmsd_row);
		}
	}

	// If number of dummy spinlabels exceeds the maximum number do clustering
	if ( radical_atom_wghts_coords.size() > max_ensemble_size_ ) {
		cluster_conformers_and_set_weights_and_coordinates(radical_atom_wghts_coords, reduced_rmsd_mat);
	}
	set_weights_and_coordinates(radical_atom_wghts_coords);
	return radical_atom_wghts_coords;
}

/// @brief filter dummy ensemble and keep those spinlabels provided in boolean mask
NMRSpinlabel::WeightCoordVector
NMRSpinlabel::filter_spinlabel_ensemble_by_mask(
	pose::Pose const & pose,
	Size const & target_resid,
	utility::vector1<bool> const & mask,
	utility::vector1<Real> const & scores
)
{
	WeightCoordVector radical_atom_wghts_coords;
	numeric::HomogeneousTransform<Real> Q(dummy_ensemble_->coordinate_transform_onto_target_site(pose, target_resid));

	if ( TR.Debug.visible() ) {
		Size nonclash(std::count(mask.begin(),mask.end(),true));
		TR.Debug << "Filter NMRDummySpinlabelConformers on target site " << target_resid << " with boolean mask" << std::endl;
		TR.Debug << "# Non-clashing conformers: " << nonclash << " / " << dummy_ensemble_->get_ensemble_size() << std::endl;
	}

	// Get dummy spinlabels and pull out their elements from the RMSD matrix
	utility::vector1< NMRDummySpinlabelConformerOP > & ndsl_all_conformers = dummy_ensemble_->get_conformer_table();
	utility::vector1<utility::vector1<Real>> const & rmsd_mat = dummy_ensemble_->get_rmsd_mat();
	runtime_assert_msg(ndsl_all_conformers.size()==rmsd_mat.size(), "Vector of NMRDummySpinlabelConformers and RMSD matrix have different length.");
	utility::vector1<utility::vector1<Real>> reduced_rmsd_mat;
	reduced_rmsd_mat.reserve(rmsd_mat.size());
	Real minscore(*std::min_element(scores.begin(), scores.end()));

	runtime_assert_msg(ndsl_all_conformers.size()==mask.size(), "Vector of NMRDummySpinlabelConformers and filter mask have different length.");
	runtime_assert_msg(ndsl_all_conformers.size()==scores.size(), "Vector of NMRDummySpinlabelConformers and scores have different length.");

	// Iterate over positions in mask. If mask[i] is true, we keep this spinlabel conformer.
	for ( Size i(1); i <= mask.size(); ++i ) {
		// Assign clash score and do Boltzmann weighting
		ndsl_all_conformers[i]->clash_score() = scores[i];
		Real wght(exp((minscore-scores[i])/boltzmann_factor_));
		ndsl_all_conformers[i]->frequency() = wght;
		if ( mask[i] ) { // keep spinlabel in ensemble
			ndsl_all_conformers[i]->clash_off();
			NMRDummySpinlabelAtomTable const & atom_table = ndsl_all_conformers[i]->get_atom_table();
			if ( atom_table.find(radical_atom_) != atom_table.end() ) {
				radical_atom_wghts_coords.push_back( std::make_pair( wght, Q * atom_table.at(radical_atom_).get_coordinates() ) );
			} else {
				utility_exit_with_message("Unable to find radical atom " + radical_atom_ + " in NMRDummySpinlabel representation");
			}
			utility::vector1<Real> rmsd_row;
			rmsd_row.reserve(rmsd_mat[i].size());
			for ( Size j(1); j <= mask.size(); ++j ) {
				if ( !(ndsl_all_conformers[j]->has_clash()) ) {
					rmsd_row.push_back(rmsd_mat[i][j]);
				}
			}
			reduced_rmsd_mat.push_back(rmsd_row);

		} else { // mark spinlabel as clashing
			ndsl_all_conformers[i]->clash_on();
		}
	}

	// If number of dummy spinlabels exceeds the maximum number do clustering
	if ( radical_atom_wghts_coords.size() > max_ensemble_size_ ) {
		cluster_conformers_and_set_weights_and_coordinates(radical_atom_wghts_coords, reduced_rmsd_mat);
	}
	set_weights_and_coordinates(radical_atom_wghts_coords);
	return radical_atom_wghts_coords;
}

void
NMRSpinlabel::cluster_conformers_and_set_weights_and_coordinates(
	WeightCoordVector & weights_coords,
	utility::vector1<utility::vector1<Real>> & rmsd_mat
)
{
	// Sanity check of input vectors
	debug_assert(weights_coords.size()==rmsd_mat.size());

	// Do average linkage agglomerative clustering with input RMSD matrix
	TR.Debug << "Running agglomerative hierarchical clustering of NMRSpinlabelEnsemble" << std::endl;
	numeric::AverageLinkClusterer alc;
	utility::vector1<numeric::ClusteringTreeNodeOP> clusters = alc.cluster(rmsd_mat, max_ensemble_size_);
	runtime_assert(max_ensemble_size_ == clusters.size());
	TR.Debug << "Found " << clusters.size() << " clusters" << std::endl;
	WeightCoordVector weights_coords_out;
	weights_coords_out.reserve(max_ensemble_size_);

	// Find the size and representative of each cluster.
	// As representative we pick the index of the row of the RMSD matrix
	// that gives the lowest average RMSD to all other cluster members
	utility::vector1<Size> idx_in(rmsd_mat.size());
	std::iota(idx_in.begin(), idx_in.end(), 1);
	Size ncluster(1);
	for ( auto cluster : clusters ) {
		utility::vector1<Size> idx_out;
		numeric::get_cluster_data(idx_in, cluster, idx_out);
		Size center = idx_out[1];
		Real weightsum(0);
		Real minrmsd = std::numeric_limits<Real>::max();
		for ( Size i(1); i<= idx_out.size(); ++i ) {
			Real rmsdsum(0);
			for ( Size j(1); j <= idx_out.size(); ++j ) {
				rmsdsum+=rmsd_mat[ idx_out[i] ][ idx_out[j] ];
			}
			if ( rmsdsum < minrmsd ) {
				center=idx_out[i];
				minrmsd=rmsdsum;
			}
			// Sum of spinlabel weights
			weightsum+=weights_coords[ idx_out[i] ].first;
		}
		weights_coords_out.push_back(std::make_pair(weightsum, weights_coords[center].second));
		TR.Debug << "Cluster: " << ncluster << ", Size: " << idx_out.size() << ", Members: [ " << utility::join(idx_out,std::string(", "))
			<< " ], Center: " << center << ", Weight: " << weightsum << std::endl;
		ncluster++;
	}
	weights_coords=weights_coords_out;
}

/// @brief register command line options
void
NMRSpinlabel::register_options() {
	using namespace basic::options;
	option.add_relevant(OptionKeys::nmr::spinlabel::max_ensemble_size);
	option.add_relevant(OptionKeys::nmr::spinlabel::highres_conformer_filter_type);
	option.add_relevant(OptionKeys::nmr::spinlabel::boltzmann_kt);
}

void
NMRSpinlabel::init_from_cml() {
	using namespace basic::options;
	max_ensemble_size_ = option[ basic::options::OptionKeys::nmr::spinlabel::max_ensemble_size ]();
	convert_string_to_conformer_filter_type(option[ basic::options::OptionKeys::nmr::spinlabel::highres_conformer_filter_type ]());
	boltzmann_factor_ = option[ basic::options::OptionKeys::nmr::spinlabel::boltzmann_kt ]();
}

void
NMRSpinlabel::convert_string_to_conformer_filter_type(std::string const & filter_type) {
	std::string str = boost::to_upper_copy<std::string>(filter_type);
	if ( str == "DISTANCE" ) {
		highres_conformer_filter_type_ = NMRSpinlabel::DISTANCE;
	} else if ( str == "BUMP_ENERGY" ) {
		highres_conformer_filter_type_ = NMRSpinlabel::BUMP_ENERGY;
	} else {
		utility_exit_with_message( "ERROR: Provided string does not match NMRSpinlabel conformer filter type. Possible options are \"DISTANCE\" and \"BUMP_ENERGY\".");
	}
}

void
NMRSpinlabel::set_highres_conformer_filter_type(std::string const & filter_type) {
	convert_string_to_conformer_filter_type(filter_type);
}

} // namespace nmr
} // namespace scoring
} // namespace core
