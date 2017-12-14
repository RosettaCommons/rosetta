// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/energy_based_clustering/EnergyBasedClusteringOptions.cc
/// @brief A container for the options used by the EnergyBasedClusteringProtocol.
/// @author Vikram K. Mulligan (vmullig@u.washington.edu)

// Header
#include <protocols/energy_based_clustering/EnergyBasedClusteringOptions.hh>

// Basic includes
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/cluster.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>

static basic::Tracer TR( "protocols.cluster.energy_based_clustering.EnergyBasedClusteringOptions" );


namespace protocols {
namespace energy_based_clustering {

EnergyBasedClusteringOptions::EnergyBasedClusteringOptions( bool const initialize_from_options ):
	prerelax_(false),
	relax_rounds_(1),
	cluster_by_( EBC_bb_cartesian ),
	use_CB_(false),
	cluster_radius_(1.0),
	//weight_by_energy_(true),
	//kbt_(0.62),
	residues_to_ignore_(),
	chains_to_ignore_(),
	limit_structures_per_cluster_(0),
	limit_clusters_(0),
	cyclic_(false),
	cyclic_symmetry_(0),
	cyclic_symmetry_mirroring_(false),
	cyclic_symmetry_threshold_(10.0),
	cyclic_symmetry_threshold_specified_(false),
	cluster_cyclic_permutations_(false),
	cyclic_permutation_offset_(1),
	mutate_to_ala_(false),
	disulfide_positions_(),
	homooligomer_swap_(false),
	silent_output_(false),
	cst_files_(),
	extra_rms_atoms_(),
	rebuild_all_in_dihedral_mode_(true),
	output_prefix_("")
{
	if ( initialize_from_options ) initialize_from_global_options();
}

EnergyBasedClusteringOptions::~EnergyBasedClusteringOptions() = default;

EnergyBasedClusteringOptionsOP
EnergyBasedClusteringOptions::clone() const {
	return EnergyBasedClusteringOptionsOP( new EnergyBasedClusteringOptions( *this ) );
}

/// @brief Initialize this option from the global options system.
/// @details Called by default constructor.
void
EnergyBasedClusteringOptions::initialize_from_global_options() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys::cluster::energy_based_clustering;

	prerelax_ = option[ prerelax ]();

	runtime_assert_string_msg( option[relax_rounds]() > 0, "Error in global options: the -cluster:energy_based_clustering:relax_rounds option must have a set value greater than zero." );
	relax_rounds_ = static_cast<core::Size>( option[relax_rounds]() );

	std::string const clusterby( option[cluster_by]() );
	if ( !clusterby.compare("bb_cartesian") ) {
		cluster_by_ = EBC_bb_cartesian;
	} else if ( !clusterby.compare("bb_dihedral") ) {
		cluster_by_ = EBC_bb_dihedral;
	} else {
		utility_exit_with_message( "Error in global options: the -cluster:energy_based_clustering:cluster_by option was set to \"" + clusterby + "\", which I can't parse." );
	}

	use_CB_ = option[use_CB]();

	cluster_radius_ = option[cluster_radius]();
	runtime_assert_string_msg( cluster_radius_ > 0, "Error in global options: the -cluster:energy_based_clustering:cluster_radius option must be set to a value greater than zero." );

	//weight_by_energy_ = option[weight_by_energy]();

	//kbt_ = option[kbt]();
	//runtime_assert_string_msg( kbt_ > 0, "Error in global options: the -cluster:energy_based_clustering:kbt option must be set to a value greater than zero." );

	if ( option[residues_to_ignore].user() ) {
		residues_to_ignore_.resize( option[residues_to_ignore]().size() );
		for ( core::Size i(1), imax(option[residues_to_ignore]().size()); i<=imax; ++i ) {
			runtime_assert_string_msg( option[residues_to_ignore]()[i] > 0, "Error in global options: all residues specified with the -cluster:energy_based_clustering:residues_to_ignore option must have indices greater than zero." );
			residues_to_ignore_[i] = static_cast<core::Size>( option[residues_to_ignore]()[i] );
		}
	} else {
		residues_to_ignore_.clear();
	}

	if ( option[chains_to_ignore].user() ) {
		chains_to_ignore_.resize( option[chains_to_ignore]().size() );
		for ( core::Size i(1), imax(option[chains_to_ignore]().size()); i<=imax; ++i ) {
			runtime_assert_string_msg( option[chains_to_ignore]()[i] > 0, "Error in global options: all chains specified with the -cluster:energy_based_clustering:chains_to_ignore option must have indices greater than zero." );
			chains_to_ignore_[i] = static_cast<core::Size>( option[chains_to_ignore]()[i] );
		}
	} else {
		chains_to_ignore_.clear();
	}

	runtime_assert_string_msg( option[limit_structures_per_cluster]() >= 0,
		"Error in global options: the -cluster:energy_based_clustering:limit_structures_per_cluster option cannot be set to a negative value." );
	limit_structures_per_cluster_ = static_cast<core::Size>( option[limit_structures_per_cluster]() );

	runtime_assert_string_msg( option[limit_clusters]() >= 0,
		"Error in global options: the -cluster:energy_based_clustering:limit_clusters option cannot be set to a negative value." );
	limit_clusters_ = static_cast<core::Size>( option[limit_clusters]() );

	cyclic_ = option[cyclic]();

	runtime_assert_string_msg( option[cyclic_symmetry]() >= 0,
		"Error in global options: the -cluster:energy_based_clustering:cyclic_symmetry option cannot be set to a negative value." );
	cyclic_symmetry_ = static_cast<core::Size>( option[cyclic_symmetry]() );

	cyclic_symmetry_mirroring_ = option[cyclic_symmetry_mirroring]();

	cyclic_symmetry_threshold_ = option[cyclic_symmetry_threshold]();
	runtime_assert_string_msg( cyclic_symmetry_threshold_ >= 0.0, "Error in global options: the -cluster:energy_based_clustering:cyclic_symmetry_threshold option cannot be set to a negative value." );
	cyclic_symmetry_threshold_specified_ = option[cyclic_symmetry_threshold].user();

	cluster_cyclic_permutations_ = option[cluster_cyclic_permutations]();

	runtime_assert_string_msg( option[cyclic_permutation_offset] > 0, "Error in global options: the -cluster:energy_based_clustering:cyclic_permutation_offset option must be set to a value greater than zero." );
	cyclic_permutation_offset_ = static_cast<core::Size>(option[cyclic_permutation_offset]());

	mutate_to_ala_ = option[mutate_to_ala]();

	if ( option[disulfide_positions].user() ) {
		disulfide_positions_.resize( option[disulfide_positions]().size() );
		runtime_assert_string_msg( disulfide_positions_.size() % 2 == 0, "Error in global options: the number of positions specified with the -cluster:energy_based_clustering:disulfide_positions flag must be even." );
		for ( core::Size i(1), imax(disulfide_positions_.size()); i<=imax; ++i ) {
			runtime_assert_string_msg( option[disulfide_positions]()[i] > 0, "Error in global options: each position specified with the -cluster:energy_based_clustering:disulfide_positions flag must have an index greater than zero." );
			disulfide_positions_[i] = option[disulfide_positions]()[i];
		}
	} else {
		disulfide_positions_.clear();
	}

	homooligomer_swap_ = option[homooligomer_swap]();

	silent_output_ = option[silent_output]();

	if ( option[cst_file].user() ) {
		cst_files_ = option[cst_file]();
	} else {
		cst_files_.clear();
	}

	if ( option[extra_rms_atoms].user() ) {
		extra_rms_atoms_ = option[extra_rms_atoms]();
	} else {
		extra_rms_atoms_.clear();
	}

	rebuild_all_in_dihedral_mode_ = option[ rebuild_all_in_dihedral_mode ]();

	if ( option[ basic::options::OptionKeys::out::prefix ].user() ) {
		output_prefix_ = option[ basic::options::OptionKeys::out::prefix ]();
	} else {
		output_prefix_ = "";
	}
}

} //energy_based_clustering
} //protocols
