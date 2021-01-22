// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/energy_methods/DEEREnergy.cc
/// @brief  Score term for data obtained with double electron-electron resonance (DEER)
/// @details This method is initiated by assigning a weight to deer_decay and providing an input
///       text file via the option -epr_deer:input_files data.txt.
///
/// @author  Diego del Alamo ( del.alamo@vanderbilt.edu )

// Input headers
#include <basic/options/option.hh>
#include <basic/options/keys/epr_deer.OptionKeys.gen.hh>

// Unit headers
#include <core/energy_methods/DEEREnergy.hh>
#include <core/energy_methods/DEEREnergyCreator.hh>

// Unit headers
#include <core/scoring/epr_deer/metrics/DEERData.hh>
#include <core/scoring/epr_deer/DEERDataCache.hh>
#include <core/scoring/epr_deer/EPRSpinLabel.hh>
#include <core/scoring/epr_deer/util.hh>

// Package headers
#include <core/chemical/AA.hh>
#include <core/pose/Pose.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <core/pose/PDBInfo.hh>
#include <core/scoring/methods/EnergyMethodCreator.hh>
#include <core/scoring/methods/EnergyMethod.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/AtomType.hh>
#include <core/select/util.hh>
#include <core/scoring/AtomVDW.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/methods/Methods.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/DenseEnergyContainer.hh>
#include <core/scoring/LREnergyContainer.hh>
#include <basic/datacache/BasicDataCache.hh>
#include <core/types.hh>

// Basic headers
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/database/open.hh>

// Utility headers
#include <utility/exit.hh>
#include <utility/vector1.hh>

// Numeric headers
#include <numeric/constants.hh>
#include <numeric/xyzVector.hh>
#include <numeric/HomogeneousTransform.hh>

// ObjexxFCL headers
// Objexx headers
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/FArray1D.hh>

// C++ headers
#include <string>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <set>
#include <stdlib.h>
#include <algorithm>
#include <time.h>
#include <limits>

namespace core {
namespace energy_methods {

/// @brief Tracer used for error messages
/// @details Global to avoid re-instantiating tracer with every new object
static basic::Tracer TR( "core.energy_methods.DEEREnergy" );

/// @brief  Energy Method creator required by framework
/// @param options: Energy method options passed to all energy methods
/// @return Pointer to new method
core::scoring::methods::EnergyMethodOP
DEEREnergyCreator::create_energy_method(
	scoring::methods::EnergyMethodOptions const & options
) const {
	return core::scoring::methods::EnergyMethodOP( new DEEREnergy( options ) );
}

/// @brief  Energy Method descriptor to identify which weights it uses
/// @return Vector of scoretype used by this method
scoring::ScoreTypes
DEEREnergyCreator::score_types_for_method() const {
	return { scoring::epr_deer_score };
}

/// @brief  Constructor
/// @param  energymethodcreator: Creator object
DEEREnergy::DEEREnergy(
	scoring::methods::EnergyMethodCreatorOP energymethodcreator
) : scoring::methods::ContextDependentLRTwoBodyEnergy( energymethodcreator ) {}

/// @brief  Constructor for obtaining command-line options
DEEREnergy::DEEREnergy(
	scoring::methods::EnergyMethodOptions const &
) :
	scoring::methods::ContextDependentLRTwoBodyEnergy(
	scoring::methods::EnergyMethodCreatorOP( new DEEREnergyCreator ) )
{}

/// @brief  Default Constructor
DEEREnergy::DEEREnergy() : scoring::methods::ContextDependentLRTwoBodyEnergy(
	scoring::methods::EnergyMethodCreatorOP( new DEEREnergyCreator ) ) {}

/// @brief   Copy constructor
DEEREnergy::DEEREnergy(
	DEEREnergy const & other
) : ContextDependentLRTwoBodyEnergy( other ) {}

/// @brief  Destructor
DEEREnergy::~DEEREnergy() {}

/// @brief Copy function
/// @return Pointer to new object that is a copy of this object
core::scoring::methods::EnergyMethodOP
DEEREnergy::clone() const {
	return scoring::methods::EnergyMethodOP( new DEEREnergy( *this ) );
}

///////////////////////////
/// INHERITED FUNCTIONS

/// @brief Inherited function specifying where the data is getting stored
/// @return Descriptor of this energy method
scoring::methods::LongRangeEnergyType
DEEREnergy::long_range_type() const {
	return scoring::methods::epr_deer_lr;
}

/// @brief Set up scoring process
/// @param pose: Pose being scored
void
DEEREnergy::setup_for_scoring(
	pose::Pose & pose,
	scoring::ScoreFunction const &
) const {

	// Cursory check that everything is in order
	initialize_energy_method( pose );

	// For clash evaluation
	pose.update_residue_neighbors();

	// The actual scoring proceure
	iter_over_labels( pose );
}

/// @brief Add energy/score for a given pair of residues to map.
/// @param rsd1: First residue
/// @param rsd2: Second residue
/// @param pose: Pose being scored
/// @param emap: Where we are adding the score
/// @detail Because some sets of data may be between more than two
/// @detail residues, this function looks at a map kept in
/// @detail DEERDataCache specifically maintained for this purpose
void
DEEREnergy::residue_pair_energy(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	pose::Pose const & pose,
	scoring::ScoreFunction const &,
	scoring::EnergyMap & emap
) const {

	// Get the list of data
	auto const & edges =
		get_const_deer_data( pose )->edge( rsd1.seqpos(), rsd2.seqpos() );

	// If this object is empty, no need to iterate through it
	if ( edges.size() == 0 ) {
		return;
	}
	for ( auto const & i_w : edges ) {
		emap[ scoring::ScoreType::epr_deer_score ] +=
			get_const_deer_data( pose )->at( i_w.first )->score() * i_w.second;
	}
}

/// @brief  A declared but unused function
void
DEEREnergy::setup_for_minimizing(
	pose::Pose &,
	scoring::ScoreFunction const &,
	kinematics::MinimizerMapBase const &
) const {}

/// @brief  Returns "false" by default
/// @return Returns false (DEER cannot evaluate anything intra-residue)
bool
DEEREnergy::defines_intrares_energy(
	scoring::EnergyMap const &
) const {
	return false;
}

/// @brief Fetches if two residues are "connected" (i.e. there is a score)
/// @param pose: Pose being scored
/// @param rsd1: First residue
/// @param rsd2: Second residue
/// @return Whether the pair of residues are linked by data
bool
DEEREnergy::defines_residue_pair_energy(
	pose::Pose const & pose,
	Size rsd1,
	Size rsd2
) const {
	return get_const_deer_data( pose )->edge( rsd1, rsd2 ).size() > 0;
}

/// @brief  Calculates the derivatives for minimization. Done numerically
/// @param pose: Pose being scored/modified
void
DEEREnergy::setup_for_derivatives(
	pose::Pose & pose,
	scoring::ScoreFunction const &
) const {

	// For clash evaluation
	pose.update_residue_neighbors();

	// Instantiate the maps where these terms are stored
	std::map< Size, numeric::xyzVector< Real > > f1_forces, f2_forces;

	// Store the scores for forward and reverse movement
	std::map< int, std::map< Size, Real > > total_scores;
	auto deer_data = get_deer_data( pose );
	total_scores[ -1 ] = iter_over_labels( pose, -1 );
	total_scores[  1 ] = iter_over_labels( pose,  1 );

	TR << "Starting derivative calculation" << std::endl;

	// Iterate over all data
	for ( Size i = 1; i <= get_const_deer_data( pose )->size(); ++i ) {
		auto const & residues = get_const_deer_data( pose )->at( i )->residues();

		// Movement magnitude
		Real mag = total_scores[ -1 ][ i ] - total_scores[ 1 ][ i ];

		// Iterate over all pairs of residues in dataset
		for ( Size i = 1; i < residues.size(); ++i ) {

			// We will be calculating derivaties for this one first
			auto const & res = residues[ i ].first;

			// Get the XYZ of the CA atom
			auto const & res1_ca = pose.residue( res ).xyz( "CA" );

			for ( Size j = i + 1; j <= residues.size(); ++j ) {

				// Get the XYZ of the CA atom
				auto const & res2_ca = pose.residue( residues[ j ].first ).xyz( "CA" );

				// Add the values to those vectors
				scoring::epr_deer::add_to_map( f2_forces, res,
					( res1_ca - res2_ca ).normalize() * mag );
				scoring::epr_deer::add_to_map( f1_forces, res,
					cross( res1_ca, res2_ca ).normalize() * mag / 2 );
			}
		}
	}
	// Now save these values to the DEERDataCache object
	for ( auto const & res_vec : f1_forces ) {
		get_deer_data( pose )->f1_force( res_vec.first, res_vec.second );
	}
	for ( auto const & res_vec : f2_forces ) {
		get_deer_data( pose )->f2_force( res_vec.first, res_vec.second );
	}
}

/// @brief  A declared but unused function
void
DEEREnergy::eval_intrares_energy(
	conformation::Residue const &,
	pose::Pose const &,
	scoring::ScoreFunction const &,
	scoring::EnergyMap &
) const {}

/// @brief  A declared but unused function
void
DEEREnergy::finalize_total_energy(
	pose::Pose &,
	scoring::ScoreFunction const &,
	scoring::EnergyMap &
) const {}

/// @brief  Apply atom derivative to atom
/// @param id: Atom ID in pose
/// @param pose: Pose being modified
/// @param emap: Energy map with scores used by derivatives
/// @param F1 vector to modify for atom
/// @param F2 vector to modify for atom
void
DEEREnergy::eval_atom_derivative(
	id::AtomID const & id,
	pose::Pose const & pose,
	kinematics::DomainMap const &,
	scoring::ScoreFunction const & /*sfxn*/,
	scoring::EnergyMap const & emap,
	numeric::xyzVector< Real > & F1,
	numeric::xyzVector< Real > & F2
) const {

	// Check if backbone atom
	if ( id.atomno() > 3 ) {
		return;
	}

	// Add F1 force
	if ( !get_const_deer_data( pose )->has_f1_force( id.rsd() ) ) {
		F1 += emap[ scoring::ScoreType::epr_deer_score ]
			* get_const_deer_data( pose )->f1_force( id.rsd() );
	}

	// Add F2 force
	if ( !get_const_deer_data( pose )->has_f2_force( id.rsd() ) ) {
		F2 += emap[ scoring::ScoreType::epr_deer_score ]
			* get_const_deer_data( pose )->f2_force( id.rsd() );
	}
}

/// @brief  A declared but unused function
void
DEEREnergy::indicate_required_context_graphs(
	utility::vector1< bool > &
) const {}

/// @brief  Version 1 (as of 3 January 2021)
/// @return Version 1
Size
DEEREnergy::version() const {
	return 1;
}

////////////////////////////////////////////
/// NON-VIRTUAL FUNCTIONS

/// @brief  Returns the DEERDataCache in the pose (const)
/// @param pose: Pose with DEERDataCache stashed in it
/// @return Pointer to const DEERDataCache object
/// @detail This will crash if the DEERDataCache object is absent!
scoring::epr_deer::DEERDataCacheCOP const
DEEREnergy::get_const_deer_data(
	pose::Pose const & pose
) const {
	return utility::pointer::dynamic_pointer_cast<
		const scoring::epr_deer::DEERDataCache >( pose.data().get_const_ptr(
		pose::datacache::CacheableDataType::EPR_DEER_DATA ) );
}

/// @brief  Returns the DEERDataCache in the pose (non-const)
/// @param pose: Pose with DEERDataCache stashed in it
/// @return Pointer to non-const DEERDataCache object
/// @detail This will crash if the DEERDataCache object is absent!
scoring::epr_deer::DEERDataCacheOP
DEEREnergy::get_deer_data(
	pose::Pose & pose
) const {
	return utility::pointer::dynamic_pointer_cast<
		scoring::epr_deer::DEERDataCache >( pose.data().get_ptr(
		pose::datacache::CacheableDataType::EPR_DEER_DATA ) );
}


/// @brief This checks if multiple sets of spin labels is being used
/// @param pose: Pose being scored
/// @param mod: Modifier to X-axis (for calculating derivatives)
/// @return Map of edges/data indices to scores for those edges/indices
/// @detail This is used when spin labels are obtained for specific
/// @detail proteins a priori, e.g., by triangulation/multilateration,
/// @detail and the solution is underdetermined, so multiple solutions
/// @detail are required.
std::map< Size, Real >
DEEREnergy::iter_over_labels(
	pose::Pose & pose,
	int const & mod // = 0
) const {

	// DEER data to modify
	auto deer_data = get_deer_data( pose );

	// Output
	std::map< Size, Real > total_scores;

	// Code in case we are using multiple custom spin labels
	if ( get_const_deer_data( pose )->labels().size() > 0 ) {

		// First the scores for each set are calculated and added
		for ( core::Size i = 1; i <= deer_data->labels().size(); ++i ) {

			// Weight and details of the spin label
			Real const w = deer_data->sl_weights()[ i ];
			scoring::epr_deer::EPRSpinLabel sl = deer_data->labels()[ i ];

			// Now score (increment score if previously recorded)
			for ( core::Size j = 1; j <= deer_data->size(); ++j ) {
				scoring::epr_deer::add_to_map(
					total_scores, j, w * get_score( pose, sl, j, mod ) );
			}
		}

		// Store score in datacache for retrieval by residue_pair_energy()
		for ( core::Size j = 1; j <= deer_data->size(); ++j ) {
			( *deer_data )[ j ]->score( total_scores[ j ] );
		}

		// Default method if no special custom spin labels are used
	} else {
		scoring::epr_deer::EPRSpinLabel sl;
		for ( core::Size i = 1; i <= deer_data->size(); ++i ) {
			total_scores[ i ] = get_score( pose, sl, i, mod );
		}
	}

	// Return
	return total_scores;
}

/// @brief Calculates the score for a specific data set
/// @param pose: Pose being scored
/// @param sl: Spin label being used to obtain score
/// @param edge_id: Index of data being scored (in DEERDataCache)
/// @param mod: Modifier to X-axis (zero when scored)
/// @return Score corresponding to the data
Real
DEEREnergy::get_score(
	pose::Pose & pose,
	scoring::epr_deer::EPRSpinLabel & sl,
	Size const & edge_id,
	int const & mod // = 0.0
) const {

	// Get data
	scoring::epr_deer::DEERDataCacheOP deer_data = get_deer_data( pose );

	// Get data for this index (metrics::DEERDataOP object)
	auto & edge = ( *deer_data )[ edge_id ];

	// Calculate DEER distance distribuion between simulated spin labels
	sl.label_everything( pose, edge->residues(), false );
	auto const histogram = sl.histogram( edge->residues(), edge->bins_per_a(),
		mod, edge->stdev(), edge->dist_map() );

	// This is just in case the data needs to be printed
	// Probably, a future version of this method can remove this
	if ( pose.pdb_info() ) {
		edge->print_histogram( histogram, pose.pdb_info()->name() );
	} else {
		edge->print_histogram( histogram );
	}

	// Pass distribution to scoring approach (encoded by DEERDataOP)
	return edge->get_score( histogram, ( mod == 0 ) );
}


/// @brief Initializes energy method
/// @param pose: Pose being evaluated/used for scoring
void
DEEREnergy::initialize_energy_method(
	pose::Pose & pose
) const {

	// Get relevant objects
	scoring::methods::LongRangeEnergyType const & lr_type(
		long_range_type() );
	scoring::Energies & energies( pose.energies() );

	// Determine whether to create new container
	bool create_new_lre_container( false );
	if ( energies.long_range_container( lr_type ) == nullptr ) {
		create_new_lre_container = true;
	} else {
		scoring::LREnergyContainerOP lrc
			= energies.nonconst_long_range_container( lr_type );
		scoring::DenseEnergyContainerOP dec(
			utility::pointer::static_pointer_cast<
			scoring::DenseEnergyContainer > ( lrc ) );
		if ( dec->size() != pose.size() ) {
			create_new_lre_container = true;
		}
	}

	// Create the container
	if ( create_new_lre_container ) {
		scoring::LREnergyContainerOP new_dec(
			new scoring::DenseEnergyContainer( pose.size(),
			scoring::epr_deer_score ) );
		energies.set_long_range_container( lr_type, new_dec );
	}

	// Initialize DEERDataCache object and stuff it in the pose
	if ( !pose.data().has(
			pose::datacache::CacheableDataType::EPR_DEER_DATA )
			) {
		scoring::epr_deer::DEERDataCacheOP ptr(
			new scoring::epr_deer::DEERDataCache() );
		ptr->fetch_and_organize_data( pose );
		pose.data().set(
			pose::datacache::CacheableDataType::EPR_DEER_DATA, ptr );
	}
}

} // namespace energy_methods
} // namespace core
