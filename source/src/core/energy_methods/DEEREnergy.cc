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
#include <core/scoring/epr_deer/DEERData.hh>
#include <core/scoring/epr_deer/DEERDataCache.hh>
#include <core/scoring/epr_deer/EPRSpinLabel.hh>

// Package headers
#include <core/chemical/AA.hh>
#include <core/pose/Pose.hh>
#include <core/pose/datacache/CacheableDataType.hh>
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

static basic::Tracer TR( "core.energy_methods.DEEREnergy" );

/// @brief  Energy Method creator required by framework (registrator in devel/init.cc)
scoring::methods::EnergyMethodOP
DEEREnergyCreator::create_energy_method(
	scoring::methods::EnergyMethodOptions const & options
) const {
	return scoring::methods::EnergyMethodOP( new DEEREnergy( options ) );
}

/// @brief  Energy Method descriptor to identify it if a weight is provided in the command line
scoring::ScoreTypes
DEEREnergyCreator::score_types_for_method() const {
	scoring::ScoreTypes sts;
	sts.push_back( scoring::epr_deer_score );
	return sts;
}

/// @brief  Constructor
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
scoring::methods::EnergyMethodOP
DEEREnergy::clone() const {
	return scoring::methods::EnergyMethodOP( new DEEREnergy() );
}

/// @brief Derived function for specifying where the data is getting stored
scoring::methods::LongRangeEnergyType
DEEREnergy::long_range_type() const {
	return scoring::methods::epr_deer_lr;
}

/// @brief Overrides virtual fxn for finding scores for a given pose
void
DEEREnergy::setup_for_scoring(
	pose::Pose & pose,
	scoring::ScoreFunction const &
) const {
	// Cursory check that everything is in order
	initialize_energy_method( pose );

	scoring::epr_deer::DEERDataCacheOP deer_data = get_deer_data( pose );

	// This checks if multiple sets of spin labels is being used for scoring
	// This is used when spin labels are obtained for specific proteins a priori,
	// e.g., by triangulation/multilateration, and the solution is underdetermined.
	if ( get_const_deer_data( pose )->labels().size() > 0 ) {
		std::map< Size, Real > total_scores;
		// First the scores for each set are calculated and added
		for ( core::Size i = 1; i <= deer_data->labels().size(); ++i ) {
			Real weight = deer_data->sl_weights()[ i ];
			scoring::epr_deer::EPRSpinLabel sl = deer_data->labels()[ i ];
			for ( core::Size j = 1; j <= deer_data->size(); ++j ) {
				if ( total_scores.find( j ) == total_scores.end() ) {
					total_scores[ j ] = weight * get_score( pose, sl, j );
				} else {
					total_scores[ j ] += weight * get_score( pose, sl, j );
				}
			}
		}
		// Then the scores are added and stored
		for ( core::Size j = 1; j <= deer_data->size(); ++j ) {
			( *deer_data )[ j ]->score( total_scores[ j ] );
		}
	} else {
		// Default method if no special custom spin labels are used
		scoring::epr_deer::EPRSpinLabel sl;
		for ( core::Size i = 1; i <= deer_data->size(); ++i ) {
			get_score( pose, sl, i, 0.0 );
		}
	}
}

/// @brief Calculates the score for a specific data set
Real
DEEREnergy::get_score(
	pose::Pose & pose,
	scoring::epr_deer::EPRSpinLabel & sl,
	Size const & edge_id,
	Real const & modifier // = 0.0
) const {
	scoring::epr_deer::DEERDataCacheOP deer_data = get_deer_data( pose );
	for ( auto const & res : deer_data->labeled_residues() ) {
		sl.label( res.first, res.second, pose, false );
	}
	auto & edge_data = ( *deer_data )[ edge_id ];
	auto histogram = sl.histogram( edge_data->residues(), edge_data->bins_per_angstrom(), modifier );
	edge_data->print_histogram( histogram );
	return edge_data->get_score( histogram, ( modifier == 0.0 ) );
}

/// @brief Initializes energy method
void
DEEREnergy::initialize_energy_method(
	pose::Pose & pose
) const {
	scoring::methods::LongRangeEnergyType const & lr_type( long_range_type() );
	scoring::Energies & energies( pose.energies() );

	bool create_new_lre_container( false );
	if ( energies.long_range_container( lr_type ) == nullptr ) {
		create_new_lre_container = true;
	} else {
		scoring::LREnergyContainerOP lrc = energies.nonconst_long_range_container( lr_type );
		scoring::DenseEnergyContainerOP dec( utility::pointer::static_pointer_cast<
			scoring::DenseEnergyContainer > ( lrc ) );
		if ( dec->size() != pose.size() ) {
			create_new_lre_container = true;
		}
	}
	if ( create_new_lre_container ) {
		scoring::LREnergyContainerOP new_dec( new scoring::DenseEnergyContainer( pose.size(), scoring::epr_deer_score ) );
		energies.set_long_range_container( lr_type, new_dec );
	}
	// Initialize DEERDataCache object and stuff it in the pose
	if ( !pose.data().has( pose::datacache::CacheableDataType::EPR_DEER_DATA ) ) {
		scoring::epr_deer::DEERDataCacheOP deer_data_ptr( new scoring::epr_deer::DEERDataCache() );
		deer_data_ptr->fetch_and_organize_data( pose );
		pose.data().set( pose::datacache::CacheableDataType::EPR_DEER_DATA, deer_data_ptr );
	}
	// Necessary for clash evaluation using the TenANeighborGraph later (EPRSpinLabel)
	pose.update_residue_neighbors();
}

/// @brief Get energy/score for a given pair of residues.
void
DEEREnergy::residue_pair_energy(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	pose::Pose const & pose,
	scoring::ScoreFunction const &,
	scoring::EnergyMap & emap
) const {
	auto const & edges = get_const_deer_data( pose )->edge( rsd1.seqpos(), rsd2.seqpos() );
	// If this object isn't empty, iterate through it
	if ( edges.size() > 0 ) {
		for ( auto const & edge_weight : edges ) {
			Size edge_id = edge_weight.first;
			Real weight_edge = edge_weight.second; // This is to avoid double-counting
			Real score = get_const_deer_data( pose )->at( edge_id )->score();
			Real weight_total = get_const_deer_data( pose )->at( edge_id )->relative_weight();
			emap[ scoring::ScoreType::epr_deer_score ] += score * weight_total * weight_edge;
		}
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
bool
DEEREnergy::defines_intrares_energy(
	scoring::EnergyMap const &
) const {
	return false;
}

/// @brief Fetches if two residues are "connected" (i.e. there is a score)
bool
DEEREnergy::defines_residue_pair_energy(
	pose::Pose const & pose,
	Size rsd1,
	Size rsd2
) const {
	return ( get_const_deer_data( pose )->edge( rsd1, rsd2 ).size() > 0 ) ? true : false;
}

/// @brief  Calculates the derivatives for minimization. Done numerically
void
DEEREnergy::setup_for_derivatives(
	pose::Pose & pose,
	scoring::ScoreFunction const &
) const {
	// Instantiate the maps where these terms are stored
	std::map< Size, numeric::xyzVector< Real > > f1_forces;
	std::map< Size, numeric::xyzVector< Real > > f2_forces;

	// Store the scores for forward and reverse movement
	std::map< int, std::map< Size, Real > > total_scores;
	total_scores[ -1 ] = std::map< Size, Real >();
	total_scores[ +1 ] = std::map< Size, Real >();

	// As with setup_for_scoring(), this is in the event that there are multiple
	// sets of custom spin labels used. See above for details on why
	if ( get_const_deer_data( pose )->labels().size() > 0 ) {
		for ( core::Size i = 1; i <= get_const_deer_data( pose )->labels().size(); ++i ) {
			Real sl_weight = get_const_deer_data( pose )->sl_weights()[ i ];
			scoring::epr_deer::EPRSpinLabel sl = get_const_deer_data( pose )->labels()[ i ];
			for ( core::Size j = 1; j <= get_const_deer_data( pose )->size(); ++j ) {
				if ( total_scores[ -1 ].find( j ) == total_scores[ -1 ].end() ) {
					total_scores[ -1 ][ j ] =  sl_weight * get_score( pose, sl, j, -1 );
					total_scores[ +1 ][ j ] =  sl_weight * get_score( pose, sl, j, +1 );
				} else {
					total_scores[ -1 ][ j ] += sl_weight * get_score( pose, sl, j, -1 );
					total_scores[ +1 ][ j ] += sl_weight * get_score( pose, sl, j, +1 );
				}
			}
		}
		// Normal case scenario
	} else {
		scoring::epr_deer::EPRSpinLabel sl;
		for ( auto const & res : get_const_deer_data( pose )->labeled_residues() ) {
			sl.label( res.first, res.second, pose, false );
		}
		for ( core::Size i = 1; i <= get_const_deer_data( pose )->size(); ++i ) {
			Real edge_weight = get_const_deer_data( pose )->at( i )->relative_weight();
			total_scores[ -1 ][ i ] = get_score( pose, sl, i, -1 ) * edge_weight;
			total_scores[ +1 ][ i ] = get_score( pose, sl, i, +1 ) * edge_weight;
		}
	}

	// Now go over all the scores and apply to residues.
	// Uses inter-residue vector to calculate the precise vector for minimization
	for ( Size i = 1; i <= get_const_deer_data( pose )->size(); ++i ) {
		auto const & residues = get_const_deer_data( pose )->at( i )->residues();
		for ( Size res1_iter = 1; res1_iter <= residues.size(); ++res1_iter ) {
			Size res1 = residues[ res1_iter ].first;
			numeric::xyzVector< Real > const & res1_ca = pose.residue( res1 ).xyz( "CA" );
			for ( Size res2_iter = 1; res2_iter <= residues.size(); ++res2_iter ) {
				if ( res1_iter == res2_iter ) {
					continue;
				}
				Size res2 = residues[ res2_iter ].first;
				numeric::xyzVector< Real > const & res2_ca = pose.residue( res2 ).xyz( "CA" );
				if ( f2_forces.find( res1 ) == f2_forces.end() ) {
					f2_forces[ res1 ] = numeric::xyzVector< Real >(0,0,0);
				}
				f2_forces[ res1 ] += ( res2_ca - res1_ca ).normalize() * ( total_scores[ +1 ][ i ] - total_scores[ -1 ][ i ] );

				if ( f1_forces.find( res1 ) == f1_forces.end() ) {
					f1_forces[ res1 ] = numeric::xyzVector< Real >(0,0,0);
				}
				f1_forces[ res1 ] += cross( res2_ca, res1_ca ) * ( total_scores[ +1 ][ i ] - total_scores[ -1 ][ i ] ) / 2.0;
			}
		}
	}
	// Now apply these values to the pose
	for ( auto const & res_vec : f1_forces ) {
		get_deer_data( pose )->f1_force( res_vec.first, res_vec.second );
	}
	for ( auto const & res_vec : f2_forces ) {
		get_deer_data( pose )->f2_force( res_vec.first, res_vec.second );
	}
}

/// @brief  Returns the parsed input file (const)
scoring::epr_deer::DEERDataCacheCOP const
DEEREnergy::get_const_deer_data(
	pose::Pose const & pose
) const {
	return utility::pointer::dynamic_pointer_cast< const scoring::epr_deer::DEERDataCache > (
		pose.data().get_const_ptr( pose::datacache::CacheableDataType::EPR_DEER_DATA ) );
}

/// @brief  Returns the parsed input file (non-const)
scoring::epr_deer::DEERDataCacheOP
DEEREnergy::get_deer_data(
	pose::Pose & pose
) const {
	return utility::pointer::dynamic_pointer_cast< scoring::epr_deer::DEERDataCache >(
		pose.data().get_ptr( pose::datacache::CacheableDataType::EPR_DEER_DATA ) );
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
	if ( !get_const_deer_data( pose )->has_f1_force( id.rsd() ) ||
			!get_const_deer_data( pose )->has_f2_force( id.rsd() )
			) {
		return;
	}
	if ( id.atomno() > 3 ) {
		return;
	}
	F1 += emap[ scoring::ScoreType::epr_deer_score ] * get_const_deer_data( pose )->f1_force( id.rsd() );
	F2 += emap[ scoring::ScoreType::epr_deer_score ] * get_const_deer_data( pose )->f2_force( id.rsd() );
}

/// @brief  A declared but unused function
void
DEEREnergy::indicate_required_context_graphs(
	utility::vector1< bool > &
) const {}

/// @brief  Version 1 (as of 20 November 2017)
Size
DEEREnergy::version() const {
	return 1;
}

} // namespace energy_methods
} // namespace core
