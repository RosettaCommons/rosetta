// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.


/// @file protocols/simple_moves/DEEROptimizeCoordsMover.cc
/// @brief
/// @author Diego del Alamo

#include <iostream>

#include <protocols/simple_moves/DEEROptimizeCoordsMover.hh>
#include <protocols/simple_moves/DEEROptimizeCoordsMoverCreator.hh>

#include <basic/datacache/CacheableData.hh>
#include <basic/datacache/DataMap.hh>
#include <basic/options/keys/epr_deer.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <basic/options/after_opts.hh>
#include <basic/Tracer.hh>

#include <core/scoring/epr_deer/metrics/DEERData.hh>
#include <core/scoring/epr_deer/metrics/DEERDecayData.hh>
#include <core/scoring/epr_deer/DEERDataCache.hh>
#include <core/scoring/epr_deer/EPRSpinLabel.hh>
#include <core/scoring/epr_deer/Simulated4PDEERTrace.hh>
#include <core/scoring/epr_deer/Simulated4PDEERTraceFactory.hh>
#include <core/scoring/epr_deer/util.hh>

#include <core/chemical/ResidueType.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/types.hh>

#include <protocols/filters/Filter.hh>
#include <protocols/moves/mover_schemas.hh>

#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <utility/vector1.hh>

#include <numeric/random/random_permutation.hh>
#include <numeric/HomogeneousTransform.hh>

#include <string>
#include <unordered_set>
#include <time.h>
#include <cstdlib>
#include <random>

#if defined(_OPENMP)
#include "omp.h"
#endif // OPENMPI

namespace protocols {
namespace simple_moves {

using namespace core::scoring::epr_deer;

static basic::Tracer TR( "protocols.simple_moves.DEEROptimizeCoordsMover" );

/*
ResPair is std::pair< PairSizeString, PairSizeString >
PairSizeString is std::pair< core::Size, std::string >
SizePair is std::pair< core::Size, core::Size >
CoordModel is std::map< PairSizeString, utility::vector1< core::Real > >
*/

/// @brief Default constructor
DEEROptimizeCoordsMover::DEEROptimizeCoordsMover() {}

/// @brief Default constructor
DEEROptimizeCoordsMover::~DEEROptimizeCoordsMover() {}

/// @brief Makes a new, identical version of itself
moves::MoverOP
DEEROptimizeCoordsMover::clone() const {
	return moves::MoverOP( new DEEROptimizeCoordsMover( *this ) );
}

/// @brief Main function
/// @param  pose: Input pose to be modified (in this case it isn't)
void
DEEROptimizeCoordsMover::apply(
	core::pose::Pose & pose
) {

	// Get the data
	core::scoring::epr_deer::DEERDataCache datacache;
	datacache.fetch_and_organize_data( pose );

	// Get the coordinates
	auto default_coords = init_coords( pose, datacache );

	// Get DEER data
	auto trace_map = init_data( default_coords, datacache );

	// Make a model randomly
	auto model = init_model( default_coords, datacache );

	// Initialize starting and ending resolutions
	core::Real const START_RESOLUTION = 1e-1;
	core::Real const END_RESOLUTION = 1e-4;
	auto resolution = START_RESOLUTION;

	// Score being tracked
	core::Real current_score = 0.0;

	// Do an initial refinement for non-3D background
	initial_search( model, resolution, trace_map, datacache );

	// Now loop through and start refining
	do {
		current_score = minimize( model, resolution, trace_map, datacache );
	} while ( resolution >= END_RESOLUTION );

	// Print everything
	print_coords( model, default_coords, current_score, datacache, trace_map, pose );

	// Coordinates applied to pose
	pose.data().set( core::pose::datacache::CacheableDataType::EPR_DEER_DATA,
		datacache.clone() );

}

/// @brief  Initializes coords used throughout method
/// @param  pose: Pose being labeled
/// @return Coords that pass clash evaluation
std::map< PairSizeString, utility::vector1< PseudoSL > >
DEEROptimizeCoordsMover::init_coords(
	core::pose::Pose & pose,
	core::scoring::epr_deer::DEERDataCache const & datacache
) {

	// Output
	std::map< PairSizeString, utility::vector1< PseudoSL > > output;

	// Critical for calculating clashes
	pose.update_residue_neighbors();

	// Initialize various other objects we will need
	core::scoring::epr_deer::EPRSpinLabel sl;

	// Iterate through all residues in datacache
	for ( auto const & res : datacache.labeled_residues() ) {

		// Label and get coordinates
		sl.label( res, pose );
		output[ res ] = utility::vector1< PseudoSL >();

		// Iterate through each PseudoSL, save in vector with weight=0
		for ( auto const & c : sl[ res ] ) {
			if ( c.second > sl.cutoff() ) {
				output[ res ].push_back( std::make_pair( c.first, 0 ) );
			}
		}
	}
	return output;
}

/// @brief Initializes the data saved in the Mover
/// @param pose: Pose being modified
/// @param  datacache: Datacache
/// @detail ResPair is std::pair< PairSizeString, PairSizeString >
/// @detail PairSizeString is std::pair< core::Size, std::string >
/// @detail SizePair is std::pair< core::Size, core::Size >
/// @detail Map used: data[ resPair ][ std::make_pair( res1_id, res2_id ) ]
std::map< DEEROptimizeCoordsMover::ResPair, std::map< DEEROptimizeCoordsMover::SizePair, utility::vector1< core::Real > > >
DEEROptimizeCoordsMover::init_data(
	std::map< PairSizeString, utility::vector1< PseudoSL > > const & default_coords,
	core::scoring::epr_deer::DEERDataCache & datacache
) {

	// Initialize output
	std::map< ResPair, std::map< SizePair, utility::vector1< core::Real > > > output;

	// Now go through each DEER trace. We will simulate a single trace for
	// each pair of PseudoSLs across both spin labeled residues for each
	// DEER trace. This allows the recomputation to be avoided

	for ( auto const & i : datacache.indices() ) {

		// TODO: Failsafe to ensure that all the data is decay data
		// Case the data to get the DEER decay
		core::scoring::epr_deer::metrics::DEERDecayDataOP data
			= utility::pointer::dynamic_pointer_cast<
			metrics::DEERDecayData
			>( datacache[ i ] );

		// Bin width for distribution, to be used below
		auto const bpa = data->bins_per_a();

		// OUTER LOOP: Go through residues, pair by pair
		// Usually, this will not be a loop at all, as two residues
		//  are labeled per DEER trace. This is in the off-chance that
		// more than two are labeled (rare).
		auto const residues = data->residues();
		for ( core::Size r1 = 1; r1 < residues.size(); ++r1 ) {
			auto const res1 = residues[ r1 ];
			for ( core::Size r2 = r1 + 1; r2 <= residues.size(); ++r2 ) {
				auto const res2 = residues[ r2 ];
				auto const key1 = std::make_pair( res1, res2 );
				if ( output.find( key1 ) == output.end() ) {
					output[ key1 ] = std::map< SizePair, utility::vector1< core::Real > >();
				}

				// Calculate the width of the distribution for the pair
				core::Real const s = 4 * data->stdev();

				// INNER LOOP: Go through PseudoSLs, pair by pair
				auto const & c1 = default_coords.at( res1 );
				auto const & c2 = default_coords.at( res2 );
				for ( core::Size e1 = 1; e1 <= c1.size(); ++e1 ) {
					for ( core::Size e2 = 1; e2 <= c2.size(); ++e2 ) {
						// Get the distance between the two coordinates
						auto const d = c1[ e1 ].first.distance( c2[ e2 ].first );

						// Get the minimum and maximum bin
						// Lowest possible bin is 1! No negative distances!
						auto const start = ( d > s ) ? round( bpa * ( d - s ) ) : 1;
						auto const end = round( bpa * ( d + s ) );

						// Set up the DEER distance distribution
						std::map< core::Size, core::Real > distr;
						core::Real total = 0.0;
						for ( auto bin = start; bin <= end; ++bin ) {

							// Calculate amplitude at that bin
							auto d_bin = bin / core::Real( bpa );
							distr[ bin ] = gauss( d_bin, d, data->stdev() );
							total += distr[ bin ];
						}

						// Normalize the distribution so it adds up to 1.0
						for ( auto & bin_amp : distr ) {
							bin_amp.second /= total;
						}

						// Stash this in the DEER trace map
						auto const trace = data->factory().kernel_mult( distr );
						auto const key2 = std::make_pair( e1, e2 );
						output[ key1 ][ key2 ] = trace;
					}
				}
			}
		}
	}
	return output;
}

/// @brief Initialize and run the optimization procedure
/// @param  start_resolution: Starting resolution for modifications
/// @param  end_resolution: End resolution for modifications
/// @return Best coordinate model obtained (by sum of squared residuals)
DEEROptimizeCoordsMover::CoordModel
DEEROptimizeCoordsMover::init_model(
	std::map< PairSizeString, utility::vector1< PseudoSL > > const & coords,
	core::scoring::epr_deer::DEERDataCache & datacache
) {

	// Initialize model
	DEEROptimizeCoordsMover::CoordModel model;

	// Randomly assign some starting coordinates
	for ( auto const & res : datacache.labeled_residues() ) {
		auto n = coords.at( res ).size();
		model[ res ] = utility::vector1< core::Real >( n, 0 );
		model[ res ][ numeric::random::random_range( 1, n ) ] = 1.0;
	}
	return model;
}

/// @brief Performs an initial search (for non-3D backgrounds)
/// @param  model: Model to refine
/// @param  resolution: Refinement resolution
/// @param  trace_map: Map of DEER traces for calculation
/// @param  datacache: Datacache with data
void
DEEROptimizeCoordsMover::initial_search(
	CoordModel & model,
	core::Real & resolution,
	std::map< ResPair, std::map< SizePair, utility::vector1< core::Real > > > const & trace_map,
	core::scoring::epr_deer::DEERDataCache & datacache
) {

	// Assemble a set of indices with non-3D backgrounds
	std::set< core::Size > non_3d_ids;
	for ( auto const & idx : datacache.indices() ) {
		core::scoring::epr_deer::metrics::DEERDecayDataOP data
			= utility::pointer::dynamic_pointer_cast<
			core::scoring::epr_deer::metrics::DEERDecayData
			>( datacache[ idx ] );

		// Add to set if background is non-3D and change temporarily
		if ( data->factory().bckg_type() == "NON_3D" ) {
			non_3d_ids.insert( idx );
			data->factory().bckg_type( "3D" );
		}
	}

	// Check in case no data meets this criterion
	if ( non_3d_ids.empty() ) {
		return;
	}

	// Now run the initial minimization
	minimize( model, resolution, trace_map, datacache );

	// Now go back and reset the backgrounds
	for ( auto const & idx : non_3d_ids ) {
		core::scoring::epr_deer::metrics::DEERDecayDataOP data
			= utility::pointer::dynamic_pointer_cast<
			core::scoring::epr_deer::metrics::DEERDecayData
			>( datacache[ idx ] );
		data->factory().bckg_type( "NON_3D" );
	}
}

/// @brief Minimizes the model to fit the data
/// @param  model: Model to refine
/// @param  resolution: Refinement resolution
/// @param  trace_map: Map of DEER traces for calculation
/// @param  datacache: Datacache with data
/// @return Score of the model
core::Real
DEEROptimizeCoordsMover::minimize(
	CoordModel & model,
	core::Real & resolution,
	std::map< ResPair, std::map< SizePair, utility::vector1< core::Real > > > const & trace_map,
	core::scoring::epr_deer::DEERDataCache & datacache
) {

	// Get residues
	auto residues = datacache.labeled_residues();

	// Initialize
	auto current_score = score( model, trace_map, datacache );
	auto best_score = current_score;
	auto best_model = model;

	// Initialize variables
	core::Size improvements = 0;
	core::Size failures = 0;
	core::Size const MAX_TRIALS = datacache.size() * 2500;
	core::Size const MAX_ATTEMPTS = 500;
	core::Real t = 1.5;

	// Go through the loop
	for ( core::Size trials = 1; trials <= MAX_TRIALS; ++trials ) {

		// Increment the variables
		++failures;
		t -= ( t / MAX_TRIALS );

		// Print periodically
		if ( trials % 100 == 0 ) {
			TR.Trace << "\tTRIAL " << trials << ":\t" << current_score << std::endl;
		}

		// Exit possibility: no consecutive improvements
		if ( failures > MAX_ATTEMPTS ) {
			break;
		}

		// Choose a random residue and a random PseudoSL
		auto res = residues[ numeric::random::random_range( 1, residues.size() ) ];
		auto e_id = numeric::random::random_range( 1, model[ res ].size() );

		// Deep copy the model
		auto model_copy = model;

		// Assign the new weight of the PseudoSL
		auto const & w = model[ res ][ e_id ];
		auto delta = std::max( -1 * w, 2 * resolution
			* ( numeric::random::uniform() - 0.5 ) );
		model_copy[ res ][ e_id ] += delta;

		// Correct any negatives
		if ( model_copy[ res ][ e_id ] < 0.0 ) {
			model_copy[ res ][ e_id ] = 0.0;
		}

		// Score the model and check if we have an improvement
		auto new_score = score( model_copy, trace_map, datacache );
		if ( boltzmann( current_score - new_score, t ) ) {

			// Normalize the coordinates so that they add up to one
			normalize( model_copy, res );

			// Also save these new coordinates and the new score
			model = model_copy;
			current_score = new_score;
			if ( new_score < best_score ) {
				best_score = new_score;
				best_model = model;
				failures = 0;
				++improvements;
			}
		}
	}

	// Save the best model
	model = best_model;

	// Reduce resolution if no improvements were found
	if ( improvements == 0 ) {
		TR.Trace << "Reducing resolution to from " << resolution;
		resolution /= sqrt( 10 );
		TR.Trace << " to " << resolution << std::endl;
	}

	// Return the best score
	return best_score;
}

bool
DEEROptimizeCoordsMover::boltzmann(
	core::Real const & diff,
	core::Real const & temp // = 1.0
) const {
	return ( numeric::random::uniform() < ( diff / temp ) );
}

/// @brief Recursive function to normalize the weights of a model
/// @param model: The model being normalized
/// @param res: Residue to normalize
/// @detail If residue is zero, the function calls itself for each residue
void
DEEROptimizeCoordsMover::normalize(
	CoordModel & model,
	PairSizeString const & res // = std::make_pair( 0, "" )
) const {

	// If no residue is specified
	if ( res.first == 0 && res.second == "" ) {
		for ( auto const & res_vec : model ) {
			normalize( model, res_vec.first );
		}

		// Otherwise...
	} else {

		// ...add up the total weights
		core::Real total = 0.0;
		for ( auto & coord : model[ res ] ) {
			if ( coord < 0.0 ) {
				coord = 0.0;
			} else {
				total += coord;
			}
		}

		// If they are equal to zero
		if ( total <= 0.0 ) {
			utility_exit_with_message( "Weights equal zero! Quitting." );
		}
		for ( auto & coord: model[ res ] ) {
			coord /= total;
		}
	}
}

core::Real
DEEROptimizeCoordsMover::score(
	CoordModel const & coords,
	std::map< ResPair, std::map< SizePair, utility::vector1< core::Real > > > const & trace_map,
	core::scoring::epr_deer::DEERDataCache & datacache,
	bool const print // = false
) {

	// Tracks total number of parameters
	core::Size n = 0;

	// Actual score
	core::Real total_score = 0.0;

	// Iterate through each dataset
	for ( auto const & idx : datacache.indices() ) {

		// Get the data and convert to proper format
		core::scoring::epr_deer::metrics::DEERDecayDataOP data
			= utility::pointer::dynamic_pointer_cast<
			core::scoring::epr_deer::metrics::DEERDecayData
			>( datacache[ idx ] );

		// Simulated DEER trace
		utility::vector1< core::Real > sim_trace(
			data->factory().trace().size(), 0.0 );

		n += sim_trace.size();

		// Now iterate through pairs of residues
		// Make aliases for the residue ID, the PseudoSL coordinates
		for ( core::Size r1 = 1; r1 < data->residues().size(); ++r1 ) {
			auto const & r1_idx = data->residues()[ r1 ];
			auto const & r1_evec = coords.at( r1_idx );
			for ( core::Size r2 = r1 + 1; r2 <= data->residues().size(); ++r2 ) {
				auto const & r2_idx = data->residues()[ r2 ];
				auto const & r2_evec = coords.at( r2_idx );

				// Alias for pair
				auto const rpair = std::make_pair( r1_idx, r2_idx );

				// Now iterate through pairs of SLs
				for ( core::Size e1 = 1; e1 <= r1_evec.size(); ++e1 ) {
					auto const & w1 = coords.at( r1_idx ).at( e1 );
					for ( core::Size e2 = 1; e2 <= r2_evec.size(); ++e2 ) {
						auto const & w2 = coords.at( r2_idx ).at( e2 );

						// Alias for pair
						auto const epair = std::make_pair( e1, e2 );

						// Alias for DEER trace between these
						auto const & pairtrace = trace_map.at( rpair ).at( epair );

						// Add it to the total trace
						for ( core::Size i = 1; i <= pairtrace.size(); ++i ) {
							sim_trace[ i ] += w1 * w2 * pairtrace.at( i );
						}
					}
				}
			}

			// Now normalize so that the max element is 1.0
			auto max_v = *std::max_element( sim_trace.begin(), sim_trace.end() );
			for ( auto & v : sim_trace ) {
				v /= max_v;
			}

			// Now add the background and score}

			if ( numeric::random::random_range( 1, 100 ) == 1 ) {
				data->factory().reset();
			}
			auto bckg_trace = data->factory().opt_bckg( sim_trace );
			total_score += data->sum_of_squares( bckg_trace, false );
			if ( print ) {
				TR.Info << "\tDataset " << idx << std::endl;
				for ( core::Size i = 1; i <= bckg_trace.size(); ++i ) {
					TR.Info << "\t\t" << data->factory().time_pts()[ i ] << "\t"
						<< data->factory().trace()[ i ] << "\t"
						<< bckg_trace[ i ] << std::endl;
				}
			}
		}
	}

	// Return negative log-likelihood
	auto nll = core::Real( n ) * 0.5 * log( total_score / std::max( 1, int( n ) ) );
	return nll;
}

void
DEEROptimizeCoordsMover::print_coords(
	CoordModel const & model,
	std::map< PairSizeString, utility::vector1< PseudoSL > > const & default_coords,
	core::Real const & current_score,
	core::scoring::epr_deer::DEERDataCache & datacache,
	std::map< ResPair, std::map< SizePair, utility::vector1< core::Real > > > const & trace_map,
	core::pose::Pose & pose
) {

	// Header
	TR.Info << "\n===================" << std::endl;;
	TR.Info << "LOWEST: " << current_score << std::endl;

	// Needed to get put coords in global frame
	core::scoring::epr_deer::EPRSpinLabel sl;

	// Iterate through residues
	for ( auto const & r : model ) {

		// Label
		sl.label( r.first, pose, true );

		// Global-to-local frame
		numeric::HomogeneousTransform< core::Real > frame(
			pose.residue( r.first.first ).xyz( "N"  ),
			pose.residue( r.first.first ).xyz( "C"  ),
			pose.residue( r.first.first ).xyz( "CA" ) );

		// Iterate through each PseudoSL
		for ( core::Size i = 1; i <= r.second.size(); ++i ) {

			// Print out all relevant information
			// Residue, weight, local frame, global frame
			TR.Info << "COORD\t" << r.first.first << "\t" << r.second[ i ] << "\t";

			// Print local data (for Rosetta)
			auto local = frame.to_local_coordinate( default_coords.at( r.first ).at( i ).first );
			TR.Info << local.x() << "\t" << local.y() << "\t" << local.z() << "\t#\t";

			// Print global XYZ data
			auto global = sl[ r.first ][ i ].first;
			TR.Info << global.x() << "\t" << global.y() << "\t" << global.z() << "\t#\t";
			TR.Info << "\n";
		}
	}
	TR.Info << std::endl;
	TR.Info << "===================" << std::endl;
	score( model, trace_map, datacache, true );
	// using CoordModel = std::map< PairSizeString, utility::vector1< core::Real > >;
	datacache.set_labels( { sl } );
	datacache.set_sl_weights( { 1.0 } );
}

std::string
DEEROptimizeCoordsMover::get_name() const {
	return mover_name();
}

std::string
DEEROptimizeCoordsMover::mover_name() const {
	return "DEEROptimizeCoordsMover";
}

std::string
DEEROptimizeCoordsMover::keyname() const {
	return mover_name();
}

moves::MoverOP
DEEROptimizeCoordsMover::create_mover() const {
	return moves::MoverOP( new DEEROptimizeCoordsMover );
}

void
DEEROptimizeCoordsMover::provide_xml_schema(
	utility::tag::XMLSchemaDefinition & xsd
) {
	using namespace utility::tag;
	AttributeList attributes;
	moves::xsd_type_definition_w_attributes( xsd, mover_name(),
		"Mover for optimizing rotamer positions for DEER data", attributes );
}

void
DEEROptimizeCoordsMover::parse_my_tag(
	utility::tag::TagCOP,
	basic::datacache::DataMap &
) {}

//////////////////////////////////

std::string
DEEROptimizeCoordsMoverCreator::keyname() const {
	return DEEROptimizeCoordsMover().mover_name();
}

protocols::moves::MoverOP
DEEROptimizeCoordsMoverCreator::create_mover() const {
	return moves::MoverOP( new DEEROptimizeCoordsMover );
}

void
DEEROptimizeCoordsMoverCreator::provide_xml_schema(
	utility::tag::XMLSchemaDefinition & xsd
) const {
	DEEROptimizeCoordsMover().provide_xml_schema( xsd );
}

}
}
