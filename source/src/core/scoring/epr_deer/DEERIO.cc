// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/epr_deer/DEERIO.cc
/// @brief  IO class for data obtained with double electron-electron resonance (DEER)
/// @details This class is only called once at the beginning of the application by the
///      energy method. It then parses the input file and the electron coordinate file.
///      The input file is then used to create a "map" (essentially a fake graph) which
///      then becomes the basis of the DEERData object.
/// @author  Diego del Alamo ( del.alamo@vanderbilt.edu )

#include <core/pose/selection.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/ResidueIndexDescription.hh>
#include <core/scoring/epr_deer/metrics/DEERData.hh>
#include <core/scoring/epr_deer/DEERIO.hh>
#include <core/scoring/epr_deer/EPRSpinLabel.hh>
#include <core/scoring/epr_deer/util.hh>
#include <core/scoring/epr_deer/metrics/DEERData.hh>
#include <core/scoring/epr_deer/metrics/DEERDistanceDistribution.hh>
#include <core/scoring/epr_deer/metrics/DEERDistanceDistribution.fwd.hh>
#include <core/scoring/epr_deer/metrics/DEERDistanceBounds.hh>
#include <core/scoring/epr_deer/metrics/DEERDistanceBounds.fwd.hh>
#include <core/scoring/epr_deer/metrics/DEERDecayData.hh>
#include <core/scoring/epr_deer/metrics/DEERMiscMethod.hh>
#include <core/scoring/epr_deer/metrics/DEERMiscMethod.fwd.hh>
#include <core/scoring/epr_deer/metrics/DEEROverlapMethod.hh>
#include <core/scoring/epr_deer/metrics/DEEROverlapMethod.fwd.hh>
#include <core/scoring/epr_deer/metrics/DEERJaccardMethod.hh>
#include <core/scoring/epr_deer/metrics/DEERJaccardMethod.fwd.hh>
#include <core/scoring/epr_deer/metrics/DEERJSMethod.hh>
#include <core/scoring/epr_deer/metrics/DEERJSMethod.fwd.hh>
#include <core/scoring/epr_deer/metrics/DEERChiSqMethod.hh>
#include <core/scoring/epr_deer/metrics/DEERChiSqMethod.fwd.hh>
#include <core/types.hh>

// Numeric headers
#include <numeric/xyzVector.hh>
#include <numeric/constants.hh>

// Basic headers
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/epr_deer.OptionKeys.gen.hh>
#include <basic/options/option_macros.hh>
#include <basic/database/open.hh>

// Utility headers
#include <utility/io/izstream.hh>
#include <utility/fixedsizearray1.hh>
#include <utility/excn/Exceptions.hh>

// ObjexxFCL headers
#include <ObjexxFCL/string.functions.hh>

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

#include <boost/algorithm/string.hpp>

namespace core {
namespace scoring {
namespace epr_deer {

/// @brief Tracer used for error messages
/// @details Global to avoid re-instantiating tracer with every new object
static basic::Tracer TR( "core.scoring.epr_deer.DEERIO" );

/// @brief Split an input line into vector of strings
/// @param  line: Line of interest
/// @return Vector of strings
/// @detail Ignores lines starting with "#"
SplitLine
DEERIO::split_line(
	std::string const & line
) const {

	// First clear out any whitespace from the beginning and end
	utility::trim( line, " \t\n" );

	// Only proceed if the line actually has something worth reading
	if ( line.size() == 0 || line[ 0 ] == '#' ) {
		return SplitLine();
	}

	// Instantiate output
	SplitLine sl, sl_out;

	// Split line into vector
	boost::algorithm::split( sl, line, boost::is_any_of( "\t " ) );

	// Remove empty strings (caused by two consecutive whitespace chars)
	for ( Size i = 1; i <= sl.size(); ++i ) {
		if ( sl[ i ].size() > 0 ) {
			sl_out.push_back( sl[ i ] );
		}
	}

	// Return the output
	return sl_out;
}

/// @brief Generates the data for use. Core function
/// @param  pose: Input pose (used for residue numbering)
/// @return Map of input data
/// @detail Will only compute this once! Otherwise, the data is saved
/// @detail  in a static/global object (epr_deer_input_data_) defined in
/// @detail  header file
std::map< Size, metrics::DEERDataOP >
DEERIO::generate_data(
	core::pose::Pose const & pose
) const {

	// Return pre-calculated data if file was read before (saves time)
	if ( !epr_deer_input_data_.empty() ) {
		return epr_deer_input_data_;
	}

	// Define output
	std::map< Size, metrics::DEERDataOP > output;

	// Iterate over lines grouped by object
	for ( auto const & data : parse_input() ) {
		output[ data.first ] = parse_data( pose, data.second );
	}

	// Return data
	return output;
}

/// @brief  Parses file
/// @param  pose: Pose (for residue numbering)
/// @return Map of strings (grouped by data index)
std::map< Size, SplitLines >
DEERIO::parse_input() const {

	// To make things less ugly
	using namespace basic::options;
	SplitLines lines;

	// Fetch input files (assuming the option has been passed)
	utility::vector1< std::string > input_files =
		option[ OptionKeys::epr_deer::input_files ]();

	// Outermost loop over input files
	for ( auto const & input_file : input_files ) {

		// Open the file
		std::ifstream infile;
		infile.open( input_file.c_str(), std::ios::in );
		if ( !infile.is_open() ) {
			throw CREATE_EXCEPTION( utility::excn::FileNotFound,
				input_file );
		}

		// Innermost loop: iterate through each line of the file
		std::string line;
		while ( std::getline( infile, line ) ) {
			auto splitline = split_line( line );
			if ( splitline.size() > 0 ) {
				lines.push_back( splitline );
			}
		}
	}

	// Actual output
	// Data for lines starting with "PAIR" are treated differently due
	//  to indexing (not mentioned in file)
	std::map< Size, SplitLines > output;
	SplitLines pair_data;

	// Iterate over lines
	for ( auto line : lines ) {

		// Set aside "PAIR" data as described above
		if ( line[ 1 ] == "PAIR" ) {
			pair_data.push_back( line );

			// Sort the other lines
		} else {
			Size const & edge_id = std::stoi( line[ 2 ] );
			append_to_map( output, edge_id, line );
		}
	}

	// Now we process the "PAIR" data
	Size new_index = output.rbegin()->first + 1;

	for ( auto const & pair_line : pair_data ) {

		// These two lines are placeholders for the lines that would be present
		//  if the data had been passed using "BOUNDS".
		auto line1 = pairline1_;
		auto line2 = pairline2_;

		// Replace the empty spaces with information from the "PAIR" object
		line1[ 2 ] = std::to_string( new_index );
		line1[ 4 ] = pair_line[ 2 ];
		line1[ 6 ] = pair_line[ 3 ];
		line2[ 2 ] = std::to_string( new_index );
		line2[ 3 ] = pair_line[ 4 ];
		line2[ 4 ] = pair_line[ 4 ];

		// Now set aside a new index for these lines and increment index
		output[ new_index ] = SplitLines{ line1, line2 };
		++new_index;
	}

	// Return output
	return output;
}

/// @brief Parse the input lines into actual data containers
/// @param  pose: Pose (used for residue numbering)
/// @param  idx: Index of datum (used for error messages)
/// @param  lines: Lines associated with this dataset
/// @return Proper data container that can be used by Rosetta
metrics::DEERDataOP
DEERIO::parse_data(
	pose::Pose const & pose,
	SplitLines const & lines
) const {

	// Start by sorting the lines by linetype (first element in splitline)
	std::map< std::string, SplitLines > sorted_data{
		{ "DESC", SplitLines() }, { "INFO", SplitLines() } };
	for ( auto const & line : lines ) {
		auto linetype = line[ 1 ];

		// Carve out a spot in the DataStringMap for this line if none exists
		append_to_map( sorted_data, linetype, line );
	}

	// Now figure out what kind of DEERData object we need
	// We will do this by looking for the actual lines used by each one
	if ( sorted_data.find( "DECAY" ) != sorted_data.end() ) {
		return parse_decay_data( pose, sorted_data );
	} else if ( sorted_data.find( "BOUNDS" ) != sorted_data.end() ) {
		return parse_bounds_data( pose, sorted_data );
	} else if ( sorted_data.find( "DIST" ) != sorted_data.end() ) {
		return parse_dist_data( pose, sorted_data );
	} else {
		throw CREATE_EXCEPTION( utility::excn::BadInput,
			"Insufficient information in dataset to deduce data type!" );
	}
}

/// @brief  Parse inputs lines into DEERDecayDataOP
/// @param  pose: Pose (used for residue renumbering)
/// @param  sorted_data: Input lines sorted by linetype
/// @return Proper data container that can be used by Rosetta
metrics::DEERDataOP
DEERIO::parse_decay_data(
	pose::Pose const & pose,
	std::map< std::string, SplitLines > const & sorted_data
) const {

	// Output type is known at this point
	metrics::DEERDataOP data( new metrics::DEERDecayData() );

	// Get residues
	data->residues( parse_desc_lines( pose, sorted_data.at( "DESC" ) ) );

	// Process the actual distance data (data passed by reference)
	parse_decay_info_lines( sorted_data.at( "INFO" ), data );

	// Initialize the factory object used for simulating the raw data
	parse_decay_lines( sorted_data.at( "DECAY" ), data );

	// Return output
	return data;
}

/// @brief  Parse inputs lines into DEERDistanceBoundsOP
/// @param  pose: Pose (used for residue renumbering)
/// @param  sorted_data: Input lines sorted by linetype
/// @return Proper data container that can be used by Rosetta
metrics::DEERDataOP
DEERIO::parse_bounds_data(
	pose::Pose const & pose,
	std::map< std::string, SplitLines > const & sorted_data
) const {

	// Output type is known at this point
	metrics::DEERDataOP data( new metrics::DEERDistanceBounds() );

	// Get residues
	data->residues( parse_desc_lines( pose, sorted_data.at( "DESC" ) ) );

	// Process the actual distance data (data passed by reference)
	parse_bounds_lines( sorted_data.at( "BOUNDS" ), data );

	// Return output
	return data;
}

/// @brief  Parse inputs lines into DEERDistanceDistributionOP
/// @param  pose: Pose (used for residue renumbering)
/// @param  sorted_data: Input lines sorted by linetype
/// @return Proper data container that can be used by Rosetta
metrics::DEERDataOP
DEERIO::parse_dist_data(
	pose::Pose const & pose,
	std::map< std::string, SplitLines > const & sorted_data
) const {

	// Output type is NOT known at this point, since the way the data are
	//  scored is determined by information in the "INFO" lines.
	// Therefore we need to parse those first
	auto data = parse_dist_info_lines( sorted_data.at( "INFO" ) );

	// Get residues
	data->residues( parse_desc_lines( pose, sorted_data.at( "DESC" ) ) );

	// Process the actual distance data (data passed by reference)
	parse_dist_lines( sorted_data.at( "DIST" ), data );

	// Return output
	return data;
}

/// @brief Parses lines labeled "INFO" into information usable by Rosetta
/// @param  info_lines: Lines with the information (as strings)
/// @return Proper data container that can be used by Rosetta
/// @detail Note that this container, upon output, is incomplete!
metrics::DEERDataOP
DEERIO::parse_dist_info_lines(
	SplitLines const & info_lines
) const {

	// Determine the datatype
	metrics::DEERDataOP data = parse_dist_datatype( info_lines );

	// Set defaults
	data->bins_per_a( BINS_PER_A_ );
	data->stdev( 1.0 );

	// Spot check in case there are no lines
	if ( info_lines.size() == 0 ) {

		// Go through and get this information, if not provided
		for ( auto const & line : info_lines ) {
			for ( Size i = 3; i <= line.size() - ( line.size() % 2 ); i += 2 ) {
				if ( line.at( i ) == "BINS_PER_ANGSTROM" ) {
					data->bins_per_a( std::stoi( line.at( i + 1 ) ) );
				} else if ( line.at( i ) == "DISTR_STDEV" ) {
					data->stdev( std::stod( line.at( i + 1 ) ) );
				}
			}
		}
	}

	// Return the output
	return data;
}

/// @brief Parses lines labeled "DIST" into information usable by Rosetta
/// @param  dist_lines: Lines with the distance distributions
/// @param  data: DEERData object (passed by reference and modified here)
void
DEERIO::parse_dist_lines(
	SplitLines const & dist_lines,
	metrics::DEERDataOP data
) const {

	// Three parts of distance distribution: best fit and confidence bands
	std::map< Real, std::tuple< Real, Real, Real > > best_vals;

	// Total amplitude of best fit (for normalization)
	Real total = 0.0;

	// Iterate over lines in input files
	for ( auto const & line : dist_lines ) {

		// Label appropriate items
		Real const dist = std::stod( line[ 3 ] );
		Real const amp_best = std::stod( line[ 4 ] );
		total += amp_best;

		// If there are enough values in the line to read conf bands, do so
		Real const amp_lo = ( line.size() >= 6 ) ? std::stod( line[ 5 ] ) : 0.0;
		Real const amp_hi = ( line.size() >= 6 ) ? std::stod( line[ 6 ] ) : 0.0;

		// Set the tuple
		best_vals[ dist ] = std::make_tuple( amp_best,
			std::min( amp_lo, amp_hi ), std::max( amp_lo, amp_hi ) );
	}

	// Pull this information to determine how distances are binned
	Size const bins_per_a = data->bins_per_a();

	// The actual output maps to be passed to DEERDataOP object
	std::map< Size, Real > distr, lower, upper;

	Size i = 1;
	for ( auto const & d_amp : best_vals ) {
		Real amp = 0.0;
		Real lo = 0.0;
		Real hi = 0.0;
		std::tie( amp, lo, hi ) = d_amp.second;

		// We need to save the exact distance bins for later
		if ( bins_per_a == 0 ) {
			data->append_dist_id( i, d_amp.first );
			distr[ i ] = amp / total;
			lower[ i ] = lo / total;
			upper[ i ] = hi / total;

			// Otherwise, just read the distance bins as normal
		} else {
			Size const bin = round( bins_per_a * d_amp.first );
			distr[ bin ] = amp / total;
			lower[ bin ] = lo / total;
			upper[ bin ] = hi / total;
		}
		++i;
	}

	// Now assign these maps to the data
	metrics::DEERDistanceDistributionOP distr_data =
		utility::pointer::dynamic_pointer_cast<
		metrics::DEERDistanceDistribution >( data );

	distr_data->best_fit( distr );
	distr_data->lower_bound( lower );
	distr_data->upper_bound( upper );
}

/// @brief  Determines function used by DEERDistanceDistribution
/// @param  info_lines: Lines labeled "INFO" in input file
/// @return Proper data container that can be used by Rosetta
/// @detail Note that this container, upon output, is incomplete!
metrics::DEERDataOP
DEERIO::parse_dist_datatype(
	SplitLines const & info_lines
) const {

	// Set a default value
	std::string mode = "CROSS_ENTROPY_BOUNDS_BB";

	// Spot check to make sure there is even any content here
	if ( info_lines.size() > 0 ) {

		// Go through lines in map and find the string
		for ( auto const & line : info_lines ) {
			Size n_elements = line.size() - ( line.size() % 2 );
			for ( Size i = 3; i <= n_elements; i += 2 ) {
				if ( line.at( i ) == "MODE" ) {
					mode = line.at( i + 1 );
				}
			}
		}
	}

	// Initialize a default one (that can be replaced later)
	metrics::DEERDataOP data = metrics::DEERDataOP(
		new metrics::DEERDistanceDistribution() );

	// Use default function (cross entropy)
	if ( mode.find( "CROSS_ENTROPY" ) != std::string::npos ) {
		data = metrics::DEERDataOP(
			new metrics::DEERDistanceDistribution() );

		// Use overlap calculation
	} else if ( mode.find( "OVERLAP" ) != std::string::npos ) {
		data = metrics::DEERDataOP( new metrics::DEEROverlapMethod() );

		// Use Kolmogorov-Smirnov
	} else if ( mode.find( "KS" ) != std::string::npos ) {
		data = metrics::DEERDataOP( new metrics::DEEROverlapMethod() );

		// Need to set this to use the maximum distance between integrals
		metrics::DEERDistanceDistributionOP distr_data =
			utility::pointer::dynamic_pointer_cast<
			metrics::DEERDistanceDistribution >( data );
		distr_data->integral( true );
		distr_data->singleval( true );

		// Use maximum discrepancy between distributions
	} else if ( mode.find( "DISCREPANCY" ) != std::string::npos ) {
		data = metrics::DEERDataOP( new metrics::DEEROverlapMethod() );
		metrics::DEERDistanceDistributionOP distr_data =
			utility::pointer::dynamic_pointer_cast<
			metrics::DEERDistanceDistribution >( data );
		distr_data->singleval( true );

		// Use area between integrals
	} else if ( mode.find( "WASSERSTEIN" ) != std::string::npos ) {
		data = metrics::DEERDataOP( new metrics::DEEROverlapMethod() );
		metrics::DEERDistanceDistributionOP distr_data =
			utility::pointer::dynamic_pointer_cast<
			metrics::DEERDistanceDistribution >( data );
		distr_data->integral( true );
	} else if ( mode.find( "JENSEN-SHANNON" ) != std::string::npos
			// Use Jensen-shannon divergence
			|| mode.find( "JS" ) != std::string::npos
			) {
		data = metrics::DEERDataOP( new metrics::DEERJSMethod() );

		// Use Chi-squared
	} else if ( mode.find( "CHI-SQUARED" ) != std::string::npos
			|| mode.find( "CHISQ" ) != std::string::npos
			) {
		data = metrics::DEERDataOP( new metrics::DEERChiSqMethod() );

		// Use Jaccard index
	} else if ( mode.find( "JACCARD" ) != std::string::npos ) {
		data = metrics::DEERDataOP( new metrics::DEERJaccardMethod() );

		// Use one of the miscellaneous methods
	} else if ( mode.find( "HELLINGER" ) != std::string::npos
			|| mode.find( "BHATTACHARRYYA" ) != std::string::npos
			|| mode.find( "CONDITIONAL" ) != std::string::npos
			) {
		data = metrics::DEERDataOP( new metrics::DEERMiscMethod() );
		metrics::DEERMiscMethodOP misc_data =
			utility::pointer::dynamic_pointer_cast<
			metrics::DEERMiscMethod >( data );
		misc_data->mode( mode );

		// Error message in case none of the options are found
	} else {
		TR.Error << "Mode " << mode << " not known! ";
		TR.Error << "Viable options are: " << std::endl;
		TR.Error << "\tCROSS_ENTROPY" << std::endl;
		TR.Error << "\tOVERLAP" << std::endl;
		TR.Error << "\tKS" << std::endl;
		TR.Error << "\tDISCREPANCY" << std::endl;
		TR.Error << "\tWASSERSTEIN" << std::endl;
		TR.Error << "\tJENSEN-SHANNON (OR JS)" << std::endl;
		TR.Error << "\tCHI-SQUARED (OR CHISQ)" << std::endl;
		TR.Error << "\tJACCARD" << std::endl;
		TR.Error << "\tHELLINGER" << std::endl;
		TR.Error << "\tBHATTACHARYYA" << std::endl;
		TR.Error << "\tCONDITIONAL" << std::endl;
		throw CREATE_EXCEPTION( utility::excn::BadInput,
			"Non-viable mode for calculating DEER distribution scores." );
	}

	// Now set the remainder of the options
	metrics::DEERDistanceDistributionOP distr_data =
		utility::pointer::dynamic_pointer_cast<
		metrics::DEERDistanceDistribution >( data );

	// If the backbone needs to be broadened (by Cross Entropy method)
	if ( mode.find( "BB" ) != std::string::npos ) {
		distr_data->bb( true );
	}

	// If confidence bands are to be used
	if ( mode.find( "BOUNDS" ) != std::string::npos ) {
		distr_data->bounds( true );
	}

	// If the reverse comparison needs to be made (for asymmetric methods)
	if ( mode.find( "REV" ) != std::string::npos ) {
		distr_data->reverse( true );
	}

	// Return the data
	return data;
}

/// @brief  Parses data for DEERDistanceBounds object
/// @param  bounds_lines: Lines labeled "BOUNDS" in input file
/// @param  data: Proper data container that can be used by Rosetta
/// @detail Note that only the last BOUNDS line will be used
void
DEERIO::parse_bounds_lines(
	SplitLines const & bounds_lines,
	metrics::DEERDataOP bounds_data
) const {

	// Cast to correct type so we can save the data
	metrics::DEERDistanceBoundsOP data =
		utility::pointer::dynamic_pointer_cast<
		metrics::DEERDistanceBounds >( bounds_data );

	// Remind the user that only the last line that was read will be used
	if ( bounds_lines.size() > 1 ) {
		TR.Warning << "Only the last BOUNDS line is used!" << std::endl;
	}

	// Save the data
	for ( auto const & line : bounds_lines ) {
		data->bounds( std::stod( line[ 3 ] ), std::stod( line[ 4 ] ) );
		data->step( std::stod( line[ 5 ] ) );
	}
}

/// @brief  Parses information data for DEERDecayData object
/// @param  info_lines: Lines labeled "INFO" in input file
/// @param  data: Proper data container that can be used by Rosetta
void
DEERIO::parse_decay_info_lines(
	SplitLines const & info_lines,
	metrics::DEERDataOP data
) const {

	// Cast to DEERDecayData object to access the setter functions
	metrics::DEERDecayDataOP decay_data =
		utility::pointer::dynamic_pointer_cast<
		metrics::DEERDecayData >( data );

	data->bins_per_a( BINS_PER_A_ );
	decay_data->max_dist( ANGSTROM_LIMIT_ * BINS_PER_A_ );
	// Iterate through each line
	for ( auto const & line : info_lines ) {

		// Go through each option
		for ( Size i = 3; i <= line.size() - ( line.size() % 2 ); i += 2 ) {

			// If the standard deviation can be floated as a parameter
			if ( line[ i ] == "FIT_STDEV"
					&& ObjexxFCL::uppercased( line.at( i + 1 ) ) != "FALSE"
					) {
				decay_data->fit_stdev( true );

				// The background type ("3D", "NON-3D", "NONE")
			} else if ( line.at( i ) == "BCKG_TYPE"
					|| line.at( i ) == "BACKGROUND"
					) {
				decay_data->bckg( line.at( i + 1 ) );

				// The noise level
			} else if ( line.at( i ) == "NOISE" ) {
				decay_data->noise( std::stod( line.at( i + 1 ) ) );

				// The granularity of the distance distribution
			} else if ( line.at( i ) == "BINS_PER_ANGSTROM" ) {
				decay_data->bins_per_a( std::stoi( line.at( i + 1 ) ) );

				// The breadth of the distribution (overridden by FIT_STDEV)
			} else if ( line.at( i ) == "DISTR_STDEV" ) {
				data->stdev( std::stod( line.at( i + 1 ) ) );

				// The maximum distance for the kernel method
			} else if ( line.at( i ) == "MAX_DIST" ) {
				decay_data->max_dist( std::stoi( line.at( i + 1 ) ) );

				// What gets printed if the string is not found
			} else {
				TR.Warning << "Couldn't read " << line.at( i ) << std::endl;
				TR.Warning << "Skipping." << std::endl;
			}
		}
	}
}

/// @brief  Parses decay data for DEERDecayData object
/// @param  decay_lines: Lines labeled "DECAY" in input file
/// @param  data: Proper data container that can be used by Rosetta
void
DEERIO::parse_decay_lines(
	SplitLines const & decay_lines,
	metrics::DEERDataOP & data
) const {

	// Cast to DEERDecayData object to access the setter functions
	metrics::DEERDecayDataOP decay_data =
		utility::pointer::dynamic_pointer_cast<
		metrics::DEERDecayData >( data );

	// Now read the data one-by-one and append to vectors
	// Note that the initial set of v = 1.0 when t = 0.0 is to correctly
	// normalize the DEER trace in case data has no t = 0.0 point
	utility::vector1< Real > decay{ 1.0 };
	utility::vector1< Real > time_pts{ 0.0 };
	for ( auto const & splitline : decay_lines ) {
		time_pts.push_back( std::stod( splitline[ 3 ] ) );
		decay.push_back( std::stod( splitline[ 4 ] ) );
	}

	// Instantiate factory object passed to decay_data
	decay_data-> init_factory( decay, time_pts );
}

/// @brief  Parse input lines for information on residues involved in data
/// @param  pose: Pose (used for residue renumbering)
/// @param  desc_lines: Input lines containing this information
/// @return List of residues describing both the ID# and spin label type
utility::vector1< PairSizeString >
DEERIO::parse_desc_lines(
	pose::Pose const & pose,
	SplitLines const & desc_lines
) const {

	// Output type
	utility::vector1< PairSizeString > output;

	// Iterate through data one by one
	for ( auto const & line : desc_lines ) {

		// Number of residues expected based on size of vector
		Size const n_res = ( line.size() - 2 ) / 2;

		// Iterate through each one
		for ( Size i = 1; i <= n_res; ++i ) {

			// Get the correct amino acid ID for each
			pose::ResidueIndexDescriptionCOP const res = pose::parse_resnum(
				line[ 2 + ( 2 * i ) ] );
			Size const resnum = res->resolve_index( pose );

			// If this residue is not in the pose, skip
			if ( resnum > pose.size() ) {
				TR.Warning << "Skipping residue # " << resnum << std::endl;
			}

			// Otherwise, add it to the list
			output.push_back( std::make_pair( resnum, line[ 1 + ( 2 * i ) ] ) );
		}
	}

	// Send it back
	return output;
}

/// @brief  Parse coordinate files for custom spin labels
/// @param  pose: Pose (used for residue renumbering)
/// @return List of spin labels and corresponding weights
utility::vector1< std::pair< EPRSpinLabel, Real > >
DEERIO::pull_coords(
	pose::Pose const & pose
) const {

	// To make things look a little less crazy
	using namespace basic::options;

	// Defining output
	utility::vector1< std::pair< EPRSpinLabel, Real > > output;

	std::map< Size, Real > temp_weights;
	std::map< Size, std::map< Size, utility::vector1< PseudoSL > > > all_coords;

	// Only proceed if the right option was passed
	if ( option[ OptionKeys::epr_deer::coords_files ].user() ) {

		// Outer loop iterate through files
		utility::vector1< std::string > coords_files =
			option[ OptionKeys::epr_deer::coords_files ]();
		for ( auto const & file : coords_files ) {

			// Open the file, or try to
			std::ifstream infile;
			std::string line;
			infile.open( file.c_str(), std::ios::in );
			if ( !infile.is_open() ) {
				TR.Error << "Cannot open file (-epr_deer:coords_files "
					<< file << ")." << std::endl;
				TR.Error << "Skipping." << std::endl;
				continue;
			}

			// Inner loop: iterate through lines
			while ( std::getline( infile, line ) ) {

				// Split the line into substrings
				auto const splitline = split_line( line );

				// Skip if the string starts with "#" or is empty
				if ( splitline.size() == 0 || line[ 0 ] == '#' ) {
					continue;
				}

				// Read either DESC file or COORD file
				// DESC lines split the custom coordinates if more than
				// one is used and assign a weight (obtained using AICc,
				// for example)
				if ( splitline[ 1 ] == "DESC" ) {
					Size const idx = std::stoi( splitline[ 2 ] );
					temp_weights[ idx ] = std::stod( splitline[ 3 ] );

					// COORD lines describe residue-specific custom coordinates
				} else if ( splitline[ 1 ] == "COORD" ) {
					Size const idx =  std::stoi( splitline[ 2 ] );

					// Get the right residue number
					pose::ResidueIndexDescriptionCOP resnum =
						pose::parse_resnum( splitline[ 3 ] );
					Size const res = resnum->resolve_index( pose );
					if ( res > pose.size() ) {
						TR.Warning << "Skipping residue number " << res <<
							" line " << splitline << std::endl;
						break;
					}

					// Coordinate information
					Real const w = std::stod( splitline[ 4 ] );
					Real const x = std::stod( splitline[ 5 ] );
					Real const y = std::stod( splitline[ 6 ] );
					Real const z = std::stod( splitline[ 7 ] );

					// Make a new slot for this specific idx if not found
					if ( all_coords.find( idx ) == all_coords.end() ) {
						all_coords[ idx ] =
							std::map< Size, utility::vector1< PseudoSL > >();
					}

					// Make a new slot for this specific residue if not found
					append_to_map( all_coords[ idx ], res, std::make_pair(
						numeric::xyzVector< Real >( x, y, z ), w ) );
				}
			}
		}

		// Now normalize the probability distributions of these coords
		Real baseline = 0.0;
		for ( auto const & weight : temp_weights ) {
			baseline += weight.second;
		}

		// Check if baseline is zero to avoid divide-by-zero error
		if ( baseline == 0.0 ) {
			throw CREATE_EXCEPTION( utility::excn::BadInput,
				"Coord weights add up to zero (-epr_deer:coords_file)!" );
		}

		// Normalize
		for ( auto & weight : temp_weights ) {
			weight.second /= baseline;
		}

		// Rearrange everything into pairs
		for ( core::Size i = 1; i <= all_coords.size(); ++i ) {
			EPRSpinLabel temp_sl;
			temp_sl.load_custom_electrons( all_coords[ i ] );
			output.push_back( std::make_pair( temp_sl, temp_weights[ i ] ) );
		}
	}

	// Return this rearranged object (or empty, if nothing was read)
	return output;
}

} // namespace epr_deer
} // namespace scoring
} // namespace core
