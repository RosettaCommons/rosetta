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
#include <core/scoring/epr_deer/DEERData.hh>
#include <core/scoring/epr_deer/DEERIO.hh>
#include <core/scoring/epr_deer/EPRSpinLabel.hh>
#include <core/types.hh>
#include <numeric/xyzVector.hh>
#include <numeric/constants.hh>
#include <utility/io/izstream.hh>
#include <utility/fixedsizearray1.hh>

#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/epr_deer.OptionKeys.gen.hh>
#include <basic/options/option_macros.hh>
#include <basic/database/open.hh>

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

#include <boost/algorithm/string.hpp>

namespace core {
namespace scoring {
namespace epr_deer {

static basic::Tracer TR( "core.scoring.epr_deer.DEERIO" );

/// @brief  Constructor
DEERIO::DEERIO() {}

/// @brief  Destructor
DEERIO::~DEERIO() {}

/// @brief Generates the data for use. Core function
std::map< Size, DEERDataOP >
DEERIO::generate_data(
	core::pose::Pose const & pose
) {
	std::map< Size, DEERDataOP > output;

	// Get data from file. Tab-separated
	utility::vector1< utility::vector1< std::string > > splitlines = get_splitlines();

	// Arrangement of data: key is type of line, value is a list of splitlines
	std::map< std::string, utility::vector1< utility::vector1< std::string > > > lines_map;

	// Iterate over lines and sort
	for ( auto const & splitline : splitlines ) {
		std::string linetype = splitline[ 1 ];
		Size edge_id = std::stoi( splitline[ 2 ] );
		if ( lines_map.find( linetype ) == lines_map.end() ) {
			lines_map[ linetype ] = utility::vector1< utility::vector1< std::string > >();
		}
		lines_map[ linetype ].push_back( splitline );

		// Make the array of objects. Since they are derived classes we need to do this first
		// Since PAIR is a one-line summary, it comes last
		if ( output.find( edge_id ) == output.end() ) {
			if ( linetype == "BOUNDS" ) {
				output[ edge_id ] = DEERDataOP( new DEERDistanceBounds() );
			} else if ( linetype == "DIST" || linetype == "GAUSS" ) {
				output[ edge_id ] = DEERDataOP( new DEERDistanceDistribution() );
			} else if ( linetype == "DECAY" ) {
				output[ edge_id ] = DEERDataOP( new DEERDecayData() );
			}
		}
	}

	// Now parse the input for real
	read_desc_lines(   output, lines_map[ "DESC"   ], pose );
	read_info_lines(   output, lines_map[ "INFO"   ]      );
	read_decay_lines(  output, lines_map[ "DECAY"  ]     );
	read_bounds_lines( output, lines_map[ "BOUNDS" ]      );
	read_dist_lines(   output, lines_map[ "DIST"   ]      );
	read_gauss_lines(  output, lines_map[ "GAUSS"  ]     );

	// Must be last
	read_pair_lines(   output, lines_map[ "PAIR"   ], pose );

	return output;
}

/// @brief Reads the input file(s) and makes an unsorted vector of their whitespace-separated contents
utility::vector1< utility::vector1< std::string > >
DEERIO::get_splitlines() {
	utility::vector1< std::string > lines;
	utility::vector1< std::string > input_files =
		basic::options::option[ basic::options::OptionKeys::epr_deer::input_files ]();
	// Outer loop: files listed as inputs
	for ( auto const & input_file : input_files ) {
		std::ifstream infile;
		std::string line;
		infile.open( input_file.c_str(), std::ios::in );
		if ( !infile.is_open() ) {
			utility_exit_with_message( "Unable to open DEER input file." );
		}

		// Inner loop: lines in each input file
		while ( std::getline( infile, line ) ) {
			utility::trim( line, " \t\n" );
			if ( ( line.size() >= 1 ) && ( line[ 0 ] != '#' ) ) {
				lines.push_back( line );
			}
		}
	}
	utility::vector1< utility::vector1< std::string > > output;

	// Iterate over lines
	for ( auto line : lines ) {
		utility::vector1< std::string > splitline;
		utility::vector1< std::string > splitline_out;
		boost::algorithm::split( splitline, line, boost::is_any_of( "\t " ) );

		// Check if any of these are empty, which can happen if there are two consecutive spaces
		for ( Size i = 1; i <= splitline.size(); ++i ) {
			if ( splitline[ i ].size() > 0 ) {
				splitline_out.push_back( splitline[ i ] );
			}
		}
		output.push_back( splitline_out );
	}
	return output;
}

/// @brief Reads lines that start with DESC, which contain residue and spin label information
void
DEERIO::read_desc_lines(
	std::map< Size, DEERDataOP > & output,
	utility::vector1< utility::vector1< std::string > > const & splitlines,
	pose::Pose const & pose
) {
	for ( auto const & splitline : splitlines ) {
		utility::vector1< std::pair< Size, std::string > > res_label_map;
		if ( ( splitline.size() - 2 ) % 2 != 0 ) {
			TR.Error << "Each residue must be listed alongside a spin label identity! Example:" << std::endl;
			TR.Error << "'DESC' Experiment_ID Rotamer_res1 Res1 Rotamer_res2 Res2 ... Rotamer_resN ResN" << std::endl;
			TR.Error << "\tDESC 1 DEFAULT 059 DEFAULT 159" << std::endl;
		}
		Size n_res = ( splitline.size() - 2 ) / 2;
		for ( Size i = 1; i <= n_res; ++i ) {
			Size res = pose::parse_resnum( splitline[ 2 + ( 2 * i ) ], pose );
			res_label_map.push_back( std::make_pair( res, splitline[ 1 + ( 2 * i ) ] ) );
		}
		Size edge_id = std::stoi( splitline[ 2 ] );
		if ( output.find( edge_id ) == output.end() ) {
			TR.Error << "Could not find edge_id " << edge_id << " in data!" << std::endl;
			TR.Error << "Skipping." << std::endl;
			continue;
		} else {
			output[ edge_id ]->bins_per_angstrom( BINS_PER_ANGSTROM_ );
			output[ edge_id ]->residues( res_label_map );
		}
	}
}

/// @brief Reads lines that start with INFO, which gives additional options specific for each bit of data
void
DEERIO::read_info_lines(
	std::map< Size, DEERDataOP > & output,
	utility::vector1< utility::vector1< std::string > > const & splitlines
) {
	for ( auto const & splitline : splitlines ) {
		Size data_n = std::stoi( splitline[ 2 ] );
		bool decay_option = false;
		Size n_options = ( splitline.size() - 2 ) / 2;

		// First check generic options:
		for ( Size option = 1; option <= n_options; ++option ) {
			std::string key = splitline[ option * 2 + 1 ];
			std::string value = splitline[ option * 2 + 2 ];
			if ( key == "RELATIVE_WEIGHT" ) {
				output[ data_n ]->relative_weight( std::stod( value ) );
				TR.Info << "Relative weight of data " << data_n << ": " << output[ data_n ]->relative_weight() << std::endl;
			} else if ( key == "BINS_PER_ANGSTROM" ) {
				output[ data_n ]->bins_per_angstrom( std::stoi( value ) );
				TR.Info << "Bins per angstrom for data " << data_n << ": " << output[ data_n ]->bins_per_angstrom() << std::endl;
			} else {
				decay_option = true;
			}
		}

		// All other options are specific for decay traces, so that's why this is here
		if ( decay_option ) {
			DEERDecayDataOP data = utility::pointer::dynamic_pointer_cast< DEERDecayData >( output[ data_n ] );
			for ( Size option = 1; option <= n_options; ++option ) {
				std::string key = splitline[ option * 2 + 1 ];
				std::string value = splitline[ option * 2 + 2 ];
				if ( key == "MOD_DEPTH_MAX" ) {
					data->mod_depth_bounds( data->mod_depth_bounds().first, std::stod( value ) );
					TR.Info << "Lower modulation depth bounds for data " << data_n << ": " << data->mod_depth_bounds().second << std::endl;
				} else if ( key == "MOD_DEPTH_MIN" ) {
					data->mod_depth_bounds( std::stod( value ), data->mod_depth_bounds().second );
					TR.Info << "Upper modulation depth bounds for data " << data_n << ": " << data->mod_depth_bounds().first << std::endl;
				} else if ( key == "FIT_STDEV" && value != "FALSE" ) {
					data->fit_stdev( true );
					TR.Info << "Setting standard deviation as a fitting parameter for data " << data_n << std::endl;
				} else if ( key == "NOISE" ) {
					data->noise( std::stod( value ) );
					TR.Info << "Setting noise value of " << data->noise() << " for data " << data_n << std::endl;
				} else if ( key == "BACKGROUND" || key == "BCKG" ) {
					TR.Info << "Reading background function of " << value << " for data " << data_n << std::endl;
					if ( value == "NON_3D" || value == "NON3D" ) {
						data->bckg( value );
					} else if ( value == "NONE" ) {
						data->bckg( value );
					} else {
						//TR.Error << "Unable to parse background value for " << data_n << std::endl;
						//TR.Error << "Skipping." << std::endl;
						continue;
					}
				} else {
					continue;
				}
			}
		}
	}
}

/// @brief Read lines for decay data (raw data)
void
DEERIO::read_decay_lines(
	std::map< Size, DEERDataOP > & output,
	utility::vector1< utility::vector1< std::string > > const & data_lines
) {
	for ( auto const & splitline : data_lines ) {
		Size data_n = std::stoi( splitline[ 2 ] );
		DEERDecayDataOP data = utility::pointer::dynamic_pointer_cast< DEERDecayData >( output[ data_n ] );
		if ( output.find( data_n ) == output.end() ) {
			TR.Error << "Could not find data indexed to " << data_n << " in DEER data file! Skipping." << std::endl;
		} else {
			Real time = Real( std::stod( splitline[ 3 ] ) );
			Real sign = Real( std::stod( splitline[ 4 ] ) );
			data->append_trace_data_and_calculate( time, sign, ANGSTROM_LIMIT_ );
		}
	}
}

/// @brief Read lines for bounded distance restraints
void
DEERIO::read_bounds_lines(
	std::map< Size, DEERDataOP > & output,
	utility::vector1< utility::vector1< std::string > > const & splitlines
) {
	for ( auto const & splitline : splitlines ) {
		Size data_n = std::stoi( splitline[ 2 ] );
		Real lower  = std::stod( splitline[ 3 ] );
		Real upper  = std::stod( splitline[ 4 ] );
		Real step  = std::stod( splitline[ 5 ] );
		DEERDistanceBoundsOP data = utility::pointer::dynamic_pointer_cast< DEERDistanceBounds >( output[ data_n ] );
		if ( data->bounds() != std::make_pair( 0.0, 0.0 ) ) {
			TR.Warning << "BOUNDS for " << data_n << " previously declared and are being overwritten!" << std::endl;
		}
		data->bounds( lower, upper );
		if ( data->step() != 1.0 ) {
			TR.Warning << "BOUNDS for " << data_n << " previously declared and are being overwritten!" << std::endl;
		}
		data->step( step );
	}
}

/// @brief Read lines for gaussian-distributed distance data
void
DEERIO::read_gauss_lines(
	std::map< Size, DEERDataOP > & output,
	utility::vector1< utility::vector1< std::string > > const & splitlines
) {
	EPRSpinLabel sl_for_gauss;
	std::map< Size, std::map< Size, Real > > gauss_map;
	for ( auto const & splitline : splitlines ) {
		Size data_n = std::stoi( splitline[ 2 ] );
		if ( gauss_map.find( data_n ) == gauss_map.end() ) gauss_map[ data_n ] = std::map< Size, Real >();
		if ( output.find( data_n )  == output.end()  ) {
			TR.Error << "Could not find data indexed to " << data_n << " in DEER data file! Skipping." << std::endl;
		} else {
			Real mean  = Real( std::stod( splitline[ 3 ] ) );
			Real stdev  = Real( std::stod( splitline[ 4 ] ) );
			Real amp   = Real( std::stod( splitline[ 5 ] ) );
			Size start_bin = std::max( 1.0, BINS_PER_ANGSTROM_ * ( mean - ( 4 * stdev ) ) );
			Size end_bin   = std::min( Real( BINS_PER_ANGSTROM_ * ANGSTROM_LIMIT_ ), BINS_PER_ANGSTROM_ * ( mean + ( 4 * stdev ) ) );
			for ( Size bin = start_bin; bin <= end_bin; ++bin ) {
				if ( gauss_map[ data_n ].find( bin ) == gauss_map[ data_n ].end() ) gauss_map[ data_n ][ bin ] = 0.0;
				gauss_map[ data_n ][ bin ] += sl_for_gauss.gauss( bin / Real( BINS_PER_ANGSTROM_ ), mean, stdev ) * amp;
			}
		}
	}
	for ( auto & distr : gauss_map ) {
		normalize_distribution( distr.second );
		DEERDistanceDistributionOP data = utility::pointer::dynamic_pointer_cast< DEERDistanceDistribution >( output[ distr.first ] );
		data->best_fit( distr.second );
		data->lower_bound( distr.second );
		data->upper_bound( distr.second );
	}
}

/// @brief Read lines for non-gaussian-distributed distance data
void
DEERIO::read_dist_lines(
	std::map< Size, DEERDataOP > & output,
	utility::vector1< utility::vector1< std::string > > const & splitlines
) {
	// Set containers: key is data ID, value is a map of distances to amplitudes
	// Distances are multiplied by DEERDataOP->bins_per_angstrom() and rounded to integers
	std::map< Size, std::map< Size, Real > > deer_map;
	std::map< Size, std::map< Size, Real > > lower_map;
	std::map< Size, std::map< Size, Real > > upper_map;
	// Iterate over lines
	for ( auto const & splitline : splitlines ) {
		Size data_n = std::stoi( splitline[ 2 ] );
		Size dist = round( BINS_PER_ANGSTROM_ * Real( std::stod( splitline[ 3 ] ) ) );
		if ( output.find( data_n ) == output.end() ) {
			TR.Error << "Could not find data indexed to " << data_n << " in DEER data file! Skipping." << std::endl;
		} else {
			if ( deer_map.find( data_n ) == deer_map.end() ) {
				deer_map[ data_n ] = std::map< Size, Real >();
			}
			deer_map[ data_n ][ dist ] = Real( std::stod( splitline[ 4 ] ) );

			// If the vector is 4 long, the confidence bands are not provided.
			// If it is six long, the lower and upper bands are assumed to be provided
			if ( splitline.size() == 6 ) {
				Real temp_min = Real( std::max( 0.0, std::stod( splitline[ 5 ] ) ) );
				Real temp_max = Real( std::max( 0.0, std::stod( splitline[ 6 ] ) ) );
				if ( lower_map.find( data_n ) == lower_map.end() ) {
					lower_map[ data_n ] = std::map< Size, Real >();
				}
				lower_map[ data_n ][ dist ] = std::min( temp_min, temp_max );
				if ( upper_map.find( data_n ) == upper_map.end() ) {
					upper_map[ data_n ] = std::map< Size, Real >();
				}
				upper_map[ data_n ][ dist ] = std::max( temp_min, temp_max );
			} else {
				lower_map[ data_n ][ dist ] = Real( std::stod( splitline[ 4 ] ) );
				upper_map[ data_n ][ dist ] = Real( std::stod( splitline[ 4 ] ) );
			}
		}
	}
	// Now organize these
	for ( auto & distribution : deer_map ) {
		Real baseline = normalize_distribution( distribution.second );
		DEERDistanceDistributionOP data = utility::pointer::dynamic_pointer_cast< DEERDistanceDistribution >( output[ distribution.first ] );
		if ( lower_map.find( distribution.first ) != lower_map.end() && upper_map.find( distribution.first ) != upper_map.end() ) {
			for ( auto & val : lower_map[ distribution.first ] ) {
				val.second /= baseline;
			}
			data->lower_bound( lower_map[ distribution.first ] );
			for ( auto & val : upper_map[ distribution.first ] ) {
				val.second /= baseline;
			}
			data->upper_bound( upper_map[ distribution.first ] );
		}
		data->best_fit( distribution.second );
	}
}

/// @brief Read the lines starting with PAIR, which are for easy input
void
DEERIO::read_pair_lines(
	std::map< Size, DEERDataOP > & output,
	utility::vector1< utility::vector1< std::string > > const & splitlines,
	pose::Pose const & pose
) {
	// These are arranged like this:
	// PAIR RES1 RES2 DISTANEC
	Size n_start = output.size();
	for ( core::Size i = 1; i <= splitlines.size(); ++i ) {
		auto res1 = std::make_pair( pose::parse_resnum( splitlines[ i ][ 2 ], pose ), "DEFAULT" );
		auto res2 = std::make_pair( pose::parse_resnum( splitlines[ i ][ 3 ], pose ), "DEFAULT" );
		Real distance = std::stod( splitlines[ i ][ 4 ] );
		output[ n_start + i ] = DEERDataOP( new DEERDistanceBounds() );
		DEERDistanceBoundsOP new_data = utility::pointer::dynamic_pointer_cast< DEERDistanceBounds >( output[ n_start + i ] );
		new_data->residues( { res1, res2 } );
		new_data->bins_per_angstrom( BINS_PER_ANGSTROM_ );
		new_data->bounds( distance, distance );
		new_data->step( 1.0 );
	}
}

/// @brief Normalize distance distributions so their integrals equal 1
Real
DEERIO::normalize_distribution( std::map< Size, Real > & in_map ) {
	Real baseline( 0.0 );
	for ( auto & item : in_map ) {
		baseline += item.second;
	}
	for ( auto & item : in_map ) {
		item.second /= baseline;
	}
	return baseline;
}

/// @brief Get custom spin label coordinates
utility::vector1< std::pair< EPRSpinLabel, Real > >
DEERIO::pull_coords() const {
	utility::vector1< std::pair< EPRSpinLabel, Real > > output;
	if ( basic::options::option[ basic::options::OptionKeys::epr_deer::coords_files ].user() ) {
		utility::vector1< std::string > coords_files = basic::options::option[ basic::options::OptionKeys::epr_deer::coords_files ]();
		std::map< Size, Real > temp_weights;
		std::map< Size, std::map< Size, utility::vector1< PseudoElectron > > > all_coords;

		// Outer loop iterate through files
		for ( auto const & file : coords_files ) {
			std::ifstream infile;
			std::string line;
			infile.open( file.c_str(), std::ios::in );
			if ( !infile.is_open() ) {
				utility_exit_with_message( "Unable to open COORD file (-epr_deer:coords_files " + file + ")." );
			}
			// Inner loop: iterate through lines
			while ( std::getline( infile, line ) ) {
				utility::trim( line, " \t\n" );
				if ( ( line.size() >= 1 ) && ( line[ 0 ] != '#' ) ) {
					std::istringstream linestream(line);
					utility::vector1< std::string > splitline_raw, splitline;
					// Split into a vector1 of strings
					boost::algorithm::split( splitline_raw, line, boost::is_any_of( "\t " ) );

					// Remove empty contents that occur from two consecutive whitespaces
					for ( Size i = 1; i <= splitline_raw.size(); ++i ) {
						if ( splitline_raw[ i ].size() > 0 ) {
							splitline.push_back( splitline_raw[ i ] );
						}
					}

					// Files are split between a DESC file which provides a weight and COORD files with coords
					if ( splitline[ 1 ] == "DESC" ) {
						temp_weights[ std::stoi( splitline[ 2 ] ) ] = std::stod( splitline[ 3 ] );
					} else if ( splitline[ 1 ] == "COORD" ) {
						Size id =  std::stoi( splitline[ 2 ] );
						Size res = std::stoi( splitline[ 3 ] );
						Real w =   std::stod( splitline[ 4 ] );
						Real x =   std::stod( splitline[ 5 ] );
						Real y =   std::stod( splitline[ 6 ] );
						Real z =   std::stod( splitline[ 7 ] );
						if ( all_coords.find( id ) == all_coords.end() ) {
							all_coords[ id ] = std::map< Size, utility::vector1< PseudoElectron > >();
						}
						if ( all_coords[ id ].find( res ) == all_coords[ id ].end() ) {
							all_coords[ id ][ res ] = utility::vector1< PseudoElectron >();
						}
						all_coords[ id ][ res ].push_back( std::make_pair( numeric::xyzVector< Real >( x, y, z ), w ) );
					}
				}
			}
		}
		// Normalize coordinates for total
		Real baseline = 0.0;
		for ( auto const & weight : temp_weights ) {
			baseline += weight.second;
		}
		for ( auto & weight : temp_weights ) {
			weight.second /= baseline;
		}
		// Return output
		for ( core::Size i = 1; i <= all_coords.size(); ++i ) {
			EPRSpinLabel temp_sl;
			temp_sl.load_custom_electrons( all_coords[ i ] );
			output.push_back( std::make_pair( temp_sl, temp_weights[ i ] ) );
		}
	}
	return output;
}

} // namespace epr_deer
} // namespace scoring
} // namespace core
