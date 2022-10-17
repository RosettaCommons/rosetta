// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/DensityZscoresStatsSetup.cc
/// @brief  Electron density zscores database class
/// @author Gabriella Reggiano

// Unit headers
#include <core/scoring/DensityZscoresStatsSetup.hh>

//Project headers
#include <basic/database/open.hh>

//Utility headers
#include <utility/io/izstream.hh>
#include <utility/vector1.hh>
#include <utility/exit.hh>

#include <basic/Tracer.hh>

static basic::Tracer TR( "core.scoring.DensityZscoresStatsSetup" );

namespace core {
namespace scoring {


/// @brief constructor
DensityZscoresStatsSetup::DensityZscoresStatsSetup()
{
	set_statistics_from_database();
}

/// @brief deconstructor
DensityZscoresStatsSetup::~DensityZscoresStatsSetup() = default;

/// @brief return window 1 or window 3 zscores given an amino acid identity and a bfactor
core::Real
DensityZscoresStatsSetup::eval_zscore_by_resn_bfactor(
	core::Real const input_val,
	std::string const & resn,
	core::Real const bfactor,
	core::Size const window_size ) const
{
	core::Real zscore = 0;
	core::Size bfactor_bin_num = 0;
	for ( core::Size i=1; i<= bfac_bins_start_end_.size(); i++ ) {
		if ( std::get<0>( bfac_bins_start_end_[i] ) <= bfactor
				&& std::get<1>( bfac_bins_start_end_[i] ) > bfactor ) {
			bfactor_bin_num = i;
			break;
		}
	}
	if ( bfactor_bin_num == 0 ) { bfactor_bin_num = bfac_bins_start_end_.size(); } // TODO - should remove and have this fail here instead?

	if ( window_size == 3 ) {
		zscore = calc_weighted_zscore( bfactor_bin_num, input_val, resn, win3res_db_means_, win3res_db_stdevs_);
	} else if ( window_size == 1 ) {
		zscore = calc_weighted_zscore( bfactor_bin_num, input_val, resn, win1res_db_means_, win1res_db_stdevs_);
	} else {
		std::stringstream message;
		message << "window size " << window_size << " not supported. only 1 or 3 are valid" << std::endl;
		utility_exit_with_message(message.str());
	}
	return zscore;

}

////////////////////////////////////////////////////////////////////////////////////////////////
/// private methods ///
//////////////////////
typedef std::map < std::string, std::map < core::Real, core::Real > > DensZScoreMap;

DensZScoreMap
DensityZscoresStatsSetup::read_database_file(std::string const & dbfile)
{
	utility::io::izstream db_stream;
	db_stream.open( basic::database::full_name(dbfile));
	std::string line;
	DensZScoreMap stats;
	db_stream.getline(line);
	std::string substr;
	while ( db_stream >> line ) {
		std::stringstream sline (line);
		utility::vector1 < std::string > parsed_line;
		while ( sline.good() ) {
			std::getline( sline, substr, ',' );
			parsed_line.push_back(substr);
		}
		// 1 is the residue 3 leter abv. 2 is first instance of stats
		for ( core::Size i=2; i<=parsed_line.size(); i++ ) {
			stats[parsed_line[1]][std::get<0>( bfac_bins_start_end_[i-1] )] = (core::Real) std::stod(parsed_line[i]);
		}
	}
	db_stream.close();
	return stats;
}

/// @brief set headers from database files
void
DensityZscoresStatsSetup::set_db_headers( std::string const & dbfile )
{
	utility::io::izstream db_stream;
	db_stream.open( basic::database::full_name(dbfile) );
	std::string line;
	db_stream.getline(line);
	std::string sub_header;
	std::stringstream sline(line);
	// TODO - don't push back empty stuff at beginning
	while ( std::getline(sline, sub_header, ',') ) {
		if ( !sub_header.empty() ) {
			core::Real n1 = (core::Real) std::stod ( sub_header.substr(0, sub_header.find(header_delim_)) );
			core::Real n2 = (core::Real) std::stod ( sub_header.substr(sub_header.find(header_delim_)+header_delim_.length(), sub_header.length()) );
			bfac_bins_start_end_.push_back( std::make_tuple( n1, n2 ) );
		}
	}
	db_stream.close();
}

///@brief read database files and store data in maps
void
DensityZscoresStatsSetup::set_statistics_from_database()
{

	set_db_headers( WIN3RES_MEANS_DBFILE_ );
	win3res_db_means_ = read_database_file( WIN3RES_MEANS_DBFILE_ );
	win1res_db_means_ = read_database_file( WIN1RES_MEANS_DBFILE_ );
	win3res_db_stdevs_ = read_database_file( WIN3RES_STDEVS_DBFILE_ );
	win1res_db_stdevs_ = read_database_file( WIN1RES_STDEVS_DBFILE_ );

}

/// @brief calculate weighted z-score
core::Real
DensityZscoresStatsSetup::calc_weighted_zscore(
	core::Size const bfbin,
	core::Real const input_val,
	std::string const & resn,
	DensZScoreMap const & means,
	DensZScoreMap const & stdevs) const
{
	// gaussian - 1 kernel radius, std = 0.5
	utility::vector1< core::Real > gauss_smooth(3);
	gauss_smooth[1] = 0.157731160;
	gauss_smooth[2] = 0.684537679;
	gauss_smooth[3] = 0.157731160;

	// bin keys for accessing from map
	core::Real current_bin_key = std::get<0>(bfac_bins_start_end_[bfbin]);
	core::Real prev_bin_key;
	core::Real next_bin_key;

	// do the gaussian weighted_mean and std
	core::Real weighted_mean;
	core::Real weighted_stdev;
	if ( bfbin == 1 ) { // 85% n, 15% n+1
		next_bin_key = std::get<0>(bfac_bins_start_end_[bfbin+1]);

		weighted_mean = means.at(resn).at(current_bin_key)*gauss_smooth[1] +
			means.at(resn).at(current_bin_key)*gauss_smooth[2] +
			means.at(resn).at(next_bin_key)*gauss_smooth[3];

		weighted_stdev = stdevs.at(resn).at(current_bin_key)*gauss_smooth[1] +
			stdevs.at(resn).at(current_bin_key)*gauss_smooth[2] +
			stdevs.at(resn).at(next_bin_key)*gauss_smooth[3];
	} else if ( bfbin == bfac_bins_start_end_.size() ) { //15% n-1, 85% n
		prev_bin_key = std::get<0>(bfac_bins_start_end_[bfbin-1]);

		weighted_mean = means.at(resn).at(prev_bin_key)*gauss_smooth[1] +
			means.at(resn).at(current_bin_key)*gauss_smooth[2] +
			means.at(resn).at(current_bin_key)*gauss_smooth[3];

		weighted_stdev = stdevs.at(resn).at(prev_bin_key)*gauss_smooth[1] +
			stdevs.at(resn).at(current_bin_key)*gauss_smooth[2] +
			stdevs.at(resn).at(current_bin_key)*gauss_smooth[3];
	} else { //15% n-1, 70% n, 15% n+1
		prev_bin_key = std::get<0>(bfac_bins_start_end_[bfbin-1]);
		next_bin_key = std::get<0>(bfac_bins_start_end_[bfbin+1]);

		weighted_mean = means.at(resn).at(prev_bin_key)*gauss_smooth[1] +
			means.at(resn).at(current_bin_key)*gauss_smooth[2] +
			means.at(resn).at(next_bin_key)*gauss_smooth[3];

		weighted_stdev = stdevs.at(resn).at(prev_bin_key)*gauss_smooth[1] +
			stdevs.at(resn).at(current_bin_key)*gauss_smooth[2] +
			stdevs.at(resn).at(next_bin_key)*gauss_smooth[3];
	}
	return ( input_val - weighted_mean )/weighted_stdev;
}

} // scoring
} // core
