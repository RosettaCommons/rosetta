// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/DensityZscoresStatsSetup.hh
/// @brief  Electron density zscores database class
/// @author Gabriella Reggiano

#ifndef INCLUDED_core_scoring_DensityZscoresStatsSetup_hh
#define INCLUDED_core_scoring_DensityZscoresStatsSetup_hh

// Unit Headers
#include <core/scoring/DensityZscoresStatsSetup.fwd.hh>

// Project headers
#include <core/types.hh>

// Utility Headers
#include <utility/VirtualBase.hh>
#include <utility/vector1.hh>

// C++ headers
#include <map>
#include <math.h>
#include <string>

namespace core {
namespace scoring {

class DensityZscoresStatsSetup : public utility::VirtualBase {

	typedef std::map < std::string, std::map < core::Real, core::Real > > DensZScoreMap;
	typedef utility::vector1< std::tuple< core::Real, core::Real > > DensZscoreBfacBins;

public:

	/// @brief default constructor
	DensityZscoresStatsSetup();

	/// @brief default destructor
	~DensityZscoresStatsSetup() override;

	core::Real
	eval_zscore_by_resn_bfactor(core::Real const  input_val,
		std::string const & resn,
		core::Real const bfactor,
		core::Size const window_size) const;

private: // private methods

	DensZScoreMap
	read_database_file( std::string const & dbfile);

	void
	set_db_headers( std::string const & dbfile );

	void
	set_statistics_from_database();

	core::Real
	calc_weighted_zscore( core::Size const bfbin,
		core::Real const input_val,
		std::string const & resn,
		DensZScoreMap const & means,
		DensZScoreMap const & stdevs) const;

private: // private member variables

	std::string const WIN3RES_MEANS_DBFILE_ = "scoring/electron_density/lgwin_resCC_mean.csv";
	std::string const WIN1RES_MEANS_DBFILE_ = "scoring/electron_density/smwin_resCC_mean.csv";
	std::string const WIN3RES_STDEVS_DBFILE_ = "scoring/electron_density/lgwin_resCC_std.csv";
	std::string const WIN1RES_STDEVS_DBFILE_ = "scoring/electron_density/smwin_resCC_std.csv";

	// containers for statistics from dbfiles
	DensZscoreBfacBins bfac_bins_start_end_;
	std::string const header_delim_ = " to ";
	DensZScoreMap win3res_db_means_;
	DensZScoreMap win1res_db_means_;
	DensZScoreMap win3res_db_stdevs_;
	DensZScoreMap win1res_db_stdevs_;


}; // DensityZscores

} // scoring
} // core


#endif
