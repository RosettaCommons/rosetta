// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/thermal_sampling/util.cc
/// @brief
/// @details
///
/// @author Andy Watkins


#include <core/types.hh>
#include <utility/io/ozstream.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/pose/Pose.hh>

// Utility headers

// ObjexxFCL headers

//// C++ headers
#include <string>
#include <map>
#include <cmath>
#include <utility/vector1.hh>


namespace protocols {
namespace thermal_sampling {

//////////////////////////////////////////////////////////////////////////////
// Histogram class for accumulating samples
class Histogram {
public:
	Histogram( core::Real const min, core::Real const max, core::Real const spacing ):
		min_( min ),
		max_( max ),
		spacing_( spacing )
	{
		runtime_assert( max > min );
		n_elem_ = static_cast<core::Size>( ( max - min ) / spacing ) + 1;
		for ( core::Size i = 1; i <= n_elem_; ++i ) hist_.push_back( 0 );
	}

	void add( float const value, core::Size const n_items ) {
		core::Size bin_index;
		if ( value <= min_ ) {
			bin_index = 1;
		} else if ( value >= max_ ) {
			bin_index = n_elem_;
		} else {
			bin_index = static_cast<core::Size>( ( value - min_ ) / spacing_ ) + 1;
		}
		hist_[bin_index] += n_items;
	}

	void clear() {
		for ( core::Size i = 0; i <= n_elem_; ++i ) hist_[i] = 0;
	}

	utility::vector1<core::Real> get_scores() const {
		utility::vector1<core::Real> scores;
		for ( core::Size i = 1; i <= n_elem_; ++i ) {
			scores.push_back( min_ + spacing_ * ( i - 0.5 ) );
		}
		return scores;
	}

	utility::vector1<core::Size> get_hist() const { return hist_; }

private:
	core::Real const min_, max_, spacing_;
	core::Size n_elem_;
	utility::vector1<core::Size> hist_;
};

// score types to be recorded
utility::vector1<core::scoring::ScoreType> const & get_scoretypes() {
	using namespace core::scoring;
	static utility::vector1< ScoreType > scoretypes;
	if ( !scoretypes.empty() ) return scoretypes;

	// List of score types to be cached
	scoretypes.push_back( fa_atr );
	scoretypes.push_back( fa_rep );
	scoretypes.push_back( fa_intra_rep );
	scoretypes.push_back( fa_stack );
	scoretypes.push_back( rna_torsion );
	scoretypes.push_back( hbond_sc );
	scoretypes.push_back( lk_nonpolar );
	scoretypes.push_back( geom_sol_fast );
	scoretypes.push_back( stack_elec );
	scoretypes.push_back( fa_elec_rna_phos_phos );
	return scoretypes;
}

////////////////////////////////////////////////////////////////////////////////
core::Size data_dim() {
	using namespace core::scoring;
	utility::vector1<ScoreType> const & score_types( get_scoretypes() );
	return score_types.size() + 2;
}

//////////////////////////////////////////////////////////////////////////////
void update_scores(
	utility::vector1<float> & scores,
	core::pose::Pose & pose,
	core::scoring::ScoreFunctionCOP scorefxn
) {
	using namespace core::scoring;
	scores.clear();
	scores.push_back( ( *scorefxn )( pose ) );
	utility::vector1< ScoreType > const & score_types = get_scoretypes();
	for ( core::Size i = 1; i<= score_types.size(); ++i ) {
		scores.push_back( scorefxn->score_by_scoretype( pose, score_types[i], false /*weighted*/ ) );
	}
}

//////////////////////////////////////////////////////////////////////////////
void fill_data(
	utility::vector1<float> & data,
	core::Size const count,
	utility::vector1<float> const & scores
) {
	data.push_back( count );
	data.insert( data.end(), scores.begin(), scores.end() );
}

//////////////////////////////////////////////////////////////////////////////
// Simple heuristic for gaussian stdev
core::Real gaussian_stdev( core::Real const n_rsd, core::Real const temp, bool const is_bp ) {
	// Negative temp is infinite
	if ( temp < 0 ) return -1;
	if ( is_bp ) return 5 * temp / n_rsd;
	return 6 * std::pow( temp / n_rsd, 0.75 );
}

} //thermal_sampling
} //protocols
