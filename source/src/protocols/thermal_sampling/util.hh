// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/thermal_sampling/util.hh
/// @brief
/// @details
///
/// @author Andy Watkins


#ifndef INCLUDED_protocols_thermal_sampling_util_HH
#define INCLUDED_protocols_thermal_sampling_util_HH

#include <core/types.hh>
#include <utility/io/ozstream.hh>

#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/ScoreType.hh>
#include <core/pose/Pose.fwd.hh>

// Utility headers

// ObjexxFCL headers

//// C++ headers
#include <string>
#include <map>
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
	core::Real min_, max_, spacing_;
	core::Size n_elem_;
	utility::vector1<core::Size> hist_;
};

//////////////////////////////////////////////////////////////////////////////
utility::vector1<core::scoring::ScoreType> const & get_scoretypes();

core::Size data_dim();

//////////////////////////////////////////////////////////////////////////////
void update_scores(
	utility::vector1<float> & scores,
	core::pose::Pose & pose,
	core::scoring::ScoreFunctionCOP scorefxn
);

//////////////////////////////////////////////////////////////////////////////
void fill_data(
	utility::vector1<float> & data,
	core::Size const count,
	utility::vector1<float> const & scores
);

core::Real gaussian_stdev( core::Real const n_rsd, core::Real const temp, bool const is_bp );

//////////////////////////////////////////////////////////////////////////////
template<typename T>
void vector2disk_in1d(
	std::string const & out_filename,
	utility::vector1<T> const & out_vector
) {
	utility::io::ozstream out( out_filename.c_str(), std::ios::out | std::ios::binary );
	out.write( (const char*) &out_vector[1], sizeof(T) * out_vector.size() );
	out.close();
}
//////////////////////////////////////////////////////////////////////////////
template<typename T>
void vector2disk_in2d(
	std::string const & out_filename,
	core::Size const dim1,
	core::Size const dim2,
	utility::vector1<T> const & out_vector
) {
	//std::cout << "dim1 " << dim1 << " dim2 " << dim2 << " out " << out_vector.size() << std::endl;
	utility::io::ozstream out( out_filename.c_str(), std::ios::out | std::ios::binary );
	runtime_assert( dim1 * dim2 == out_vector.size() );
	out.write( (const char*) &dim1, sizeof(core::Size) );
	out.write( (const char*) &dim2, sizeof(core::Size) );
	if ( out_vector.size() == 0 ) {
		std::cout << "Warning: no data available for the requested condition. Output file may be malformed.\n";
		// noop
	} else {
		out.write( (const char*) &out_vector[1], sizeof(T) * out_vector.size() );
	}
	out.close();
}

} //thermal_sampling
} //protocols

#endif
