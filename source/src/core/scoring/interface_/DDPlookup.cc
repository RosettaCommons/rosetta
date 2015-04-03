// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   src/core/scoring/Interface_/DDPlookup.cc
/// @author Hermann Zellner (hermann1.zellner@biologie.uni-regensburg.de)
#include <core/scoring/interface_/DDPlookup.hh>

#include <basic/database/open.hh>
#include <basic/Tracer.hh>
#include <utility/io/izstream.hh>

#include <utility/vector1.hh>


namespace core {
namespace scoring {
namespace interface_ {

static thread_local basic::Tracer TR( "core.scoring.Interface_.DDPlookup" );

DDPlookup::DDPlookup(std::string filename)
{
	core::Size const n_bins ( 50 );
	core::Real const n_max (20.);
	core::Size const max_aa( 20 );

	std::string tag, line;
	core::chemical::AA aa1, aa2;

	utility::vector1< core::Real > myvec;
	utility::vector1< core::Real > myvec_2;
	utility::vector1< utility::vector1< core::Real > > mymyvec;
	utility::vector1< utility::vector1< utility::vector1< core::Real > > > known_values;

	for (core::Size i=1; i<=n_bins; i++)
	{
		myvec.push_back(0.);
	}
	for (core::Size i=1; i<=max_aa; i++)
	{
		myvec_2.push_back(0.);
	}
	for (core::Size i=1; i<=max_aa; i++)
	{
		left_.push_back(myvec_2);
		right_.push_back(myvec_2);
		mymyvec.push_back(myvec);
	}
	for (core::Size i=1; i<=max_aa; i++)
	{
		known_values.push_back(mymyvec);
	}
	utility::vector1< core::Real > distance_bin;


	core::Real stepsize = n_max / (core::Real) n_bins;
	core::Real e = 0 + stepsize/2;

	for (core::Size i=1; i<=n_bins; i++)
	{
		distance_bin.push_back(e);
		e += stepsize;
	}

	utility::io::izstream stream;
	basic::database::open( stream, filename );
	while ( getline( stream, line ) ) {
		std::istringstream l(line);
		l >> tag >> aa1 >> aa2;
		for ( Size i=1; i<=n_bins; ++i ) {
			l >> known_values[aa1][aa2][i];
		}
		if ( l.fail() || tag != "DDIPP:"  )
		{
			TR << "bad format for interface_dd_score.txt" << std::endl;
			TR << "fail: " << l.fail() << "  tag: " << tag << std::endl;
			exit(1);
		}

		// #######################################
		// Identify borders of the valey. Left border ist the smallest distance with score < 0. of a valey with an integral larger than .5
		// The right border is either the smallest distance with score >= 0 or a maximum of the score (if you have 2 separated minima)
		bool valey = false;
		core::Size leftidx=1;
		while( !valey && leftidx < n_bins )
		{
			core::Real auc=0;
			for (core::Size i=leftidx; i<=n_bins; i++)
			{
				if ( known_values[aa1][aa2][i] < -1e-10 )
				{
					left_[aa1][aa2] = distance_bin[i-1];
					leftidx = i;
					break;
				}
				leftidx = n_bins;
			}
			bool up=false;

			for (core::Size i=leftidx; i<=n_bins; i++)
			{
				auc += known_values[aa1][aa2][i] * (n_max/n_bins);
				if ( !up && i>1 && known_values[aa1][aa2][i-1] < known_values[aa1][aa2][i])
				{
					up=true;
				}
				if ( known_values[aa1][aa2][i] > -1e-10 )
				{
					right_[aa1][aa2] = distance_bin[i];
					break;
				}
				if ( up && known_values[aa1][aa2][i-1] > known_values[aa1][aa2][i] )
				{
					right_[aa1][aa2] = distance_bin[i-3];
					break;
				}
				right_[aa1][aa2] = distance_bin[1];
			}
			if ( auc < .5 )
			{
				valey=true;
			}
			else
			{
				left_[aa1][aa2] = distance_bin[n_bins];
				right_[aa1][aa2] = distance_bin[1];
			}
			leftidx++; // to avoid an endless loop
		}
	}
	stream.close();

	//generate spline for every pair of amino acid
	lookup_table_ = new numeric::interpolation::spline::SplineGenerator**[max_aa+1];
	for (core::Size i=1; i<=max_aa; i++)
	{
		lookup_table_[i] = new numeric::interpolation::spline::SplineGenerator*[max_aa+1];
		for (core::Size j=1; j<=max_aa; j++)
		{
			utility::vector1<core::Real>::iterator lower_bound_y_iterator(std::min_element(known_values[i][j].begin(),known_values[i][j].end()));

			core::Real lower_bound_x = 0 - stepsize/2;
			core::Real upper_bound_x = distance_bin[n_bins] + stepsize/2;
		    core::Real lower_bound_y = *lower_bound_y_iterator;
		    core::Real upper_bound_y = 0.;

	        lower_bound_y -= 0.1; //Bounds are exclusive, so we pad out the min and max a bit
	        upper_bound_y += 0.1;

	        lookup_table_[i][j] = new numeric::interpolation::spline::SplineGenerator(lower_bound_x, lower_bound_y, 0, upper_bound_x, upper_bound_y, 0);
			for (core::Size energy_index = 1; energy_index <= n_bins; energy_index++)
			{
				if (distance_bin[energy_index] > left_[i][j] && distance_bin[energy_index] < right_[i][j])
				{
					lookup_table_[i][j]->add_known_value( distance_bin[energy_index], known_values[i][j][energy_index]);
				}
				else
				{
					lookup_table_[i][j]->add_known_value( distance_bin[energy_index], -1e-5);
				}
			}
		}
	}
}


core::Real
DDPlookup::get_potentials( const core::chemical::AA & aa1, const core::chemical::AA & aa2, core::Real distance ) const
{
	if ( aa1 < 1 || aa1 > 20 || aa2 < 1 || aa2 > 20 )
	{
		TR << "bad amino acid for DDPlookup" << std::endl;
		exit(1);
	}

	core::Real delta_potential_energy = 0;
	core::Real potential_energy = 0;
	core::Real result;


	numeric::interpolation::spline::InterpolatorOP interpolator =
			lookup_table_[aa1][aa2]->get_interpolator();
	interpolator->interpolate(distance, potential_energy, delta_potential_energy);
	if ( potential_energy > 0. || distance < 2. || distance < left_[aa1][aa2] || distance > right_[aa1][aa2] )
	{
		result = -1e-5;
	}
	else
	{
		result = potential_energy;
	}


	return result;
}

} //NVscore
} //scoring
} //core
