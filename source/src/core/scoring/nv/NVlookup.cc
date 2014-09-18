// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/nv/NVlookup.cc
/// @brief  Neighbor Vector algorithm lookup table processing class
/// @detailed The neighbor vector algorithm relies on a table of knowledge based potentials based on statistical analysis of the PDB
/// This class loads the lookup table and provides an interface for accessing it.  In the future this should and will probably be implemented as a spline rather than a table
/// @author Sam DeLuca (samuel.l.deluca@vanderbilt.edu)

#include <core/scoring/nv/NVlookup.hh>

#include <numeric/interpolation/spline/SplineGenerator.hh>

#include <utility/vector1.hh>
#include <utility/io/izstream.hh>
#include <utility/json_spirit/json_spirit_reader.h>

#include <basic/Tracer.hh>

#include <algorithm>
#include <fstream>
#include <sstream>
#include <vector>

using namespace std;

static thread_local basic::Tracer TR( "core.scoring.NV" );

namespace core {
namespace scoring {
namespace nv {

/// @details Auto-generated virtual destructor
NVlookup::~NVlookup() {}

using namespace numeric::interpolation::spline;

void NVlookup::set_up_spline_from_data(core::chemical::AA const & aa_type, utility::vector1<core::Real> const & bin_centers, utility::vector1<core::Real> const & data)
{

	assert(bin_centers.size() == data.size()); //You need to have the same number of bin centers and spline points

	Real lower_bound_x = 0;
	Real upper_bound_x = 1.1;

	core::Real lower_bound_y = *std::min_element(data.begin(),data.end());
	core::Real upper_bound_y = *std::max_element(data.begin(),data.end());

	lower_bound_y -= 0.1; //Bounds are exclusive, so we pad out the min and max a bit
	upper_bound_y += 0.1;

	SplineGenerator spline(lower_bound_x,lower_bound_y,0,upper_bound_x,upper_bound_y,0);

	for(Size potential_index = 1; potential_index <= data.size(); ++potential_index)
	{
		Real x_value = bin_centers[potential_index];
		Real y_value = data[potential_index];
		spline.add_known_value(x_value,y_value);
	}

	//InterpolatorOP interpolator = spline.get_interpolator();

	lookup_table_[aa_type] = spline.get_interpolator();


}

NVlookup::NVlookup(std::string filename) : lookup_table_(core::chemical::num_canonical_aas,0)
{
	//Look at rosetta_database/scoring/score_function/NV/neighbor_vector_score.histogram to see what formatting should look like
	utility::io::izstream infile;
	infile.open(filename.c_str(),ifstream::in);
	utility::json_spirit::mValue histogram_object;
	utility::json_spirit::read(infile,histogram_object);

	//The bin centers are in an array with the key "bin_centers"
	utility::json_spirit::mArray const & bin_array_data = histogram_object.get_obj()["bin_centers"].get_array();
	utility::vector1<core::Real> bin_centers;

	//the rest of the code expects a vector of reals, so we convert the array to one
	for(core::Size i =0; i < bin_array_data.size();++i)
	{
		bin_centers.push_back(bin_array_data[i].get_real());
	}

	//the score information is a dictionary, key is amino acid type, data is an array of reals
	utility::json_spirit::mObject const & spline_objects = histogram_object.get_obj()["scores"].get_obj();

	utility::json_spirit::mObject::const_iterator spline_object_iterator = spline_objects.begin();
	for(; spline_object_iterator != spline_objects.end();++spline_object_iterator)
	{
		std::string aa_name_string = spline_object_iterator->first;
		assert(aa_name_string.size() == 1); //The key should be a single letter amino acid
		char aa_name = aa_name_string[0];

		core::chemical::AA aa_type(core::chemical::aa_from_oneletter_code(aa_name));

		utility::json_spirit::mArray const & spline_array_data = spline_object_iterator->second.get_array();
		utility::vector1<core::Real> spline_values;
		for(core::Size i = 0; i < spline_array_data.size(); ++i)
		{
			spline_values.push_back(spline_array_data[i].get_real());
		}
		set_up_spline_from_data(aa_type,bin_centers,spline_values);

	}

	infile.close();

}
//Given a neighbor vector score and a single letter AA abbrev, get the potential from the lookup table
core::Real NVlookup::get_potentials(core::chemical::AA const & aa_type, core::Real const & score) const
{

	if(aa_type <= core::chemical::num_canonical_aas)
	{
		InterpolatorOP row(lookup_table_[aa_type]);
		Real potential_energy;
		Real delta_potential_energy;

		row->interpolate(score,potential_energy,delta_potential_energy);

		return potential_energy;
	}

	return 0;


}

} //NVscore
} //scoring
} //core
