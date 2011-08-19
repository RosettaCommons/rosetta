// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/NV/NVlookup.cc
/// @brief  Neighbor Vector algorithm lookup table processing class
/// @detailed The neighbor vector algorithm relies on a table of knowledge based potentials based on statistical analysis of the PDB
/// This class loads the lookup table and provides an interface for accessing it.  In the future this should and will probably be implemented as a spline rather than a table
/// @author Sam DeLuca (samuel.l.deluca@vanderbilt.edu)

#include <utility/vector1.hh>
#include <core/scoring/nv/NVlookup.hh>
// AUTO-REMOVED #include <core/pose/Pose.fwd.hh>
#include <numeric/interpolation/spline/SplineGenerator.hh>
// AUTO-REMOVED #include <numeric/interpolation/spline/SimpleInterpolator.hh>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
// AUTO-REMOVED #include <string>
#include <vector>
#include <basic/Tracer.hh>

using namespace std;

static basic::Tracer TR("core.scoring.NV");

namespace core {
namespace scoring {
namespace nv {

using namespace numeric::interpolation::spline;

NVlookup::NVlookup(std::string filename)
{
	ifstream infile;
	infile.open(filename.c_str(),ifstream::in);
	string line;
	do
	{
		//get every line in the db file
		getline(infile,line);
		if(line.size() > 0)
		{
			istringstream linebuffer(line,istringstream::in);

			char residue;
			//the first field of the line is a single letter aa abbrev.
			linebuffer >> residue;
			//assert(residue.size() == 1); //Something is wrong with your input file

			//0<x<1
			Real lower_bound_x = 0;
			Real upper_bound_x = 1.1;

			utility::vector1<core::Real> potential_vect;

			core::Real potential;
			//read the bins into a vector
			while(linebuffer >> potential)
			{
				potential_vect.push_back(potential);
			}

			utility::vector1<core::Real>::iterator upper_bound_y_iterator(std::max_element(potential_vect.begin(),potential_vect.end()));
			utility::vector1<core::Real>::iterator lower_bound_y_iterator(std::min_element(potential_vect.begin(),potential_vect.end()));

			core::Real lower_bound_y = *lower_bound_y_iterator;
			core::Real upper_bound_y = *upper_bound_y_iterator;

			lower_bound_y -= 0.1; //Bounds are exclusive, so we pad out the min and max a bit
			upper_bound_y += 0.1;

			SplineGenerator common_spline(lower_bound_x,lower_bound_y,0,upper_bound_x,upper_bound_y,0);

			for(Size potential_index = 1; potential_index <= potential_vect.size(); ++potential_index)
			{
				Real x_value = potential_index*0.05;
				Real y_value = potential_vect[potential_index];
				common_spline.add_known_value(x_value,y_value);
			}

			std::pair<char,SplineGenerator > residueRow(residue,common_spline);

			lookup_table_.insert(residueRow);
		}
	}while(!infile.eof());

	infile.close();

}
//Given a neighbor vector score and a single letter AA abbrev, get the potential from the lookup table
core::Real NVlookup::get_potentials(char &name, core::Real &score) const
{

	map<char, SplineGenerator >::const_iterator it = lookup_table_.find(name);
	if(it != lookup_table_.end())
	{
		SplineGenerator row(it->second);
		//table is binned in increments of 0.05
		Real potential_energy;
		Real delta_potential_energy;

		InterpolatorOP row_interpolator(row.get_interpolator());
		row_interpolator->interpolate(score,potential_energy,delta_potential_energy);

		return potential_energy;
	}

	return 0;


}

} //NVscore
} //scoring
} //core
