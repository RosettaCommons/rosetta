// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @author James Thompson

// Unit headers
#include <protocols/simple_filters/PredictedBurialEvaluator.hh>

// Package headers
#include <core/types.hh>
#include <ObjexxFCL/string.functions.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <protocols/evaluation/PoseEvaluator.hh>
#include <utility/io/izstream.hh>
#include <string>

#include <utility/vector1.hh>


namespace protocols {
namespace simple_filters {

PredictedBurialEvaluator::PredictedBurialEvaluator(
	std::string const & fn
)  : evaluation::SingleValuePoseEvaluator< core::Real >( "burial" )
{
	init_from_file(fn);
}

PredictedBurialEvaluator::~PredictedBurialEvaluator() {}

void PredictedBurialEvaluator::apply(
	core::pose::Pose & pose,
	std::string,
	core::io::silent::SilentStruct & ss
) const {
	using core::Real;
	using core::Size;

	Real score(0.0);
	for ( Size ii = 1; ii <= pose.total_residue(); ++ii ) {
		core::conformation::Residue const & resi( pose.residue(ii) );
		Size countN(0);
		for ( Size jj = 1; jj <= pose.total_residue(); ++jj ) {
			core::conformation::Residue const & resj( pose.residue(jj) );
			core::Real const distance(
				resi.xyz(resi.nbr_atom()).distance( resj.xyz(resj.nbr_atom()) )
			);

			if ( distance < 10 ) countN++;
		}

		Real const burial_prediction( predicted_burial_[ii]);
		Real const score_ii( -1 * burial_prediction * countN );
		score += score_ii;
	}

	ss.add_energy( "burial", score );
}

void PredictedBurialEvaluator::init_from_file(
	std::string const &
) {
	utility::io::izstream input;
	std::string line;

	predicted_burial_.clear();
	while ( getline(input,line) ) {
		if ( line.substr(9,1) == "[" ) {
			std::istringstream line_input(line);
			std::string token;
			line_input >> token;
			while ( !line_input.fail() ) line_input >> token;

			float pred = ObjexxFCL::float_of( token.substr(2,4) );
			if ( token.substr(1,1) == "1" ) pred *= -1;
			std::cout << "line = " << line << std::endl;
			std::cout << "pred = " << pred << std::endl;
			predicted_burial_.push_back(pred);
		}
	}
}

} // simple_filter
} // protocols
