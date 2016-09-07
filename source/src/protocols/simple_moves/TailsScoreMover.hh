// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   TailsScoreMover.hh
/// @brief
/// @author Monica Berrondo

#ifndef INCLUDED_protocols_simple_moves_TailsScoreMover_hh
#define INCLUDED_protocols_simple_moves_TailsScoreMover_hh

// Unit headers
#include <protocols/simple_moves/ScoreMover.hh>

#include <core/types.hh>

#include <core/scoring/ScoreFunction.hh>

#include <protocols/moves/Mover.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace simple_moves {

class TailsScoreMover : public ScoreMover
{
public:
	/// @brief
	///  empty constructor fills values with the values
	///  read in from the commandline
	TailsScoreMover() :
		ScoreMover()
	{}

	/// @brief
	TailsScoreMover( core::scoring::ScoreFunctionOP scorefxn_in ) :
		ScoreMover(scorefxn_in)
	{}

	void apply( core::pose::Pose & pose ) override;
	std::string get_name() const override;
private:
	double score_mode1(int& out_min_ltail_length, int& out_min_rtail_length,std::ofstream & in_tail_output, core::pose::Pose & pose);
	double score_mode2(int& out_min_ltail_length, int& out_min_rtail_length,std::ofstream & in_tail_output, core::pose::Pose & pose);
	double score_mode3(int& out_min_ltail_length, int& out_min_rtail_length,std::ofstream & in_tail_output, core::pose::Pose & pose);

	double visit(double in_current_min, int in_current_min_ltail, int in_current_min_rtail,
		int in_ltail, int in_rtail, int in_array_of_visits [][200], int &out_min_ltail,
		int &out_min_rtail, int in_sequence_length, utility::vector1< core::Size > & tail, core::pose::Pose & pose,std::ofstream& area_file);
	void make_tail(utility::vector1< core::Size > & tail,int in_ltaillength, int in_rtaillength, int in_sequence_length);
	core::Real m_hill_size;
	int m_number_of_hill_points;
	bool m_done_all;

};

} // simple_moves
} // protocols

#endif //INCLUDED_protocols_simple_moves_TailsScoreMover_HH
