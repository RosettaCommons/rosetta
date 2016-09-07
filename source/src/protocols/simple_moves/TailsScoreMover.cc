// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   TailsScoreMover.cc
///
/// @brief
/// @author Monica Berrondo

// unit headers
#include <protocols/simple_moves/TailsScoreMover.hh>

// type headers
#include <core/types.hh>

// project headers
#include <basic/options/option.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <basic/Tracer.hh>


// utility headers
#include <utility/file/FileName.hh>

#include <basic/options/keys/krassk.OptionKeys.gen.hh>

// C++ headers
#include <fstream>
#include <iostream>
#include <string>

#include <utility/vector1.hh>

using basic::T;
using basic::Error;
using basic::Warning;

static THREAD_LOCAL basic::Tracer TR( "protocols.simple_moves.TailsScoreMover" );
using namespace core;
using namespace std;

namespace protocols {
namespace simple_moves {

const int VISITED = 1;
const int PREVIOUS_HILL = 2;
const int HILL = 3;
const int VISITED_HILL = 4;
const double ED = 2;
const double KT = 15;

double TailsScoreMover::visit(
	double in_current_min,
	int in_current_min_ltail,
	int in_current_min_rtail,
	int in_ltail,
	int in_rtail,
	int in_array_of_visits[][200],
	int &out_min_ltail,
	int &out_min_rtail,
	int in_sequence_length,
	utility::vector1< core::Size > & tail,
	pose::Pose & pose,std::ofstream& area_file
)
{
	using namespace scoring;
	using namespace moves;
	using namespace basic::options;
	using namespace core::pose;

	if ( in_ltail > in_sequence_length/2 || in_ltail >= 200 || in_rtail > in_sequence_length/2 || in_rtail>= 200 ) {
		m_done_all = true;
		TR<< "done" << std::endl;
	}

	if ( in_ltail < 0 || in_ltail > in_sequence_length/2 || in_ltail >= 200 || in_rtail < 0 || in_rtail > in_sequence_length/2 || in_rtail>= 200 ) { // check that we are not outside of boundaries
		out_min_ltail = in_current_min_ltail;
		out_min_rtail = in_current_min_rtail;
		return in_current_min;
	}
	if ( in_array_of_visits[in_ltail][in_rtail] == VISITED || in_array_of_visits[in_ltail][in_rtail] == VISITED_HILL || in_array_of_visits[in_ltail][in_rtail] ==HILL ) { //if this spot was already visited or it is a hill
		out_min_ltail = in_current_min_ltail;
		out_min_rtail = in_current_min_rtail;
		return in_current_min;
	}
	make_tail(tail,in_ltail,in_rtail, in_sequence_length);
	Real tail_score = score_function()->get_sub_score_exclude_res(pose, tail);
	Real updated_tail_score = tail_score - ED*(in_ltail+in_rtail);

	//If you are not at the edge check if you are on the hill
	if ( in_ltail!=0 && in_array_of_visits[in_ltail][in_rtail] != PREVIOUS_HILL && in_array_of_visits[in_ltail][in_rtail] != HILL && in_array_of_visits[in_ltail][in_rtail] != VISITED_HILL ) {
		make_tail(tail,in_ltail+1,in_rtail, in_sequence_length);
		Real up_updated_tail_score = score_function()->get_sub_score_exclude_res(pose, tail) - ED*(in_ltail+1+in_rtail);
		make_tail(tail,in_ltail-1,in_rtail, in_sequence_length);
		Real down_updated_tail_score = score_function()->get_sub_score_exclude_res(pose, tail) - ED*(in_ltail-1+in_rtail);
		if ( up_updated_tail_score < updated_tail_score && down_updated_tail_score < updated_tail_score ) {
			// TR<< "Hill" << "    " << in_ltail<< "   " << in_rtail << "   " << updated_tail_score << "    "<< up_updated_tail_score << "     " << down_updated_tail_score << "    "<< in_array_of_visits[in_ltail][in_rtail] << std::endl;
			m_hill_size = m_hill_size + updated_tail_score;
			m_number_of_hill_points+=1;
			// TR<< "Updating hill size     " <<m_hill_size<<"    "<<m_number_of_hill_points <<  std::endl;
			in_array_of_visits[in_ltail][in_rtail] = HILL;
			out_min_ltail = in_current_min_ltail;
			out_min_rtail = in_current_min_rtail;
			return in_current_min;
		}
	}
	if ( in_rtail!=0 && in_array_of_visits[in_ltail][in_rtail] != PREVIOUS_HILL && in_array_of_visits[in_ltail][in_rtail] != HILL && in_array_of_visits[in_ltail][in_rtail] != VISITED_HILL ) {
		make_tail(tail,in_ltail,in_rtail+1, in_sequence_length);
		Real right_updated_tail_score = score_function()->get_sub_score_exclude_res(pose, tail) - ED*(in_ltail+in_rtail+1);
		make_tail(tail,in_ltail,in_rtail-1, in_sequence_length);
		Real left_updated_tail_score = score_function()->get_sub_score_exclude_res(pose, tail) - ED*(in_ltail+in_rtail-1);
		if ( left_updated_tail_score < updated_tail_score && right_updated_tail_score < updated_tail_score ) {
			// TR<< "Hill" << "    " << in_ltail<< "   " << in_rtail << "   " << updated_tail_score << "    "<< left_updated_tail_score << "     " << right_updated_tail_score << "    "<< in_array_of_visits[in_ltail][in_rtail] << std::endl;
			m_hill_size = m_hill_size + updated_tail_score;
			m_number_of_hill_points+=1;
			// TR<< "Updating hill size     " <<m_hill_size<<"    "<<m_number_of_hill_points <<  std::endl;
			in_array_of_visits[in_ltail][in_rtail] = HILL;
			out_min_ltail = in_current_min_ltail;
			out_min_rtail = in_current_min_rtail;
			return in_current_min;
		}
	}

	//Mark this as visited
	if ( in_array_of_visits[in_ltail][in_rtail] != PREVIOUS_HILL && in_array_of_visits[in_ltail][in_rtail] != HILL && in_array_of_visits[in_ltail][in_rtail] != VISITED_HILL ) {
		in_array_of_visits[in_ltail][in_rtail] = VISITED;
	}
	if ( in_array_of_visits[in_ltail][in_rtail] == PREVIOUS_HILL ) {
		in_array_of_visits[in_ltail][in_rtail] = VISITED_HILL;
	}

	area_file<< in_ltail<< "   " << in_rtail << "   " << updated_tail_score << std::endl;
	Real min_updated_tail_score = updated_tail_score;
	out_min_ltail = in_ltail;
	out_min_rtail = in_rtail;

	int t_ltail;
	int t_rtail;
	updated_tail_score = visit(in_current_min,in_current_min_ltail, in_current_min_rtail, in_ltail-1,in_rtail, in_array_of_visits,t_ltail, t_rtail,in_sequence_length,tail,  pose,area_file);
	if ( updated_tail_score < min_updated_tail_score ) {
		out_min_ltail = t_ltail;
		out_min_rtail = t_rtail;
		min_updated_tail_score = updated_tail_score;
	}
	updated_tail_score = visit(in_current_min, in_current_min_ltail, in_current_min_rtail, in_ltail,in_rtail+1, in_array_of_visits,t_ltail, t_rtail,in_sequence_length,tail,  pose,area_file);
	if ( updated_tail_score < min_updated_tail_score ) {
		out_min_ltail = t_ltail;
		out_min_rtail = t_rtail;
		min_updated_tail_score = updated_tail_score;
	}
	updated_tail_score = visit(in_current_min,in_current_min_ltail, in_current_min_rtail, in_ltail+1,in_rtail, in_array_of_visits,t_ltail, t_rtail,in_sequence_length,tail,  pose,area_file);
	if ( updated_tail_score < min_updated_tail_score ) {
		out_min_ltail = t_ltail;
		out_min_rtail = t_rtail;
		min_updated_tail_score = updated_tail_score;
	}
	updated_tail_score = visit(in_current_min,in_current_min_ltail, in_current_min_rtail, in_ltail,in_rtail-1, in_array_of_visits,t_ltail, t_rtail,in_sequence_length,tail,  pose,area_file);
	if ( updated_tail_score < min_updated_tail_score ) {
		out_min_ltail = t_ltail;
		out_min_rtail = t_rtail;
		min_updated_tail_score = updated_tail_score;
	}
	if ( min_updated_tail_score < in_current_min ) {
		return min_updated_tail_score;
	} else {
		return in_current_min;
		//out_min_ltail = in_current_min_ltail;
		//out_min_rtail = in_current_min_rtail;
	}
}
void TailsScoreMover::make_tail(utility::vector1< core::Size > & tail,int in_ltaillength, int in_rtaillength, int in_sequence_length)
{
	tail.clear();
	//populating tail with left tail residues
	for ( int rcount = 1; rcount<=in_ltaillength; rcount++ ) {
		tail.push_back(rcount);
	}
	//Putting in right tail
	for ( int rcount = 1; rcount <  in_rtaillength+1; rcount++ ) {
		tail.push_back(in_sequence_length - in_rtaillength+rcount);
	}
}
//////////////////////////////////////////////////
//Searches through all combinations of tails    //
//////////////////////////////////////////////////
double TailsScoreMover::score_mode1(int& out_min_ltail_length, int& out_min_rtail_length,std::ofstream & in_tail_output, pose::Pose & pose)
{
	TR<< "mode 1" << std::endl;
	int sequence_length = pose.total_residue(); // Need to get protein sequence length
	utility::vector1< core::Size > tail;
	Real min_updated_tail_score = score_function()->get_sub_score_exclude_res(pose, tail);
	Real updated_tail_score = score_function()->get_sub_score_exclude_res(pose, tail);
	Real min_ltail_length = 0;
	Real min_rtail_length = 0;
	int number_of_ltail_steps = sequence_length/2;
	int number_of_rtail_steps = sequence_length/2;
	for ( int ltaillength = 0; ltaillength < number_of_ltail_steps; ltaillength ++ ) {
		tail.clear();
		for ( int rcount = 1; rcount<=ltaillength; rcount++ ) { //populating tail with left tail residues
			tail.push_back(rcount);
		}
		for ( int rtaillength = 0; rtaillength < number_of_rtail_steps; rtaillength ++ ) {
			//removing previous residues
			for ( int rcount = 1; rcount < rtaillength; rcount++ ) {
				if ( !tail.empty() ) {
					tail.pop_back();
				}
			}
			//adding residues from the right tail
			for ( int rcount = 1; rcount <  rtaillength+1; rcount++ ) {
				tail.push_back(sequence_length - rtaillength+rcount);
			}
			Real tail_score = score_function()->get_sub_score_exclude_res(pose, tail);
			updated_tail_score = tail_score - ED*(ltaillength+rtaillength);
			in_tail_output <<  updated_tail_score << "    " << ltaillength <<  "     " << rtaillength <<"   "<<tail_score << std::endl;
			if ( updated_tail_score < min_updated_tail_score ) {
				min_updated_tail_score = updated_tail_score;
				min_ltail_length = ltaillength;
				min_rtail_length = rtaillength;
			}
		}
	}
	in_tail_output << "       Min updated score       " <<  min_updated_tail_score << "    " << min_ltail_length <<  "     " << min_rtail_length << std::endl;
	in_tail_output.close();
	out_min_ltail_length = (int)min_ltail_length;
	out_min_rtail_length = (int)min_rtail_length;
	return min_updated_tail_score;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//finds left tail first then right tail                                                                              //
//to stay in local minimum we are going to figure out the number of tail steps as following:                         //
//cut tail until the energy doesn't increase by KT from current min value or reach the middle of protein sequence       //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double TailsScoreMover::score_mode2(int& out_min_ltail_length, int& out_min_rtail_length,std::ofstream & in_tail_output, pose::Pose & pose)
{
	TR<< "mode 2" << std::endl;
	int sequence_length = pose.total_residue(); // Need to get protein sequence length
	//search for left tail first
	int number_of_tail_steps = sequence_length/2;
	utility::vector1< core::Size > tail;
	Real min_updated_tail_score = score_function()->get_sub_score_exclude_res(pose, tail);
	//Real original_tail_score = score_function()->get_sub_score_exclude_res(pose, tail);
	Real updated_tail_score = score_function()->get_sub_score_exclude_res(pose, tail);
	Real min_tail_length = 0;
	Real min_ltail_length = 0;
	Real min_rtail_length = 0;
	int taillength = 0;

	while ( updated_tail_score < min_updated_tail_score + KT && taillength < number_of_tail_steps )
			{
		tail.push_back(taillength);
		Real tail_score = score_function()->get_sub_score_exclude_res(pose, tail);
		updated_tail_score = tail_score - ED*taillength;
		if ( updated_tail_score < min_updated_tail_score ) {
			min_updated_tail_score = updated_tail_score;
			min_tail_length = taillength;
		}
		taillength++;
	}
	min_ltail_length = min_tail_length;

	//Clean the tail vector first
	tail.clear();
	//Put in left tail that was found
	for ( int rcount = 1; rcount<=min_ltail_length; rcount++ ) { //populating tail with left tail residues
		tail.push_back(rcount);
	}
	min_updated_tail_score = score_function()->get_sub_score_exclude_res(pose, tail);
	//original_tail_score = score_function()->get_sub_score_exclude_res(pose, tail);  // set but never used ~Labonte
	updated_tail_score = score_function()->get_sub_score_exclude_res(pose, tail);
	min_tail_length = 0;
	//search for right tail now
	taillength = 0;
	while ( updated_tail_score < min_updated_tail_score + KT && taillength < number_of_tail_steps )
			{
		//removing previous residues
		for ( int rcount = 1; rcount < taillength; rcount++ ) {
			if ( !tail.empty() ) {
				tail.pop_back();
			}
		}
		//adding residues from the right tail
		for ( int rcount = 1; rcount <  taillength+1; rcount++ ) {
			tail.push_back(sequence_length - taillength+rcount);
		}
		Real tail_score = score_function()->get_sub_score_exclude_res(pose, tail);
		updated_tail_score = tail_score - ED*(min_ltail_length + taillength);
		if ( updated_tail_score < min_updated_tail_score ) {
			min_updated_tail_score = updated_tail_score;
			min_tail_length = taillength;
		}
		taillength++;
	}
	min_rtail_length = min_tail_length;
	in_tail_output <<  min_updated_tail_score << "    " << min_ltail_length <<  "     " << min_rtail_length << std::endl;
	in_tail_output.close();

	out_min_ltail_length = (int)min_ltail_length;
	out_min_rtail_length = (int)min_rtail_length;
	return min_updated_tail_score;
}

//Looks for local min
double TailsScoreMover::score_mode3(int& out_min_ltail_length, int& out_min_rtail_length,std::ofstream & in_tail_output, pose::Pose & pose)
{
	TR<< "mode 3" << std::endl;
	int array_of_visits[200][200] = {{0}};
	int sequence_length = pose.total_residue(); // Need to get protein sequence length
	utility::vector1< core::Size > tail;

	int min_ltail_length = 0;
	int min_rtail_length = 0;
	Real min_updated_tail_score = 0;

	double hill_size = 0;
	m_done_all = false;
	while ( hill_size < KT && !m_done_all )
			{
		min_ltail_length = 0;
		min_rtail_length = 0;
		core::Real m_hill_size = 0;
		core::Real  m_number_of_hill_points = 0;
		make_tail(tail, 0,0,sequence_length);
		double min_current_tail_score = score_function()->get_sub_score_exclude_res(pose, tail);
		//Real barrier = min_current_tail_score + 2;

		min_updated_tail_score = visit(min_current_tail_score, 0, 0,
			0, 0, array_of_visits, min_ltail_length,
			min_rtail_length, sequence_length, tail, pose,in_tail_output);

		assert( m_number_of_hill_points != 0 );
		hill_size = (m_hill_size - min_updated_tail_score*m_number_of_hill_points)/m_number_of_hill_points;
		// TR<< "Hill size is    " << hill_size <<"    "<<m_hill_size<<"     "<< min_updated_tail_score <<"     "<<m_number_of_hill_points<< std::endl;

		// tail_output_file << "       Min updated score       " <<  min_updated_tail_score << "    " << min_ltail_length <<  "     " << min_rtail_length << std::endl;
		//Clean visits array
		// TR<< "Updating array of visits..." <<  std::endl;
		for (auto & array_of_visit : array_of_visits) {
			for ( int j=0; j < 200; j++ ) {
				if ( array_of_visit[j]==HILL ) {
					// TR<< "Hill at" << "    " << i<< "   " << j << std::endl;
					array_of_visit[j] = PREVIOUS_HILL;
				} else if ( array_of_visit[j]==PREVIOUS_HILL || array_of_visit[j]==VISITED_HILL ) {
					//TR<< "Hill at" << "    " << i<< "   " << j << std::endl;
					array_of_visit[j] = PREVIOUS_HILL;
				} else {
					array_of_visit[j] = 0;
				}
			}
		}
	}
	in_tail_output << "      Final Min updated score       " <<  min_updated_tail_score << "    " << min_ltail_length <<  "     " << min_rtail_length << std::endl;
	in_tail_output.close();

	out_min_ltail_length = min_ltail_length;
	out_min_rtail_length = min_rtail_length;
	return min_updated_tail_score;
}

////////////////////////////////////////////////////////////////////////////////
///
/// @brief Simply score a pdb
void
TailsScoreMover::apply( pose::Pose & pose )
{
	using namespace scoring;
	using namespace moves;
	using namespace basic::options;
	using namespace core::pose;
	//using core::pose::datacache::CacheableDataType::SCORE_MAP;

	(*score_function())(pose);


	int tail_mode_name = basic::options::option[ basic::options::OptionKeys::krassk::tail_mode_name ](); // What algorithm to use for tail prediction
	if ( tail_mode_name < 1 || tail_mode_name > 3 ) {
		Error() << "No such tail algorithm supported. " << std::endl;
	}

	string tail_output_file_name = basic::options::option[ basic::options::OptionKeys::krassk::tail_output_file_name ](); // Name of the tail output file
	const char * tt = tail_output_file_name.c_str();

	std::ofstream  tail_output_file(tt, std::ios::app );
	if ( !tail_output_file.good() ) {
		Error() << "Unable to open tail output file for writing. " << std::endl;
	}

	double free_energy = 0;
	int out_min_ltail_length = 0;
	int out_min_rtail_length = 0;

	switch(tail_mode_name)
			{
			case 1 :
				{
				free_energy = score_mode1(out_min_ltail_length, out_min_rtail_length,tail_output_file,pose);
				break;
			}
			case 2 :
				{
				free_energy = score_mode2(out_min_ltail_length, out_min_rtail_length,tail_output_file,pose);
				break;
			}
			case 3 :
				{
				free_energy = score_mode3( out_min_ltail_length, out_min_rtail_length,tail_output_file,pose);
				break;
			}
			default :
				{
				free_energy = score_mode1(out_min_ltail_length,out_min_rtail_length,tail_output_file,pose);
			}
			}

	setPoseExtraScore( pose, "free_energy", free_energy);
	setPoseExtraScore( pose, "left_tail",   out_min_ltail_length);
	setPoseExtraScore( pose, "right_tail",  out_min_rtail_length);
}//apply

std::string
TailsScoreMover::get_name() const {
	return "TailsScoreMover";
}

} // moves
} // protocols
