// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file apps/public/analysis/ResidueDisorder
/// @brief Application to run ResidueDisorder (predicting or measuring disordered regions of proteins).
/// @author Justin Seffernick (seffernick.9@osu.edu)

#include <basic/options/option.hh>
#include <devel/init.hh>
#include <basic/Tracer.hh>
#include <core/conformation/Conformation.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/Energies.hh>
#include <core/import_pose/pose_stream/MetaPoseInputStream.hh>
#include <core/import_pose/pose_stream/util.hh>
#include <utility/io/ozstream.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/corrections.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/option_macros.hh>

static basic::Tracer TR( "apps.public.analysis.ResidueDisorder" );

using namespace basic::options;
using namespace basic::options::OptionKeys;
//local options
OPT_KEY( Boolean, measure_disorder_from_structure )
OPT_KEY( Boolean, predict_events )
OPT_KEY( Boolean, AF )

void check_pose_length( core::Size const num_res, core::Size const min_seq_length) {
	if ( num_res < min_seq_length*2 ) utility_exit_with_message( "Pose too small!! All poses must contain at least " + std::to_string(min_seq_length*2) + " residues." );
}

void check_pose_size( core::Size const num_res_current, core::Size const num_res) {
	if ( num_res_current != num_res ) utility_exit_with_message("Pose size mismatch!! All poses must be of equal length.");
}

void check_single_chain( core::pose::PoseOP const &mypose ) {
	core::Size num_chains = mypose->conformation().num_chains();
	if ( num_chains != 1 ) utility_exit_with_message("Multiple chains!! All poses must be single chains.");
}

void check_sequence( core::pose::PoseOP const &pose, core::pose::PoseOP const &pose_current ) {
	if ( pose_current->sequence() != pose->sequence() ) utility_exit_with_message("Sequence mismatch!! All poses must have the same sequence.");
}

void check_number_of_poses ( core::Size const n_poses ) {
	if ( n_poses != 100 ) {
		TR.Warning << n_poses << " pdb's entered. Not recommended!! (100 recommended)" << std::endl;
	}
}

void calculate_average_residue_scores ( utility::vector1<core::pose::PoseOP> const &poses, utility::vector1<core::Real> &average_res_scores, std::string &sequence, core::Size const num_res, core::Size const n_poses ) {
	core::scoring::ScoreFunctionOP sfxn = core::scoring::get_score_function();
	//Loop through each pose
	for ( core::Size i=1; i<=n_poses; i++ ) {
		core::pose::PoseOP mypose = poses[i];
		//Score of the entire pose
		core::Real score = (*sfxn)(*mypose);
		TR << "SCORE: " <<  score << std::endl;
		//Make sure sequences are the same
		check_sequence( mypose, poses[1] );
		//Loop through each residue in the pose
		for ( core::Size j=1; j <= num_res; j++ ) {
			if ( i == 1 ) {
				sequence = mypose->sequence();
			}
			//Score of the current residue
			core::Real res_score = mypose->energies().residue_total_energy(j);
			TR << "residue " << j << " score: " << res_score << " REU" << std::endl;
			//Append residue scores to the average score vector
			if ( i == 1 ) average_res_scores.push_back(res_score);
			else average_res_scores[j]  = average_res_scores[j] + res_score;
		}
	}
	//Calculate average scores
	TR << "Averages: " << std::endl;
	//Loop through each residue
	for ( core::Size i=1; i <= num_res; i++ ) {
		//Average scores for each residue
		average_res_scores[i] = average_res_scores[i] / n_poses;
		TR << "avgresidue " << i << " score: " << average_res_scores[i] << " REU" << std::endl;
	}
}

void calculate_order_scores( utility::vector1<core::Real> &ORDER_RES_Score, utility::vector1<core::Real> const &average_res_scores, core::Size const num_res, core::Size const ORDER ) {
	core::Real Divide;
	//Loop through each residue
	for ( core::Size i=1; i<=num_res; i++ ) {
		core::Real TOTAL_Sum = 0.0;
		if ( i<=ORDER ) {
			//Case 1: within first 5 residues
			Divide = 0;
			for ( core::Size j=1; j<=(i+ORDER); j++ ) {
				TOTAL_Sum = TOTAL_Sum + average_res_scores[j];
				Divide++;
			}
		} else if ( (i+ORDER) > num_res ) {
			//Case 2: within 5 residues of end
			Divide = 0;
			for ( core::Size j=(i-ORDER); j<=num_res; j++ ) {
				TOTAL_Sum = TOTAL_Sum + average_res_scores[j];
				Divide++;
			}
		} else {
			//Case 3: not within 5 residues of either terminus
			Divide = 0;
			for ( core::Size j=(i-ORDER); j<=(i+ORDER); j++ ) {
				TOTAL_Sum = TOTAL_Sum + average_res_scores[j];
				Divide++;
			}
		}
		core::Real AVG = TOTAL_Sum / Divide;
		ORDER_RES_Score.push_back(AVG);
	}
}

void predict_disorder( utility::vector1<std::string> &PREDICTION, utility::vector1<core::Real> const &ORDER_RES_Score, std::string const &sequence, core::Size const num_res, utility::vector1<core::Real> const &cutoffs, bool const terminal ) {
	if ( !terminal ) TR << "Predictions:" << std::endl;
	else TR << "New Predictions:" << std::endl;
	//Loop through each residue
	for ( core::Size i=1; i<= num_res; i++ ) {
		if ( ORDER_RES_Score[i] >= cutoffs[i] ) {
			PREDICTION.push_back("D");
			TR << "residue " << i << " " << sequence[i-1] <<  " Disordered" << std::endl << "order score: " << ORDER_RES_Score[i] << " REU" << std::endl;
			if ( terminal ) TR << "New Cutoff: " << cutoffs[i] << std::endl;
		} else {
			PREDICTION.push_back("O");
			TR << "residue " << i << " " << sequence[i-1] << " Ordered" << std::endl << "order score: " << ORDER_RES_Score[i] << " REU" << std::endl;
			if ( terminal ) TR << "New Cutoff: " << cutoffs[i] << std::endl;
		}
	}
}

core::Real calculate_percent_disorder( utility::vector1<std::string> const &PREDICTION, core::Size const num_res ) {
	core::Real num_disorder;
	num_disorder = std::count(PREDICTION.begin(), PREDICTION.end(), "D");
	core::Real percent_disorder = num_disorder/num_res*100.0;
	TR << "Percent of disordered residues: " << percent_disorder << "%" << std::endl;
	return percent_disorder;
}

void do_terminal_opt( utility::vector1<core::Real> &cutoffs_new, core::Size const num_res, core::Real const CUTOFF, core::Real const term_cutoff, core::Real const endpoint ) {
	core::Real TerminalCount = (core::Real)num_res*term_cutoff;
	//Loop through each residue
	for ( core::Size i=1; i<=num_res; i++ ) {
		if ( i<=TerminalCount ) {
			//Case 1: within TerminalCount of beginning
			cutoffs_new[i] = (CUTOFF-endpoint)/(TerminalCount-1)*(i-1) + endpoint;
		} else if ( i>ceil(num_res-TerminalCount) ) {
			//Case 2: within TerminalCount of end
			cutoffs_new[i] = (endpoint-CUTOFF)/(TerminalCount-1)*(i+TerminalCount-num_res-1) + CUTOFF;
		}
		//Case 3: not within TerminalCount of either end => Do not change cutoff
	}
}

void output_results ( core::Size const num_res, std::string const &sequence, utility::vector1<std::string> const &PREDICTION, utility::vector1<core::Real> const &ORDER_RES_Score, utility::vector1<core::Real> const average_res_scores, core::Size const num_output) {
	std::ostringstream out;
	out << "Residue_number\tResidue_type\tPrediction\tOrder_score\tAverage_raw_score" << std::endl;
	for ( core::Size i=1; i<=num_res; i++ ) {
		out << i << "\t\t" << sequence.substr(i-1,1) << "\t\t" << PREDICTION[i] << "\t\t" << ORDER_RES_Score[i] << "\t" << average_res_scores[i] << std::endl;
	}
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	std::string outfile;
	if ( option[ measure_disorder_from_structure ] ) {
		outfile = "ResidueDisorder_default_" + ObjexxFCL::right_string_of( num_output, 6, '0' ) + ".out";
		if ( option[ out::prefix ].user() ) outfile = option[ out::prefix ]() + ObjexxFCL::right_string_of( num_output, 6, '0' ) + ".out";
	} else {
		outfile = "ResidueDisorder_default.out";
		if ( option[ out::file::o ].user() ) outfile = option[ out::file::o ]();
	}
	utility::io::ozstream outz( outfile.c_str() );
	outz << out.str();
	outz.close();
	outz.clear();
}

void calculate_disorder ( utility::vector1<core::pose::PoseOP> const &poses, core::Size const num_res, core::Size const n_poses, core::Size const n_output, utility::vector1<core::Real> &vector_percent_disorder, core::Size const ORDER, core::Real const CUTOFF, core::Real const term_cutoff, core::Real const endpoint, core::Real const perc_cut ) {

	utility::vector1<core::Real> average_res_scores;
	std::string sequence;

	//Calculate average residue scores
	calculate_average_residue_scores ( poses, average_res_scores, sequence, num_res, n_poses );

	//Calculate the order scores (i.e. window average score of each residue)
	utility::vector1<core::Real> ORDER_RES_Score;
	calculate_order_scores( ORDER_RES_Score, average_res_scores, num_res, ORDER );

	//Use order scores to predict disorder
	utility::vector1<std::string> PREDICTION;
	utility::vector1<core::Real> cutoffs ( num_res, CUTOFF);
	predict_disorder( PREDICTION, ORDER_RES_Score, sequence, num_res, cutoffs, false );

	//Calculate percent disorder
	core::Real percent_disorder = calculate_percent_disorder( PREDICTION, num_res );

	//If less than perc_cut% disordered, run the optimized terminal residue prediction
	if ( percent_disorder<=perc_cut ) {
		TR << "Optimizing terminal residues:" << std::endl;
		//Do terminal optimization->recalculate cutoffs
		utility::vector1<core::Real> cutoffs_new ( num_res, CUTOFF);
		do_terminal_opt( cutoffs_new, num_res, CUTOFF, term_cutoff, endpoint );
		utility::vector1<std::string> PREDICTION_new;
		//Re-predict disorder with new cutoffs
		predict_disorder( PREDICTION_new, ORDER_RES_Score, sequence, num_res, cutoffs_new, true );
		//Calculate percent disorder
		percent_disorder = calculate_percent_disorder( PREDICTION_new, num_res );
		//Output results of terminal prediction
		output_results( num_res, sequence, PREDICTION_new, ORDER_RES_Score, average_res_scores, n_output);
	} else {
		TR << "Not optimizing terminal residues!" << std::endl;
		//Output results of prediction
		output_results( num_res, sequence, PREDICTION, ORDER_RES_Score, average_res_scores, n_output);
	}
	vector_percent_disorder.push_back(percent_disorder);
}

void do_disorder_prediction ( utility::vector1<core::Real> &vector_percent_disorder, core::Size const ORDER, core::Real const CUTOFF, core::Real const term_cutoff, core::Real const endpoint, core::Real const perc_cut ) {
	core::Size n_poses = 0;
	core::Size num_res = 0;
	//import poses and calculate disorder at the end
	utility::vector1<core::pose::PoseOP> poses;
	core::import_pose::pose_stream::MetaPoseInputStream input = core::import_pose::pose_stream::streams_from_cmd_line();
	//Import each pose
	while ( input.has_another_pose() ) {
		n_poses++;
		core::pose::PoseOP mypose(utility::pointer::make_shared< core::pose::Pose >() );
		input.fill_pose( *mypose );
		if ( n_poses == 1 ) {
			num_res = mypose->total_residue();
			//Make sure poses have at least 10 residues
			check_pose_length( num_res, ORDER );
		} else {
			core::Size num_res_current = mypose->total_residue();
			//Make sure pose sizes match
			check_pose_size( num_res_current, num_res );
		}
		//Make sure they are all a single chain
		check_single_chain(mypose);
		poses.push_back(mypose);
	}
	//If 100 poses are not given as input, give warning.
	check_number_of_poses( n_poses );

	calculate_disorder( poses, num_res, n_poses, 1, vector_percent_disorder, ORDER, CUTOFF, term_cutoff, endpoint, perc_cut );

}

void measure_disorder ( utility::vector1<core::Real> &vector_percent_disorder, core::Size const ORDER, core::Real const CUTOFF, core::Real const term_cutoff, core::Real const endpoint, core::Real const perc_cut ) {
	core::Size n_poses = 0;
	core::Size num_res = 0;
	//import poses and calculate disorder for each pose
	core::import_pose::pose_stream::MetaPoseInputStream input = core::import_pose::pose_stream::streams_from_cmd_line();
	//Import each pose
	core::Size count = 0;
	while ( input.has_another_pose() ) {
		n_poses++;
		count ++;
		core::pose::PoseOP mypose(utility::pointer::make_shared< core::pose::Pose >() );
		input.fill_pose( *mypose );
		if ( n_poses == 1 ) {
			num_res = mypose->total_residue();
			//Make sure poses have at least 10 residues
			check_pose_length( num_res, ORDER );
		} else {
			core::Size num_res_current = mypose->total_residue();
			//Make sure pose sizes match
			check_pose_size( num_res_current, num_res );
		}
		//Make sure they are all a single chain
		check_single_chain(mypose);

		//calculate disorder for current pose
		utility::vector1<core::pose::PoseOP> pose_current;
		pose_current.push_back(mypose);
		calculate_disorder( pose_current, num_res, 1, count, vector_percent_disorder, ORDER, CUTOFF, term_cutoff, endpoint, perc_cut );

	}

}

void do_event_prediction ( utility::vector1<core::Real> const &vector_percent_disorder ) {
	core::Size num_frames = vector_percent_disorder.size();
	//make sure at least 50 frames input
	if ( num_frames < 50 ) utility_exit_with_message("Not enough frames to do event prediction. Must input at least 50 frames.");

	utility::vector1<core::Real> difference;
	core::Size window = 25;

	//calculate differences before and after frame
	for ( core::Size i=1; i<=num_frames; i++ ) {
		core::Real before = 0.0;
		core::Real after = 0.0;

		//calculate average percent disorder of 25 steps before
		if ( i <= window+1 ) {
			for ( core::Size j=1; j<i+1; j++ ) {
				before += vector_percent_disorder[j];
			}
			before = before / (core::Real)(i);
		} else {
			for ( core::Size j=i-window; j<=i; j++ ) {
				before += vector_percent_disorder[j];
			}
			before = before / (core::Real)(window+1);
		}

		//calculate average percent disorder of 25 steps after
		if ( i + window >= num_frames ) {
			for ( core::Size j=i; j<=num_frames; j++ ) {
				after += vector_percent_disorder[j];
			}
			after = after / (core::Real)(num_frames-i+1);
		} else {
			for ( core::Size j=i; j<=i+window; j++ ) {
				after += vector_percent_disorder[j];
			}
			after = after / (core::Real)(window+1);
		}

		//calculate disorder difference
		core::Real difference_current = after - before;

		difference.push_back(difference_current);
	}

	//calculate max difference
	core::Real max_abs_diff = 0.0;
	for ( core::Size i=1; i<=difference.size(); i++ ) {
		core::Real abs_diff_current = std::abs(difference[i]);
		if ( abs_diff_current > max_abs_diff ) {
			max_abs_diff = abs_diff_current;
		}
	}

	//predict folding and unfolding events
	core::Real cutoff = max_abs_diff*0.6;
	utility::vector1<core::Size> folding_events;
	utility::vector1<core::Size> unfolding_events;

	for ( core::Size i=1; i<=difference.size(); i++ ) {
		if ( difference[i] >= cutoff ) {
			unfolding_events.push_back(i);
		} else if ( difference[i] <= (-1.0)*cutoff ) {
			folding_events.push_back(i);
		}
	}

	//group events
	//folding
	utility::vector1<utility::vector1<core::Size>> folding_events_grouped;
	core::Size num_prev = 1;
	core::Size count = 0;
	for ( core::Size i=1; i<=folding_events.size(); i++ ) {
		if ( i==1 || num_prev != folding_events[i] - 1 ) {
			utility::vector1<core::Size> temp;
			temp.push_back( folding_events[i] );
			folding_events_grouped.push_back(temp);
			count = count +1;
		} else {
			folding_events_grouped[count].push_back( folding_events[i] );
		}
		num_prev = folding_events[i];
	}
	//unfolding
	utility::vector1<utility::vector1<core::Size>> unfolding_events_grouped;
	num_prev = 1;
	count = 0;
	for ( core::Size i=1; i<=unfolding_events.size(); i++ ) {
		if ( i==1 || num_prev != unfolding_events[i] - 1 ) {
			utility::vector1<core::Size> temp;
			temp.push_back( unfolding_events[i] );
			unfolding_events_grouped.push_back(temp);
			count = count +1;
		} else {
			unfolding_events_grouped[count].push_back( unfolding_events[i] );
		}
		num_prev = unfolding_events[i];
	}

	//combine events
	//folding
	utility::vector1<core::Size> folding_events_combined;
	for ( core::Size i=1; i<=folding_events_grouped.size(); i++ ) {
		core::Real average_event = 0.0;
		for ( core::Size j=1; j<=folding_events_grouped[i].size(); j++ ) {
			average_event += folding_events_grouped[i][j];
		}
		average_event = average_event / (core::Real)(folding_events_grouped[i].size());
		folding_events_combined.push_back(floor(average_event));
	}
	//unfolding
	utility::vector1<core::Size> unfolding_events_combined;
	for ( core::Size i=1; i<=unfolding_events_grouped.size(); i++ ) {
		core::Real average_event = 0.0;
		for ( core::Size j=1; j<=unfolding_events_grouped[i].size(); j++ ) {
			average_event += unfolding_events_grouped[i][j];
		}
		average_event = average_event / (core::Real)(unfolding_events_grouped[i].size());
		unfolding_events_combined.push_back(floor(average_event));
	}

	//output results
	std::ostringstream out;
	out << "Frame\tEvent" << std::endl;
	std::string outfile = "ResidueDisorder_events_default.out";
	if ( option[ out::file::o ].user() ) outfile = option[ out::file::o ]();
	for ( core::Size i=1; i <= num_frames; i++ ) {
		//Folding events
		if ( std::find(folding_events_combined.begin(), folding_events_combined.end(), i) != folding_events_combined.end() ) {
			out << ObjexxFCL::right_string_of( i, 6, '0' ) << "\tFolding" << std::endl;
		} else if ( std::find(unfolding_events_combined.begin(), unfolding_events_combined.end(), i) != unfolding_events_combined.end() ) { //Unfolding events
			out << ObjexxFCL::right_string_of( i, 6, '0' ) << "\tUnfolding" << std::endl;
		}
	}
	utility::io::ozstream outz( outfile.c_str() );
	outz << out.str();
	outz.close();
	outz.clear();
}

void define_variables ( core::Size &win_size, core::Real &cutoff, core::Real &term_cutoff, core::Real &term_endpoint, core::Real &perc_cut ) {
	if ( option[ corrections::restore_talaris_behavior ].user() ) {
		win_size = 5;
		term_cutoff = 0.13;
		perc_cut = 60.0;
		if ( option[ AF ].user() ) {
			TR << "Using AF with Talaris options! Make sure pdb's were relaxed in Rosetta before input!" << std::endl;
			cutoff = -1.2;
			term_endpoint = -1.2;
		} else {
			TR << "Using Rosetta with Talaris options!  Make sure pdb's were relaxed in Rosetta before input!" << std::endl;
			cutoff = -1.0;
			term_endpoint = -0.3;
		}
	} else {
		win_size = 10;
		term_cutoff = 0.34;
		perc_cut = 40.0;
		if ( option[ AF ].user() ) {
			TR << "Using AF with REF15 options!" << std::endl;
			cutoff = -1.8;
			term_endpoint = -1.8;
		} else {
			TR << "Using Rosetta with REF15 options!" << std::endl;
			cutoff = -1.5;
			term_endpoint = -0.8;
		}
	}
}

int
main( int argc, char * argv [] )
{

	try{
		using namespace basic::options;
		using namespace basic::options::OptionKeys;
		NEW_OPT( measure_disorder_from_structure, "Measure disorder separately for each structure from input", false );
		NEW_OPT( predict_events, "Predict folding/unfolding events from trajectory", false );
		NEW_OPT( AF, "Use the parameters for AlphaFold structure input", false );

		using namespace core::scoring;
		devel::init( argc, argv );

		//Setup parameters based on input flags
		core::Size win_size;
		core::Real cutoff;
		core::Real term_cutoff;
		core::Real term_endpoint;
		core::Real perc_cut;
		define_variables( win_size, cutoff, term_cutoff, term_endpoint, perc_cut );
		utility::vector1<core::Real> vec_percent_disorder;

		if ( option[ measure_disorder_from_structure ].user() ) {
			//Measure disorder (for each pose)
			TR.Warning << "Make sure pdb's were relaxed before input!" << std::endl;
			measure_disorder ( vec_percent_disorder, win_size, cutoff, term_cutoff, term_endpoint, perc_cut );

			//Predict events if flag given
			if ( option[ predict_events ].user() ) {
				if ( option[ corrections::restore_talaris_behavior ].user() ) {
					if ( option[ AF ].user() ) {
						utility_exit_with_message("Flag -predict_events given. Must be run using Rosetta talaris14 parameters. Make sure input structures were relaxed in talaris14 and remove the option -AF");
					} else {
						do_event_prediction( vec_percent_disorder );
					}
				} else {
					utility_exit_with_message("Flag -predict_events given. Must be run using talaris14. Make sure input structures were relaxed in talaris14 and give the option -restore_talaris_behavior");
				}
			}
		} else {
			if ( option[ predict_events ].user() ) utility_exit_with_message("Flag -predict_events given. Must be run in -measure_disorder_from_structure mode.");

			//Predict disorder (average over all poses)
			do_disorder_prediction ( vec_percent_disorder, win_size, cutoff, term_cutoff, term_endpoint, perc_cut );
		}

	} catch (utility::excn::Exception const & e ) {
		e.display();
		return -1;
	}

	return 0;
}
