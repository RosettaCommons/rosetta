// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file apps/public/analysis/ResidueDisorder
/// @brief Application to run ResidueDisorder (predicting disordered regions of proteins).
/// @author Justin Seffernick (seffernick.9@osu.edu)

#include <basic/options/option.hh>
#include <devel/init.hh>
#include <basic/Tracer.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <utility/io/izstream.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/util.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/Energies.hh>
#include <core/import_pose/pose_stream/MetaPoseInputStream.hh>
#include <core/import_pose/pose_stream/util.hh>
#include <core/scoring/ScoreType.hh>
#include <utility/io/ozstream.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/corrections.OptionKeys.gen.hh>

static basic::Tracer TR( "apps.public.analysis.ResidueDisorder" );

void check_pose_size( core::Size num_res_current, core::Size num_res) {
	if ( num_res_current != num_res ) {
		std::stringstream err_msg;
		err_msg << "Pose size mismatch!! All poses must be of equal length." << std::endl;
		utility_exit_with_message(err_msg.str());
	}
}

void check_single_chain( const core::pose::PoseOP &mypose ) {
	core::Size num_chains = mypose->conformation().num_chains();
	if ( num_chains != 1 ) {
		std::stringstream err_msg;
		err_msg << "Multiple chains!! All poses must be single chains." << std::endl;
		utility_exit_with_message(err_msg.str());
	}
}

void check_sequence( const core::pose::PoseOP &pose, const core::pose::PoseOP &pose_current ) {
	if ( pose_current->sequence() != pose->sequence() ) {
		std::stringstream err_msg;
		err_msg << "Sequence mismatch!! All poses must have the same sequence." << std::endl;
		utility_exit_with_message(err_msg.str());
	}
}

void check_number_of_poses ( core::Size n_poses ) {
	if ( n_poses != 100 ) {
		TR.Warning << n_poses << " pdb's entered. Not recommended!! (100 recommended)" << std::endl;
	}
}

void check_score_function ( core::scoring::ScoreFunctionOP sfxn ) {
	std::string sfxn_name = sfxn->get_name();
	if ( sfxn_name != "talaris2014.wts" ) {
		TR.Warning << "talaris2014 not used! Flag -restore_talaris_behavior should be given!" << std::endl;
	}
}

void read_in_pdbs( utility::vector1<core::pose::PoseOP> &poses, core::Size &n_poses, core::Size &num_res ) {
	core::import_pose::pose_stream::MetaPoseInputStream input = core::import_pose::pose_stream::streams_from_cmd_line();
	//Import each pose
	while ( input.has_another_pose() ) {
		n_poses++;
		core::pose::PoseOP mypose( new core::pose::Pose );
		input.fill_pose( *mypose );
		if ( n_poses == 1 ) {
			num_res = mypose->total_residue();
		} else {
			core::Size num_res_current = mypose->total_residue();
			//Make sure pose sizes match
			check_pose_size( num_res_current, num_res);
		}
		//Make sure they are all a single chain
		check_single_chain(mypose);
		poses.push_back(mypose);
	}
	//If 100 poses are not given as input, give warning.
	check_number_of_poses( n_poses );
}

void calculate_average_residue_scores ( const utility::vector1<core::pose::PoseOP> &poses, utility::vector1<core::Real> &average_res_scores, std::string &sequence, core::Size num_res, core::Size n_poses ) {
	core::scoring::ScoreFunctionOP sfxn = core::scoring::get_score_function();
	check_score_function( sfxn );
	//Loop through each residue in the pose
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

void caculate_order_scores( utility::vector1<core::Real> &ORDER_RES_Score, const utility::vector1<core::Real> &average_res_scores, core::Size num_res ) {
	core::Size ORDER = 5;
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

void predict_disorder( utility::vector1<std::string> &PREDICTION, const utility::vector1<core::Real> &ORDER_RES_Score, const std::string &sequence, core::Size num_res, const utility::vector1<core::Real> &cutoffs, bool terminal ) {
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

core::Real calculate_percent_disorder( const utility::vector1<std::string> &PREDICTION, core::Size num_res ) {
	core::Real num_disorder;
	num_disorder = std::count(PREDICTION.begin(), PREDICTION.end(), "D");
	core::Real percent_disorder = num_disorder/num_res*100.0;
	TR << "Percent of disordered residues: " << percent_disorder << "%" << std::endl;
	return percent_disorder;
}

void do_terminal_opt( utility::vector1<core::Real> &cutoffs_new, core::Size num_res ) {
	core::Size TerminalCount = num_res*0.13;
	core::Real endpoint = -0.3;
	core::Real CUTOFF = -1.0;
	//Loop through each residue
	for ( core::Size i=1; i<=num_res; i++ ) {
		if ( i<=TerminalCount ) {
			//Case 1: within TerminalCount of beginning
			cutoffs_new[i] = (CUTOFF-endpoint)/(TerminalCount-1)*(i) + endpoint -(CUTOFF-endpoint)/(TerminalCount-1);
		} else if ( i>num_res-TerminalCount ) {
			//Case 2: within TerminalCount of end
			cutoffs_new[i] = (endpoint-CUTOFF)/(TerminalCount-1)*(i+TerminalCount-num_res) + CUTOFF - (endpoint-CUTOFF)/(TerminalCount-1);
		}
		//Case 3: not within TerminalCount of either end => Do not change cutoff
	}
}

void output_results ( core::Size num_res, const std::string &sequence, const utility::vector1<std::string> &PREDICTION, const utility::vector1<core::Real> &ORDER_RES_Score, utility::vector1<core::Real> average_res_scores) {
	std::ostringstream out;
	out << "Residue_number\tResidue_type\tPrediction\tOrder_score\tAverage_raw_score" << std::endl;
	for ( core::Size i=1; i<=num_res; i++ ) {
		out << i << "\t\t" << sequence.substr(i-1,1) << "\t\t" << PREDICTION[i] << "\t\t" << ORDER_RES_Score[i] << "\t" << average_res_scores[i] << std::endl;
	}
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	std::string outfile = "ResidueDisorder_default.out";
	if ( option[ out::file::o ].user() ) outfile = option[ out::file::o ]();
	utility::io::ozstream outz( outfile.c_str() );
	outz << out.str();
	outz.close();
	outz.clear();
}



int
main( int argc, char * argv [] )
{

	try{
		using namespace basic::options;
		using namespace basic::options::OptionKeys;
		using namespace core::scoring;
		devel::init( argc, argv );

		//Import poses
		core::Size n_poses = 0;
		core::Size num_res;
		utility::vector1<core::pose::PoseOP> poses;
		read_in_pdbs( poses, n_poses, num_res );

		utility::vector1<core::Real> average_res_scores;
		std::string sequence;

		//Calculate average residue scores
		calculate_average_residue_scores ( poses, average_res_scores, sequence, num_res, n_poses );

		//Calculate the order scores (i.e. window average score of each residue)
		utility::vector1<core::Real> ORDER_RES_Score;
		caculate_order_scores( ORDER_RES_Score, average_res_scores, num_res );

		//Use order scores to predict disorder
		utility::vector1<std::string> PREDICTION;
		utility::vector1<core::Real> cutoffs ( num_res, -1.0);
		predict_disorder( PREDICTION, ORDER_RES_Score, sequence, num_res, cutoffs, false );

		//Calculate percent disorder
		core::Real percent_disorder = calculate_percent_disorder( PREDICTION, num_res );

		//If less than 60% disordered, run the optimized terminal residue prediction
		if ( percent_disorder<=60.0 ) {
			TR << "Protein less than 60% disordered, optimizing terminal residues:" << std::endl;
			//Do terminal optimization->recalculate cutoffs
			utility::vector1<core::Real> cutoffs_new ( num_res, -1.0);
			do_terminal_opt( cutoffs_new, num_res );
			utility::vector1<std::string> PREDICTION_new;
			//Re-predict disorder with new cutoffs
			predict_disorder( PREDICTION_new, ORDER_RES_Score, sequence, num_res, cutoffs_new, true );
			//Calculate percent disorder
			calculate_percent_disorder( PREDICTION_new, num_res );
			//Output results of terminal prediction
			output_results( num_res, sequence, PREDICTION_new, ORDER_RES_Score, average_res_scores);
		} else {
			TR << "Protein more than 60% disordered, not optimizing terminal residues!" << std::endl;
			//Output results of prediction
			output_results( num_res, sequence, PREDICTION, ORDER_RES_Score, average_res_scores);
		}
	} catch (utility::excn::Exception const & e ) {
		TR << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;
}
