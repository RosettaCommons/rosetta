// -*- mode:c++;tab-width:2;incdent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file apps/public/analysis/SID_ERMS_rescore
/// @brief Application to rescore a structure based on ERMS data
/// @author Robert Bolz (bolz.13@osu.edu)

#include <iostream>
#include <fstream>
#include <sstream> 
#include <string>
#include <iterator>
#include <tuple>
#include <devel/init.hh>
#include <basic/options/option.hh>
#include <utility/pointer/owning_ptr.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/Pose.hh>
#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.hh>

#include <basic/Tracer.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/chains_util.hh>

#include <core/scoring/Energies.hh>
#include <core/import_pose/pose_stream/MetaPoseInputStream.hh>
#include <core/import_pose/pose_stream/util.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/corrections.OptionKeys.gen.hh>
#include <numeric/random/random.hh>
#include <numeric/random/uniform.hh>
#include <protocols/analysis/InterfaceAnalyzerMover.hh>
#include <cmath>

#include <core/pose/metrics/PoseMetricCalculatorBase.hh>
#include <protocols/pose_metric_calculators/SaltBridgeCalculator.hh>
#include <core/pose/metrics/CalculatorFactory.hh>
#include <core/pose/metrics/PoseMetricContainer.hh>

#include <boost/config.hpp>
#include <vector>
#include <boost/graph/connected_components.hpp>
#include <boost/graph/adjacency_list.hpp>

#include <utility/vectorL.hh>
#include <utility/options/OptionCollection.hh>

#include <protocols/sid_erms_prediction/sid_erms_simulate.hh>
#include <basic/options/keys/sid_erms_simulate.OptionKeys.gen.hh>

static basic::Tracer TR( "src.apps.public.analysis.robert_bolz.SID_ERMS_rescore" );

using namespace basic::options;
using namespace basic::options::OptionKeys;

//calculate the rescored model
core::Real Rescore_models( const core::Real nscore_, const core::Real RMSE_, const core::Real max_, const core::Real min_ ){
	//score structure
	core::Real rescore, midpnt, e_exp;
	midpnt = (max_ + min_ )/2;
	e_exp = (50 * RMSE_) - (0.9 * 50 * midpnt);
	rescore = 1000 * ( 1 - ( 1 /(0.1 + ( exp( e_exp ) ) ) ) );
	if (rescore < 0) {
		rescore = 0;
	}

	rescore += nscore_;
	return rescore;
}

int
main( int argc, char * argv [] )
{

	try{
		using namespace basic::options;
		using namespace basic::options::OptionKeys;
		using namespace basic::options::OptionKeys::sid_erms_simulate;
		using namespace core::scoring;
		using namespace utility;
		using namespace core;
		using namespace core::pose;
		using namespace protocols;
		using namespace protocols::sid_erms_prediction;

		devel::init( argc, argv );

	//checking and reading all inputs

		core::Real A_(0.0025);

		//breakage cutoff
		core::Real breakage_cutoff_(0.0);

		//complex type (nodes and edges)
		core::Size n_chains_;
		utility::vector1<char> nodes_;
		utility::vector1<utility::vector1<std::string>> edges_;
		if ( option[ sid_erms_simulate::complex_type ].user() ) {
			//read_complex_type(n_chains_, nodes_, edges_);
			std::string complex_filename_;
			complex_filename_ = option[ sid_erms_simulate::complex_type ]();
			tie(n_chains_, nodes_, edges_) = protocols::sid_erms_prediction::read_complex_type(complex_filename_, n_chains_, nodes_, edges_);
			TR << "check outputs " << n_chains_ << std::endl;
		//} else {
			//std::stringstream err_msg;
			//err_msg << "Must input the complex type using the complex_type option." << std::endl;
			//utility_exit_with_message(err_msg.str());
		}
		TR << "Complex type read in" << std::endl;

		//create score function object
		core::scoring::ScoreFunctionOP sfxn = core::scoring::get_score_function();

		std::ostringstream out;

		out << "Structure" << "\t" << "RMSE" << "\t" << "Ref15_Score" << "\t" << "ERMS Rescore" << std::endl;

		utility::vector1<core::Real> rmse_list;
		utility::vector1<core::Real> score_list;

		//pose specific input checks

		utility::vector1<std::string> files;
		utility::vector1<file::FileName> list;
		if ( option[in::file::l].user() ) {
			TR << "using -l option " << std::endl;
			list = basic::options::option[ in::file::l ]();
			for ( auto const & listfile : list ) {
				utility::io::izstream pdbs( listfile );
				std::string fname;
				while ( pdbs >> fname ) {
					files.push_back(fname);
				}
			}
			for ( auto const & file : files ) {
				TR << "creating pose for " << file << std::endl;
				pose::Pose pose_in;
				core::import_pose::pose_from_file(pose_in, file, core::import_pose::PDB_file);
				TR << "created pose for " << file << std::endl;
				//std::string output = pose_in->pdb_info()->name();

		//set up vector of B values. Length is number of different interface types
				TR << "edges size is: " << edges_.size() << std::endl;
				utility::vector1<core::Real> B_(edges_.size(), 0.0);
				if ( option[ basic::options::OptionKeys::sid_erms_simulate::B_vals ].user() ) {
				//read in B values if given by file input
					B_ = protocols::sid_erms_prediction::read_B_vals(B_);
					TR << "B values read in from file" << std::endl;
					if ( option[ basic::options::OptionKeys::in::file::l ].user() ) {
						TR.Warning << "PDB file input not used. B values read from file instead. To calculate B from structure, remove the flag (B_vals)." << std::endl;
					}
				} else {
					if ( !option[ basic::options::OptionKeys::in::file::l ].user() ) {
						std::stringstream err_msg;
						err_msg << "Must input either file containing B values (using flag B_vals) or PDB file to calculate B values for simulation." << std::endl;
						utility_exit_with_message(err_msg.str());
					}
					B_ = protocols::sid_erms_prediction::calc_B_values(B_, edges_, A_, n_chains_, pose_in); //read in PDB and calculate B values for each the first in each interface type
				}
				TR << "B values read in from structure" << std::endl;

		//read in ERMS data: acceleration energies and ERMS data (if applicable)
			utility::vector1<core::Real> ACE_; //Acceleration energies, read from file or input automatically
			utility::vector1<utility::vector1<core::Real>> ERMS_;
			if ( option[ basic::options::OptionKeys::sid_erms_simulate::ERMS ].user() ) {
				tie(ACE_, ERMS_) = protocols::sid_erms_prediction::read_ERMS(ACE_, ERMS_, n_chains_);
				TR << "Acceleration energies and ERMS read in from file" << std::endl;
			} else {
			//if no acceleration energies are input, the range is set from 0 to 2x max B, with 10 steps in between
				core::Real step_size = (*max_element(B_.begin(), B_.end()))*2/10;
				ACE_.push_back(0.0);
				for ( core::Size i=1; i<=10; i++ ) {
					ACE_.push_back(i*step_size);
				}
				TR << "Acceleration energies set up based on max B" << std::endl;
			}

			//simulate ERMS
			utility::vector1<utility::vector1<core::Real>> ERMS_prediction_(n_chains_);
			ERMS_prediction_ = protocols::sid_erms_prediction::simulate_ERMS(ERMS_prediction_, ACE_, B_, n_chains_, edges_, A_, breakage_cutoff_, nodes_);
			TR << "ERMS simulation complete" << std::endl;

			//calculte RMSE (if applicable)
			if ( option[ basic::options::OptionKeys::sid_erms_simulate::RMSE ].user() && option[ basic::options::OptionKeys::sid_erms_simulate::ERMS ].user() ) {
				protocols::sid_erms_prediction::calc_RMSE(ERMS_prediction_, ERMS_, n_chains_, ERMS_.size());
			}
			if ( option[ basic::options::OptionKeys::sid_erms_simulate::RMSE ].user() && !option[ basic::options::OptionKeys::sid_erms_simulate::ERMS ].user() ) {
				TR.Warning << "RMSE cannot be calculated without ERMS data input." << std::endl;
			}
	
			//get score of structure
			core::Real score_ = sfxn->score( pose_in );

			//output results to file
			core::Real RMSE_value_(0.0); 
			RMSE_value_ = protocols::sid_erms_prediction::calc_RMSE(ERMS_prediction_, ERMS_, n_chains_, ERMS_.size());
			rmse_list.push_back(RMSE_value_);
			score_list.push_back(score_);
			}
		}
		//Determine maximum and minimum RMSE
		core::Real max_rmse(0.0);
		core::Real min_rmse(1.0);

		for ( core::Size i = 1; i <= rmse_list.size(); i++ ) {
        	if ( max_rmse < rmse_list[i] ) {
            	max_rmse = rmse_list[i];
        	}
        }
		for ( core::Size i = 1; i <= rmse_list.size(); i++ ) {
        	if ( min_rmse > rmse_list[i] ) {
            	min_rmse = rmse_list[i];
        	}        
		}  

        //Calculate rescore for each structure
		for ( core::Size i=1; i <= rmse_list.size(); i++ ) {
			core::Real Rescored_value_(0.0); 
			Rescored_value_ = Rescore_models(score_list[i], rmse_list[i], max_rmse, min_rmse);
			out << files[i] << "\t" << rmse_list[i] << "\t" << score_list[i] << "\t" << Rescored_value_ << std::endl;
		}

		//output rescores as a file
		using namespace basic::options;
		using namespace basic::options::OptionKeys;
		std::string outfile("output_rescore.tsv");
		if ( option[ out::file::o ].user() ) outfile = option[ out::file::o ]();
		utility::io::ozstream outz( outfile.c_str() );
		outz << out.str();
		outz.close();
		outz.clear();
		TR << "Output results complete." << std::endl;
	}

	catch (utility::excn::Exception const & e ) {
		e.display();
		return -1;
	}

	return 0;
}
