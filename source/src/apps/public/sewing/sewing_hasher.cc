// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet;
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file SewingHasher.cc
///
/// @brief An MPI-enabled application that reads in a Hasher hashtable, scores each model in that table against all others, and
/// generates a SewGraph file.
///
/// @author Tim Jacobs

//Package headers
#include <devel/init.hh>
#include <protocols/sewing/hashing/Hasher.hh>
#include <protocols/sewing/conformation/Model.hh>
#include <protocols/sewing/util/io.hh>

//Protocol headers
#include <core/types.hh>
#include <core/id/AtomID.hh>
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>

//Utility headers
#include <utility/io/izstream.hh>
#include <utility/mpi_util.hh>
#include <utility/string_util.hh>
#include <utility/vector1.hh>
#include <utility/sql_database/DatabaseSessionManager.hh>
#include <utility/file/FileName.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/sewing.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/inout.OptionKeys.gen.hh>
#include <basic/Tracer.hh>
#include <basic/database/sql_utils.hh>

//C++ headers
#include <map>
#include <set>
#include <iterator>

static basic::Tracer TR("SewingHasher");

int
main( int argc, char * argv [] ) {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace protocols::sewing;

	try {
		// initialize core and read options
		devel::init(argc, argv);

		if ( !option[sewing::mode].user() ) {
			std::stringstream err;
			err << "You must provide a mode for SewingHasher to run in using the -sewing:mode flag. Valid options are" << std::endl;
			err << "  -generate: generates a model file from an sqlite database" << std::endl;
			err << "  -hash: score all models against each other and create a plain text score file (MPI required)" << std::endl;
			err << "  -convert: convert a plain text score file to a binary score file. This is required by the SEWING movers" << std::endl;
			utility_exit_with_message(err.str());
		}

		//Check for model file (either for reading or writing)
		std::map< int, Model > models;
		std::string model_filename = option[sewing::model_file_name];
		if ( !option[sewing::model_file_name].user() ) {
			std::stringstream err;
			err << "You must provide a model file name to the SewingHasher using the -model_file_name flag. To generate a new model file use the "
				<< "-generate_models_from_db flag and a new model file with that name will be written. Otherwise, the model file will be read";
			utility_exit_with_message(err.str());
		}

		//////////////////////////// MODEL GENERATION ///////////////////////////////////

		//Are we just generating a model file?
		if ( option[sewing::mode].value() == "generate" ) {

			//Create comments stream for model file and add the date
			std::stringstream comments;
			time_t t = time(0); // get time now
			struct tm * now = localtime( & t );
			comments << "#Model file created on " << (now->tm_year + 1900) << '-'
				<< (now->tm_mon + 1) << '-'
				<< now->tm_mday
				<< std::endl;

			//Generate model file from a list of PDBs. All models will have 1 segment
			if ( option[ in::file::l ].user() ) {
				comments << "#Models generated from PDB input (-l flag)" << std::endl;
				utility::vector1<utility::file::FileName> input_lists( option[ in::file::l ]() );
				utility::vector1<utility::file::FileName> pdb_library;
				for ( core::Size i = 1; i <= input_lists.size(); i++ ) {
					utility::io::izstream current_input_list( input_lists[i] );
					if ( !current_input_list.good() ) {
						utility_exit_with_message("unable to open input file file: "+input_lists[i].name()+"\n");
					}
					while ( current_input_list.good() ) {
						std::string name;
						current_input_list.getline(name);
						if ( current_input_list.good() ) pdb_library.push_back( utility::file::FileName(name) );
					}
				}

				for ( core::Size i=1; i<=pdb_library.size(); ++i ) {
					core::pose::Pose pose;
					core::import_pose::pose_from_pdb(pose, pdb_library[i]);
					utility::vector1< std::pair<core::Size,core::Size> > segments;
					segments.push_back(std::make_pair(1, pose.total_residue()));
					Model pdb_model = create_model_from_pose(pose, segments, (int)i);
					models.insert(std::make_pair(i, pdb_model));
				}

			} else {
				//Generate models from a features database. Each segment is a single piece of secondary structure
				if ( !option[ sewing::assembly_type ].user() ) {
					std::stringstream err;
					err << "You must provide an assembly_type (continuous or discontinuous) with the -sewing:assembly_type flag in order to extract models";
					utility_exit_with_message(err.str());
				}
				if ( option[ sewing::assembly_type ].value() == "discontinuous" ) {
					comments << "#Discontinuous models generated from sqlite database " << option[inout::dbms::database_name].value() << std::endl;
					models = get_discontinuous_models_from_db();
				} else if ( option[ sewing::assembly_type ].value() == "continuous" ) {
					comments << "#Continuous models generated from sqlite database " << option[inout::dbms::database_name].value() << std::endl;
					models = get_continuous_models_from_db();
				}
			}

			write_model_file(comments.str(), models, model_filename);
			TR << "New model file " << model_filename << " successfully written." << std::endl;
			std::exit(0);
		}


		//If we aren't generating models, then we need to read them
		models = read_model_file(model_filename);

		///////////////////////// BINARY FILE TESTING //////////////////////////////

		if ( option[sewing::mode].value() == "test" ) {
			std::string binary_filename = option[sewing::score_file_name].value();
			SewGraphOP graph( new SewGraph(models, 1) );
			graph->report_binary_stats(models, binary_filename);
		}

		///////////////////////// MODEL CONVERSION TO BINARY //////////////////////////////

		//If we are generating a binary file then do that and exit
		if ( option[sewing::mode].value() == "convert" ) {
			if ( !option[sewing::score_file_name].user() ) {
				std::stringstream err;
				err << "You must provide a score file name to the SewingHasher for binary conversion.";
				utility_exit_with_message(err.str());

			}
			std::string binary_filename = utility::to_string(option[sewing::score_file_name].value())+utility::to_string(".bin");

			SewGraphOP graph;
			if ( option[ sewing::assembly_type ].value() == "discontinuous" ) {
				graph = SewGraphOP( new SewGraph(models, 2) );
			} else if ( option[ sewing::assembly_type ].value() == "continuous" ) {
				graph = SewGraphOP( new SewGraph(models, 1) );
			}
			graph->generate_binary_score_file(option[sewing::score_file_name].value(), binary_filename);
			std::exit(0);
		}


		//////////// MODEL COMPARISON USING GEOMETRIC HASHING /////////////////

		if ( option[sewing::mode].value() == "hash" ) {

			if ( !option[sewing::score_file_name].user() ) {
				utility_exit_with_message("You must provide a graph file name to the SewingHasher using -score_file_name");
			}
			std::string score_file_name = option[sewing::score_file_name];

			core::Size min_score = option[sewing::min_hash_score].def(10);
			core::Size max_clash_score = option[sewing::max_clash_score].def(0);
			core::Size num_segments_to_match = option[sewing::num_segments_to_match].def(0);

			TR << "Bundle Hasher options:" << std::endl;
			TR << "\tMinimum Score: " << min_score << std::endl;
			TR << "\tMaximum Clash Score: " << max_clash_score << std::endl;
			TR << "\tNumber of segments to match: " << num_segments_to_match << std::endl;
			TR << "\tScore file name: " << score_file_name << std::endl;


#ifdef USEMPI
	    	//Now, score models. Split this work between multiple processors if we have them
	    	core::Size n_models;
	    	if ( option[sewing::max_models].user() ) {
	    		n_models = std::min( (core::Size)option[sewing::max_models], models.size() );
	    	} else {
	    		n_models = models.size();
	    	}

	    	int starting_model = 0;
	    	if ( option[sewing::starting_model].user() ) {
	    		starting_model = option[sewing::starting_model].value();
	    	}

	    	core::Size rank = utility::mpi_rank();
	    	if(rank == 0) {
	    		core::Size num_procs = utility::mpi_nprocs();

	    		TR << "Master node has " << n_models << " jobs to submit" << std::endl;

	    		//send out initial jobs
	    		std::map< int, Model >::const_iterator it = models.begin();
	    		++it; //Increment the iterator (no need to score the first model against just itself)

	    		std::map< int, Model >::const_iterator it_end = models.begin();
	    		std::advance(it_end, n_models);

	    		if(num_procs > n_models) {
	    			utility_exit_with_message("You have more processors than number of models to hash. Reduce your number of processors for efficiency");
	    		}

	    		core::Size curr_proc = 1;
	    		for(; it != it_end; ++it) {
	    			//send *it to curr_proc
	    			if(it->first >= starting_model) {
	    				utility::send_integer_to_node(curr_proc,it->first);
	    				TR << "Master node sent a job for model " << it->first << " to processor " << curr_proc << std::endl;
	    				++curr_proc;
	    				if(curr_proc == num_procs) { break; }
	    			}
	    		}

	    		++it;//Increment iterator once more to account for the last model sent in the above loop
	    		while( it != it_end ){
	    			//wait for a processor to send you a message
	    			//send more work
	    			core::Size received_node = utility::receive_integer_from_anyone();
	    			TR << "Master node received a message from processor " << received_node << std::endl;
	    			utility::send_integer_to_node(received_node,it->first);
	    			TR << "Master node sent a new job for model " << it->first << " to processor " << received_node << std::endl;
	    			++it;
	    		}

	    		//Once we're done with all the work, tell each processor to spin down
	    		core::Size counter = 1;
	    		while(counter != num_procs){
	    			core::Size received_node = utility::receive_integer_from_anyone();
	    			utility::send_integer_to_node(received_node,0);
	    			++counter;
	    		}
	    		TR << "Master node finished sending jobs" << std::endl;
	    	}

	    	else {
	    		while(true) {
	    			//Recieve message from master node
	    			int model_id = utility::receive_integer_from_node(0);
	    			if(model_id == 0){ break; }

	    			TR << "Processor " << rank << " received job for model number " << model_id << std::endl;

	    			//Breakup the models into groups of 1000 for scoring. This make the memory demand *much* less
	    			//with a minor hit to scoring speed.
	    			std::map< int, Model >::const_iterator it = models.begin();
	    			std::map< int, Model >::const_iterator it_end = models.find(model_id);
	    			while(true) {
	    				core::Size counter = 0;
	    				Hasher hasher;
	    				for(; it != it_end; ++it) {
	    					hasher.insert(it->second);
	    					++counter;
	    					if(counter > 100) {
	    						++it;
	    						break;
	    					}
	    				}
	    				ScoreResults group_scores = hasher.score(models[model_id], num_segments_to_match, min_score, max_clash_score, false);
	    				hasher.remove_connection_inconsistencies(models, group_scores);
	    				std::string node_score_file_name = score_file_name + "." + utility::to_string(rank);
	    				write_hashing_scores_to_file(group_scores, node_score_file_name);
	    				TR << "Processor " << rank << " done scoring. Found " << group_scores.size() << " valid comparisons" << std::endl;
	    				if(it == it_end) {
	    					break;
	    				}
	    			}

//	    			//Begin work
//	    			Hasher hasher;
//	    			ScoreResults scores;
//	    			for(; it != it_end; ++it) {
//	    				if(it->first == model_id) {
//	    					TR << "Processor " << rank << " begin scoring " << model_id << std::endl;
//	    					scores = hasher.score(it->second, 1, min_score, 0, false);
//	    					break;
//	    				}
//	    				hasher.insert(it->second);
//	    			}
//	    			hasher.remove_connection_inconsistencies(models, scores);

	    			//std::string node_score_file_name = score_file_name + "." + utility::to_string(rank);
	    			//write_hashing_scores_to_file(group_scores, node_score_file_name);

	    			//Send message back to master node
	    			utility::send_integer_to_node(0, rank);
	    		}
	    	TR << "Processor " << rank << " was told to stop working" << std::endl;
	    	}
	    	MPI_Barrier( MPI_COMM_WORLD ); // make all nodes reach this point together.
	    	MPI_Finalize();
#else
			Hasher hasher;
			std::map< int, Model >::const_iterator it1 = models.begin();
			std::map< int, Model >::const_iterator it_end = models.end();
			std::map< int, Model >::const_iterator it2 = models.begin();
			for ( ; it1 != it_end; ++it1 ) {
				ScoreResults scores;
				for ( ; it2 != it1; ++it2 ) {
					hasher.insert(it2->second);
				}
				scores = hasher.score(it1->second, num_segments_to_match, min_score, max_clash_score, true);
				hasher.remove_connection_inconsistencies(models, scores);
				TR << "Done scoring " << it2->first << " found " << scores.size() << " valid comparisons" << std::endl;
				if ( scores.size() > 0 && TR.Debug ) {
					BasisPair bp = scores.begin()->first;
					std::map< SegmentPair, core::Size > segment_matches = scores.begin()->second.segment_match_counts;
					std::map< SegmentPair, AtomMap > segment_matches2 = scores.begin()->second.segment_matches;
					TR.Debug << "After scoring." << std::endl;
					TR.Debug << "\tModels: " << bp.first.model_id << " " << bp.second.model_id << std::endl;
					TR.Debug << "\tBasis Residues: " << bp.first.resnum << " " << bp.second.resnum << std::endl;
					TR.Debug << "\tNumber of matched segments: " << segment_matches.size() << std::endl;
					std::map< SegmentPair, core::Size >::const_iterator it = segment_matches.begin();
					std::map< SegmentPair, core::Size >::const_iterator it_end = segment_matches.end();
					std::map< SegmentPair, AtomMap >::const_iterator it2 = segment_matches2.begin();
					for ( ; it != it_end; ++it ) {
						TR.Debug << "\tSegments " << it->first.first << " and " << it->first.second << " have " << it->second << " overlapping atoms." << std::endl;
						TR.Debug << "\tSegments " << it2->first.first << " and " << it2->first.second << " have " << it2->second.size() << " overlapping atoms." << std::endl;
						++it2;
					}
				}
				write_hashing_scores_to_file(scores, score_file_name);
			}
#endif
		}
	} catch ( utility::excn::EXCN_Base& excn ) {
		std::cerr << "Exception : " << std::endl;
		excn.show( std::cerr );
		return -1;
	}
	return 0;
}
