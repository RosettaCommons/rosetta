// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet;
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file SewAppend.cc
///
/// @brief Application wrapper for the SEWING protocol. Takes a pre-generated score file and creates backbone assemblies from it
///
/// @author Tim Jacobs

//Package headers
#include <devel/init.hh>
#include <protocols/sewing/sampling/SewGraph.hh>

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

#include <basic/options/option.hh>
#include <basic/Tracer.hh>
#include <basic/database/sql_utils.hh>

//C++ headers
#include <map>
#include <set>
#include <iterator>

static basic::Tracer TR("SewAppend");

namespace SewAppend {
	basic::options::FileOptionKey const model_file( "model_file" );
	basic::options::FileOptionKey const starting_pdb( "starting_pdb" );
	basic::options::IntegerOptionKey const num_assemblies( "num_assemblies" );
	basic::options::IntegerVectorOptionKey const segment_starts( "segment_starts" );
	basic::options::IntegerVectorOptionKey const segment_ends( "segment_ends" );
	basic::options::IntegerVectorOptionKey const match_segments( "match_segments" );
	basic::options::BooleanOptionKey const generate_pdb_model_only( "generate_pdb_model_only" );
}

int
main( int argc, char * argv [] ) {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace protocols::sewing;

	option.add( SewAppend::model_file, "Model file name");
	option.add( SewAppend::num_assemblies, "Number of assemblies to generate and print");
	option.add( SewAppend::starting_pdb, "starting_pdb");
	option.add( SewAppend::segment_starts, "starting residue numbers of segments");
	option.add( SewAppend::segment_ends, "ending residue numbers of segments");
	option.add( SewAppend::match_segments, "which segments should we try to match?");
	option.add( SewAppend::generate_pdb_model_only, "Generate the model pdb and quit");

	// initialize core and read options
	devel::init(argc, argv);

	//Check for valid starting pdb
	if(!option[SewAppend::starting_pdb].user()) {
		utility_exit_with_message("You must supply a starting pdb file to the SewAppend app using the -starting_pdb flag");
	}
	std::string starting_pdb = option[SewAppend::starting_pdb];

	//Get all the models from the file
	if(!option[SewAppend::model_file].user()) {
		utility_exit_with_message("You must supply a model file to the SewAppend app using the -model_file flag.");
	}
	std::string model_filename = option[SewAppend::model_file];

	if(!option[SewAppend::segment_starts].user() || !option[SewAppend::segment_ends].user()) {
		utility_exit_with_message("You must supply segment starts and ends to the SewAppend app using the -segments_starts and -segment_ends flags");
	}
	utility::vector1<int> segment_starts = option[SewAppend::segment_starts];
	utility::vector1<int> segment_ends = option[SewAppend::segment_ends];
	if(segment_starts.size() != segment_ends.size()) {
		utility_exit_with_message("You must supply the same number of segment starts and ends");
	}

	//Figure out which segments from the starting PDB we need to match
	std::set<core::Size> segments_to_match;
	if(!option[SewAppend::match_segments].user()) {
		for(core::Size i=1; i<=segment_starts.size(); ++i){
			segments_to_match.insert(i);
		}
	}
	else {
		utility::vector1<int> temp = option[SewAppend::match_segments];
		for(core::Size i=1; i<=temp.size(); ++i){
			segments_to_match.insert(core::Size(temp[i]));
		}
	}

	//Generate the PDB model and add it to the model map
	std::map< int, Model > models;

	utility::vector1< std::pair<core::Size, core::Size> > segments;
	for(core::Size i=1; i<=segment_starts.size(); ++i) {
		segments.push_back(std::make_pair(segment_starts[i], segment_ends[i]));
	}
	core::pose::Pose pose;
	core::import_pose::pose_from_pdb( pose, starting_pdb );
	Model pdb_model = create_model_from_pose(pose, segments, models);
	models.insert(std::make_pair(pdb_model.model_id_, pdb_model));

	//If we are just dumping the model pdb, do it here and exit
	if(!option[SewAppend::generate_pdb_model_only].user() || option[SewAppend::generate_pdb_model_only]) {
		write_model_file(models, model_filename);
		TR << "Finished writing model file, Exiting!" << std::endl;
		std::exit(0);
	}


	std::map< int, Model > models_from_file = read_model_file(model_filename);
	models.insert(models_from_file.begin(), models_from_file.end());


	//Insert all the models into the hash table
	Hasher hasher;
	for(std::map< int, Model >::const_iterator it = models.begin(); it != models.end(); ++it) {
		hasher.insert(it->second);
	}
	TR << "Hashing complete" << std::endl;

	//Now, score the pdb model against all others
	TR << "Begin scoring" << std::endl;
	ScoreResults scores = hasher.score(pdb_model, 2, 10, 0, segments_to_match);
	TR << "Scoring complete" << std::endl;

	//Finally, generate the assemblies
	SewGraph sew_graph(models);
	sew_graph.generate_from_scores(scores);
	Assembly assembly = sew_graph.assemble_from_starting_model(pdb_model);
	utility::vector1< utility::vector1< SewSegment > > test = assembly.find_possible_orders(20.0);

	TR << "Found " << test.size() << " possible orders" << std::endl;
	for(core::Size i=1; i<=test.size(); ++i) {
		TR << "Order # " << i << std::endl;
		assembly.reorder(test[i]);
		for(core::Size j=1; j<=test[i].size(); ++j) {
			TR << test[i][j]->residues_[1].resnum_ << "-" << test[i][j]->residues_[test[i][j]->residues_.size()].resnum_ << " ";
		}
		TR << std::endl;
		core::pose::Pose assembly_pose = assembly.to_pose();
		assembly_pose.dump_pdb("assembly_pose_" + utility::to_string(i) + ".pdb");
	}

}
