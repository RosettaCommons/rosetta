// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file DatabaseEntryWorkUnit.cc
///
/// @brief A work unit that runs a database query, processes the results, and returns a string (presumably a database insert statement)

/// @author Tim Jacobs

//Package headers
#include <devel/init.hh>
#include <protocols/legacy_sewing/hashing/Hasher.hh>
#include <protocols/legacy_sewing/conformation/Model.hh>
#include <protocols/legacy_sewing/conformation/ContinuousAssembly.hh>
#include <protocols/legacy_sewing/conformation/DisembodiedAssembly.hh>
#include <protocols/legacy_sewing/util/io.hh>

//Protocol headers
#include <core/pose/util.hh>
#include <core/types.hh>

//Utility headers
#include <basic/options/option.hh>
#include <basic/Tracer.hh>

#include <numeric/random/random.hh>

#include <utility/vector1.hh>
#include <utility/string_util.hh>

//C++ headers
#include <map>

static basic::Tracer TR("ModelDumper");

int
main( int argc, char * argv [] ) {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace protocols::legacy_sewing;

	// initialize core and read options
	devel::init(argc, argv);
	std::string model_file_name = option[legacy_sewing::model_file_name];
	core::Size num_models_to_dump = option[legacy_sewing::num_models_to_dump].def(10);
	std::map< int, Model > models = read_model_file(model_file_name).second;
	TR << "Done reading " << models.size() << " models from file: " << model_file_name << std::endl;

	if ( !option[legacy_sewing::score_file_name].user() ) {
		if ( option[legacy_sewing::models_to_dump].user() ) {
			utility::vector1<core::Size> model_ids = option[legacy_sewing::models_to_dump].value();
			for ( core::Size i=1; i<=model_ids.size(); ++i ) {
				if ( models.find(model_ids[i]) == models.end() ) {
					utility_exit_with_message("Couldn't find model with id: " + utility::to_string(model_ids[i]));
				}
				AssemblyOP model_assembly(new ContinuousAssembly());
				model_assembly->add_model(0, models.find(model_ids[i])->second, false);
				model_assembly->to_pose(core::chemical::FA_STANDARD).dump_pdb("model_" + utility::to_string(model_ids[i]) + ".pdb");
			}
		} else {
			for ( core::Size i=1; i<=num_models_to_dump; ++i ) {
				core::Size rand = numeric::random::random_range(0, models.size()-1);
				std::map< int, Model >::const_iterator it = models.begin();
				std::advance(it, rand);
				AssemblyOP model_assembly(new ContinuousAssembly());
				model_assembly->add_model(0, it->second, false);
				model_assembly->to_pose(core::chemical::FA_STANDARD).dump_pdb("model_" + utility::to_string(it->second.model_id_) + ".pdb");
			}
		}
	} else {
		std::string score_file_name = option[legacy_sewing::score_file_name];
		TR << "Done reading scores from file: " << score_file_name << std::endl;

		SewGraphOP graph;
		AssemblyOP two_model_assembly;
		if ( option[ legacy_sewing::assembly_type ].value() == "discontinuous" ) {
			graph = SewGraphOP(new SewGraph(models, 2));
			two_model_assembly = AssemblyOP(new DisembodiedAssembly());
		} else if ( option[ legacy_sewing::assembly_type ].value() == "continuous" ) {
			graph = SewGraphOP(new SewGraph(models, 1));
			two_model_assembly = AssemblyOP(new ContinuousAssembly());
		}

		//pick a random node
		ModelNode const * reference_node = graph->get_random_node();
		graph->add_edges_from_binary(score_file_name, reference_node->get_node_index());
		two_model_assembly->add_model(graph, reference_node->model());
		core::Size num_edges = reference_node->num_edges();
		TR << "Found " << num_edges << " for node " << reference_node->get_node_index() << " " << reference_node->model().model_id_ << std::endl;
		if ( num_edges > 0 ) {
			utility::graph::EdgeListConstIterator edge_it = reference_node->const_edge_list_begin();
			HashEdge const * const cur_edge = static_cast< HashEdge const * >(*edge_it);
			TR << "Following edge " << *cur_edge << std::endl;
			two_model_assembly->follow_edge(graph, cur_edge, reference_node->get_node_index());
			two_model_assembly->to_multichain_pose(core::chemical::FA_STANDARD).dump_pdb("align_" + utility::to_string(reference_node->model().model_id_) + "_multi.pdb");
			two_model_assembly->to_pose(core::chemical::FA_STANDARD).dump_pdb("align_" + utility::to_string(reference_node->model().model_id_) + ".pdb");
		}
	}
}
