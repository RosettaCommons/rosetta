// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/sewing/sampling/EnumerateAssemblyMover.cc
///
/// @brief Enumerate every possible assembly exhaustively
/// @author Doonam Kim
/// @author Tim Jacobs

// Unit Headers
#include <protocols/sewing/sampling/EnumerateAssemblyMover.hh>
#include <protocols/sewing/sampling/EnumerateAssemblyMoverCreator.hh>

// Package Headers
#include <protocols/sewing/conformation/AssemblyFactory.hh>
#include <protocols/sewing/sampling/requirements/RequirementSet.hh>
#include <protocols/sewing/scoring/AssemblyScorer.hh>

//Utility headers
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/sewing.OptionKeys.gen.hh>

#include <numeric/random/random.hh>
#include <numeric/random/random_permutation.hh>

#include <utility/tag/Tag.hh>

#include <utility/vector1.hh>

#include <algorithm> // for min
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

namespace protocols {
namespace sewing  {

static THREAD_LOCAL basic::Tracer TR( "protocols.sewing.sampling.EnumerateAssemblyMover" );

////////////////////////////////////////////////////////////////////////////////////
///////////////////////////  Boiler Plate Code   ///////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
// XRW TEMP protocols::moves::MoverOP
// XRW TEMP EnumerateAssemblyMoverCreator::create_mover() const
// XRW TEMP {
// XRW TEMP  return protocols::moves::MoverOP( new EnumerateAssemblyMover );
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP EnumerateAssemblyMoverCreator::keyname() const
// XRW TEMP {
// XRW TEMP  return EnumerateAssemblyMover::mover_name();
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP EnumerateAssemblyMover::mover_name()
// XRW TEMP {
// XRW TEMP  return "EnumerateAssemblyMover";
// XRW TEMP }


protocols::moves::MoverOP
EnumerateAssemblyMover::clone() const {
	return( protocols::moves::MoverOP( new EnumerateAssemblyMover( *this ) ) );
}

protocols::moves::MoverOP
EnumerateAssemblyMover::fresh_instance() const {
	return protocols::moves::MoverOP( new EnumerateAssemblyMover );
}

// XRW TEMP std::string
// XRW TEMP EnumerateAssemblyMover::get_name() const {
// XRW TEMP  return "EnumerateAssemblyMover";
// XRW TEMP }

/////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////  EnumerateAssemblyMover function   //////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

EnumerateAssemblyMover::EnumerateAssemblyMover():
	bogus_var_for_constructor_(10),
	min_assembly_score_(-2.0),
	models()
{}


// @brief: add one N-term edge and one C-term edge to the current node, and check all requirements
AssemblyOP
EnumerateAssemblyMover::generate_assembly(){
	TR << "EnumerateAssemblyMover::generate_assembly()" << std::endl;
	using namespace basic::options;

	/* pseudo code for EnumerateAssemblyMover::follow_edge_from_node

	Load all models from specified model_file
	For each model {
	nodes = graph->get_model_nodes(model_id)
	(all) n-term_edges = graph->get_edges_for_node (nodes[0])
	(all) c-term_edges = graph->get_edges_for_node (nodes[2/4])

	assembly = new assembly(model_id)
	for (all n_term_edges) {
	assembly -> follow_edge(each n-term_edge)

	for (all c_term_edges) {
	assembly -> follow_edge(each c-term_edge)

	if(!requirement_set_->violates(assembly)) {
	dump_pdb_of_assembly (using existing fn of dumping pdb from assembly)
	}
	}
	}*/


	/////// <begin> get model_file, edge_file
	std::map< int, Model > models;
	if ( option[basic::options::OptionKeys::sewing::score_file_name].user() &&
			option[basic::options::OptionKeys::sewing::model_file_name].user() ) {
		edge_file_ = option[basic::options::OptionKeys::sewing::score_file_name].value();
		std::string model_file = option[basic::options::OptionKeys::sewing::model_file_name].value();
		models = read_model_file(model_file);
	} else {
		utility_exit_with_message("You must give a model file and score file to an EnumerateAssemblyMover either through options or tags");
	}
	/////// <end> get model_file, edge_file

	core::Size max_ss_num = option[basic::options::OptionKeys::sewing::max_ss_num].def(5);
	// so either 3 (smotif) or 5 (5-ss_models)

	/////////////////// with each model, try to add every possible combination of n-terminal edges and c-terminal edges

	std::map< int, Model >::iterator it = models.begin();
	std::map< int, Model >::iterator it_end = models.end();

	for ( ; it != it_end; ++it ) { // 5,420 models with top8k // 22,356 models with 17k pdbs
		Model const & cur_model = it->second;

		if ( TR.Debug.visible() ) {
			TR << "[analysis] Current model_id: " << utility::to_string(cur_model.model_id_) << std::endl;
		}
		if ( max_ss_num == 5 ) {
			//Want to start with b-a-b motif (this is not generalizable to things outside of what Doo Nam is doing!)
			if ( cur_model.segments_.size() != 5 ) {
				continue;
			}
		}

		//Initialize an empty starting Assembly, but this alone can't assign XYZs of current model to assembly
		AssemblyOP assembly = AssemblyFactory::create_assembly("continuous");

		std::set<core::Size> node_indices = graph_->get_node_indices_from_model_id(cur_model.model_id_);

		ModelNode const * n_term_node = graph_->get_model_node(*node_indices.begin());

		std::set<core::Size> n_node_segments = n_term_node->segment_ids();
		runtime_assert(n_node_segments.size() == 1);
		TR << "Segments in n-terminal node: " << *(n_node_segments.begin()) << std::endl;
		core::Size test_1 = *n_node_segments.begin();

		//Tim's'
		//ModelNode const * c_term_node = graph_->get_model_node(*node_indices.begin()+4);
		// error: use of undeclared identifier 'c_node_segments'

		ModelNode const * c_term_node;
		if ( max_ss_num == 5 ) {
			//ModelNode const * c_term_node = graph_->get_model_node(*node_indices.begin()+4);
			c_term_node = graph_->get_model_node(*node_indices.begin()+4);
			std::set<core::Size> c_node_segments = c_term_node->segment_ids();
			runtime_assert(c_node_segments.size() == 1);
			TR << "Segments in c-terminal node: " << *(c_node_segments.begin()) << std::endl;
			core::Size test_2 = *c_node_segments.begin();

			runtime_assert(test_2 - test_1 == 4);
		} else {
			//ModelNode const * c_term_node = graph_->get_model_node(*node_indices.begin()+2);
			c_term_node = graph_->get_model_node(*node_indices.begin()+2);
			std::set<core::Size> c_node_segments = c_term_node->segment_ids();
			runtime_assert(c_node_segments.size() == 1);
			TR << "Segments in c-terminal node: " << *(c_node_segments.begin()) << std::endl;
			core::Size test_2 = *c_node_segments.begin();

			runtime_assert(test_2 - test_1 == 2);
		}

		//ModelNode const * end_node = graph_->get_model_node(*node_indices.end());
		// but using this end_node resulted in "Didn't match on the first or last segment! OH NO!"

		if ( TR.Debug.visible() ) {
			//TR << "[analysis] Current model_id with 5 ss: " << utility::to_string(cur_model.model_id_) << std::endl;
			TR << "n_term_node->get_node_index(): " << n_term_node->get_node_index() << std::endl;
			TR << "c_term_node->get_node_index(): " << c_term_node->get_node_index() << std::endl;
			//TR << "end_node->get_node_index(): " << end_node->get_node_index() << std::endl;
		}

		assembly->add_model(graph_, n_term_node->model());
		// now current model is made as a starting assembly and dumping as a pdb shows appropriate xyzs well

		if ( option[basic::options::OptionKeys::sewing::dump_every_model].user() ) {
			core::pose::Pose to_be_dumped_pose = get_fullatom_pose(assembly);
			to_be_dumped_pose.dump_pdb( "model_id_" + utility::to_string(cur_model.model_id_) + ".pdb" );
			continue;
		}


		//////////////// (0) just for analysis purpose

		graph_->add_edges_from_binary(edge_file_, n_term_node->get_node_index());
		core::Size num_edges_from_n_term_node = n_term_node->num_edges();

		graph_->add_edges_from_binary(edge_file_, c_term_node->get_node_index());
		core::Size num_edges_from_c_term_node = c_term_node->num_edges();

		if ( TR.Debug.visible() ) {
			TR << "[analysis] num_edges_from_n_term_node: " << num_edges_from_n_term_node << std::endl;
			TR << "[analysis] num_edges_from_c_term_node: " << num_edges_from_c_term_node << std::endl;
			TR << "[analysis] possible combination (=num_edges_from_n_term_node*num_edges_from_c_term_node) with this model_id ("<< cur_model.model_id_ << "): " << num_edges_from_n_term_node*num_edges_from_c_term_node << std::endl;
		}
		if ( !num_edges_from_n_term_node || !num_edges_from_c_term_node ) {
			continue;
		}


		//////////////// (1) the first try to add n-terminal edges
		graph_->add_edges_from_binary(edge_file_, n_term_node->get_node_index());
		num_edges_from_n_term_node = n_term_node->num_edges();

		if ( TR.Debug.visible() ) {
			TR << "============ now this model is a candidate to add both n-term and c-term edge ============" << std::endl;
		}

		////////////TIM MODS/////////////
		for ( core::Size n_edge_index=1; n_edge_index<=num_edges_from_n_term_node; ++n_edge_index ) {
			if ( TR.Debug.visible() ) {
				TR << std::endl << "~~~~~~~~ n_edge_index: " << n_edge_index << " ~~~~~~~~" << std::endl;
			}
			graph_->add_edges_from_binary(edge_file_, n_term_node->get_node_index());

			utility::graph::EdgeListConstIterator n_edge_it = n_term_node->const_edge_list_begin();
			for ( core::Size foo=1; foo< n_edge_index; ++foo ) {
				++n_edge_it;
			}

			HashEdge const * const cur_edge_n = static_cast< HashEdge const * >(*n_edge_it);

			AssemblyOP original_assembly = assembly->clone();

			core::Size mobile_node_index_n = cur_edge_n->get_other_ind(n_term_node->get_node_index());
			core::Size mobile_model_id_n( graph_->get_model_node(mobile_node_index_n)->model().model_id_ );
			if ( TR.Debug.visible() ) {
				TR.Debug << "mobile_node_index_n for n-terminal addition: " << mobile_node_index_n << std::endl;
				TR.Debug << "mobile_model_id_n for n-terminal addition: " << mobile_model_id_n << std::endl;
			}

			assembly->follow_edge(graph_, cur_edge_n, n_term_node->get_node_index());
			AssemblyOP one_edge_assembly = assembly->clone();

			//Dump temporary pose for debugging
			// core::pose::Pose to_be_dumped_full_pose_after_n_term_add = get_fullatom_pose(assembly);
			// std::string pdb_name_after_n_term_add = utility::to_string(cur_model.model_id_) + "_" + utility::to_string(cur_edge_ind_n);
			// to_be_dumped_full_pose_after_n_term_add.dump_pdb ( pdb_name_after_n_term_add + ".pdb" );

			for ( core::Size c_edge_index=1; c_edge_index<=num_edges_from_c_term_node; ++c_edge_index ) {
				if ( TR.Debug.visible() ) {
					TR << std::endl << "~~~~~~~~ c_edge_index: " << c_edge_index << " ~~~~~~~~" << std::endl;
				}

				/////////////////// (2) the second try to add c-terminal edges
				graph_->add_edges_from_binary(edge_file_, c_term_node->get_node_index());
				num_edges_from_c_term_node = c_term_node->num_edges();

				utility::graph::EdgeListConstIterator c_edge_it = c_term_node->const_edge_list_begin();
				for ( core::Size bar=1; bar< c_edge_index; ++bar ) {
					++c_edge_it;
				}

				HashEdge const * const cur_edge_c = static_cast< HashEdge const * >(*c_edge_it);

				core::Size mobile_node_index_c = cur_edge_c->get_other_ind(c_term_node->get_node_index());
				core::Size mobile_model_id_c( graph_->get_model_node(mobile_node_index_c)->model().model_id_ );
				if ( TR.Debug.visible() ) {
					TR.Debug << "mobile_node_index_c for c-terminal addition: " << mobile_node_index_c << std::endl;
					TR.Debug << "mobile_model_id_c for c-terminal addition: " << mobile_model_id_c << std::endl;
					TR.Debug << "no checking of previously added model_id" << std::endl;
				}


				//////////////////// <begin> don't add already added model
				// (Tim Jacobs) If we've already added this model, don't add it again. This should theoretically
				//not be a problem, but due to the way model regeneration is handled in the Assembly
				//class, you get an error. This should be fixed in Assembly soon.
				// (Doonam) So far, adding already added model caused no problem for me, but use this for just in case
				std::set<core::Size> model_ids = assembly->model_ids();
				if ( model_ids.find(mobile_model_id_c) != model_ids.end() ) {
					assembly = one_edge_assembly;
					TR.Debug << "rejecting add, model previously added" << std::endl;
					continue;
				}
				//////////////////// <end> don't add already added model


				if ( TR.Debug.visible() ) {
					TR << "before cloning (one_edge_assembly = assembly->clone())" << std::endl;
					TR << "original_assembly->segments().size(): " << original_assembly->segments().size() << std::endl;
					TR << "one_edge_assembly->segments().size(): " << one_edge_assembly->segments().size() << std::endl;
					TR << "assembly->segments().size(): " << assembly->segments().size() << std::endl;
				}

				one_edge_assembly = assembly->clone();
				// It is super-weird but cloning should happen again "JUST BEFORE" 'c_term edge following' !
				// If the cloning happens before 'c-term edge' 'for loop' ONLY, then this cloning works as a shallow copy sometimes!

				if ( TR.Debug.visible() ) {
					TR << "before c_term follow_edge" << std::endl;
					TR << "original_assembly->segments().size(): " << original_assembly->segments().size() << std::endl;
					TR << "one_edge_assembly->segments().size(): " << one_edge_assembly->segments().size() << std::endl;
					TR << "assembly->segments().size(): " << assembly->segments().size() << std::endl;
				}

				// adding c-term edge to assembly
				assembly->follow_edge(graph_, cur_edge_c, c_term_node->get_node_index());

				if ( TR.Debug.visible() ) {
					TR << "after c_term follow_edge" << std::endl;
					TR << "original_assembly->segments().size(): " << original_assembly->segments().size() << std::endl;
					TR << "one_edge_assembly->segments().size(): " << one_edge_assembly->segments().size() << std::endl;
					TR << "assembly->segments().size(): " << assembly->segments().size() << std::endl;
				}

				if ( requirement_set_->violates(assembly) ) {
					if ( TR.Debug.visible() ) {
						TR.Debug << "rejecting add, violation of requirement_set" << std::endl;
					}
					assembly = one_edge_assembly;
					continue;
				}

				core::Real const score = assembly_scorefxn_->score(assembly);

				if ( score > min_assembly_score_ ) {
					if ( TR.Debug.visible() ) {
						TR.Debug << "rejecting add" << std::endl;
						TR.Debug << "score of assembly ( " << score << " ) > min_assembly_score ( " << min_assembly_score_ << " ) " << std::endl;
					}
					assembly = one_edge_assembly;
					continue;
				}

				//Dump the full pose
				std::string pdb_name = utility::to_string(cur_model.model_id_) + "_" + utility::to_string(n_edge_index) + "_" + utility::to_string(c_edge_index);
				TR << "============================ now dump " << pdb_name << " ===========================" << std::endl;


				core::pose::Pose to_be_dumped_full_pose = get_fullatom_pose(assembly);
				to_be_dumped_full_pose.dump_pdb( pdb_name + ".pdb" );

				output_stats(assembly, to_be_dumped_full_pose, pdb_name);

				assembly = one_edge_assembly;

			}//c-terminal add
			assembly = original_assembly;
		}//n-terminal add

	}// all models' for loop

	// intentionally makes this "0" return so that AssemblyMover doesn't crash
	AssemblyOP bogus_assembly;

	//AssemblyOP bogus_assembly = AssemblyFactory::create_assembly("continuous");
	return bogus_assembly;

} //generate_assembly



void
EnumerateAssemblyMover::parse_my_tag(
	TagCOP const tag,
	basic::datacache::DataMap & data,
	protocols::filters::Filters_map const & filters,
	protocols::moves::Movers_map const & movers,
	core::pose::Pose const & pose
){
	parent::parse_my_tag(tag, data, filters, movers, pose); // parent is AssemblyMover
	// indeed this "parent::parse_my_tag" is essential. For example, user can specify max_segments in parser.xml
	// inheriting classes can't recognize TR << "max_segments: " << max_segments << std::endl;

	if ( tag->hasOption("min_assembly_score") ) {
		min_assembly_score_ = tag->getOption<core::Real>("min_assembly_score"); // for Doonam's a/b design, -2.0 is recommended
	}

} //parse_my_tag

std::string EnumerateAssemblyMover::get_name() const {
	return mover_name();
}

std::string EnumerateAssemblyMover::mover_name() {
	return "EnumerateAssemblyMover";
}

void EnumerateAssemblyMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	// TO DO!
	using namespace utility::tag;
	AttributeList attlist; // TO DO: add attributes to this list
	// TO DO: perhaps this is not the right function to call? -- also, delete this comment
	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "XRW TO DO", attlist );
}

std::string EnumerateAssemblyMoverCreator::keyname() const {
	return EnumerateAssemblyMover::mover_name();
}

protocols::moves::MoverOP
EnumerateAssemblyMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new EnumerateAssemblyMover );
}

void EnumerateAssemblyMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	EnumerateAssemblyMover::provide_xml_schema( xsd );
}


} //sewing
} //protocols
