// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet;
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file GivenPathAssemblyMover.cc
///
/// @brief
/// @author Tim Jacobs

// Unit Headers
#include <protocols/sewing/sampling/GivenPathAssemblyMover.hh>
#include <protocols/sewing/sampling/GivenPathAssemblyMoverCreator.hh>

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

namespace protocols {
namespace sewing  {

static basic::Tracer TR( "protocols.sewing.sampling.GivenPathAssemblyMover" );

////////////////////////////////////////////////////////////////////////////////////
///////////////////////////  Boiler Plate Code   ///////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
protocols::moves::MoverOP
GivenPathAssemblyMoverCreator::create_mover() const
{
	return protocols::moves::MoverOP( new GivenPathAssemblyMover );
}

std::string
GivenPathAssemblyMoverCreator::keyname() const
{
	return GivenPathAssemblyMoverCreator::mover_name();
}

std::string
GivenPathAssemblyMoverCreator::mover_name()
{
	return "GivenPathAssemblyMover";
}

protocols::moves::MoverOP
GivenPathAssemblyMover::clone() const {
	return( protocols::moves::MoverOP( new GivenPathAssemblyMover( *this ) ) );
}
protocols::moves::MoverOP
GivenPathAssemblyMover::fresh_instance() const {
	return protocols::moves::MoverOP( new GivenPathAssemblyMover );
}

std::string
GivenPathAssemblyMover::get_name() const {
	return "GivenPathAssemblyMover";
}

/////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////  GivenPathAssemblyMover function   //////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

GivenPathAssemblyMover::GivenPathAssemblyMover(){}

AssemblyOP
GivenPathAssemblyMover::generate_assembly(){
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	//Add first assembly
	AssemblyOP assembly;
	if ( option[ basic::options::OptionKeys::sewing::assembly_type ].value() == "discontinuous" ) {
		assembly = AssemblyFactory::create_assembly("discontinuous");
	} else if ( option[ basic::options::OptionKeys::sewing::assembly_type ].value() == "continuous" ) {
		assembly = AssemblyFactory::create_assembly("continuous");
	}

	//build from list of nodes
	ModelNode const * reference_node = nodes_[1];
	assembly->add_model(graph_, reference_node->model());
	for ( core::Size i=2; i<=nodes_.size(); i+=2 ) {
		reference_node = nodes_[i];
		runtime_assert(nodes_[i]->model().model_id_ == nodes_[i-1]->model().model_id_);
		runtime_assert(nodes_.size() >= i+1);
		graph_->add_edges_from_binary(edge_file_, reference_node->get_node_index());
		HashEdge const * const cur_edge =
			static_cast< HashEdge const * >(graph_->find_edge(nodes_[i]->get_node_index(), nodes_[i+1]->get_node_index()));
		assembly->follow_edge(graph_, cur_edge, reference_node->get_node_index());
	}
	return assembly;

	// std::set<core::Size> reference_node_ids = graph_->get_node_indices_from_model_id(model_ids_[1]);
	// ModelNode const * const first_node_id = graph_->get_model_node(*reference_node_ids.begin());
	// assembly->add_model(graph_, first_node_id->model());
	//
	// for(core::Size i=2; i<=model_ids_.size(); ++i) {
	//  bool edge_added = false;
	//  std::set<core::Size> target_node_ids = graph_->get_node_indices_from_model_id(model_ids_[i]);
	//
	//  std::set<core::Size>::const_iterator ref_it = reference_node_ids.begin();
	//  std::set<core::Size>::const_iterator ref_it_end = reference_node_ids.end();
	//  for(; ref_it != ref_it_end; ++ref_it) {
	//   ModelNode const * reference_node = graph_->get_model_node(*ref_it);
	//   graph_->add_edges_from_binary(edge_file_, reference_node->get_node_index());
	//   core::graph::EdgeListConstIterator edge_it = reference_node->const_edge_list_begin();
	//   core::graph::EdgeListConstIterator edge_it_end = reference_node->const_edge_list_end();
	//   for(; edge_it != edge_it_end; ++edge_it) {
	//    std::set<core::Size>::const_iterator target_it = target_node_ids.begin();
	//    std::set<core::Size>::const_iterator target_it_end = target_node_ids.end();
	//    for(; target_it != target_it_end; ++target_it) {
	//     HashEdge const * const cur_edge = static_cast< HashEdge const * >(*edge_it);
	//     TR << "Checking edge " << *ref_it << " " << *target_it << std::endl;
	//     if(cur_edge->get_other_ind(reference_node->get_node_index()) == *target_it) {
	//      cur_edge->get_second_node_ind();
	//      assembly->follow_edge(graph_, cur_edge, reference_node->get_node_index());
	//      edge_added = true;
	//     }
	//    }
	//   }
	//   if(!edge_added) {
	//    utility_exit_with_message("No edge found to model id " + utility::to_string(model_ids_[i]));
	//   }
	//  }
	//  reference_node_ids = target_node_ids;
	// }
	// return assembly;
}

void
GivenPathAssemblyMover::parse_my_tag(
	TagCOP const tag,
	basic::datacache::DataMap & data,
	protocols::filters::Filters_map const & filters,
	protocols::moves::Movers_map const & movers,
	core::pose::Pose const & pose
){
	parent::parse_my_tag(tag, data, filters, movers, pose);

	if ( tag->hasOption("path") ) {
		std::string const path_list( tag->getOption< std::string >( "path" ) );
		utility::vector1< std::string > const path_strings( utility::string_split( path_list,',' ) );
		for ( core::Size i=1; i<=path_strings.size(); ++i ) {
			model_ids_.push_back( utility::string2int(path_strings[i]) );
		}
	} else {
		utility_exit_with_message("You must specify a path in the GivenPathAssemblyMover tag");
	}
}

} //sewing
} //protocols
