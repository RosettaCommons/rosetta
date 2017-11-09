// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file LegacyGivenPathAssemblyMover.cc
///
/// @brief
/// @author Tim Jacobs

// Unit Headers
#include <protocols/legacy_sewing/sampling/LegacyGivenPathAssemblyMover.hh>
#include <protocols/legacy_sewing/sampling/LegacyGivenPathAssemblyMoverCreator.hh>

// Package Headers
#include <protocols/legacy_sewing/conformation/AssemblyFactory.hh>
#include <protocols/legacy_sewing/sampling/requirements/LegacyRequirementSet.hh>
#include <protocols/legacy_sewing/scoring/LegacyAssemblyScorer.hh>

//Utility headers
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/legacy_sewing.OptionKeys.gen.hh>

#include <numeric/random/random.hh>
#include <numeric/random/random_permutation.hh>

#include <utility/tag/Tag.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

namespace protocols {
namespace legacy_sewing  {

static THREAD_LOCAL basic::Tracer TR( "protocols.legacy_sewing.sampling.LegacyGivenPathAssemblyMover" );

////////////////////////////////////////////////////////////////////////////////////
///////////////////////////  Boiler Plate Code   ///////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
// XRW TEMP protocols::moves::MoverOP
// XRW TEMP LegacyGivenPathAssemblyMoverCreator::create_mover() const
// XRW TEMP {
// XRW TEMP  return protocols::moves::MoverOP( new LegacyGivenPathAssemblyMover );
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP LegacyGivenPathAssemblyMoverCreator::keyname() const
// XRW TEMP {
// XRW TEMP  return LegacyGivenPathAssemblyMover::mover_name();
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP LegacyGivenPathAssemblyMover::mover_name()
// XRW TEMP {
// XRW TEMP  return "LegacyGivenPathAssemblyMover";
// XRW TEMP }

protocols::moves::MoverOP
LegacyGivenPathAssemblyMover::clone() const {
	return( protocols::moves::MoverOP( new LegacyGivenPathAssemblyMover( *this ) ) );
}
protocols::moves::MoverOP
LegacyGivenPathAssemblyMover::fresh_instance() const {
	return protocols::moves::MoverOP( new LegacyGivenPathAssemblyMover );
}

// XRW TEMP std::string
// XRW TEMP LegacyGivenPathAssemblyMover::get_name() const {
// XRW TEMP  return "LegacyGivenPathAssemblyMover";
// XRW TEMP }

/////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////  LegacyGivenPathAssemblyMover function   //////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

LegacyGivenPathAssemblyMover::LegacyGivenPathAssemblyMover(){}

AssemblyOP
LegacyGivenPathAssemblyMover::generate_assembly(){
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	//Add first assembly
	AssemblyOP assembly;
	if ( option[ basic::options::OptionKeys::legacy_sewing::assembly_type ].value() == "discontinuous" ) {
		assembly = AssemblyFactory::create_assembly("discontinuous");
	} else if ( option[ basic::options::OptionKeys::legacy_sewing::assembly_type ].value() == "continuous" ) {
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
	//   utility::graph::EdgeListConstIterator edge_it = reference_node->const_edge_list_begin();
	//   utility::graph::EdgeListConstIterator edge_it_end = reference_node->const_edge_list_end();
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
LegacyGivenPathAssemblyMover::parse_my_tag(
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
		utility_exit_with_message("You must specify a path in the LegacyGivenPathAssemblyMover tag");
	}
}

std::string LegacyGivenPathAssemblyMover::get_name() const {
	return mover_name();
}

std::string LegacyGivenPathAssemblyMover::mover_name() {
	return "LegacyGivenPathAssemblyMover";
}

void LegacyGivenPathAssemblyMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	//Extra attributes
	AttributeList extra_ats;
	extra_ats
		+ XMLSchemaAttribute::required_attribute( "path", xsct_int_cslist, "Comma-separated list of model IDs (in order) to incorporate into the assembly" );

	XMLSchemaComplexTypeGeneratorOP ct_gen = LegacyAssemblyMover::define_assembly_mover_ct_gen( xsd );
	ct_gen->element_name( mover_name() );
	ct_gen->description( "Creates an assembly from the given model IDs in the given order." );
	ct_gen->add_attributes( extra_ats );
	ct_gen->write_complex_type_to_schema( xsd );
}

std::string LegacyGivenPathAssemblyMoverCreator::keyname() const {
	return LegacyGivenPathAssemblyMover::mover_name();
}

protocols::moves::MoverOP
LegacyGivenPathAssemblyMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new LegacyGivenPathAssemblyMover );
}

void LegacyGivenPathAssemblyMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	LegacyGivenPathAssemblyMover::provide_xml_schema( xsd );
}


} //legacy_sewing
} //protocols
