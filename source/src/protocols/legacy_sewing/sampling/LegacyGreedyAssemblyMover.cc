// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file LegacyGreedyAssemblyMover.cc
///
/// @brief
/// @author Tim Jacobs

// Unit Headers
#include <protocols/legacy_sewing/sampling/LegacyGreedyAssemblyMover.hh>
#include <protocols/legacy_sewing/sampling/LegacyGreedyAssemblyMoverCreator.hh>

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

static basic::Tracer TR( "protocols.legacy_sewing.sampling.LegacyGreedyAssemblyMover" );

////////////////////////////////////////////////////////////////////////////////////
///////////////////////////  Boiler Plate Code   ///////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
// XRW TEMP protocols::moves::MoverOP
// XRW TEMP LegacyGreedyAssemblyMoverCreator::create_mover() const
// XRW TEMP {
// XRW TEMP  return protocols::moves::MoverOP(new LegacyGreedyAssemblyMover);
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP LegacyGreedyAssemblyMoverCreator::keyname() const
// XRW TEMP {
// XRW TEMP  return LegacyGreedyAssemblyMover::mover_name();
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP LegacyGreedyAssemblyMover::mover_name()
// XRW TEMP {
// XRW TEMP  return "LegacyGreedyAssemblyMover";
// XRW TEMP }

protocols::moves::MoverOP
LegacyGreedyAssemblyMover::clone() const {
	return( protocols::moves::MoverOP( new LegacyGreedyAssemblyMover( *this ) ) );
}
protocols::moves::MoverOP
LegacyGreedyAssemblyMover::fresh_instance() const {
	return protocols::moves::MoverOP( new LegacyGreedyAssemblyMover );
}

// XRW TEMP std::string
// XRW TEMP LegacyGreedyAssemblyMover::get_name() const {
// XRW TEMP  return "LegacyGreedyAssemblyMover";
// XRW TEMP }

/////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////  LegacyGreedyAssemblyMover function   //////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

LegacyGreedyAssemblyMover::LegacyGreedyAssemblyMover():
	best_complete_assembly_(0),
	best_score_(10000),
	cycles_(1000),
	max_edges_per_node_(300)
{}

AssemblyOP
LegacyGreedyAssemblyMover::generate_assembly(){
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	//Main loop
	core::Size starttime = time(NULL);
	core::Size cur_cycle=1;

	for ( cur_cycle=1; cur_cycle <= cycles_; ++cur_cycle ) {

		//Initialize starting Assembly for the current cycle
		AssemblyOP assembly;
		if ( option[ basic::options::OptionKeys::legacy_sewing::assembly_type ].value() == "discontinuous" ) {
			assembly = AssemblyFactory::create_assembly("discontinuous");
		} else if ( option[ basic::options::OptionKeys::legacy_sewing::assembly_type ].value() == "continuous" ) {
			assembly = AssemblyFactory::create_assembly("continuous");
		}
		add_starting_model(assembly);

		TR << "Cycle " << cur_cycle << std::endl;
		while ( requirement_set_->can_be_added_to(assembly) ) {
			TR << "Looking to add edge " << assembly->path().size()+1 << std::endl;
			core::Real best_edge_score_ = 10000;
			AssemblyOP best_edge_assembly_ = 0;

			ModelNode const * reference_node = graph_->get_model_node(assembly->get_next_reference_node(graph_));
			graph_->add_edges_from_binary(edge_file_, reference_node->get_node_index());
			core::Size num_edges = reference_node->num_edges();

			//If there are no edges from the selected node, abort this trajectory and go to the next cycle
			if ( num_edges == 0 ) {
				TR << "No edges for model " << reference_node->model().model_id_ << std::endl;
				break;
			}

			core::Size max_edge_attempts = std::min(num_edges, max_edges_per_node_);

			//Randomize the edge order
			utility::vector1<core::Size> edge_order(num_edges);
			for ( core::Size i = 0; i < num_edges; ++i ) {
				edge_order[i+1]=i;
			}
			numeric::random::random_permutation(edge_order, numeric::random::rg());

			for ( core::Size cur_edge_ind=1; cur_edge_ind<=max_edge_attempts; ++cur_edge_ind ) {

				AssemblyOP pre_edge_assembly = assembly->clone();

				utility::graph::EdgeListConstIterator edge_it = reference_node->const_edge_list_begin();
				for ( core::Size j=0; j<edge_order[cur_edge_ind]; ++j ) {
					++edge_it;
				}

				//Cast the edge to a proper HashEdge, check if the new model satisfies requirements, and follow it
				HashEdge const * const cur_edge = static_cast< HashEdge const * >(*edge_it);
				assembly->follow_edge(graph_, cur_edge, reference_node->get_node_index());
				if ( !requirement_set_->violates(assembly) ) {
					core::Real edge_score = assembly_scorefxn_->score(assembly);
					if ( edge_score < best_edge_score_ ) {
						best_edge_score_ = edge_score;
						best_edge_assembly_ = assembly->clone();
					}
				}
				assembly = pre_edge_assembly;
			}

			//If we couldn't find an edge that satisfies node requirements, go to the next cycle
			if ( best_edge_assembly_ == 0 ) { break; }

			//Revert to the best scoring Assembly for the most recent edge addition, check to see if
			//this assembly is complete. If so, check to see if it's the best one and continue on
			assembly = best_edge_assembly_;
			if ( requirement_set_->satisfies(assembly) ) {
				core::Real complete_score = assembly_scorefxn_->score(assembly);
				if ( complete_score < best_score_ ) {
					best_score_ = complete_score;
					best_complete_assembly_ = assembly->clone();
					TR << "SAVING BEST " << best_score_ << std::endl;
				}
			}
		}// While can be added to
	}//Cycles

	core::Size endtime = time(NULL);
	TR << "Completed " << cur_cycle << " cycles in " << endtime - starttime << " seconds" << std::endl;
	return best_complete_assembly_;
}

void
LegacyGreedyAssemblyMover::parse_my_tag(
	TagCOP const tag,
	basic::datacache::DataMap & data,
	protocols::filters::Filters_map const & filters,
	protocols::moves::Movers_map const & movers,
	core::pose::Pose const & pose
){

	parent::parse_my_tag(tag, data, filters, movers, pose);

	if ( tag->hasOption("cycles") ) {
		cycles_ = tag->getOption<core::Size>("cycles");
	}

	if ( tag->hasOption("max_edges_per_node") ) {
		max_edges_per_node_ = tag->getOption<core::Size>("max_edges_per_node");
	}
}

std::string LegacyGreedyAssemblyMover::get_name() const {
	return mover_name();
}

std::string LegacyGreedyAssemblyMover::mover_name() {
	return "LegacyGreedyAssemblyMover";
}

void LegacyGreedyAssemblyMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	//Extra attributes
	AttributeList extra_ats;
	extra_ats
		+ XMLSchemaAttribute( "cycles", xsct_non_negative_integer, "Number of rounds during which we will add to the assembly" )
		+ XMLSchemaAttribute( "max_edges_per_node", xsct_non_negative_integer, "Number of edges to try per cycle" );

	XMLSchemaComplexTypeGeneratorOP ct_gen = LegacyAssemblyMover::define_assembly_mover_ct_gen( xsd );
	ct_gen->element_name( mover_name() );
	ct_gen->description( "At each cycle, it tries a given number of edges and chooses the best scoring assembly for that round." );
	ct_gen->add_attributes( extra_ats );
	ct_gen->write_complex_type_to_schema( xsd );

}

std::string LegacyGreedyAssemblyMoverCreator::keyname() const {
	return LegacyGreedyAssemblyMover::mover_name();
}

protocols::moves::MoverOP
LegacyGreedyAssemblyMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new LegacyGreedyAssemblyMover );
}

void LegacyGreedyAssemblyMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	LegacyGreedyAssemblyMover::provide_xml_schema( xsd );
}


} //legacy_sewing
} //protocols
