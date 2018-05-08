// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file StartingNodeMotifScorer.cc
///
/// @brief
/// @author Tim Jacobs

//Unit headers
#include <protocols/sewing/scoring/StartingNodeMotifScorer.hh>
#include <protocols/sewing/scoring/StartingNodeMotifScorerCreator.hh>

//Package headers
#include <protocols/sewing/data_storage/SmartAssembly.hh>
#include <protocols/sewing/data_storage/SmartSegment.hh>
#include <protocols/sewing/data_storage/SmartSewingResidue.hh>
#include <protocols/sewing/scoring/AssemblyScorerFactory.hh>
//Core headers
#include <core/types.hh>

#include <core/chemical/ChemicalManager.hh>

#include <core/pose/motif/reference_frames.hh>

#include <core/scoring/dssp/Dssp.hh>

#include <core/scoring/motif/util.hh>

//Utility headers
#include <basic/Tracer.hh>
#include <numeric/xyzTransform.hh>
#include <utility/io/ozstream.hh>
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
namespace protocols {
namespace sewing  {
namespace scoring {

static basic::Tracer TR("protocols.sewing.scoring.StartingNodeMotifScorer");

StartingNodeMotifScorer::StartingNodeMotifScorer():
	MotifScorer()
{}

StartingNodeMotifScorer::StartingNodeMotifScorer( StartingNodeMotifScorer const & src ):
	MotifScorer( src ),
	weight_( src.weight_ )
{
}

///@details Use the negative normalized motif score
core::Real
StartingNodeMotifScorer::score(
	data_storage::SmartAssemblyCOP assembly
) {
	this->set_old_last_score( this->get_last_score() );
	this->set_last_score( -1.0 * full_motif_score( assembly ) );
	return this->get_last_score();
}


//Loop through segments
//For each segment, check if it is vital or if it is a chimera and its parent is vital
//For each of those, loop through and check them against all segments that aren't vital
//No need to check for double counting in this case





///@details use Will's Motif score to calculate the motif score for interactions between
///a given segment and segments from other models. Divide by total number of segments
core::Real
StartingNodeMotifScorer::full_motif_score(
	data_storage::SmartAssemblyCOP assembly
) {
	TR << "Scoring StartingNodeMotifScorer"<< std::endl;
	core::Real score = 0.0;
	core::Size counter = 0; //This keeps track of how many interactions we have scored
	core::Size vital_seg_res = 0; //This keeps track of the number of residues in vital segments (for normalization)
	data_storage::SmartSegmentOP current_segment = assembly->get_n_terminal_segment();
	//Note that this segment is pulled from the assembly's set of segments, not the master set
	while ( current_segment != nullptr ) {
		TR << "Current segment ID " << current_segment->get_segment_id() << std::endl;
		TR << "Current segment dssp code " << current_segment->get_dssp_code() << std::endl;
		if ( !current_segment->is_vital() ) {
			current_segment = current_segment->get_c_terminal_neighbor();
			TR.Debug <<"Current segment is not vital" <<std::endl;
			continue; //Only proceed for vital segments
		}

		TR.Debug <<"Current segment is vital" <<std::endl;
		char ss1 = current_segment->get_dssp_code();
		//Iterate through the residues in this segment
		for ( data_storage::SmartSewingResidueOP current_res : current_segment->get_residue_vector() ) {
			//TR <<"Scoring a residue" << std::endl;
			//New version of get_stub just takes the SmartSewingResidueOP
			++vital_seg_res;
			numeric::xyzTransform< core::Real > stub1 = get_stub( current_res );
			//We already have the DSSP, just get the aa type
			char aa1 = res_type_set_->name_map( current_res->get_amino_acid_type() ).name1();

			//For our second loop, we'll score all the non-vital segments
			data_storage::SmartSegmentOP second_segment = assembly->get_n_terminal_segment();
			while ( second_segment != nullptr ) {
				if ( second_segment->is_vital() ) {
					second_segment = second_segment->get_c_terminal_neighbor();
					continue;
				}
				//Now we can find the motif score between them like normal
				char ss2 = second_segment->get_dssp_code();
				for ( data_storage::SmartSewingResidueOP second_res : second_segment->get_residue_vector() ) {
					//TR << "Comparing two residues" << std::endl;
					numeric::xyzTransform< core::Real > stub2 = get_stub( second_res );
					char aa2 = res_type_set_->name_map( second_res->get_amino_acid_type() ).name1();
					score += get_score(stub1, ss1, aa1, stub2, ss2, aa2);
					++counter;
				}//for second_res
				second_segment = second_segment->get_c_terminal_neighbor();
			}//while second_segment
		} //for current_res
		current_segment = current_segment->get_c_terminal_neighbor();
	} //while current_segment
	if ( counter == 0 ) { return score; }
	//return score / counter;
	//To bring the scale of this in line with other scorers, we should instead normalize by the length of the vital section
	return score / vital_seg_res;
}


void
StartingNodeMotifScorer::set_options_from_tag(
	utility::tag::TagCOP scorer,
	basic::datacache::DataMap& ,
	protocols::filters::Filters_map const & ,
	protocols::moves::Movers_map const & ,
	core::pose::Pose const & )
{
	//Data members? It can figure out starting node from is_vital
	if ( !scorer->hasOption( "weight" ) ) {
		TR << "No weight provided for scorer " << scorer->getName() << " Using default value" << std::endl;
	}

	weight_ = scorer->getOption<core::Real>( "weight", 1.0 );
	TR << "Created StartingNodeMotifScorer from tag." << std::endl;

}
std::string
StartingNodeMotifScorer::get_name() const {
	return "StartingNodeMotifScorer";
}

core::Real
StartingNodeMotifScorer::get_weight() const{
	return weight_;
}

void
StartingNodeMotifScorer::set_weight( core::Real weight ){
	weight_ = weight;
}

void
StartingNodeMotifScorer::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ){
	using namespace utility::tag;
	AttributeList attributes;
	attributes
		+ XMLSchemaAttribute::attribute_w_default( "weight", xsct_real, "How heavily will this term be weighted during scoring?", "1.0"  );

	XMLSchemaComplexTypeGenerator ct_gen;
	ct_gen.complex_type_naming_func( & AssemblyScorerFactory::assembly_scorer_ct_namer )
		.element_name( StartingNodeMotifScorer::type_name() )
		.description( "Specifically scores packing against the starting node")
		.add_attributes( attributes )
		.add_optional_name_attribute()
		.write_complex_type_to_schema( xsd ); //We won't have regular score functions/task operations
}


//////Creator methods//////////

void
StartingNodeMotifScorerCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const{
	StartingNodeMotifScorer::provide_xml_schema( xsd );

}


AssemblyScorerOP
StartingNodeMotifScorerCreator::create_assembly_scorer() const{
	return AssemblyScorerOP( new StartingNodeMotifScorer() );
}

std::string
StartingNodeMotifScorer::type_name(){
	return "StartingNodeMotifScorer";
}


std::string
StartingNodeMotifScorerCreator::keyname() const{
	return StartingNodeMotifScorer::type_name();
}

/////End Creator methods///////






} //scoring namespace
} //sewing namespace
} //protocols namespace
