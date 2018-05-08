// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file SubsetPartnerMotifScorer.cc
///
/// @brief
/// @author Frank Teets

//Unit headers
#include <protocols/sewing/scoring/SubsetPartnerMotifScorer.hh>
#include <protocols/sewing/scoring/SubsetPartnerMotifScorerCreator.hh>
#include <protocols/sewing/scoring/AssemblyScorerFactory.hh>
//Core headers
#include <core/types.hh>

#include <core/chemical/ChemicalManager.hh>

#include <core/pose/motif/reference_frames.hh>

#include <core/scoring/motif/util.hh>

#include <core/scoring/dssp/Dssp.hh>

#include <core/pose/Pose.hh>

//Utility headers
#include <basic/Tracer.hh>
#include <numeric/xyzTransform.hh>
#include <utility/io/ozstream.hh>
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
namespace protocols {
namespace sewing  {
namespace scoring {

static basic::Tracer TR("protocols.sewing.scoring.SubsetPartnerMotifScorer");

SubsetPartnerMotifScorer::SubsetPartnerMotifScorer():
	MotifScorer()
{}

SubsetPartnerMotifScorer::SubsetPartnerMotifScorer( SubsetPartnerMotifScorer const & src):
	MotifScorer( src ),
	weight_( src.weight_ ),
	starts_and_ends_( src.starts_and_ends_ ),
	region_start_( src.region_start_ ),
	region_end_( src.region_end_ )
{
}

///@details Use the negative normalized motif score
core::Real
SubsetPartnerMotifScorer::score(
	data_storage::SmartAssemblyCOP assembly
) {
	this->set_old_last_score( this->get_last_score() );
	this->set_last_score( -1.0 * interface_motif_score(assembly) );
	return this->get_last_score();
}

///@details use Will's Motif score to calculate the motif score for the entire Assembly
core::Real
SubsetPartnerMotifScorer::interface_motif_score(
	data_storage::SmartAssemblyCOP assembly
) {
	core::Real score = 0.0;

	// utility::vector1< data_storage::SmartSegmentOP > const segments = assembly->get_segment_vector();
	data_storage::SmartSegmentOP current_segment = assembly->get_n_terminal_segment();
	core::pose::PoseOP partner = assembly->get_partner();
	if ( !partner ) {
		return score;
	}
	if ( partner->secstruct() == "" ) {
		core::scoring::dssp::Dssp dssp(*partner);
		dssp.insert_ss_into_pose(*partner);
	}


	//Iterate over segments
	//for ( core::Size i=1; i<=segments.size(); ++i ) {
	while ( current_segment != nullptr ) {
		//secondary structure is stored in the segment
		char ss1 = current_segment->get_dssp_code();
		for ( core::Size res_i = 1; res_i <= current_segment->get_length(); ++res_i ) {
			//get_stub now just takes a SmartSewingResidueOP
			numeric::xyzTransform<core::Real> stub1 = get_stub( current_segment->get_residue( res_i ) );
			char aa1 = res_type_set_->name_map( current_segment->get_residue( res_i )->get_amino_acid_type() ).name1();
			//Iterate over partner residues
			for ( core::Size j=region_start_; j<=region_end_; ++j ) {

				numeric::xyzTransform<core::Real> stub2 = core::pose::motif::get_backbone_reference_frame(*partner, j);
				char ss2 = partner->secstruct(j);
				char aa2 = partner->residue(j).name1();
				score += get_score(stub1, ss1, aa1, stub2, ss2, aa2);
			} //end for partner residue
		} //end for res_i
		current_segment = current_segment->get_n_terminal_neighbor();
	} //end while
	return score / ( assembly->get_length() );
}

void
SubsetPartnerMotifScorer::set_options_from_tag(
	utility::tag::TagCOP scorer,
	basic::datacache::DataMap& ,
	protocols::filters::Filters_map const & ,
	protocols::moves::Movers_map const & ,
	core::pose::Pose const &)
{
	//No data members to set!

	if ( !scorer->hasOption( "weight" ) ) {
		TR << "No weight provided for scorer " << scorer->getName() << " Using default value" << std::endl;
	}

	weight_ = scorer->getOption<core::Real>( "weight", 1.0 );

	region_start_ = scorer->getOption<core::Size>("region_start", 1);
	region_end_ = scorer->getOption<core::Size>("region_end", 1);

	TR << "Created SubsetPartnerMotifScorer from tag!" << std::endl;
}


void
SubsetPartnerMotifScorer::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ){
	using namespace utility::tag;
	AttributeList attributes;
	attributes
		+ XMLSchemaAttribute::attribute_w_default( "weight", xsct_real, "How heavily will this term be weighted during scoring?", "1.0"  )
		+ XMLSchemaAttribute::attribute_w_default( "region_start", xsct_non_negative_integer, "What is the first residue of the scored subset?","1"  )
		+ XMLSchemaAttribute::attribute_w_default( "region_end", xsct_non_negative_integer, "What is the last residue of the scored subset?","2"  );
	XMLSchemaComplexTypeGenerator ct_gen;
	ct_gen.complex_type_naming_func( & AssemblyScorerFactory::assembly_scorer_ct_namer )
		.element_name( SubsetPartnerMotifScorer::type_name() )
		.add_optional_name_attribute()
		.description( "Motif score to measure packing of assembly against partner PDB" )
		.add_attributes( attributes )
		.write_complex_type_to_schema( xsd ); //We won't have regular score functions/task operations

}

core::Real
SubsetPartnerMotifScorer::get_weight() const {
	return weight_;
}

void
SubsetPartnerMotifScorer::set_weight( core::Real weight ){
	weight_ = weight;
}

std::string
SubsetPartnerMotifScorer::get_starts_and_ends() const{
	return starts_and_ends_;
}

void
SubsetPartnerMotifScorer::set_starts_and_ends(std::string starts_and_ends){
	TR << "This variable doesn't do anything!" << std::endl;
	starts_and_ends_ = starts_and_ends;
}


core::Size
SubsetPartnerMotifScorer::get_region_start() const{
	return region_start_;
}

core::Size
SubsetPartnerMotifScorer::get_region_end() const{
	return region_end_;
}

void
SubsetPartnerMotifScorer::set_region_start( core::Size setting ){
	region_start_ = setting;
}

void
SubsetPartnerMotifScorer::set_region_end( core::Size setting ){
	region_end_ = setting;
}



//////Creator methods//////////
void
SubsetPartnerMotifScorerCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const{
	SubsetPartnerMotifScorer::provide_xml_schema( xsd );

}



AssemblyScorerOP
SubsetPartnerMotifScorerCreator::create_assembly_scorer() const{
	return AssemblyScorerOP( new SubsetPartnerMotifScorer() );
}

std::string
SubsetPartnerMotifScorer::type_name(){
	return "SubsetPartnerMotifScorer";
}


std::string
SubsetPartnerMotifScorerCreator::keyname() const{
	return SubsetPartnerMotifScorer::type_name();
}

/////End Creator methods///////











} //scoring namespace
} //sewing namespace
} //protocols namespace
