// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.


/// @file protocols/seeded_abinitio/SwapSegment.cc
/// @brief
/// @author Eva-Maria Strauch (evas01@u.washington.edu)

#include <protocols/seeded_abinitio/SwapSegment.hh>
#include <protocols/seeded_abinitio/SwapSegmentCreator.hh>

#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/pose/variant_util.hh>
#include <core/import_pose/import_pose.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

//other
#include <protocols/rosetta_scripts/util.hh>
#include <core/pose/selection.hh>

// C++ headers
#include <string>

#include <basic/Tracer.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>

//parser
#include <utility/tag/Tag.hh>
#include <utility/tag/Tag.fwd.hh>
#include <basic/datacache/DataMap.hh>

//loops
#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.fwd.hh>

//util
#include <utility/vector1.hh>
#include <utility/vector0.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>


using namespace core;
using namespace protocols::seeded_abinitio;
static basic::Tracer TR( "protocols.seeded_abinitio.SwapSegment" );


namespace protocols {
namespace seeded_abinitio {

using namespace protocols::moves;
using namespace core;

// XRW TEMP std::string
// XRW TEMP SwapSegmentCreator::keyname() const
// XRW TEMP {
// XRW TEMP  return SwapSegment::mover_name();
// XRW TEMP }

// XRW TEMP protocols::moves::MoverOP
// XRW TEMP SwapSegmentCreator::create_mover() const {
// XRW TEMP  return protocols::moves::MoverOP( new SwapSegment() );
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP SwapSegment::mover_name(){
// XRW TEMP  return "SwapSegment";
// XRW TEMP }


SwapSegment::~SwapSegment() = default;

SwapSegment::SwapSegment():
	protocols::moves::Mover( SwapSegment::mover_name() ){
	previously_grown_ = false;
}

protocols::moves::MoverOP
SwapSegment::clone() const {
	return( protocols::moves::MoverOP( new SwapSegment( *this ) ) );
}

protocols::moves::MoverOP
SwapSegment::fresh_instance() const {
	return protocols::moves::MoverOP( new SwapSegment );
}

void
SwapSegment::copying_side_chains(
	core::pose::Pose & pose,
	core::pose::PoseOP & swap_segment,
	protocols::loops::Loops & seeds
){

	//counter for seeds
	Size offsetres = 0;

	for ( Size pos = 1; pos <= pose.size(); ++pos ) {
		TR.Debug<<"iterating through pose: " <<pos << std::endl;
		if ( seeds.is_loop_residue( pos ) ) {
			TR.Debug<<"pos "<< pos <<std::endl;
			++offsetres;
			TR.Debug<<"offsetres " << offsetres <<std::endl;
			pose.replace_residue(pos, swap_segment->residue( offsetres ), true);
		}
	}

}

//for this method, use only the actual fragments that are supposed to be swapped under the swap_segment
void
SwapSegment::swap_segment(
	core::pose::Pose & pose,
	core::pose::PoseOP & swap_segment,
	protocols::loops::Loops & seeds
){

	using namespace core::conformation;
	Size offsetres = 0;

	// preparing the segments for swap
	core::pose::remove_lower_terminus_type_from_pose_residue(*swap_segment, 1 );
	core::pose::remove_upper_terminus_type_from_pose_residue(*swap_segment, swap_segment->size());

	for ( Size pos = 1; pos <= pose.size(); ++pos ) {
		TR.Debug<<"iterating through pose: " <<pos << std::endl;
		if ( seeds.is_loop_residue( pos ) ) {
			TR.Debug <<"pos "<< pos <<std::endl;
			++offsetres ;
			TR.Debug <<"offsetres " << offsetres <<std::endl;

			pose.replace_residue( pos, swap_segment->residue( offsetres ), true);
			TR.Debug << "After Swap loop Residue  " << offsetres << std::endl;

		}
	}

}
void
SwapSegment::swap_chain(
	core::pose::Pose & pose,
	core::pose::PoseOP & target_chain,
	core::Size chain_to_swap
){

	TR<<"replacing residues from chain " <<chain_to_swap << " with appropriate rotamers" <<std::endl;
	for ( core::Size pos = pose.conformation().chain_begin( chain_to_swap ); pos <= pose.conformation().chain_end( chain_to_swap ); ++pos ) {
		pose.replace_residue( pos, target_chain->residue( pos ), true);
	}
}


void
SwapSegment::apply( core::pose::Pose & pose )
{

	core::pose::PoseOP segment = seeds_pdb_->split_by_chain( from_chain_ );
	core::pose::PoseOP target_c = seeds_pdb_->split_by_chain( swap_chain_ );
	//assert sizes of target chain and input pose chain

	TR.Debug<<"chains: "<< pose.conformation().num_chains() <<" total residues: " <<pose.size() <<std::endl;
	//need to put in some assertions for stupid things like number bigger than there are chains and
	//something that automatically adjusts the residue number and hook in the to_chain option
	//for now hardcoding for first chain
	//if( pose.conformation().num_chains() > 1

	//core::Size adjust_numbering = 0;

	//if( to_chain > 1)
	core::Size adjust_numbering = pose.conformation().chain_begin( to_chain_ ) - 1;

	//if( pose.conformation().num_chains() > 1 )
	// adjust_numbering = pose.conformation().chain_end( 1 );

	//need to assert that the seed elements are the same residue numbers as in the input pdb!!
	TR.Debug<<"all_seeds_.loop_size() " << all_seeds_.loop_size() <<" total residues in seeds pdb: " << segment->size() <<std::endl;
	TR<<"adjusting for chain numbering by: " <<adjust_numbering <<std::endl;

	all_seeds_.make_sequence_shift( adjust_numbering );
	TR.Debug <<"new loops: " <<all_seeds_ <<std::endl;

	if ( swap_chain_ > 0 ) {
		TR<<"swapping chain: " << swap_chain_ << std::endl;
		swap_chain( pose, target_c , swap_chain_ );
	}

	if ( all_seeds_.loop_size() != segment->size() ) {
		utility_exit_with_message("residues specified under the seeds does not agree with the number of residues provided as segment");
	}

	if ( copy_sidechains_ ) {
		copying_side_chains( pose, segment, all_seeds_);
		(*scorefxn_)(pose);
		TR<<"adopting side chains" << std::endl;
	}

	if ( swap_segment_ ) {
		swap_segment( pose, segment, all_seeds_);
		TR<<"swapping segment" <<std::endl;
		(*scorefxn_)(pose);
	}

	TR.flush();
}


// XRW TEMP std::string
// XRW TEMP SwapSegment::get_name() const {
// XRW TEMP  return SwapSegment::mover_name();
// XRW TEMP }

void
SwapSegment::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data ,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	core::pose::Pose const & pose ){

	all_seeds_.clear();//just in case

	//need scorefxn to score after swapping
	scorefxn_ = protocols::rosetta_scripts::parse_score_function( tag, data )->clone();
	TR<<"scoring with following scorefunction: " << *scorefxn_ <<std::endl;

	TR<<"SwapSegment mover has been instantiated" <<std::endl;

	copy_sidechains_ = tag->getOption< bool > ("copy_sidechains", true );
	swap_segment_ = tag->getOption< bool > ("swap_segment", false );
	swap_chain_ = tag->getOption< core::Size > ("swap_chain", 0 );
	from_chain_ = tag->getOption< core::Size > ("from_chain", 1 );
	to_chain_ = tag->getOption< core::Size > ("to_chain", 1 );

	TR<<"swap sidechains: " << copy_sidechains_ << ", swapping segment: "<< swap_segment_ <<", taking segments from chain: " << from_chain_ << std::endl;

	if ( copy_sidechains_ && swap_segment_ ) {
		TR<<"not necessary to swap segments AND to copy side chains" << std::endl;
	}

	if ( tag->hasOption( "seeds_pdb" ) || tag->hasOption( "template_pdb" ) ) {
		std::string const template_pdb_fname( tag->getOption< std::string >( "seeds_pdb" ));
		seeds_pdb_ = core::pose::PoseOP( new core::pose::Pose ) ;
		core::import_pose::pose_from_file( *seeds_pdb_, template_pdb_fname , core::import_pose::PDB_file);
		TR<<"read in a template pdb with " <<seeds_pdb_->size() <<"residues"<<std::endl;
		seeds_presence_ = true;
	}

	if ( !seeds_presence_ ) {
		utility_exit_with_message("need to specify a template to swap from!!!!");
	}
	//use the input pdb too ( todo )

	if ( tag->hasOption( "previously_grown") ) {
		previously_grown_ = tag->getOption< bool > ("previously_grown", false );
		TR<<"decoy was previously changed in its length, aka was \"grown\" "<<std::endl;
	}

	//parsing branch tags
	utility::vector0< TagCOP > const & branch_tags( tag->getTags() );
	for ( TagCOP btag : branch_tags ) {

		//there is a problem with the parsing of the seeds if the pose was previously grown. This is due to
		//the difference between parse time and computing time. Since the seeds are parsed before the pose has been
		//grown to its full length, a different numbering needs to be used

		if ( btag->getName() == "Seeds" ) { //need an assertion for the presence of these or at least for the option file

			core::Size begin;
			core::Size end;

			if ( !previously_grown_ ) {
				std::string const beginS( btag->getOption<std::string>( "begin" ) );
				std::string const endS( btag->getOption<std::string>( "end" ) );
				begin =( core::pose::parse_resnum( beginS, pose ) );
				end   =( core::pose::parse_resnum( endS, pose ) );
				TR.Debug<<"string seeds: \n"<< beginS <<" and " << endS <<std::endl;
			} else {
				begin = ( btag->getOption< core::Size >( "begin" ) );
				end = ( btag->getOption< core::Size >( "end" ) );
			}

			all_seeds_.add_loop( begin , end , 0, 0, false );

			TR.Debug<<"parsing seeds: \n"<< begin <<" and " << end <<std::endl;
			TR.Debug<<"seeds: "<< all_seeds_ <<std::endl;
		}//end seed tags
	}//end branch tags

	if ( all_seeds_.size() == 0 ) {
		utility_exit_with_message("NEED TO SPECIFY A SEGMENT TO REPLACE!!!, whole pose swap currently not supported");
	}

}

std::string SwapSegment::get_name() const {
	return mover_name();
}

std::string SwapSegment::mover_name() {
	return "SwapSegment";
}

void SwapSegment::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;

	protocols::rosetta_scripts::attributes_for_parse_score_function( attlist );
	attlist
		+ XMLSchemaAttribute::attribute_w_default("copy_sidechains", xsct_rosetta_bool, "XRW TO DO", "1")
		+ XMLSchemaAttribute::attribute_w_default("swap_segment", xsct_rosetta_bool, "XRW TO DO", "0")
		+ XMLSchemaAttribute::attribute_w_default("swap_chain", xsct_non_negative_integer, "XRW TO DO", "0")
		+ XMLSchemaAttribute::attribute_w_default("from_chain", xsct_non_negative_integer, "XRW TO DO", "1")
		+ XMLSchemaAttribute::attribute_w_default("to_chain", xsct_non_negative_integer, "XRW TO DO", "1")
		+ XMLSchemaAttribute("seeds_pdb", xs_string, "XRW TO DO, required if template_pdb is not specified.")
		+ XMLSchemaAttribute("template_pdb", xs_string, "XRW TO DO, required if seeds_pdb is not specified.")
		+ XMLSchemaAttribute::attribute_w_default("previously_grown", xsct_rosetta_bool, "XRW TO DO", "0");

	AttributeList subelement_attributes;
	subelement_attributes
		+ XMLSchemaAttribute::required_attribute("begin", xs_string, "XRW TO DO")
		+ XMLSchemaAttribute::required_attribute("end", xs_string, "XRW TO DO");
	XMLSchemaSimpleSubelementList subelement_list;
	subelement_list.add_simple_subelement("Seeds", subelement_attributes, "XRW TO DO");

	protocols::moves::xsd_type_definition_w_attributes_and_repeatable_subelements( xsd, mover_name(), "XRW TO DO", attlist, subelement_list );
}

std::string SwapSegmentCreator::keyname() const {
	return SwapSegment::mover_name();
}

protocols::moves::MoverOP
SwapSegmentCreator::create_mover() const {
	return protocols::moves::MoverOP( new SwapSegment );
}

void SwapSegmentCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	SwapSegment::provide_xml_schema( xsd );
}

}
}
