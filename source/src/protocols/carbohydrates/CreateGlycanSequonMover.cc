// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/carbohydrates/CreateGlycanSequonMover.cc
/// @brief Mutates residues to create a potential glycosylation site using known sequence motifs of N- or C- linked glycans.  Includes options for Enhanced Sequons for N-linked glycans that have been shown to have higher rates of glycosylation as well as other positions that have been shown to influence the glycosylation chemistry.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

// Unit headers
#include <protocols/carbohydrates/CreateGlycanSequonMover.hh>
#include <protocols/carbohydrates/CreateGlycanSequonMoverCreator.hh>

#include <protocols/simple_moves/CreateSequenceMotifMover.hh>

// Core headers
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

// Basic/Utility headers
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>
#include <utility/string_util.hh>
#include <core/types.hh>

#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/select/residue_selector/ReturnResidueSubsetSelector.hh>
#include <core/select/residue_selector/util.hh>


// XSD Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>
#include <protocols/rosetta_scripts/util.hh>

static basic::Tracer TR( "protocols.carbohydrates.CreateGlycanSequonMover" );

namespace protocols {
namespace carbohydrates {
using namespace core::select::residue_selector;
using namespace protocols::simple_moves;
using namespace protocols::rosetta_scripts;

/////////////////////
/// Constructors  ///
/////////////////////

/// @brief Default constructor
CreateGlycanSequonMover::CreateGlycanSequonMover():
	protocols::moves::Mover( CreateGlycanSequonMover::mover_name() )
{
	sequons_ = create_sequons();
	motif_mover_ = CreateSequenceMotifMoverOP( new CreateSequenceMotifMover() );
}

CreateGlycanSequonMover::CreateGlycanSequonMover( ResidueSelectorCOP selector ):
	protocols::moves::Mover( CreateGlycanSequonMover::mover_name() )
{
	set_residue_selector(selector);
	sequons_ = create_sequons();
	motif_mover_ = CreateSequenceMotifMoverOP( new CreateSequenceMotifMover() );
}


////////////////////////////////////////////////////////////////////////////////
/// @brief Copy constructor
CreateGlycanSequonMover::CreateGlycanSequonMover( CreateGlycanSequonMover const & src ):
	protocols::moves::Mover( src ),
	sequon_type_( src.sequon_type_ ),
	positions_as_start_of_sequon_( src.positions_as_start_of_sequon_ ),
	design_x_positions_( src.design_x_positions_),
	pack_neighbors_( src.pack_neighbors_ ),
	design_neighbors_( src.design_neighbors_ ),
	pack_distance_( src.pack_distance_ ),
	pack_rounds_( src.pack_rounds_)

{
	if ( src.selector_ ) {
		selector_ = src.selector_->clone();
	}
	if ( src.scorefxn_ ) {
		scorefxn_ = src.scorefxn_->clone();
	}
	motif_mover_ = CreateSequenceMotifMoverOP( new CreateSequenceMotifMover( *src.motif_mover_ ) );
	sequons_ = create_sequons();
}

////////////////////////////////////////////////////////////////////////////////
/// @brief Destructor (important for properly forward-declaring smart-pointer members)
CreateGlycanSequonMover::~CreateGlycanSequonMover(){}


////////////////////////////////////////////////////////////////////////////////
/// Mover Methods ///
/////////////////////

void
CreateGlycanSequonMover::set_glycosylation_position(core::Size position, core::pose::Pose const & pose ){
	utility::vector1<bool> subset(pose.size(), false);
	subset[position] = true;

	selector_ = ReturnResidueSubsetSelectorOP( new ReturnResidueSubsetSelector( subset ));
}

void
CreateGlycanSequonMover::set_glycan_sequon_type(protocols::carbohydrates::GlycanSequon sequon){
	sequon_type_ = sequon;
}


void
CreateGlycanSequonMover::set_use_basic_enhanced_n_linked_sequon(bool enhanced){
	if ( enhanced ) {
		sequon_type_ = n_linked_basic_enhanced;
	} else {
		sequon_type_ = n_linked_typical;
	}
}


void
CreateGlycanSequonMover::set_positions_as_start_of_sequon(bool positions_as_start){
	positions_as_start_of_sequon_ = positions_as_start;
}

void
CreateGlycanSequonMover::set_residue_selector(core::select::residue_selector::ResidueSelectorCOP selector){
	selector_ = selector->clone();
}

void
CreateGlycanSequonMover::set_pack_rounds( core::Size pack_rounds ){
	pack_rounds_ = pack_rounds;
}

void
CreateGlycanSequonMover::set_pack_neighbors(const bool pack_neighbors){
	pack_neighbors_ = pack_neighbors;
}

void
CreateGlycanSequonMover::set_design_neighbors(const bool design_neighbors){
	design_neighbors_ = design_neighbors;
}

void
CreateGlycanSequonMover::set_pack_distance(const core::Real pack_distance){
	pack_distance_ = pack_distance;
}

void
CreateGlycanSequonMover::set_design_x_positions(const bool design_x_positions){
	design_x_positions_ = design_x_positions;
}

////////////////////////////////////////////////////////////////////////////////
/// @brief Show the contents of the Mover
void
CreateGlycanSequonMover::show(std::ostream & output) const
{
	protocols::moves::Mover::show(output);
}

////////////////////////////////////////////////////////////////////////////////
/// Rosetta Scripts Support ///
///////////////////////////////

/// @brief parse XML tag (to use this Mover in Rosetta Scripts)
void
CreateGlycanSequonMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap& data,
	protocols::filters::Filters_map const & ,
	protocols::moves::Movers_map const & ,
	core::pose::Pose const & )
{

	if ( tag->hasOption("residue_selector") ) {
		set_residue_selector(parse_residue_selector(tag, data));
	} else {
		utility_exit_with_status("Must pass a residue_selector!");
	}

	set_design_x_positions(tag->getOption< bool >( "design_x_positions", false));

	//Scorefunction.
	if ( tag->hasOption("scorefxn") ) {
		set_score_function(protocols::rosetta_scripts::parse_score_function( tag, data ));
	}


	//Glycan Sequon Type.
	if ( tag->getOption<bool>("basic_enhanced_n_sequon", false) ) {
		sequon_type_ = n_linked_basic_enhanced;
	} else if ( tag->getOption<bool>("c_linked_NxC", false) ) {
		sequon_type_ = c_linked_NxC;
	} else if ( tag->getOption<bool>("c_linked_WxxW", false) ) {
		sequon_type_ = c_linked_WxxW;
	} else if ( tag->getOption<bool>("c_linked_WSTxC", false) ) {
		sequon_type_ = c_linked_WSTxC;
	} else {
		sequon_type_ = n_linked_typical;
	}

	set_pack_neighbors( tag->getOption<bool>("pack_neighbors", pack_neighbors_ ) );
	set_design_neighbors( tag->getOption<bool>("design_neighbors", design_neighbors_ ) );
	set_pack_distance( tag->getOption<core::Real>("neighbor_distance", pack_distance_ ) );
	set_pack_rounds( tag->getOption< core::Size >("pack_rounds", pack_rounds_ ));

}

void CreateGlycanSequonMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{

	using namespace utility::tag;
	AttributeList attlist;

	//here you should write code to describe the XML Schema for the class.  If it has only attributes, simply fill the probided AttributeList.

	std::string enhanced_references =
		"\tEnhanced Aromatic Sequons Increase Oligosaccharyltransferase Glycosylation Efficiency and Glycan Homogeneity\n"
		"    Murray et al., 2015, Chemistry and Biology 22, 1052 1062 http://dx.doi.org/10.1016/j.chembiol.2015.06.017\n"
		"\n"
		"   Residues Comprising the Enhanced Aromatic Sequon Influence Protein N-Glycosylation Efficiency\n"
		"    Yen-Wen Huang, J. Am. Chem. Soc. 2017, 139, 12947 12955 DOI: 10.1021/jacs.7b03868\n"
		"\n"
		" If this is in beginning or end of protein and could not be created, will set fail do not retry mover status. \n";

	attlist + XMLSchemaAttribute::attribute_w_default( "basic_enhanced_n_sequon", xsct_rosetta_bool, "Set to use the enhanced sequon which has been shown to result in higher rates of glycosylation.\n REFERENCES: \n" +enhanced_references, "false");

	attlist + XMLSchemaAttribute::attribute_w_default("c_linked_NxC", xsct_rosetta_bool,
		"Design in a c-linked sequon using NxC as the motif (Glycosylation at the N)", "false");

	attlist + XMLSchemaAttribute::attribute_w_default("c_linked_WxxW", xsct_rosetta_bool,
		"Design in a c-linked sequon using WxxW as the motif (Glycosylation at the first W)", "false");

	attlist + XMLSchemaAttribute::attribute_w_default("c_linked_WSTxC", xsct_rosetta_bool,
		"Design in a c-linked sequon using WSTxC as the motif (Glycosylation at W)", "false");

	attlist + XMLSchemaAttribute::attribute_w_default("positions_as_start_of_sequon", xsct_rosetta_bool, "Instead of having each position be the glycosylation position, have each position be the start of the sequon. This just makes using ResidueSelectors with this class a bit easier if needed.  Only really applies to enhanced motifs", "false");

	attlist + XMLSchemaAttribute::attribute_w_default("design_x_positions", xsct_rosetta_bool, "Design any X position (which can be ANY reside, instead of keeping it as the input residue type in the pose.", "false");

	attlist +XMLSchemaAttribute::attribute_w_default("pack_neighbors", xsct_rosetta_bool, "Should we pack the neighboring residues to the motif during design?", "true");
	//here you should write code to describe the XML Schema for the class.  If it has only attributes, simply fill the probided AttributeList.

	attlist +XMLSchemaAttribute::attribute_w_default("design_neighbors", xsct_rosetta_bool, "Should we design the neighboring residues to the motif during design?", "false");

	attlist +XMLSchemaAttribute::attribute_w_default("neighbor_distance", xsct_real, "Packing distance to neighbors?", "6.0");

	attlist + XMLSchemaAttribute::attribute_w_default("pack_rounds", xsct_non_negative_integer, "Number of rounds to run packing/design. ", "5");

	attributes_for_parse_residue_selector_default_option_name(attlist, "Set a number of positions using a residue selector.");
	attributes_for_get_score_function_name( attlist );

	std::string doc_string =
		" @brief Mutates residues to create a potential glycosylation site using known sequence motifs of N- or C- linked glycans.\n"
		"  Includes options for Enhanced Sequons for N-linked glycans that have been shown to have higher rates of glycosylation\n"
		"  as well as other positions that have been shown to influence the glycosylation chemistry.\n"
		"\n"
		" @details\n"
		"  Creates the glycan sequence motif around (and including) the potential glycosylation site.\n"
		"  If the site could not be created due to the position being too close to the beginning or end of the sequence,\n"
		"  Will set the mover status to fail, do not retry.\n"
		"\n"
		" Creates an N-Linked Sequence Motif by default using the non-enhanced motif.  This can be changed in options.\n";


	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), doc_string, attlist );
}


////////////////////////////////////////////////////////////////////////////////
/// @brief required in the context of the parser/scripting scheme
protocols::moves::MoverOP
CreateGlycanSequonMover::fresh_instance() const
{
	return protocols::moves::MoverOP( new CreateGlycanSequonMover );
}

/// @brief required in the context of the parser/scripting scheme
protocols::moves::MoverOP
CreateGlycanSequonMover::clone() const
{
	return protocols::moves::MoverOP( new CreateGlycanSequonMover( *this ) );
}

std::string CreateGlycanSequonMover::get_name() const {
	return mover_name();
}

std::string CreateGlycanSequonMover::mover_name() {
	return "CreateGlycanSequonMover";
}



/////////////// Creator ///////////////

protocols::moves::MoverOP
CreateGlycanSequonMoverCreator::create_mover() const
{
	return protocols::moves::MoverOP( new CreateGlycanSequonMover );
}

std::string
CreateGlycanSequonMoverCreator::keyname() const
{
	return CreateGlycanSequonMover::mover_name();
}

void CreateGlycanSequonMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	CreateGlycanSequonMover::provide_xml_schema( xsd );
}

void
CreateGlycanSequonMover::set_score_function(core::scoring::ScoreFunctionCOP scorefxn){
	scorefxn_ = scorefxn;
}

/// @brief Apply the mover
void
CreateGlycanSequonMover::apply( core::pose::Pose& pose ){

	if ( selector_ == nullptr ) {
		utility_exit_with_message("CreateGlycanSequonMover requires a residue selector!");
	}

	TR << "Applying... " << std::endl;
	std::string motif = sequons_[ sequon_type_ ];

	if ( design_x_positions_ ) {
		motif = utility::replace_in( motif, "-", "X");
	}

	//Make sure numbering is OK.  For basic_enhanced, we not at 1!
	if ( (sequon_type_ == n_linked_basic_enhanced || sequon_type_ == n_linked_best_enhanced ) && (! positions_as_start_of_sequon_) ) {

		utility::vector1< bool > subset = selector_->apply( pose );
		utility::vector1< bool > new_subset( pose.total_residue(), false );
		for ( core::Size i = 1; i <= subset.size(); ++i ) {
			if ( subset[i] ) {
				core::Size new_resnum = i - 2;
				if ( new_resnum >= 1 ) {
					new_subset[new_resnum] = true;
				} else {
					TR.Warning <<"Skipping glycosylation position "+utility::to_string( i )+" as the start of the motif is before the start of the protein!";
				}
			}
		}

		ReturnResidueSubsetSelectorOP subset_selector = ReturnResidueSubsetSelectorOP( new ReturnResidueSubsetSelector( new_subset ) );
		motif_mover_->set_residue_selector( subset_selector );
	} else {
		motif_mover_->set_residue_selector( selector_ );
	}

	motif_mover_->set_pack_neighbors( pack_neighbors_ );
	motif_mover_->set_neighbor_distance( pack_distance_ );
	motif_mover_->set_design_neighbors( design_neighbors_ );

	motif_mover_->set_motif( motif );
	motif_mover_->set_pack_rounds( pack_rounds_ );
	motif_mover_->set_score_function( scorefxn_ );

	motif_mover_->apply( pose );
}


///@brief Create a map of the name and full sequon
std::map< GlycanSequon, std::string >
create_sequons(){

	std::map< GlycanSequon, std::string > sequons;

	sequons[n_linked_typical] = "N[^P][ST]";
	sequons[n_linked_basic_enhanced] ="[%AROMATIC]-N[^P][ST]";
	sequons[c_linked_NxC] = "N-C";
	sequons[c_linked_WxxW] ="W--W";
	sequons[c_linked_WSTxC] = "W[ST]-C";

	TR << "Sequons created " << std::endl;
	return sequons;
}






////////////////////////////////////////////////////////////////////////////////
/// private methods ///
///////////////////////


std::ostream &
operator<<( std::ostream & os, CreateGlycanSequonMover const & mover )
{
	mover.show(os);
	return os;
}




} //protocols
} //carbohydrates
