// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/simple_moves/CreateSequenceMotifMover.cc
/// @brief Create a sequence motif in a region of protein using the SequenceMotifTaskOperation.  Uses psueo-regular expressions to define the motif.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

// Unit headers
#include <protocols/simple_moves/CreateSequenceMotifMover.hh>
#include <protocols/simple_moves/CreateSequenceMotifMoverCreator.hh>
#include <protocols/toolbox/task_operations/SequenceMotifTaskOperation.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>

// Core headers
#include <core/pose/Pose.hh>
#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/select/residue_selector/util.hh>
#include <core/sequence/sequence_motif.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/task/operation/ResLvlTaskOperations.hh>
#include <core/pack/task/operation/OperateOnResidueSubset.hh>
#include <core/select/residue_selector/NeighborhoodResidueSelector.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

// Basic/Utility headers
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>
#include <basic/datacache/DataMap.hh>
#include <utility/string_util.hh>

// XSD Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>
#include <protocols/rosetta_scripts/util.hh>

static basic::Tracer TR( "protocols.simple_moves.CreateSequenceMotifMover" );

namespace protocols {
namespace simple_moves {

using namespace core::select::residue_selector;
using namespace protocols::simple_moves;
using namespace protocols::toolbox::task_operations;
using namespace core::pack::task::operation;
using namespace core::pack::task;
using namespace protocols::rosetta_scripts;


/////////////////////
/// Constructors  ///
/////////////////////

/// @brief Default constructor
CreateSequenceMotifMover::CreateSequenceMotifMover():
	protocols::moves::Mover( CreateSequenceMotifMover::mover_name() )
{

}

CreateSequenceMotifMover::CreateSequenceMotifMover( ResidueSelectorCOP selector ):
	protocols::moves::Mover( CreateSequenceMotifMover::mover_name() )
{
	set_residue_selector( selector );
}

////////////////////////////////////////////////////////////////////////////////
/// @brief Copy constructor
CreateSequenceMotifMover::CreateSequenceMotifMover( CreateSequenceMotifMover const & src ):
	protocols::moves::Mover( src.mover_name() ),
	motif_( src.motif_ ),
	pack_neighbors_( src.pack_neighbors_ ),
	design_neighbors_( src.design_neighbors_ ),
	pack_distance_( src.pack_distance_ ),
	pack_rounds_( src.pack_rounds_ )
{
	if ( src.selector_ ) {
		selector_ = src.selector_->clone();
	}
	if ( src.scorefxn_ ) {
		scorefxn_ = src.scorefxn_->clone();
	}

	neighbor_operations_.clear();
}

////////////////////////////////////////////////////////////////////////////////
/// @brief Destructor (important for properly forward-declaring smart-pointer members)
CreateSequenceMotifMover::~CreateSequenceMotifMover(){}

////////////////////////////////////////////////////////////////////////////////
/// Mover Methods ///
/////////////////////

void
CreateSequenceMotifMover::set_residue_selector(core::select::residue_selector::ResidueSelectorCOP selector){
	selector_ = selector->clone();
	neighbor_operations_.clear();
}

void
CreateSequenceMotifMover::set_score_function( core::scoring::ScoreFunctionCOP scorefxn ){
	scorefxn_ = scorefxn;
}

void
CreateSequenceMotifMover::set_design_neighbors(bool design_neighbors){
	design_neighbors_ = design_neighbors;
	neighbor_operations_.clear();
}

void
CreateSequenceMotifMover::set_pack_neighbors(bool pack_neighbors){
	pack_neighbors_ = pack_neighbors;
	neighbor_operations_.clear();
}

void
CreateSequenceMotifMover::set_motif(std::string motif){
	motif_ = motif;
	neighbor_operations_.clear();
}

void
CreateSequenceMotifMover::set_neighbor_distance(core::Real neighbor_distance){
	pack_distance_ = neighbor_distance;
	neighbor_operations_.clear();
}

void
CreateSequenceMotifMover::set_pack_rounds( core::Size pack_rounds ){
	pack_rounds_ = pack_rounds;
}

void
CreateSequenceMotifMover::initialize_objects( core::pose::Pose & pose){

	if ( ! selector_ ) {
		utility_exit_with_message("CreateSequenceMotifMover requires a residue selector! ");
	}


	//create TF if not present.
	if ( ! tf_ ) {
		tf_ = TaskFactoryOP( new TaskFactory());
	} else {
		tf_->clear();
	}

	//Create operations if not present.

	if ( ! cmd_line_operation_ ) {
		cmd_line_operation_ = InitializeFromCommandlineOP( new InitializeFromCommandline() );
	}

	if ( ! motif_operation_ ) {
		motif_operation_ = SequenceMotifTaskOperationOP( new SequenceMotifTaskOperation() );
	}

	if ( ! scorefxn_ ) {
		scorefxn_ = core::scoring::get_score_function();
	}


	tf_->push_back( cmd_line_operation_);

	//Restrict to our motif at each position of the selector, for the stretch.
	//Restrict restypes using the new task op.

	if ( neighbor_operations_.size() == 0 ) {
		OperateOnResidueSubsetOP operate_on_neighbors = OperateOnResidueSubsetOP( new OperateOnResidueSubset() );
		OperateOnResidueSubsetOP operate_on_others = OperateOnResidueSubsetOP( new OperateOnResidueSubset() );

		PreventRepackingRLTOP prevent_repacking = PreventRepackingRLTOP( new PreventRepackingRLT() );
		RestrictToRepackingRLTOP restrict_to_repacking = RestrictToRepackingRLTOP( new RestrictToRepackingRLT() );
		utility::vector1< bool > subset = selector_->apply( pose );
		utility::vector1< bool > focus( pose.total_residue(), false);

		//TR << "Subset " << utility::to_string( core::select::residue_selector::selection_positions(subset) );
		//Create the focus based on the motif and the residue selector.
		core::Size motif_length = core::sequence::get_motif_length( motif_ );
		for ( core::Size i = 1; i <= subset.size(); ++i ) {
			if ( subset[ i ] ) {
				for ( core::Size resnum = i; resnum <= i + motif_length - 1; ++resnum ) {
					focus[ resnum ] = true;
				}
			}
		}

		TR << "Focus " << utility::to_string( core::select::residue_selector::selection_positions(focus) );

		NeighborhoodResidueSelector nbr_selector = NeighborhoodResidueSelector();
		scorefxn_->score( pose ); //Needed for nbr selector to build neighbor graph.

		nbr_selector.set_focus( focus );
		nbr_selector.set_distance( pack_distance_ );

		nbr_selector.set_include_focus_in_subset(true);
		utility::vector1< bool > focus_and_neighbors = nbr_selector.apply( pose );

		nbr_selector.set_include_focus_in_subset( false );
		utility::vector1< bool > neighbors = nbr_selector.apply( pose );


		//Add subset operations to pack/design neighbors.

		//Close design and packing to focus and neighbors
		operate_on_others->subset( focus_and_neighbors );
		operate_on_others->flip_subset( true );
		operate_on_others->op( prevent_repacking  );


		operate_on_neighbors->subset( neighbors );

		if ( pack_neighbors_ ) {
			operate_on_neighbors->op( restrict_to_repacking );
		} else if ( ! design_neighbors_ ) {

			//Here we are not designing and not packing - so we turn it off.
			operate_on_neighbors->op( prevent_repacking );
		}

		//Add our operations to the TF!
		tf_->push_back( operate_on_others );
		tf_->push_back( operate_on_neighbors );

	}

	motif_operation_->set_residue_selector( selector_ );
	motif_operation_->set_motif( motif_ );
	tf_->push_back( motif_operation_ );

}

/// @brief Apply the mover
void
CreateSequenceMotifMover::apply( core::pose::Pose & pose){

	initialize_objects(pose);

	TR << "Applying... " << std::endl;

	PackRotamersMover packer = PackRotamersMover( scorefxn_ );
	packer.nloop( pack_rounds_ );
	packer.task_factory( tf_ );
	packer.apply( pose );

}

////////////////////////////////////////////////////////////////////////////////
/// @brief Show the contents of the Mover
void
CreateSequenceMotifMover::show(std::ostream & output) const
{
	protocols::moves::Mover::show(output);
}

////////////////////////////////////////////////////////////////////////////////
/// Rosetta Scripts Support ///
///////////////////////////////

/// @brief parse XML tag (to use this Mover in Rosetta Scripts)
void
CreateSequenceMotifMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap& data,
	protocols::filters::Filters_map const & ,
	protocols::moves::Movers_map const & ,
	core::pose::Pose const & )
{
	if ( tag->hasOption("residue_selector") ) {
		selector_ = parse_residue_selector(tag, data);
	} else {
		utility_exit_with_status("residue_selector tag Required");
	}

	//Scorefunction.
	if ( tag->hasOption("scorefxn") ) {
		scorefxn_ = protocols::rosetta_scripts::parse_score_function( tag, data );
	}

	motif_ = tag->getOption<std::string>("motif");
	pack_neighbors_ = tag->getOption<bool>("pack_neighbors", pack_neighbors_);
	design_neighbors_ = tag->getOption<bool>("design_neighbors", design_neighbors_);
	pack_distance_ = tag->getOption<core::Real>("neighbor_distance", pack_distance_);
	pack_rounds_ = tag->getOption< core::Size >("pack_rounds", pack_rounds_);

}

////////////////////////////////////////////////////////////////////////////////
/// @brief required in the context of the parser/scripting scheme
protocols::moves::MoverOP
CreateSequenceMotifMover::fresh_instance() const
{
	return protocols::moves::MoverOP( new CreateSequenceMotifMover );
}

/// @brief required in the context of the parser/scripting scheme
protocols::moves::MoverOP
CreateSequenceMotifMover::clone() const
{
	return protocols::moves::MoverOP( new CreateSequenceMotifMover( *this ) );
}

std::string CreateSequenceMotifMover::get_name() const {
	return mover_name();
}

std::string CreateSequenceMotifMover::mover_name() {
	return "CreateSequenceMotifMover";
}

void CreateSequenceMotifMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{

	using namespace utility::tag;
	AttributeList attlist;

	std::string motif_str = core::sequence::get_design_sequence_motif_syntax();
	attlist + XMLSchemaAttribute("motif", xs_string, motif_str);

	attlist +XMLSchemaAttribute::attribute_w_default("pack_neighbors", xsct_rosetta_bool, "Should we pack the neighboring residues to the motif during design?", "true");
	//here you should write code to describe the XML Schema for the class.  If it has only attributes, simply fill the probided AttributeList.

	attlist +XMLSchemaAttribute::attribute_w_default("design_neighbors", xsct_rosetta_bool, "Should we design the neighboring residues to the motif during design?", "false");

	attlist +XMLSchemaAttribute::attribute_w_default("neighbor_distance", xsct_real, "Packing distance to neighbors?", "6.0");

	attlist + XMLSchemaAttribute::attribute_w_default("pack_rounds", xsct_non_negative_integer, "Number of rounds to run packing/design. ", "5");

	std::string doc_string = "@brief Simple mover to Create a sequence motif in a region of protein using the\n" "SequenceMotifTaskOperation.  Uses psueo-regular expressions to define the motif.\n"
		"\n"
		"@details"
		" Simply calls the packer using the operation, optionally packing neighbor residues as it does so.\n"
		" If you need something more complex, use the SequenceMotifTaskOperation directly.\n";

	attributes_for_parse_residue_selector( attlist );
	attributes_for_get_score_function_name( attlist );

	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), doc_string, attlist );
}

/////////////// Creator ///////////////

protocols::moves::MoverOP
CreateSequenceMotifMoverCreator::create_mover() const
{
	return protocols::moves::MoverOP( new CreateSequenceMotifMover );
}

std::string
CreateSequenceMotifMoverCreator::keyname() const
{
	return CreateSequenceMotifMover::mover_name();
}

void CreateSequenceMotifMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	CreateSequenceMotifMover::provide_xml_schema( xsd );
}

////////////////////////////////////////////////////////////////////////////////
/// private methods ///
///////////////////////


std::ostream &
operator<<( std::ostream & os, CreateSequenceMotifMover const & mover )
{
	mover.show(os);
	return os;
}

} //protocols
} //simple_moves
