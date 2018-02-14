// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief
/// @author Carl Walkey ( cwalkey at uw dot edu )

// Unit Headers
#include <protocols/simple_moves/AddResidueLabelMover.hh>
#include <protocols/simple_moves/AddResidueLabelMoverCreator.hh>


// Core Headers
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/PDBInfo.hh>
#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/select/residue_selector/ResidueSelectorFactory.hh>
#include <core/select/residue_selector/TrueResidueSelector.hh>
#include <core/select/residue_selector/util.hh>

// Utility Headers
#include <basic/Tracer.hh>
#include <utility>
#include <utility/tag/Tag.hh>

// Protocols Headers
#include <protocols/moves/MoverContainer.hh>
#include <protocols/moves/mover_schemas.hh>
#include <utility/tag/XMLSchemaGeneration.hh>


static basic::Tracer TR( "protocols.simple_moves.AddResidueLabelMover" );

namespace protocols {
namespace simple_moves {

// @brief default constructor
AddResidueLabelMover::AddResidueLabelMover():
	protocols::moves::Mover(),
	selector_( new core::select::residue_selector::TrueResidueSelector ),
	label_( "" )
{}

AddResidueLabelMover::AddResidueLabelMover(
	core::select::residue_selector::ResidueSelectorCOP selector,
	std::string label ):
	protocols::moves::Mover(),
	selector_( selector ),
	label_( label )
{}

// @brief destructor
AddResidueLabelMover::~AddResidueLabelMover() = default;

void AddResidueLabelMover::apply(core::pose::Pose & pose) {

	runtime_assert( selector_ );
	if ( label_.empty() ) {
		utility_exit_with_message( "No label specified. A label must be specified" );
	}

	TR.Info << "Executing AddResidueLabelMover..." << std::endl;

	// Input: ResidueSelector defining residues to add labels to

	// Input: Label to add


	//core::select::residue_selector::ResidueSubsetCOP const subset( new core::select::residue_selector::ResidueSubset( selector_->apply( pose ) ) );
	core::select::residue_selector::ResidueSubset const subset = selector_->apply( pose );

	core::Size count=0;
	for ( core::Size resid=1; resid<=pose.size(); ++resid ) {
		if ( subset[resid] ) { // Add label if residue is in the selector
			TR.Info << "Adding to residue #" << std::to_string(resid) << ", label: " << label_ << std::endl;
			pose.pdb_info()->add_reslabel(resid, label_);
			++count;
		}
	}
	TR.Info << "Number of residues labeled: " << count << std::endl;

	//pymol selection output
	TR.Info << "select " << label_ << ", resi ";
	for ( core::Size resid=1; resid<=pose.size(); ++resid ) {
		if ( subset[resid] ) {
			TR.Info << std::to_string(resid) << "+";
		}
	} TR.Info << std::endl;

	//core::pose::PoseOP & target_;
	//std::string label;
	// i is the index of the ResidueSelector
	//pose->pdb_info()->add_reslabel(i, label)
	//
}

void AddResidueLabelMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	core::pose::Pose const & ) {

	TR.Info << "Reading XML parameters..." << std::endl;

	selector_ = core::select::residue_selector::parse_residue_selector( tag, data );
	label_ = tag->getOption< std::string >( "label");

	/*
	if( tag->hasOption("label") ) {
	} else {
	throw utility::excn::EXCN_RosettaScriptsOption("Must specify a residue-level label");
	}
	if ( selector ) {
	}
	if ( !selector_ ) {
	throw utility::excn::EXCN_RosettaScriptsOption( "CoordinateConstraintGenerator::parse_tag(): Error obtaining ResidueSelector from tag\n" );
	}*/
	/*
	if( tag->hasOption("selector") ) {
	selector_name = tag->getOption< std::string >( "selector" );
	try {
	selector = datamap.get_ptr< ResidueSelector const >( "ResidueSelector", selector_name );
	} catch ( utility::excn::EXCN_Msg_Exception & e ) {
	std::string error_message = "Failed to find ResidueSelector named '" + selector_name + "' from the Datamap from OperateOnResidueSubset::parse_tag\n" + e.msg();
	throw utility::excn::EXCN_Msg_Exception( error_message );
	}
	}*/
}

protocols::moves::MoverOP
AddResidueLabelMover::clone() const {
	return protocols::moves::MoverOP( new AddResidueLabelMover( *this ) );
}

protocols::moves::MoverOP
AddResidueLabelMover::fresh_instance() const {
	return protocols::moves::MoverOP( new AddResidueLabelMover );
}

std::string
AddResidueLabelMover::get_name() const {
	return "AddResidueLabelMover";
}

// @brief Identification
std::string
AddResidueLabelMoverCreator::keyname() const {
	return AddResidueLabelMoverCreator::mover_name();
}

std::string
AddResidueLabelMoverCreator::mover_name() {
	return "AddResidueLabel";
}

protocols::moves::MoverOP
AddResidueLabelMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new AddResidueLabelMover );
}

void
AddResidueLabelMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const {
	using namespace utility::tag;
	AttributeList attlist;
	core::select::residue_selector::attributes_for_parse_residue_selector(attlist);
	attlist + XMLSchemaAttribute::attribute_w_default( "label" , xs_string , "Label to add to residues defined in residue selector" , "" );
	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "Adds PDB label to residues", attlist );
}


}
}
