// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/fold_from_loops/LabelPoseFromResidueSelectorMover.cc
/// @brief Adds labels to residues selected through a ResidueSelector
/// @author Jaume Bonet (jaume.bonet@gmail.com)

#include <protocols/fold_from_loops/LabelPoseFromResidueSelectorMover.hh>
#include <protocols/fold_from_loops/LabelPoseFromResidueSelectorMoverCreator.hh>

// Protocol headers
#include <protocols/rosetta_scripts/util.hh>

// Core headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/select/residue_selector/TrueResidueSelector.hh>
#include <core/select/residue_selector/NotResidueSelector.hh>
#include <core/select/residue_selector/util.hh>
#include <core/scoring/constraints/Constraint.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/pose/selection.hh>
#include <core/import_pose/import_pose.hh>

// Basic/Utility headers
#include <basic/datacache/DataMap.hh>
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <utility>
#include <utility/tag/Tag.hh>
#include <utility/vector1.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

static basic::Tracer TR( "protocols.fold_from_loops.LabelPoseFromResidueSelectorMover" );

namespace protocols {
namespace fold_from_loops {

LabelPoseFromResidueSelectorMover::LabelPoseFromResidueSelectorMover():
	protocols::moves::Mover( mover_name() ),
	selector_( default_selector() ),
	label_( default_label() ),
	reverse_( default_reverse() )
{}

LabelPoseFromResidueSelectorMover::LabelPoseFromResidueSelectorMover( core::select::residue_selector::ResidueSelectorCOP selector ):
	protocols::moves::Mover( mover_name() ),
	selector_(std::move( selector )),
	label_( default_label() ),
	reverse_( default_reverse() )
{}

LabelPoseFromResidueSelectorMover::LabelPoseFromResidueSelectorMover( core::select::residue_selector::ResidueSelectorCOP selector, std::string const & label ):
	protocols::moves::Mover( mover_name() ),
	selector_(std::move( selector )),
	label_( label ),
	reverse_( default_reverse() )
{}

LabelPoseFromResidueSelectorMover::~LabelPoseFromResidueSelectorMover()= default;

void
LabelPoseFromResidueSelectorMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data,
	protocols::filters::Filters_map const & ,
	protocols::moves::Movers_map const & ,
	core::pose::Pose const & )
{

	core::select::residue_selector::ResidueSelectorCOP selector = core::select::residue_selector::parse_residue_selector( tag, data );
	if ( selector ) residue_selector( selector );

	label( tag->getOption< std::string >( "property", default_label() ) );
	reverse( tag->getOption< bool >( "reverse", default_reverse() ) );

}

protocols::moves::MoverOP
LabelPoseFromResidueSelectorMover::clone() const
{
	return protocols::moves::MoverOP( new LabelPoseFromResidueSelectorMover( *this ) );
}

protocols::moves::MoverOP
LabelPoseFromResidueSelectorMover::fresh_instance() const
{
	return protocols::moves::MoverOP( new LabelPoseFromResidueSelectorMover );
}

void
LabelPoseFromResidueSelectorMover::apply( core::pose::Pose & pose )
{
	using namespace core::select::residue_selector;
	ResidueSubset const subset = selector_->apply( pose );
	if ( all_false_selection( subset ) ) return;
	runtime_assert_msg( not label_.empty(), "A label needs to be provided to insert!" );

	for ( core::Size i=1; i<=subset.size(); ++i ) {
		if ( subset[i] ) {
			if ( not reverse_ ) {
				pose.pdb_info()->add_reslabel(i, label_ );
			} else {
				if ( pose.pdb_info()->res_haslabel(i, label_ ) ) {
					utility::vector1< std::string> labels = pose.pdb_info()->get_reslabels( i );
					pose.pdb_info()->clear_reslabel( i );
					for ( auto const & label : labels ) {
						if ( label != label_ ) {
							pose.pdb_info()->add_reslabel(i, label );
						}
					}
				}
			}
		}
	}
}

void
LabelPoseFromResidueSelectorMover::label( std::string const & label )
{
	label_ = label;
}

void
LabelPoseFromResidueSelectorMover::reverse( bool pick )
{
	reverse_ = pick;
}

void
LabelPoseFromResidueSelectorMover::residue_selector( core::select::residue_selector::ResidueSelectorCOP const & selector )
{
	selector_ = selector;
}

void
LabelPoseFromResidueSelectorMover::residue_selector( core::select::residue_selector::ResidueSelector const & selector )
{
	selector_ = selector.clone();
}

std::string LabelPoseFromResidueSelectorMover::get_name() const {
	return mover_name();
}

std::string LabelPoseFromResidueSelectorMover::mover_name() {
	return "LabelPoseFromResidueSelector";
}

void LabelPoseFromResidueSelectorMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist
		+ XMLSchemaAttribute::required_attribute( "property", xs_string, "New label to add to the residues." )
		+ XMLSchemaAttribute::attribute_w_default( "reverse", xsct_rosetta_bool, "Remove the label (if any) instead of adding it", std::to_string( default_reverse() ) );
	core::select::residue_selector::attributes_for_parse_residue_selector_when_required( attlist, "residue_selector", "Selector specifying residues to be labeled." );
	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "Add labels to residues defined by ResidueSelector.", attlist );
}

core::select::residue_selector::ResidueSelectorCOP
LabelPoseFromResidueSelectorMover::default_selector()
{
	using namespace core::select::residue_selector;
	ResidueSelectorCOP true_selector( new TrueResidueSelector );
	ResidueSelectorCOP false_selector( new NotResidueSelector( true_selector ) );
	return false_selector;
}

std::string
LabelPoseFromResidueSelectorMover::default_label()
{
	return std::string();
}

bool
LabelPoseFromResidueSelectorMover::default_reverse()
{
	return false;
}

std::string LabelPoseFromResidueSelectorMoverCreator::keyname() const {
	return LabelPoseFromResidueSelectorMover::mover_name();
}

protocols::moves::MoverOP
LabelPoseFromResidueSelectorMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new LabelPoseFromResidueSelectorMover );
}

void LabelPoseFromResidueSelectorMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	LabelPoseFromResidueSelectorMover::provide_xml_schema( xsd );
}


} //protocols
} //fold_from_loops
